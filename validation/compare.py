#!/usr/bin/env python
"""Cross-validation: PhyloSlicer vs treesliceR reference values.

Reads the exact inputs exported by ``validation/r/generate_reference.R`` (Newick
tree + presence matrix + queen adjacency per scenario) and the treesliceR
reference CSVs, recomputes everything with phyloslicer and reports equivalence.

Deterministic quantities (per-slice PD/PE/PB profiles, whole-tree index values)
use rtol 1e-8; optimiser outputs (rates) are compared at rtol 1e-4 because
treesliceR's single-start ``nls`` and phyloslicer's bounded multi-start least
squares are different algorithms (see validation/NOTES.md).

The beta families (``CpB``, ``CpB_RW``, ``r_phylo(index = "PB")``) are scored
only for scenarios that ship an adjacency reference, because a beta neighbourhood
set cannot be reconstructed from a tree and a matrix alone.  Where the shipped
rate and the reference disagree, the row also carries a *control* fit that uses
the reference implementation's own curve normalisation
(``cumsum(series) / total``, no re-scaling of the profile to sum 1); a row whose
shipped value fails while the control passes is reported as
``documented divergence`` - the disagreement is attributable to the profile
renormalisation recorded in validation/NOTES.md, which is a pinned design decision, rather than
to the arithmetic.

Provenance is a precondition, not a courtesy (see validation/NOTES.md).  Every
reference file must be registered in ``reference/manifest.txt`` with its SHA256,
its row/column counts, its column-name pattern and its orientation:

* the run refuses to start when the manifest is missing or does not cover a file
  it is about to read (or when the bytes differ from the registered hash);
* a reference whose shape differs from the registered one is a hard error.  The
  comparator used to transpose such a file silently, which is exactly how a
  slice-by-site CSV and a site-by-slice CSV became interchangeable;
* ``pDO`` / ``pEO`` / ``pBO`` NA *patterns* are compared element by element, and
  the reference ``PD`` / ``PE`` / ``PB`` columns are compared numerically - the
  old code only counted non-finite values on each side, which any two files pass.

Files the manifest marks ``status = unverified`` (the committed per-slice PD/PE
profiles, whose layout contradicts the generator that was supposed to have
written them - see the manifest header) are hash- and shape-checked but *not*
numerically cross-checked; the skipped rows are reported as such.  With
``--strict`` they fail the run instead, along with every ``documented divergence``
and every family for which the reference implementation produces no file.
"""

# The thread pin below is a statement, so every import after it is "not at the
# top of the file" by ruff's counting.  That ordering is the point: it has to
# happen before numpy binds its BLAS back end, so it cannot move below the
# imports it is there to condition.
# ruff: noqa: E402

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

# Pin the BLAS/thread-pool environment before numpy loads its back ends.  With
# the default thread counts this run is *not* byte-reproducible: measured on the
# shipped references, 7 of the 162 report rows wobble in their last significant
# digits between two consecutive runs (a per-slice profile error moving
# 9.161e-14 -> 9.141e-14, a fitted rate moving 1.47158404e-06 ->
# 1.47158405e-06), because multi-threaded reductions sum in a
# scheduling-dependent order.  No status changes, and every wobble is four
# orders of magnitude inside tolerance, but a validation report that is quoted
# digit by digit in the documentation has to be a function of its inputs alone -
# and anyone re-running it must land on the same bytes, not on "close enough".
# Set, not forced: an outer harness that already pinned the threads wins.
for _var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_var, "1")

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

from phyloslicer.cli.main import (
    MANIFEST_NAME,
    ReferenceProvenanceError,
    check_reference_columns,
    check_reference_shape,
    load_reference_manifest,
    reference_provenance,
    sha256_file,
)
from phyloslicer.indices import pd_per_slice, pe_per_slice, r_compatible_profile
from phyloslicer.io.matrices import align_tree_matrix, as_matrix
from phyloslicer.io.newick import read_newick
from phyloslicer.rates.api import cpb_rate, cpd_rate, cpe_rate
from phyloslicer.rates.fitting import adaptive_bounds
from phyloslicer.slicing import slice_pieces

TOL_PROFILE = 1e-8
TOL_RATE = 1e-4

#: Magnitude below which a reference entry counts as "literally zero".  A
#: per-slice beta profile *is* 0 in the oldest slices (no lineage pair has
#: diverged yet), and a relative error against an exact zero is undefined;
#: scoring those cells as failures would mark a bit-for-bit agreement as a
#: mismatch.  At such entries the comparison falls back to absolute error.
ZERO_FLOOR = 1e-12

BETA_COMPONENTS = ("sorensen", "turnover", "nestedness")
BETA_APPROACHES = ("multisite", "pairwise")


def _beta_rate_families():
    """(file suffix, component, approach, weighted, rate column, total column).

    ``CpB`` is spelled per component x approach; ``CpB_RW`` has no component
    argument in the reference (``CpB_RW.R:50``), so its files carry the approach
    only and are scored as the Sorensen variant it computes.
    """
    families = [
        (f"CpB__{comp}_{appr}", comp, appr, False, "CpB", "PB")
        for appr in BETA_APPROACHES
        for comp in BETA_COMPONENTS
    ]
    families += [
        (f"CpB_RW__{appr}", "sorensen", appr, True, "CpB_RW", "PB_RW") for appr in BETA_APPROACHES
    ]
    return tuple(families)


BETA_RATE_FAMILIES = _beta_rate_families()
BETA_PROFILE_FAMILIES = tuple(
    (f"rPB__{comp}_{appr}", comp, appr) for appr in BETA_APPROACHES for comp in BETA_COMPONENTS
)

#: Families whose *fitted curve* in treesliceR is not the one ``method`` names.
#: For ``comp = "turnover"`` the pairwise branch (CpB.R:605,
#: ``mean(apply(pw_piece - t(combn(piece, 2)), 1, min)) / turn`` with
#: ``turn <- mean(min_bc)``) reduces algebraically to
#: ``sum_pairs(min_slice) / sum_pairs(min_full)`` - word for word the multisite
#: branch of CpB.R:609.  The reference's turnover *rate* is therefore independent
#: of ``method`` (measured on the archived scenarios: its turnover_multisite and
#: turnover_pairwise ``CpB`` columns agree to 1.2e-09) and only its reported
#: ``PB`` differs.  Phyloslicer averages the per-pair ratios instead, which is the
#: convention its own ``CpB_RW`` branch uses, so the control for this family
#: refits the reference's curve rather than ours.
CURVE_CONTROL = {"CpB__turnover_pairwise": "multisite"}

#: every beta row label this comparator can emit, used to bucket the report
BETA_LABELS = tuple([f[0] for f in BETA_RATE_FAMILIES] + [f[0] for f in BETA_PROFILE_FAMILIES])

#: Families for which the *reference implementation* produces no file, with the
#: failure measured by ``validation/r/generate_reference.R`` (R 4.5.3 / ape 5.8.1
#: / treesliceR 1.1.0, 2026-09-21).  They are reported as ``no reference`` and
#: must never be counted towards an equivalence claim.
NO_REFERENCE = {
    "CpB__sorensen_pairwise": (
        "treesliceR::CpB(comp='sorensen', method='pairwise') aborts in nls:"
        " \"'qr' and 'y' must have the same number of rows\" (its pairwise branch"
        " recycles a length-P per-pair Sorensen vector against a length-k slice"
        " vector, CpB.R:479 / CpB.R:525)"
    ),
    "rPB__sorensen_pairwise": (
        "treesliceR::r_phylo(index='PB', method='pairwise') returns a ragged list"
        " (element lengths 300/750/1800 for 50 slices x 6-36 pairs), not one"
        " vector per assemblage, so there is no site-by-slice reference"
    ),
    "rPB__turnover_pairwise": (
        "treesliceR::r_phylo(index='PB', method='turnover') returns a ragged list"
        " (element lengths 300/750/1800 for 50 slices x 6-36 pairs), not one"
        " vector per assemblage, so there is no site-by-slice reference"
    ),
}


def rel_err(mine: np.ndarray, ref: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        out = np.abs(mine - ref) / np.abs(ref)
    return np.where(np.isfinite(out), out, np.nan)


def count_identical(mine: np.ndarray, ref: np.ndarray) -> int:
    """Cells that agree *exactly*, not merely within tolerance.

    ``max_rel_err`` says how far apart the two implementations are and
    ``pass_rate`` says how much of the array is inside tolerance; neither says
    how much of it is bit-for-bit the same number, which is the part a reader
    uses to judge how much of the agreement is floating-point luck.  Counted on
    the intersection of the finite cells, the same set ``compare_pair`` scores.
    """
    finite = np.isfinite(mine) & np.isfinite(ref)
    if not finite.any():
        return 0
    return int(np.count_nonzero(mine[finite] == ref[finite]))


def compare_pair(
    a: np.ndarray, b: np.ndarray, tol: float, *, zero_floor: float = 0.0
) -> tuple[float, float]:
    """Max relative error and fraction of entries inside ``tol``.

    ``zero_floor`` switches the comparison to absolute error where the reference
    is numerically zero (see :data:`ZERO_FLOOR`); without it those entries yield
    ``NaN`` and a NaN never satisfies ``<= tol``, so an exactly-reproduced zero
    would be counted as a mismatch.
    """
    finite = np.isfinite(a) & np.isfinite(b)
    if not finite.any():
        return float("nan"), 0.0
    a, b = a[finite], b[finite]
    if zero_floor:
        at_zero = np.abs(b) <= zero_floor
        with np.errstate(divide="ignore", invalid="ignore"):
            errs = np.where(at_zero, np.abs(a - b), np.abs(a - b) / np.abs(b))
    else:
        errs = rel_err(a, b)
    return float(np.nanmax(errs)), float(np.mean(errs <= tol))


def fit_reference_convention(series: np.ndarray, total: np.ndarray, ages: np.ndarray) -> np.ndarray:
    """Re-fit each row the way treesliceR's ``nls`` does, on its own curve.

    treesliceR fits ``cumsum(series) / total ~ exp(-r * age)`` (CpB.R:392), while
    :func:`phyloslicer.rates.cpb_rate` fits the *profile* normalisation
    ``series / series.sum()``.  The model has no amplitude parameter, so the two
    curves - which differ only by the constant factor ``total / series.sum()`` -
    have different least-squares optima whenever that factor is not 1.  This is
    the A2 divergence, and fitting the reference's own objective is what
    separates "the implementations disagree" from "the implementations disagree
    *because of* the documented normalisation".
    """
    series = np.atleast_2d(np.asarray(series, dtype=float))
    total = np.asarray(total, dtype=float)
    ages = np.asarray(ages, dtype=float)
    lo, hi = adaptive_bounds(ages)
    out = np.full(series.shape[0], np.nan)
    for i in range(series.shape[0]):
        row, tot = series[i], total[i]
        if not np.isfinite(row).all() or not np.isfinite(tot) or tot == 0:
            continue
        target = np.cumsum(row) / tot
        if not np.isfinite(target).all():
            continue

        def sse(r: float, target=target) -> float:
            resid = target - np.exp(-r * ages)
            return float(resid @ resid)

        best = None
        # three disjoint starts across the band, so a bounded search cannot be
        # trapped in a local basin of this one-parameter problem
        for start in (lo + (hi - lo) * f for f in (0.05, 0.35, 0.8)):
            cand = minimize_scalar(
                sse,
                bounds=(max(lo, start * 0.2), min(hi, start * 5.0)),
                method="bounded",
                options={"xatol": 1e-12},
            )
            if cand.success and (best is None or cand.fun < best.fun):
                best = cand
        if best is not None:
            out[i] = best.x
    return out


def na_pattern_report(a: np.ndarray, b: np.ndarray, *, label: str) -> tuple[str, bool]:
    """Element-wise NA-pattern agreement: the count-equality check it replaced
    passed for any pair of files, because both counts are 0 in particular when
    one side is NA everywhere and the other nowhere."""
    na_mine = ~np.isfinite(a)
    na_ref = ~np.isfinite(b)
    diff = np.flatnonzero(na_mine != na_ref)
    if diff.size:
        shown = ", ".join(f"{i}:{'mine' if na_mine[i] else 'ref'}" for i in diff[:8])
        return f"{label}: NA pattern differs at {diff.size} row(s) [{shown}]", False
    return f"{label}: NA patterns identical ({int(na_ref.sum())} NA rows)", True


def check_registration(ref_dir: Path, records: dict[str, dict[str, str]]) -> list[str]:
    """Return provenance problems for every file this run is about to read."""
    problems: list[str] = []
    for path in sorted(ref_dir.iterdir()):
        if path.name == MANIFEST_NAME or path.name == "equivalence_report.csv":
            continue
        if path.suffix not in (".csv", ".tre"):
            continue
        entry = records.get(path.name)
        if entry is None:
            problems.append(f"{path.name}: no entry in {MANIFEST_NAME}")
            continue
        if not entry.get("sha256"):
            problems.append(f"{path.name}: manifest entry carries no sha256")
            continue
        actual = sha256_file(path)
        if actual != entry["sha256"]:
            problems.append(
                f"{path.name}: bytes {actual[:12]}… differ from registered {entry['sha256'][:12]}…"
            )
    return problems


def read_profile_reference(
    name: str, ref_dir: Path, entry: dict[str, str]
) -> tuple[np.ndarray, str]:
    """Load a per-slice profile reference, returning (array, status note).

    Never transposes.  A shape that disagrees with the manifest (or with what
    phyloslicer computes) raises: that is the only way a mis-oriented reference
    stays visible instead of becoming a spurious 1e-14 agreement.
    """
    path = ref_dir / name
    frame = pd.read_csv(path, index_col=0)
    check_reference_shape(entry, frame.shape, name=name)
    check_reference_columns(entry, frame.columns, name=name)
    if entry.get("status") == "unverified":
        return (
            np.full((0, 0), np.nan),
            "declared unverified in the manifest (orientation="
            + entry["orientation"]
            + ", shape "
            + f"{entry['rows']}x{entry['cols']}): "
            + entry.get("reason", "no reason recorded"),
        )
    return frame.to_numpy(dtype=float), ""


def main(ref_dir: str = "validation/reference", strict: bool = False, out=None) -> int:
    ref_dir = Path(ref_dir)
    manifest_path = ref_dir / MANIFEST_NAME
    if not manifest_path.is_file():
        print(
            f"error: no provenance manifest at {manifest_path}; refusing to run.\n"
            "Without it a 'reference' cannot be told apart from the output of the "
            "code under test.  Register every reference file by "
            "SHA256, shape, column pattern and orientation before comparing.",
            file=sys.stderr,
        )
        return 2
    try:
        _header, records, _used = load_reference_manifest(ref_dir, manifest_path)
    except ReferenceProvenanceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    problems = check_registration(ref_dir, records)
    if problems:
        print(
            "error: reference files are not properly registered, refusing to run:\n  "
            + "\n  ".join(problems),
            file=sys.stderr,
        )
        return 2

    rows = []
    for tree_path in sorted(ref_dir.glob("*__tree.tre")):
        tag = tree_path.name.replace("__tree.tre", "")
        tree = read_newick(tree_path)
        for mat_path in sorted(ref_dir.glob(f"{tag}__*__mat.csv")):
            scenario = mat_path.name.replace("__mat.csv", "")
            mat = pd.read_csv(mat_path, index_col=0)
            mat.columns = mat.columns.astype(str)
            mat = (mat > 0).astype(float)

            tree_a, mat_a = align_tree_matrix(tree, as_matrix(mat))
            slices = slice_pieces(tree_a, n=50, criterion="time")

            # --- deterministic profiles (rtol 1e-8) ---------------------- #
            for index, ref_suffix, fn in (
                ("PD", "rPD", pd_per_slice),
                ("PE", "rPE", pe_per_slice),
            ):
                name = f"{scenario}__{ref_suffix}.csv"
                entry = records[name]
                mine = fn(tree_a, mat_a, slices)
                try:
                    ref_arr, note = read_profile_reference(name, ref_dir, entry)
                except ReferenceProvenanceError as exc:
                    rows.append(
                        {
                            "scenario": f"{scenario}:{index}",
                            "n": mine.size,
                            "max_rel_err": np.nan,
                            "pass_rate": 0.0,
                            "tol": TOL_PROFILE,
                            "status": "provenance failure",
                            "note": str(exc),
                        }
                    )
                    continue
                if ref_arr.size == 0:  # declared unverified
                    rows.append(
                        {
                            "scenario": f"{scenario}:{index}",
                            "n": mine.size,
                            "max_rel_err": np.nan,
                            "pass_rate": np.nan,
                            "tol": TOL_PROFILE,
                            "status": "not evaluated",
                            "ref_basis": f"{name} (declared unverified in the manifest)",
                            "note": "hash + shape verified; numeric cross-check skipped - " + note,
                        }
                    )
                    continue
                if mine.shape != ref_arr.shape:
                    # no transpose: report it as the orientation error it is
                    rows.append(
                        {
                            "scenario": f"{scenario}:{index}",
                            "n": mine.size,
                            "max_rel_err": np.nan,
                            "pass_rate": 0.0,
                            "tol": TOL_PROFILE,
                            "status": "orientation failure",
                            "note": f"phyloslicer produced {mine.shape[0]}x{mine.shape[1]} "
                            f"but the registered reference is "
                            f"{ref_arr.shape[0]}x{ref_arr.shape[1]}",
                        }
                    )
                    continue
                worst, rate = compare_pair(mine, ref_arr, TOL_PROFILE)
                rows.append(
                    {
                        "scenario": f"{scenario}:{index}",
                        "n": mine.size,
                        "max_rel_err": worst,
                        "pass_rate": rate,
                        "tol": TOL_PROFILE,
                        "status": "verified",
                        "ref_basis": name,
                        "identical_cells": count_identical(mine, ref_arr),
                        "note": "",
                    }
                )

            # --- optimiser outputs (rtol 1e-4) --------------------------- #
            for kind, runner, rate_col, total_col, origin_col in (
                ("CpD", cpd_rate, "CpD", "PD", "pDO"),
                ("CpE", cpe_rate, "CpE", "PE", "pEO"),
            ):
                name = f"{scenario}__{kind}.csv"
                try:
                    entry = reference_provenance(ref_dir / name, manifest_path)
                    ref = pd.read_csv(ref_dir / name, index_col=0)
                    check_reference_shape(entry, ref.shape, name=name)
                    check_reference_columns(entry, ref.columns, name=name)
                    if entry.get("rate_column", rate_col) != rate_col:
                        raise ReferenceProvenanceError(
                            f"{name}: manifest registers rate_column="
                            f"{entry.get('rate_column')!r}, this comparison needs "
                            f"{rate_col!r}"
                        )
                except ReferenceProvenanceError as exc:
                    rows.append(
                        {
                            "scenario": f"{scenario}:{kind}",
                            "n": 0,
                            "max_rel_err": np.nan,
                            "pass_rate": 0.0,
                            "tol": TOL_RATE,
                            "status": "provenance failure",
                            "note": str(exc),
                        }
                    )
                    continue
                missing = [c for c in (rate_col, total_col, origin_col) if c not in ref]
                if missing:
                    rows.append(
                        {
                            "scenario": f"{scenario}:{kind}",
                            "n": len(ref),
                            "max_rel_err": np.nan,
                            "pass_rate": 0.0,
                            "tol": TOL_RATE,
                            "status": "provenance failure",
                            "note": f"reference lacks column(s) {missing}",
                        }
                    )
                    continue
                mine = runner(tree_a, mat_a, n_slices=50)
                a_r = mine[rate_col].to_numpy(dtype=float)
                b_r = ref[rate_col].to_numpy(dtype=float)
                a_t = mine[total_col].to_numpy(dtype=float)
                b_t = ref[total_col].to_numpy(dtype=float)
                a_o = mine[origin_col].to_numpy(dtype=float)
                b_o = ref[origin_col].to_numpy(dtype=float)
                worst, rate = compare_pair(a_r, b_r, TOL_RATE)
                worst_t, rate_t = compare_pair(a_t, b_t, TOL_PROFILE)
                worst_o, rate_o = compare_pair(a_o, b_o, TOL_RATE)
                notes = []
                passed_all = rate >= 0.95 and rate_t >= 0.99 and rate_o >= 0.95
                for label, arr_mine, arr_ref in (
                    (rate_col, a_r, b_r),
                    (total_col, a_t, b_t),
                    (origin_col, a_o, b_o),
                ):
                    note, agree = na_pattern_report(arr_mine, arr_ref, label=label)
                    notes.append(note)
                    passed_all = passed_all and agree
                notes.append(f"{total_col} max rel err {worst_t:.3g} (tol {TOL_PROFILE:g})")
                notes.append(f"{origin_col} max rel err {worst_o:.3g} (tol {TOL_RATE:g})")
                rows.append(
                    {
                        "scenario": f"{scenario}:{kind}",
                        "n": a_r.size,
                        "max_rel_err": worst,
                        "pass_rate": rate if passed_all else min(rate, rate_t, rate_o),
                        "tol": TOL_RATE,
                        "status": "verified" if passed_all else "FAILED",
                        "ref_basis": f"{name}:{rate_col}",
                        "note": "; ".join(notes),
                    }
                )
                # The whole-tree totals that treesliceR prints alongside its rates are a
                # *separately registered* reference, and they are the only certified basis on
                # which this build can check the deterministic per-site diversity value.  The
                # unverified r-profile files are skipped, so per-slice
                # agreement has no certified reference here and must not be reported as if it
                # did.  Recording the totals as their own rows keeps that distinction in the
                # machine-readable report instead of hiding it inside a free-text note.
                rows.append(
                    {
                        "scenario": f"{scenario}:{total_col}_total",
                        "n": a_t.size,
                        "max_rel_err": worst_t,
                        "pass_rate": rate_t,
                        "tol": TOL_PROFILE,
                        "status": "verified" if rate_t >= 0.99 else "FAILED",
                        "ref_basis": f"{name}:{total_col}",
                        "note": "per-site whole-tree value reported by the reference "
                        f"alongside its {rate_col} fit; not a per-slice profile",
                    }
                )

            # --- beta diversity: CpB / CpB_RW / r_phylo(index = "PB") ------- #
            # A beta reference is only meaningful next to the adjacency matrix
            # the neighbourhoods were read from, so a scenario without one simply
            # contributes no beta rows (rather than a vacuous agreement).
            adj_name = f"{scenario}__adj.csv"
            if not (ref_dir / adj_name).is_file():
                continue
            adj = None
            try:
                # an adjacency matrix is an input, not a fitted output, so it is
                # registered without a rate_column
                adj_entry = reference_provenance(
                    ref_dir / adj_name, manifest_path, require_rate_column=False
                )
                adj_frame = pd.read_csv(ref_dir / adj_name, index_col=0)
                check_reference_shape(adj_entry, adj_frame.shape, name=adj_name)
                check_reference_columns(adj_entry, adj_frame.columns, name=adj_name)
                adj = adj_frame.to_numpy(dtype=float)
            except ReferenceProvenanceError as exc:
                rows.append(
                    {
                        "scenario": f"{scenario}:adjacency",
                        "n": 0,
                        "max_rel_err": np.nan,
                        "pass_rate": 0.0,
                        "tol": TOL_PROFILE,
                        "status": "provenance failure",
                        "ref_basis": "none",
                        "note": f"cannot load the adjacency reference: {exc}",
                    }
                )
            if adj is None:
                continue

            for suffix, comp, appr, weighted, rate_col, total_col in BETA_RATE_FAMILIES:
                name = f"{scenario}__{suffix}.csv"
                label = f"{scenario}:{suffix}"
                if not (ref_dir / name).is_file():
                    rows.append(
                        {
                            "scenario": label,
                            "n": 0,
                            "max_rel_err": np.nan,
                            "pass_rate": np.nan,
                            "tol": TOL_RATE,
                            "status": "no reference",
                            "ref_basis": "none",
                            "note": NO_REFERENCE.get(
                                suffix,
                                "the reference generator produced no file for this family",
                            ),
                        }
                    )
                    continue
                try:
                    entry = reference_provenance(ref_dir / name, manifest_path)
                    ref = pd.read_csv(ref_dir / name, index_col=0)
                    check_reference_shape(entry, ref.shape, name=name)
                    check_reference_columns(entry, ref.columns, name=name)
                    if entry.get("rate_column", rate_col) != rate_col:
                        raise ReferenceProvenanceError(
                            f"{name}: manifest registers rate_column="
                            f"{entry.get('rate_column')!r}, this comparison needs "
                            f"{rate_col!r}"
                        )
                except ReferenceProvenanceError as exc:
                    rows.append(
                        {
                            "scenario": label,
                            "n": 0,
                            "max_rel_err": np.nan,
                            "pass_rate": 0.0,
                            "tol": TOL_RATE,
                            "status": "provenance failure",
                            "ref_basis": "none",
                            "note": str(exc),
                        }
                    )
                    continue
                missing = [c for c in (rate_col, total_col, "pBO") if c not in ref]
                if missing:
                    rows.append(
                        {
                            "scenario": label,
                            "n": len(ref),
                            "max_rel_err": np.nan,
                            "pass_rate": 0.0,
                            "tol": TOL_RATE,
                            "status": "provenance failure",
                            "ref_basis": "none",
                            "note": f"reference lacks column(s) {missing}",
                        }
                    )
                    continue
                mine = cpb_rate(
                    tree_a,
                    mat_a,
                    adj,
                    n_slices=50,
                    component=comp,
                    approach=appr,
                    weighted=weighted,
                )
                profile = r_compatible_profile(
                    tree_a,
                    mat_a,
                    adj,
                    slices,
                    component=comp,
                    approach=appr,
                    weighted=weighted,
                )
                a_r = mine[rate_col].to_numpy(dtype=float)
                b_r = ref[rate_col].to_numpy(dtype=float)
                a_t = mine[total_col].to_numpy(dtype=float)
                b_t = ref[total_col].to_numpy(dtype=float)
                a_o = mine["pBO"].to_numpy(dtype=float)
                b_o = ref["pBO"].to_numpy(dtype=float)
                worst, rate = compare_pair(a_r, b_r, TOL_RATE)
                worst_t, rate_t = compare_pair(a_t, b_t, TOL_PROFILE)
                worst_o, rate_o = compare_pair(a_o, b_o, TOL_RATE)
                # the control: our own arithmetic, the reference's own curve
                control_approach = CURVE_CONTROL.get(suffix, appr)
                control_profile = (
                    profile
                    if control_approach == appr
                    else r_compatible_profile(
                        tree_a,
                        mat_a,
                        adj,
                        slices,
                        component=comp,
                        approach=control_approach,
                        weighted=weighted,
                    )
                )
                control = fit_reference_convention(
                    control_profile["series"], control_profile["total"], slices.ages
                )
                worst_c, rate_c = compare_pair(control, b_r, TOL_RATE)
                control_how = (
                    "same curve, the reference's normalisation"
                    if control_approach == appr
                    else f"the reference's method-independent '{control_approach}' "
                    "curve (CpB.R:605 reduces to CpB.R:609), its normalisation"
                )
                with np.errstate(divide="ignore", invalid="ignore"):
                    amplitude = np.abs(profile["series"].sum(axis=1) / profile["total"] - 1.0)
                amp_max = float(np.nanmax(amplitude)) if np.isfinite(amplitude).any() else np.nan

                notes = []
                na_ok = True
                for col, arr_mine, arr_ref in (
                    (rate_col, a_r, b_r),
                    (total_col, a_t, b_t),
                    ("pBO", a_o, b_o),
                ):
                    note, agree = na_pattern_report(arr_mine, arr_ref, label=col)
                    notes.append(note)
                    na_ok = na_ok and agree
                notes.append(f"{total_col} max rel err {worst_t:.3g} (tol {TOL_PROFILE:g})")
                notes.append(f"pBO max rel err {worst_o:.3g} (tol {TOL_RATE:g})")
                notes.append(
                    f"control fit ({control_how}): max rel err {worst_c:.3g}, "
                    f"pass rate {rate_c:.3g}"
                )
                notes.append(
                    f"shipped profile sums to 1 instead of to the whole-tree index: "
                    f"max |sum(series)/total - 1| = {amp_max:.3g}"
                )
                if na_ok and rate_t >= 0.99 and rate >= 0.95:
                    status = "verified"
                elif na_ok and rate_t >= 0.99 and rate_c >= 0.95:
                    # the beta arithmetic and the deterministic index agree; only
                    # the rate differs, because the two implementations fit
                    # differently-constructed curves of it (the renormalised profile
                    # renormalisation, or the pairwise aggregation named in
                    # CURVE_CONTROL).  Attributed, not assumed: the control row
                    # above is the evidence.
                    status = "documented divergence"
                    notes.append(
                        "the rate is an accepted design divergence, not a mismatch of "
                        "the beta arithmetic: under the reference's own curve the same "
                        "arithmetic reproduces it, while the shipped profile "
                        "normalisation / aggregation differs by construction"
                    )
                else:
                    status = "FAILED"
                    notes.append(
                        "no reference-convention control closes the gap, so this is an "
                        "unexplained difference"
                    )
                rows.append(
                    {
                        "scenario": label,
                        "n": a_r.size,
                        "max_rel_err": worst,
                        "pass_rate": rate,
                        "tol": TOL_RATE,
                        "status": status,
                        "ref_basis": f"{name}:{rate_col}",
                        "note": "; ".join(notes),
                    }
                )
                # The whole-tree beta index is a deterministic quantity, so it is
                # scored on its own row at 1e-8 and stays certified even on the
                # families whose *rate* is an accepted divergence.
                rows.append(
                    {
                        "scenario": f"{label}:{total_col}",
                        "n": a_t.size,
                        "max_rel_err": worst_t,
                        "pass_rate": rate_t,
                        "tol": TOL_PROFILE,
                        "status": "verified" if rate_t >= 0.99 and na_ok else "FAILED",
                        "ref_basis": f"{name}:{total_col}",
                        "note": "per-neighbourhood whole-tree beta index reported by "
                        f"the reference alongside its {rate_col} fit; not a per-slice "
                        "profile",
                    }
                )

            for suffix, comp, appr in BETA_PROFILE_FAMILIES:
                name = f"{scenario}__{suffix}.csv"
                label = f"{scenario}:{suffix}"
                if not (ref_dir / name).is_file():
                    rows.append(
                        {
                            "scenario": label,
                            "n": 0,
                            "max_rel_err": np.nan,
                            "pass_rate": np.nan,
                            "tol": TOL_PROFILE,
                            "status": "no reference",
                            "ref_basis": "none",
                            "note": NO_REFERENCE.get(
                                suffix,
                                "the reference generator produced no file for this family",
                            ),
                        }
                    )
                    continue
                entry = records[name]
                try:
                    ref_arr, note = read_profile_reference(name, ref_dir, entry)
                except ReferenceProvenanceError as exc:
                    rows.append(
                        {
                            "scenario": label,
                            "n": 0,
                            "max_rel_err": np.nan,
                            "pass_rate": 0.0,
                            "tol": TOL_PROFILE,
                            "status": "provenance failure",
                            "ref_basis": "none",
                            "note": str(exc),
                        }
                    )
                    continue
                if ref_arr.size == 0:
                    rows.append(
                        {
                            "scenario": label,
                            "n": 0,
                            "max_rel_err": np.nan,
                            "pass_rate": np.nan,
                            "tol": TOL_PROFILE,
                            "status": "not evaluated",
                            "ref_basis": f"{name} (declared unverified in the manifest)",
                            "note": "hash + shape verified; numeric cross-check skipped - " + note,
                        }
                    )
                    continue
                mine = r_compatible_profile(
                    tree_a, mat_a, adj, slices, component=comp, approach=appr
                )["series"]
                if mine.shape != ref_arr.shape:
                    rows.append(
                        {
                            "scenario": label,
                            "n": mine.size,
                            "max_rel_err": np.nan,
                            "pass_rate": 0.0,
                            "tol": TOL_PROFILE,
                            "status": "orientation failure",
                            "note": f"phyloslicer produced {mine.shape[0]}x{mine.shape[1]} "
                            f"but the registered reference is "
                            f"{ref_arr.shape[0]}x{ref_arr.shape[1]}",
                        }
                    )
                    continue
                worst, rate = compare_pair(mine, ref_arr, TOL_PROFILE, zero_floor=ZERO_FLOOR)
                at_zero = np.abs(ref_arr) <= ZERO_FLOOR
                zero_note = (
                    f"{int(at_zero.sum())} reference cells are numerically 0 (the oldest "
                    f"slices, where no pair has diverged); compared absolutely there, "
                    f"max abs err {np.abs(mine - ref_arr)[at_zero].max():.3g}"
                    if at_zero.any()
                    else "no numerically-zero reference cells"
                )
                rows.append(
                    {
                        "scenario": label,
                        "n": mine.size,
                        "max_rel_err": worst,
                        "pass_rate": rate,
                        "tol": TOL_PROFILE,
                        "status": "verified" if rate >= 0.99 else "FAILED",
                        "ref_basis": name,
                        "identical_cells": count_identical(mine, ref_arr),
                        "note": "undivided per-slice beta series, i.e. treesliceR's "
                        "r_phylo(index='PB') curve before the whole-tree division; " + zero_note,
                    }
                )

    report = pd.DataFrame(rows)
    if "identical_cells" in report:
        # an integer column, not the float pandas would infer from the rows that
        # have no per-cell array to count (a 2.0 in a CSV reads like a rounding)
        report["identical_cells"] = report["identical_cells"].astype("Int64")
    # the report is a build artefact of this run; writing it next to the pinned
    # references dirties a tracked file every time a reader follows the
    # documented command, so the destination is selectable with --out
    out = Path(out) if out else ref_dir / "equivalence_report.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    report.to_csv(out, index=False)
    pd.set_option("display.max_colwidth", 60, "display.width", 200)
    print(report.to_string(index=False))

    checked = report[report["status"] == "verified"]
    skipped = report[report["status"] == "not evaluated"]
    gaps = report[report["status"] == "no reference"]
    diverged = report[report["status"] == "documented divergence"]
    tolerated = ("verified", "not evaluated", "no reference", "documented divergence")
    failed = report[~report["status"].isin(tolerated)]
    rates = checked["pass_rate"].dropna()
    # bucketed on the tolerance each row was scored at, so a new index family
    # cannot land in neither (or both) bucket by accident
    deterministic = checked[checked["tol"] == TOL_PROFILE]
    optimiser = checked[checked["tol"] == TOL_RATE]
    if len(deterministic):
        print(
            f"\ndeterministic profiles: worst rel err = "
            f"{deterministic['max_rel_err'].max():.3g} (tol {TOL_PROFILE:g}) over "
            f"{len(deterministic)} rows / {int(deterministic['n'].sum())} values"
        )
    else:
        print(
            "\ndeterministic profiles: NONE cross-checked - every per-slice "
            "reference is declared unverified in the manifest"
        )
    if len(optimiser):
        print(
            f"optimiser outputs:        worst rel err = "
            f"{optimiser['max_rel_err'].max():.3g} (tol {TOL_RATE:g}) over "
            f"{len(optimiser)} scenario x index rows (rate, total and origin columns)"
        )
    family = report["scenario"].str.split(":", n=1).str[1].fillna("")
    is_beta = family.str.startswith(BETA_LABELS)
    if bool(is_beta.any()):
        parts = []
        bver = report[is_beta & (report["status"] == "verified")]
        if len(bver):
            parts.append(
                f"{len(bver)} verified (worst rel err {bver['max_rel_err'].max():.3g}, "
                f"{int(bver['n'].sum())} values)"
            )
        if len(diverged):
            per = (
                diverged.assign(family=diverged["scenario"].str.split(":", n=1).str[1])
                .groupby("family")["max_rel_err"]
                .max()
                .sort_values()
            )
            parts.append(
                f"{len(diverged)} documented divergences in the *rate* only - the "
                "whole-tree index of every one of them agrees at "
                f"{TOL_PROFILE:g}, and refitting the reference's own curve closes the "
                "gap (the renormalisation recorded in validation/NOTES.md, or the pairwise "
                "aggregation named in "
                "CURVE_CONTROL); worst shipped error "
                + ", ".join(f"{k} {v:.3g}" for k, v in per.items())
                + ")"
            )
        if len(gaps):
            parts.append(
                f"{len(gaps)} families treesliceR 1.1.0 itself cannot produce "
                f"({', '.join(sorted(set(gaps['scenario'].str.split(':', n=1).str[1])))})"
            )
        print(f"\nbeta diversity (CpB / CpB_RW / r_phylo(index='PB')): {'; '.join(parts)}")
    if len(skipped):
        print(
            f"\nNOT EVALUATED: {len(skipped)} reference file(s) are registered but "
            "declared unverified - see the manifest provenance section."
        )
    if len(failed):
        kinds = ", ".join(f"{k} x{n}" for k, n in failed["status"].value_counts().items())
        print(f"\nFAILED ROWS: {len(failed)} ({kinds})")
        for _, row in failed.iterrows():
            print(f"  {row['scenario']}: rel err {row['max_rel_err']:.3g} - {row['note'][:120]}")

    ok = bool(len(rates)) and bool((rates >= 0.99).all()) and not len(failed)
    if strict and (len(skipped) or len(gaps) or len(diverged)):
        ok = False
    print(f"report -> {out}")
    print("verdict:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("ref_dir", nargs="?", default="validation/reference")
    parser.add_argument(
        "--out",
        default=None,
        help="where to write equivalence_report.csv (default: inside ref_dir, "
        "which is a tracked file; point it at a build directory in CI)",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="fail if any reference is declared unverified, any family has no "
        "reference, or any scored row is a documented divergence",
    )
    args = parser.parse_args()
    sys.exit(main(args.ref_dir, strict=args.strict, out=args.out))
