"""Command-line interface: phyloslicer slice | rates | posterior | validate."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------- #
# reference provenance (see validation/NOTES.md)
#
# ``validate`` used to take ``iloc[:, 0]`` of any CSV and call it "the R
# reference", which meant a reference could be phyloslicer's own output.  The
# fix is a manifest: every reference file is registered by SHA256 together with
# its shape, its column-name pattern and its orientation, and a file with no
# entry - or whose bytes differ from the registered hash - is refused.
#
# The reader/writer lives here rather than in ``validation/compare.py`` because
# the installed package cannot import a repository script; ``compare.py`` imports
# it back from this module so that exactly one implementation exists.
# ---------------------------------------------------------------------- #
MANIFEST_FORMAT = "phyloslicer-reference-manifest/1"
MANIFEST_NAME = "manifest.txt"
#: default places to look for a manifest when ``--manifest`` is not given
MANIFEST_SEARCH_NAMES = (Path("validation") / "reference" / MANIFEST_NAME,)

#: fields every record must carry to be usable as a reference; ``rows`` /
#: ``cols`` are added by the consumers that need a tabular shape (a Newick file
#: is registered with ``rows = -``, since "rows" has no meaning there)
REQUIRED_RECORD_FIELDS = ("sha256", "orientation")


class ReferenceProvenanceError(RuntimeError):
    """Raised when a reference file is not backed by registered provenance."""


def sha256_file(path) -> str:
    """SHA256 of a file's bytes (chunked; the reference files are small)."""
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def make_reference_record(
    path,
    *,
    orientation: str,
    columns_pattern: str,
    rows: int | None = None,
    cols: int | None = None,
    **extra: str,
) -> dict[str, str]:
    """Build one manifest record for ``path`` (shape auto-detected for CSVs)."""
    p = Path(path)
    rec: dict[str, str] = {
        "sha256": sha256_file(p),
        "orientation": orientation,
        "columns_pattern": columns_pattern,
    }
    if rows is None or cols is None:
        if p.suffix != ".csv":
            raise ValueError(
                f"shape must be given explicitly for non-CSV file {p.name} "
                "(pass rows='-' / cols='-' for a layout that has none)"
            )
        frame = pd.read_csv(p, index_col=0)
        rows = frame.shape[0] if rows is None else rows
        cols = frame.shape[1] if cols is None else cols
    rec["rows"] = str(rows) if isinstance(rows, str) else str(int(rows))
    rec["cols"] = str(cols) if isinstance(cols, str) else str(int(cols))
    for key, value in extra.items():
        rec[key] = str(value)
    return rec


def format_reference_manifest(header: dict[str, str], records: dict[str, dict[str, str]]) -> str:
    """Serialise a manifest (``key = value`` header plus ``[file NAME]`` blocks)."""
    lines = [
        "# PhyloSlicer cross-validation reference manifest",
        f"format = {MANIFEST_FORMAT}",
    ]
    for key, value in header.items():
        if key == "format":
            continue
        text = str(value)
        if "\n" in text:
            lines.append(f"{key} =")
            lines.extend(f"# {row}" for row in text.splitlines())
        else:
            lines.append(f"{key} = {text}")
    for name in sorted(records):
        lines.append("")
        lines.append(f"[file {name}]")
        for key, value in records[name].items():
            lines.append(f"{key} = {value}")
    return "\n".join(lines) + "\n"


def parse_reference_manifest(text: str) -> tuple[dict[str, str], dict[str, dict[str, str]]]:
    """Parse manifest text into ``(header_fields, records_by_file_name)``.

    Accepted syntax: ``#`` comment lines, ``key = value`` header lines and
    ``[file NAME]`` blocks whose lines are ``field = value``.  A multi-line value
    is written as ``key =`` followed by indented ``# `` continuation lines and is
    read back joined with a space, which keeps the long provenance notes both
    readable and greppable.
    """
    header: dict[str, str] = {}
    records: dict[str, dict[str, str]] = {}
    current: dict[str, str] | None = None
    pending: str | None = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            if pending is not None:
                body = line[1:].strip()
                if body:
                    target = header if current is None else current
                    target[pending] = f"{target[pending]} {body}".strip()
            continue
        pending = None
        block = re.fullmatch(r"\[file (.+)\]", line)
        if block:
            current = records.setdefault(block.group(1).strip(), {})
            continue
        if "=" not in line:
            raise ReferenceProvenanceError(f"unparseable manifest line: {line!r}")
        key, value = (s.strip() for s in line.split("=", 1))
        target = header if current is None else current
        if value:
            target[key] = value
        else:
            target[key] = ""
            pending = key
    if header.get("format") != MANIFEST_FORMAT:
        raise ReferenceProvenanceError(
            f"unsupported manifest format {header.get('format')!r}; expected {MANIFEST_FORMAT!r}"
        )
    return header, records


def default_manifest_paths(reference) -> list[Path]:
    """Where to look for a manifest: ``$PHYLOSLICER_REFERENCE_MANIFEST``, the
    reference file's own directory, then the walking-up
    ``validation/reference/manifest.txt`` convention."""
    out: list[Path] = []
    env = os.environ.get("PHYLOSLICER_REFERENCE_MANIFEST")
    if env:
        out.append(Path(env).expanduser())
    ref = Path(reference).expanduser().resolve()
    out.append(ref.parent / MANIFEST_NAME)
    for parent in ref.parents:
        for name in MANIFEST_SEARCH_NAMES:
            candidate = parent / name
            if candidate not in out:
                out.append(candidate)
    return out


def load_reference_manifest(reference, manifest=None):
    """Locate and parse the manifest governing ``reference``.

    Hard-fails (rather than degrading) when no manifest is found: without it a
    comparison cannot tell an external reference from the file under test.
    """
    candidates = [Path(manifest).expanduser()] if manifest else default_manifest_paths(reference)
    tried: list[Path] = []
    for path in candidates:
        tried.append(path)
        if path.is_dir():
            path = path / MANIFEST_NAME
            tried.append(path)
        if path.is_file():
            header, records = parse_reference_manifest(path.read_text(encoding="utf-8"))
            return header, records, path
    raise ReferenceProvenanceError(
        "no provenance manifest found (looked in: "
        + ", ".join(str(p) for p in tried)
        + "); refusing to treat an unregistered CSV as an R reference.  Add a "
        "``[file <name>]`` block recording sha256, rows, cols, orientation, "
        "columns_pattern and rate_column - see ``phyloslicer.cli.main."
        "make_reference_record()`` for the exact field set, or "
        "``validation/reference/manifest.txt`` for a worked example - then pass it "
        "with --manifest PATH (or $PHYLOSLICER_REFERENCE_MANIFEST)."
    )


def reference_provenance(
    reference, manifest=None, *, require_rate_column: bool = True
) -> dict[str, str]:
    """Return the manifest record for ``reference``, verifying its hash.

    Raises :class:`ReferenceProvenanceError` when the file is not registered,
    when its bytes differ from the registered hash, or when it lacks the fields
    needed to interpret it (rate column / orientation / shape).
    """
    path = Path(reference).expanduser()
    if not path.is_file():
        raise ReferenceProvenanceError(f"reference file not found: {path}")
    _header, records, manifest_path = load_reference_manifest(path, manifest)
    name = path.name
    if name not in records:
        raise ReferenceProvenanceError(
            f"{name} has no entry in the provenance manifest {manifest_path}, so"
            " nothing distinguishes it from phyloslicer's own output.  Register it"
            " there - a ``[file "
            + name
            + "]`` block with sha256/rows/cols/orientation/columns_pattern"
            + ("/rate_column" if require_rate_column else "")
            + " - or pass the manifest that already carries it with --manifest."
        )
    entry = dict(records[name])
    missing = [f for f in REQUIRED_RECORD_FIELDS if not entry.get(f)]
    if require_rate_column and not entry.get("rate_column"):
        missing.append("rate_column")
    if missing:
        raise ReferenceProvenanceError(
            f"manifest entry for {name} is incomplete, missing field(s): " + ", ".join(missing)
        )
    actual = sha256_file(path)
    if actual != entry["sha256"]:
        raise ReferenceProvenanceError(
            f"{name} does not match the hash registered in {manifest_path}: "
            f"manifest says {entry['sha256']}, file is {actual}.  The reference"
            " has been edited, truncated or regenerated since it was registered -"
            " re-run the generator and refresh the manifest."
        )
    return entry


def check_reference_shape(entry: dict[str, str], shape, *, name: str = "reference") -> None:
    """Refuse a reference whose shape differs from the registered one.

    This replaces the silent auto-transpose the comparator used to apply: a
    transposed file used to pass, which is what made a site-by-slice reference
    and a slice-by-site reference indistinguishable.
    """
    rows, cols = entry.get("rows", ""), entry.get("cols", "")
    if rows in ("", "-", "none") or cols in ("", "-", "none"):
        return  # a layout without a tabular shape (Newick) cannot be transposed
    rows, cols = int(rows), int(cols)
    if (rows, cols) != tuple(shape):
        raise ReferenceProvenanceError(
            f"{name} is {int(shape[0])}x{int(shape[1])} but the manifest registers "
            f"{rows}x{cols} (orientation={entry['orientation']}, "
            f"columns_pattern={entry.get('columns_pattern')}).  A transposed or "
            "otherwise re-shaped reference is a hard error; regenerate it and "
            "refresh the manifest."
        )


def check_reference_columns(entry: dict[str, str], columns, *, name: str = "reference") -> None:
    """Check column names against the registered pattern."""
    pattern = entry.get("columns_pattern")
    if not pattern or pattern in ("-", "none"):
        return
    bad = [str(c) for c in columns if not re.fullmatch(pattern, str(c))]
    if bad:
        raise ReferenceProvenanceError(
            f"{name} has column(s) {bad[:5]} that violate the registered pattern "
            f"{pattern!r} (manifest orientation={entry.get('orientation')})"
        )


# ---------------------------------------------------------------------- #
# argument parsing
# ---------------------------------------------------------------------- #
def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="phyloslicer",
        description="Temporal slicing of phylogenies and cumulative-diversity rates",
    )
    p.add_argument("--version", action="store_true", help="print version and exit")
    sub = p.add_subparsers(dest="command")

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--log", type=str, default=None, help="write a JSON run summary")

    sp_slice = sub.add_parser("slice", parents=[common], help="slice a tree")
    sp_slice.add_argument("tree", help="input Newick tree")
    sp_slice.add_argument("--rootward", type=float, default=None, help="keep youngest X Ma")
    sp_slice.add_argument("--tipward", type=float, default=None, help="drop youngest X Ma")
    sp_slice.add_argument("--n-slices", type=int, default=None, help="cut into N equal slices")
    sp_slice.add_argument("--criterion", choices=["time", "pd"], default="time")
    sp_slice.add_argument(
        "--out", required=True, help="output Newick file (single slice) or directory"
    )
    sp_slice.add_argument(
        "--force",
        action="store_true",
        help="overwrite an existing --out (a slice directory is cleared of its"
        " slice_*.tre files first, so a later 'posterior' run cannot mix two"
        " generations of slices)",
    )

    sp_rates = sub.add_parser("rates", parents=[common], help="cumulative rates")
    sp_rates.add_argument("tree")
    sp_rates.add_argument("sites", help="presence/absence CSV (sites x species)")
    sp_rates.add_argument("--rate", choices=["cpd", "cpe"], default="cpd")
    sp_rates.add_argument("--n-slices", type=int, default=100)
    sp_rates.add_argument("--criterion", choices=["time", "pd"], default="time")
    sp_rates.add_argument("--percentile", type=float, default=0.95)
    sp_rates.add_argument("--out", required=True, help="output CSV")

    sp_post = sub.add_parser("posterior", parents=[common], help="posterior tree sets")
    sp_post.add_argument("trees", help="directory or multi-tree Newick file")
    sp_post.add_argument("sites", help="presence/absence CSV")
    sp_post.add_argument("--rate", choices=["cpd", "cpe"], default="cpd")
    sp_post.add_argument("--n-slices", type=int, default=100)
    sp_post.add_argument("--hdi", type=float, default=0.95)
    sp_post.add_argument("--out", required=True)

    sp_val = sub.add_parser(
        "validate", parents=[common], help="cross-validate against an R reference CSV"
    )
    sp_val.add_argument("tree")
    sp_val.add_argument("sites")
    sp_val.add_argument(
        "--reference",
        required=True,
        help="reference CSV produced by treesliceR (must be registered in a manifest)",
    )
    sp_val.add_argument(
        "--manifest",
        default=None,
        help=f"path of the provenance {MANIFEST_NAME} governing --reference "
        "(default: $PHYLOSLICER_REFERENCE_MANIFEST, the reference's own directory, "
        "then validation/reference/manifest.txt)",
    )
    sp_val.add_argument(
        "--rate-column",
        default=None,
        help="assert which column the manifest names as the reference rate "
        "(any other value is a hard error; the manifest remains the authority)",
    )
    sp_val.add_argument("--n-slices", type=int, default=100)
    sp_val.add_argument(
        "--tolerance",
        type=float,
        default=1e-4,
        help="maximum relative error for a site to count as equivalent "
        "(default 1e-4: the optimiser agreement the shipped references actually "
        "reach is ~2e-6, so the previous 1e-8 default reported 0.0000 on a correct run)",
    )
    sp_val.add_argument(
        "--no-strict",
        dest="strict",
        action="store_false",
        help="exit 0 even when the equivalence rate is below 1 (default: a "
        "mismatch is a nonzero exit, so this can gate CI)",
    )
    sp_val.add_argument("--out", required=True)

    return p


def _load_tree(path: str):
    from ..core.tree import TreeArray

    return TreeArray.from_newick(path)


def _load_mat(path: str) -> pd.DataFrame:
    mat = pd.read_csv(path, index_col=0)
    mat.columns = mat.columns.astype(str)
    return (mat > 0).astype(float)


def _write_log(args, started: float, extra: dict) -> None:
    if not getattr(args, "log", None):
        return
    from .. import __version__

    payload = {
        "version": __version__,
        "command": args.command,
        # both clocks come from the same source: `started` is a perf_counter
        # reading (see main()), and subtracting it from time.time() used to put
        # the Unix epoch into every run summary (elapsed_s ~= 1.79e9)
        "elapsed_s": round(time.perf_counter() - started, 3),
        "started_at": _STARTED_AT,
        **extra,
    }
    Path(args.log).write_text(json.dumps(payload, indent=2))


def _refuse_dirty_out(path: Path, force: bool, what: str) -> int | None:
    """Return an exit code when ``path`` already holds an earlier generation."""
    if force:
        return None
    exists = path.is_dir() and any(path.glob("slice_*.tre")) if path.is_dir() else path.exists()
    if exists:
        print(
            f"error: {what} {path} already exists and is not empty; passing --force "
            "overwrites it.  Refusing by default because a slice directory is fed "
            "straight back into `posterior`, and stale slices from an earlier run "
            "silently join the pool.",
            file=sys.stderr,
        )
        return 2
    return None


def cmd_slice(args) -> int:
    from ..slicing.api import slice_pieces, slice_rootward, slice_tipward

    tree = _load_tree(args.tree)
    out = Path(args.out)
    n_given = sum(x is not None for x in (args.rootward, args.tipward, args.n_slices))
    if n_given != 1:
        print("error: pass exactly one of --rootward / --tipward / --n-slices", file=sys.stderr)
        return 2
    if (rc := _refuse_dirty_out(out, args.force, "output path")) is not None:
        return rc
    if args.n_slices is not None:
        stack = slice_pieces(tree, n=args.n_slices, criterion=args.criterion)
        out.mkdir(parents=True, exist_ok=True)
        if args.force:
            for stale in out.glob("slice_*.tre"):
                stale.unlink()
        for j, sliced in enumerate(stack.to_sliced_trees()):
            (out / f"slice_{j + 1:04d}.tre").write_text(sliced.tree.to_newick() + "\n")
        _write_log(args, _START["t0"], {"n_slices": stack.n_slices})
        print(f"wrote {stack.n_slices} slices to {out}")
        return 0
    if args.rootward is not None:
        sliced = slice_rootward(tree, args.rootward, criterion=args.criterion)
    else:
        sliced = slice_tipward(tree, args.tipward, criterion=args.criterion)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(sliced.tree.to_newick() + "\n")
    _write_log(args, _START["t0"], {"pd": sliced.pd})
    print(f"slice PD = {sliced.pd:.6g} -> {out}")
    return 0


def cmd_rates(args) -> int:
    from ..rates.api import cpd_rate, cpe_rate

    tree, mat = _load_tree(args.tree), _load_mat(args.sites)
    runner = cpd_rate if args.rate == "cpd" else cpe_rate
    res = runner(
        tree, mat, n_slices=args.n_slices, criterion=args.criterion, percentile=args.percentile
    )
    res.to_csv(args.out, index=False)
    _write_log(
        args,
        _START["t0"],
        {
            "n_sites": len(res),
            "n_converged": int(res["converged"].sum()),
            "median_rate": float(res.iloc[:, 1].median()),
        },
    )
    print(f"wrote {len(res)} rows to {args.out}")
    return 0


def cmd_posterior(args) -> int:
    from ..uncertainty import posterior_rate

    mat = _load_mat(args.sites)
    res = posterior_rate(args.trees, mat, rate=args.rate, n_slices=args.n_slices, hdi=args.hdi)
    res.to_csv(args.out, index=False)
    _write_log(args, _START["t0"], {"n_sites": len(res)})
    print(f"wrote {len(res)} rows to {args.out}")
    return 0


def cmd_validate(args) -> int:
    # 1. provenance first: an unregistered CSV is not a reference (A13)
    try:
        entry = reference_provenance(args.reference, args.manifest, require_rate_column=True)
    except ReferenceProvenanceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    orientation = entry["orientation"]
    if orientation != "site-by-index":
        print(
            f"error: validate compares one reference row per site, but the manifest "
            f"registers {Path(args.reference).name} as orientation={orientation!r} "
            "(expected 'site-by-index').  Pass the file that holds the rates, not a "
            "per-slice profile.",
            file=sys.stderr,
        )
        return 2
    rate_col = entry["rate_column"]
    if args.rate_column and args.rate_column != rate_col:
        print(
            f"error: --rate-column {args.rate_column!r} disagrees with the manifest, "
            f"which registers the reference rate as {rate_col!r}; fix the manifest "
            "entry rather than overriding it on the command line.",
            file=sys.stderr,
        )
        return 2
    # Which index this reference actually holds decides which analysis to run.
    # Comparing a CpB/CpE reference against CpD numbers used to "work" and print
    # a 0 % rate under a `reference_CpD` column that contained CpB values.  The
    # index is read off the manifest's own rate_column, tolerating suffixed
    # spellings such as "CpD_ref" but refusing anything that is not a CpD/CpE.
    index_kind = re.match(r"^(CpD|CpE)", rate_col)
    runner_name = {"CpD": "cpd", "CpE": "cpe"}[index_kind.group(1)] if index_kind else None
    if runner_name is None:
        print(
            f"error: {Path(args.reference).name} is registered as the {rate_col!r} "
            "index, but `validate` compares per-site cumulative rates for CpD and "
            "CpE only.  Beta references (CpB/CpB_RW) also need an adjacency matrix, "
            "which this subcommand does not take; score them with "
            "`python validation/compare.py`.",
            file=sys.stderr,
        )
        return 2
    # The manifest records how the reference was produced; comparing a run with
    # a different window count is not a comparison of the same quantity.
    try:
        header, _records, _manifest_path = load_reference_manifest(args.reference, args.manifest)
    except ReferenceProvenanceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    registered_slices = header.get("n_slices")
    if registered_slices and int(registered_slices) != args.n_slices:
        print(
            f"error: the manifest produced {Path(args.reference).name} with "
            f"n_slices = {registered_slices}, but --n-slices {args.n_slices} was "
            "requested; the two curves are not comparable.  Pass "
            f"--n-slices {registered_slices}.",
            file=sys.stderr,
        )
        return 2

    from ..rates.api import cpd_rate, cpe_rate

    tree, mat = _load_tree(args.tree), _load_mat(args.sites)
    res = (cpd_rate if runner_name == "cpd" else cpe_rate)(tree, mat, n_slices=args.n_slices)
    ref = pd.read_csv(args.reference, index_col=0)
    name = Path(args.reference).name
    try:
        check_reference_shape(entry, ref.shape, name=name)
        check_reference_columns(entry, ref.columns, name=name)
        if rate_col not in ref.columns:
            raise ReferenceProvenanceError(
                f"{name} carries no column {rate_col!r} (it has {list(ref.columns)}), "
                "contrary to its manifest entry"
            )
    except ReferenceProvenanceError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    if len(ref) != len(res):
        print(
            f"error: reference has {len(ref)} rows but analysis produced {len(res)}",
            file=sys.stderr,
        )
        return 2
    our_col = "CpD" if runner_name == "cpd" else "CpE"
    mine = res[our_col].to_numpy(dtype=float)
    theirs = ref[rate_col].to_numpy(dtype=float)
    ok = np.isfinite(mine) & np.isfinite(theirs)
    both_na = ~np.isfinite(mine) & ~np.isfinite(theirs)
    one_sided_na = np.isfinite(mine) != np.isfinite(theirs)
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(ok, np.abs(mine - theirs) / np.maximum(np.abs(theirs), 1e-300), np.nan)
    # agreement rules (matching validation/NOTES.md): a finite pair compares
    # by relative error, a pair that is NA on BOTH sides counts as agreement,
    # and a one-sided NA counts as disagreement
    passed = np.where(
        both_na, True, np.where(ok, np.nan_to_num(rel, nan=1.0) <= args.tolerance, False)
    )
    report = pd.DataFrame(
        {
            "site_id": res["site_id"],
            f"phyloslicer_{our_col}": mine,
            f"reference_{rate_col}": theirs,
            "relative_error": rel,
            "equivalent": passed,
        }
    )
    report.to_csv(args.out, index=False)
    rate = float(np.mean(passed)) if passed.size else float("nan")
    observed = float(np.nanmax(rel)) if np.isfinite(rel).any() else float("nan")
    _write_log(
        args,
        _START["t0"],
        {
            "equivalence_rate": rate,
            "max_relative_error": observed,
            "tolerance": args.tolerance,
            "n_slices": args.n_slices,
            "reference": name,
            "reference_column": rate_col,
            "reference_sha256": entry["sha256"],
            "reference_orientation": orientation,
            "n_one_sided_na": int(one_sided_na.sum()),
        },
    )
    print(
        f"reference {name} column {rate_col!r} "
        f"({ref.shape[0]}x{ref.shape[1]}, orientation={orientation}, "
        f"sha256 {entry['sha256'][:12]}…) registered in the provenance manifest"
    )
    # the observed error is printed next to the pass rate: a bare "0.0000" told
    # a user the comparison failed when in fact it agreed to 2e-6
    print(
        f"equivalence rate: {rate:.4f} (tolerance {args.tolerance:g}"
        + (f", observed max relative error {observed:.3g})" if np.isfinite(observed) else ")")
    )
    if one_sided_na.any():
        print(f"note: {int(one_sided_na.sum())} site(s) are NA on one side only")
    if args.strict and not (rate >= 1.0):
        print(
            f"error: {int((~passed).sum())} of {passed.size} site(s) disagree beyond "
            f"{args.tolerance:g}; exiting nonzero (pass --no-strict for exploratory use)",
            file=sys.stderr,
        )
        return 1
    return 0


_START: dict[str, float] = {"t0": 0.0}
#: ISO wall-clock stamp of the current invocation, recorded alongside `_START["t0"]`
_STARTED_AT: str = ""


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    if args.version:
        from .. import __version__

        print(f"phyloslicer {__version__}")
        return 0
    if not getattr(args, "command", None):
        parser.print_help()
        return 2
    global _STARTED_AT
    _START["t0"] = time.perf_counter()
    _STARTED_AT = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    handlers = {
        "slice": cmd_slice,
        "rates": cmd_rates,
        "posterior": cmd_posterior,
        "validate": cmd_validate,
    }
    try:
        return handlers[args.command](args)
    except Exception as exc:  # CLI boundary: report, do not traceback-dump
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
