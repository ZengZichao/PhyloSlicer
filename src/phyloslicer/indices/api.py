"""Indices layer: PD / PE / PB(+RW) totals and per-slice profiles.

Every "for each site, for each slice" loop of the R implementation is
expressed here as sparse matrix algebra:

    * ``A`` : (n_sites, E) site-edge membership (0/1);
    * ``C`` : (E, k) edge-slice contribution matrix (see ``slicing.windows``);
    * per-site-per-slice PD = ``A @ C``;
    * PE uses the same product with range-weighted contributions
      ``C * (1 / range_size)``.

Beta diversity is computed per focal neighbourhood (focal site + the sites
flagged in its ``adj`` row) with vectorised pair enumeration - a single
code path for every (component x approach x weighting) combination.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as _pd
import scipy.sparse as sp

from ..core.tree import TreeArray
from ..io.matrices import align_tree_matrix, as_matrix
from ..slicing.api import SliceStack

__all__ = [
    "site_edge_membership",
    "edge_range_sizes",
    "pd_per_slice",
    "pe_per_slice",
    "pb_per_slice",
    "r_compatible_profile",
    "pd",
    "pe",
    "pb",
    "MULTISITE_DOMAINS",
    "MULTISITE_DOMAIN_ALIASES",
]

#: Summation domains the multisite beta components can be evaluated on.
#:
#: ``"mixed"`` reproduces the treesliceR arithmetic exactly: the shared term
#: ``a`` is summed over the *m* sites of the neighbourhood while the turnover
#: terms ``b + c`` and ``min(b, c)`` are summed over the *C(m, 2)* pairs.  The
#: two sums therefore grow differently with neighbourhood size (``a`` linearly,
#: ``b + c`` quadratically) and the resulting index drifts upward: on the
#: reference construction of one shared skeleton branch plus one private
#: terminal branch per site, all of length 1, Sørensen equals ``m / (m + 2)``
#: (0.5 at m = 2, 0.6 at 3, 0.714 at 5, 0.833 at 10, 0.909 at 20) while every
#: individual pair is unchanged at 0.5.  Values computed on different-sized
#: neighbourhoods are therefore **not comparable**.
#:
#: ``"paired"`` puts every term in the pair-summed domain (``a`` becomes
#: ``sum_pairs a_ij``), which removes the drift and leaves ``m = 2`` untouched;
#: for ``m > 2`` the resulting index matches what the pairwise aggregation of
#: the same neighbourhood gives.  It does **not** fix the separate
#: ``min(b, c)``-non-additivity issue - that is handled by the
#: profile-derived denominator in :func:`pb_per_slice`.
#:
#: ``"mixed"`` stays the default because it is what the reference implementation
#: publishes; switching the default would silently redefine a metric.  Opt in
#: with ``multisite_domain="paired"`` and see ``docs/USER_GUIDE.md`` section 6.
MULTISITE_DOMAINS = ("mixed", "paired")

#: Spellings accepted anywhere a domain is expected.  ``"reference"`` is published
#: as a synonym of the default ``"mixed"`` so that callers can name the domain
#: after what it is (the treesliceR arithmetic) instead of after its defect.
#: ``"mixed"`` remains the canonical spelling reported back by the code.
MULTISITE_DOMAIN_ALIASES = {"reference": "mixed"}


def _normalise_domain(domain, *, caller: str = "multisite_domain") -> str:
    """Map an accepted spelling onto a canonical :data:`MULTISITE_DOMAINS` name."""
    key = domain if isinstance(domain, str) else str(domain)
    key = MULTISITE_DOMAIN_ALIASES.get(key, key)
    if key not in MULTISITE_DOMAINS:
        raise ValueError(
            f"{caller} must be one of {MULTISITE_DOMAINS + tuple(MULTISITE_DOMAIN_ALIASES)}"
        )
    return key


# ---------------------------------------------------------------------- #
# membership machinery
# ---------------------------------------------------------------------- #
def _match_columns_to_tree(tree: TreeArray, mat, *, quantity: str = "PD"):
    """Reorder/complete ``mat``'s columns to ``tree.tip_labels``, warning if the
    label *sets* differ.

    Silent zero-filling used to turn "species name misspelled in the matrix"
    into "that species occurs at no site", while ``pd()``/``pe()`` already
    warned about exactly that situation - so the same mismatch produced a
    different contract depending on which function was called (review, item
    4/"other interface issues").  A pure reordering is still silent, because it
    changes nothing about the result.
    """
    mat = as_matrix(mat)
    if list(mat.columns) == tree.tip_labels:
        return mat
    missing = set(tree.tip_labels) - set(mat.columns)
    extra = set(mat.columns) - set(tree.tip_labels)
    if missing and extra:
        # nothing matches at all: reindexing would fill every column with zero
        # and `pd()` would answer "this assemblage holds no phylogenetic
        # diversity" for what is really a naming mistake.  align_tree_matrix
        # already refuses this; the index functions must not be more forgiving.
        raise ValueError(
            f"no matrix column matches any tree tip label ({len(missing)} tree "
            f"tips, {len(extra)} matrix species, nothing in common): "
            f"tree starts {sorted(tree.tip_labels)[:3]}, matrix starts "
            f"{sorted(mat.columns)[:3]}.  Fix the names, or align first with "
            f"phyloslicer.io.align_tree_matrix(tree, mat).  Passing a bare "
            f"numpy array gives integer column names and never matches labels."
        )
    if missing or extra:
        warnings.warn(
            f"matrix and tree tip labels differ: "
            f"{len(missing)} tree tips absent from matrix, "
            f"{len(extra)} matrix species absent from tree; "
            f"unmatched tips contribute zero {quantity}",
            stacklevel=3,
        )
    return mat.reindex(columns=tree.tip_labels, fill_value=0.0)


def site_edge_membership(tree: TreeArray, mat) -> sp.csr_matrix:
    """(n_sites, E) binary matrix: edge is on the path of >= 1 present tip.

    Built from the tip-interval representation with one ``searchsorted`` per
    site - O(E log S) per site, no (sites x tips) intermediate.
    """
    mat_df = _match_columns_to_tree(tree, mat)
    intervals = tree.tip_intervals()
    lo, hi = intervals[:, 0], intervals[:, 1]
    order = np.argsort(lo)
    lo_s, hi_s = lo[order], hi[order]

    indptr = [0]
    indices: list[int] = []
    data: list[float] = []
    present_rows = np.asarray(mat_df.values) > 0
    for row in present_rows:
        tips = np.flatnonzero(row)
        if tips.size == 0:
            indptr.append(indptr[-1])  # empty row: no covered edges
            continue
        tips.sort()
        # edge covered iff its interval [lo, hi) contains some present tip
        left = np.searchsorted(tips, lo_s, side="left")
        right = np.searchsorted(tips, hi_s, side="left")
        covered = left < right
        cols = order[covered]
        indices.append(cols)
        data.append(np.ones(cols.size))
        indptr.append(indptr[-1] + cols.size)
    indices_arr = np.concatenate(indices) if indices else np.empty(0, dtype=np.int64)
    data_arr = np.concatenate(data) if data else np.empty(0)
    return sp.csr_matrix(
        (data_arr, indices_arr, np.asarray(indptr, dtype=np.int64)),
        shape=(len(mat_df), tree.n_edges),
    )


def edge_range_sizes(tree: TreeArray, mat) -> np.ndarray:
    """(E,) range size of each edge = number of sites covering it."""
    A = site_edge_membership(tree, mat)
    return np.asarray(A.sum(axis=0)).ravel()


# ---------------------------------------------------------------------- #
# PD & PE
# ---------------------------------------------------------------------- #
def _assert_slices_match_tree(tree: TreeArray, slices: SliceStack) -> None:
    """Reject slice stacks that were not built from exactly this tree."""
    same = (
        slices.tree.n_edges == tree.n_edges
        and np.array_equal(slices.tree.edges, tree.edges)
        and np.allclose(np.sort(slices.tree.lengths), np.sort(tree.lengths))
    )
    if not same:
        raise ValueError(
            "slices were not built from the (aligned) tree passed in; "
            "align tree and matrix first, then slice"
        )


def _prepare(tree: TreeArray, mat, slices: SliceStack, weighted: bool):
    """Membership + (optionally range-weighted) contribution matrix.

    The caller must pass a tree aligned with ``mat`` and slices built from
    that *same aligned* tree (see ``phyloslicer.io.align_tree_matrix``).
    """
    mat = _match_columns_to_tree(tree, mat, quantity="PD" if not weighted else "PE")
    _assert_slices_match_tree(tree, slices)
    A = site_edge_membership(tree, mat)
    if weighted:
        ranges = np.asarray(A.sum(axis=0)).ravel()
        if (ranges == 0).any():
            # edges under zero-occurrence tips cannot happen after alignment;
            # guard numerically anyway
            ranges = np.where(ranges == 0, 1.0, ranges)
        weights = 1.0 / ranges
    else:
        weights = np.ones(tree.n_edges)
    C = slices.contribution_matrix() * weights[:, None]
    return tree, mat, A, C


def pd_per_slice(tree: TreeArray, mat, slices: SliceStack) -> np.ndarray:
    """(n_sites, k) PD stored in each slice for each site."""
    _, _, A, C = _prepare(tree, mat, slices, weighted=False)
    return np.asarray(A @ C)


def pe_per_slice(tree: TreeArray, mat, slices: SliceStack) -> np.ndarray:
    """(n_sites, k) range-weighted PD (phylogenetic endemism) per slice per site."""
    _, _, A, C = _prepare(tree, mat, slices, weighted=True)
    return np.asarray(A @ C)


def pd(tree: TreeArray, mat) -> np.ndarray:
    """Total Faith's PD per site: the branch length each assemblage stores.

    The sum runs over the edges of :class:`TreeArray`, whose convention (like
    ``ape::keep.tip(tree, ..., root.edge = 0)``) carries **no root stem**: the
    deepest nodes are the root's children, so the returned amount excludes any
    "path to the root" above the tree's own branching.  Matrix columns not
    present as tree tips are dropped and tree tips absent from the matrix
    contribute zero; either way a warning is issued.
    """
    mat_df = _match_columns_to_tree(tree, mat, quantity="PD")
    A = site_edge_membership(tree, mat_df)
    return np.asarray(A @ tree.lengths).ravel()


def pe(tree: TreeArray, mat) -> np.ndarray:
    """Total phylogenetic endemism (range-weighted PD) per site.

    Matrix columns not present as tree tips are dropped and tree tips absent
    from the matrix contribute zero; either way a warning is issued (see
    :func:`pd` for the root-stem convention).
    """
    mat_df = _match_columns_to_tree(tree, mat, quantity="PE")
    A = site_edge_membership(tree, mat_df)
    ranges = np.asarray(A.sum(axis=0)).ravel()
    ranges = np.where(ranges == 0, 1.0, ranges)
    return np.asarray(A @ (tree.lengths / ranges)).ravel()


# ---------------------------------------------------------------------- #
# beta diversity
# ---------------------------------------------------------------------- #
def _pair_vars(pd_piece: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Given per-site slice PDs (m, k) return per-pair sums and maxima ((P,k) each)."""
    m = pd_piece.shape[0]
    ii, jj = np.triu_indices(m, k=1)
    sums = pd_piece[ii] + pd_piece[jj]
    maxes = np.maximum(pd_piece[ii], pd_piece[jj])
    return sums, maxes


def _pair_indicator(ii: np.ndarray, jj: np.ndarray, m: int) -> sp.csr_matrix:
    """(P, m) 0/1 matrix whose row p selects the two members of pair p.

    Columns must interleave as [ii0, jj0, ii1, jj1, ...] so that pair p
    selects members ii[p] and jj[p] (a plain concatenate([ii, jj]) would
    silently mispair sites in neighbourhoods with >2 members).
    """
    P = len(ii)
    cols = np.empty(2 * P, dtype=np.int64)
    cols[0::2] = ii
    cols[1::2] = jj
    rows = np.repeat(np.arange(P), 2)
    return sp.csr_matrix((np.ones(2 * P), (rows, cols)), shape=(P, m))


def _neighbourhood_profiles(
    tree: TreeArray,
    mat: _pd.DataFrame,
    adj: np.ndarray | _pd.DataFrame,
    slices: SliceStack,
    weighted: bool,
):
    """Compute the raw per-neighbourhood quantities used by every PB variant.

    Returns a list (one entry per row of ``adj``) of dicts with full-tree and
    per-slice PDs for single sites and site pairs.
    """
    adj = np.asarray(adj)
    if adj.shape[0] != mat.shape[0]:
        raise ValueError("adj must have one row per site of mat")
    A_all = site_edge_membership(tree, mat)
    C = slices.contribution_matrix()
    if weighted:
        ranges = np.asarray(A_all.sum(axis=0)).ravel()
        ranges = np.where(ranges == 0, 1.0, ranges)
        vec = tree.lengths / ranges
        C = C / ranges[:, None]
    else:
        vec = tree.lengths

    results = []
    for i in range(adj.shape[0]):
        members = np.flatnonzero(adj[i] == 1)
        if members.size < 2:
            results.append({"status": "fewer_than_two_sites"})
            continue
        A_loc = A_all[members]
        # per-site quantities
        pd_full = np.asarray(A_loc @ vec).ravel()
        pd_piece = np.asarray(A_loc @ C)  # (m, k)
        # pair unions via indicator x membership
        ii, jj = np.triu_indices(members.size, k=1)
        indicator = _pair_indicator(ii, jj, members.size)
        A_pairs = (indicator @ A_loc) > 0
        A_pairs = A_pairs.astype(np.float64)
        pd_pairs_full = np.asarray(A_pairs @ vec).ravel()
        pd_pairs_piece = np.asarray(A_pairs @ C)  # (P, k)
        # union of ALL members (multisite totals)
        union_mask = (np.asarray(A_loc.sum(axis=0)).ravel() > 0).astype(np.float64)
        pd_union_full = float(union_mask @ vec)
        occupancy = float(np.asarray(A_loc.sum()).ravel()[0])
        dense_size = A_loc.shape[0] * A_loc.shape[1]

        results.append(
            {
                "status": "ok",
                "members": members,
                "n_members": members.size,
                "all_present": bool(occupancy == dense_size),
                "all_absent": bool(occupancy == 0.0),
                "n_species": int(mat.shape[1]),
                "pd_full": pd_full,  # (m,)
                "pd_piece": pd_piece,  # (m, k)
                "pairs_i": ii,
                "pairs_j": jj,  # (P,)
                "pd_pairs_full": pd_pairs_full,  # (P,)
                "pd_pairs_piece": pd_pairs_piece,  # (P, k)
                "pd_union_full": pd_union_full,
            }
        )
    return results


def _pb_series(
    res: dict, component: str, approach: str, domain: str = "mixed"
) -> tuple[np.ndarray, float]:
    """Per-slice beta series and total beta for one neighbourhood.

    Single code path for every (component x approach) combination.  The series
    holds the *raw* per-slice index values in **absolute** units (root -> tip);
    the caller turns them into a relative profile.  ``series / total`` is the
    quantity treesliceR feeds into ``nls`` (:func:`r_compatible_profile`),
    ``series / series.sum()`` is the profile that actually integrates to 1
    (:func:`pb_per_slice`).  The two differ whenever ``min(b, c)`` is involved,
    because ``min`` does not commute with the sum over slices - measured on the
    committed ``bd20__random`` scenario at 30 slices, the R-compatible turnover
    rows sum to 0.9236-1.0 and the nestedness rows to 1.0-1.5406 (review, item
    A2), while the :func:`pb_per_slice` rows sum to 1 within 4.5e-16.

    For the pairwise approach, per-pair ratios are averaged within each slice -
    the treesliceR ``CpB`` pairwise branch is dimensionally inconsistent for
    sorensen (it mixes k x P matrices with k-vectors inside ``nls``; measured
    under treesliceR 1.1.0, that branch aborts with "'qr' and 'y' must have the
    same number of rows", so no reference exists for it - see
    ``validation/NOTES.md``), while its ``CpB_RW`` pairwise branch does use the
    pairwise mean adopted here.

    ``domain`` selects the summation domain of the shared term ``a`` in the
    multisite approach: ``"mixed"`` - also spelled ``"reference"`` - keeps the
    treesliceR arithmetic, ``"paired"`` sums every term over pairs; see
    :data:`MULTISITE_DOMAINS`.  It has no effect on the pairwise approach, which
    is pair-summed by construction.
    """
    domain = _normalise_domain(domain, caller="domain")
    ii, jj = res["pairs_i"], res["pairs_j"]
    pd_full = res["pd_full"]
    pd_piece = res["pd_piece"]  # (m, k)
    pd_pairs_full = res["pd_pairs_full"]
    pd_pairs_piece = res["pd_pairs_piece"]  # (P, k)

    # ---- full-tree quantities ----------------------------------------- #
    shared_full = pd_full[ii] + pd_full[jj] - pd_pairs_full  # a_ij, (P,)
    bc_full = pd_pairs_full - shared_full  # b + c, (P,)
    min_full = pd_pairs_full - np.maximum(pd_full[ii], pd_full[jj])  # min(b, c)
    # multisite sums follow the R arithmetic: sum over the recycled
    # (P x 2) matrix gives sum_pairs (2*pw - pd_i - pd_j) = sum (pw - shared)
    bc_tot_ms = float(bc_full.sum())  # == sum_pairs (pw - shared): pair domain
    min_tot_ms = float(min_full.sum())  # pair domain
    if domain == "paired":
        # put `a` in the same pair-summed domain as b + c and min(b, c), so the
        # multisite index stops growing with the neighbourhood size
        a_tot_ms = float(shared_full.sum())
    else:  # "mixed" - treesliceR's site-summed `a` = sum_edges L_e (n_e - 1),
        # which grows linearly in m while the pair sums grow quadratically
        a_tot_ms = float(pd_full.sum() - res["pd_union_full"])

    # ---- per-slice quantities ----------------------------------------- #
    pair_sums, pair_maxes = _pair_vars(pd_piece)  # (P, k)
    shared_slice = pair_sums - pd_pairs_piece  # a_ij within slice
    min_slice = pd_pairs_piece - pair_maxes  # min(b, c) within slice
    # multisite per-slice numerator: sum_pairs (pw_slice - shared_slice)
    bc_slice_ms = (pd_pairs_piece - shared_slice).sum(axis=0)

    def _div(num: np.ndarray, den: np.ndarray) -> np.ndarray:
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.asarray(num, dtype=float) / np.asarray(den, dtype=float)

    if component not in ("sorensen", "turnover", "nestedness"):
        raise ValueError("component must be 'sorensen', 'turnover' or 'nestedness'")

    if component == "sorensen":
        if approach == "multisite":
            den = 2.0 * a_tot_ms + bc_tot_ms
            r_slice = _div(bc_slice_ms, den)
            total = _div(np.asarray(bc_tot_ms), np.asarray(den))
        else:
            den_pair = 2.0 * shared_full + bc_full  # (2a + b + c), (P,)
            r_slice = _div(pd_pairs_piece - shared_slice, den_pair[:, None]).mean(axis=0)
            total = np.nanmean(_div(bc_full, den_pair))
    elif component == "turnover":
        if approach == "multisite":
            turn_tot = min_tot_ms
            den = a_tot_ms + turn_tot
            r_slice = _div(min_slice.sum(axis=0), den)
            total = _div(np.asarray(turn_tot), np.asarray(den))
        else:
            den_pair = shared_full + min_full  # a + min(b, c)
            g_pair = _div(min_full, den_pair)
            r_slice = _div(min_slice, den_pair[:, None]).mean(axis=0)
            total = np.nanmean(g_pair)
    else:  # nestedness = sorensen - turnover
        if approach == "multisite":
            turn_tot = min_tot_ms
            sor = _div(np.asarray(bc_tot_ms), np.asarray(2.0 * a_tot_ms + bc_tot_ms))
            tur = _div(np.asarray(turn_tot), np.asarray(a_tot_ms + turn_tot))
            total = float(sor - tur)
            r_sor = _div(bc_slice_ms, 2.0 * a_tot_ms + bc_tot_ms)
            r_tur = _div(min_slice.sum(axis=0), a_tot_ms + turn_tot)
            r_slice = r_sor - r_tur
        else:
            sor_pair = _div(bc_full, 2.0 * shared_full + bc_full)
            tur_pair = _div(min_full, shared_full + min_full)
            total = float(np.nanmean(sor_pair) - np.nanmean(tur_pair))
            # per-pair (b + c) within each slice, averaged over pairs (the
            # treesliceR numerator `r_bc_sum`); NOT -shared, which would be
            # the wrong sign and magnitude
            bc_slice = pd_pairs_piece - shared_slice
            r_sor = _div(bc_slice, (2.0 * shared_full + bc_full)[:, None])
            r_tur = _div(min_slice, (shared_full + min_full)[:, None])
            r_slice = (r_sor - r_tur).mean(axis=0)

    return r_slice, float(total)


def r_compatible_profile(
    tree: TreeArray,
    mat,
    adj,
    slices: SliceStack,
    component: str = "sorensen",
    approach: str = "multisite",
    weighted: bool = False,
    multisite_domain: str = "mixed",
) -> dict[str, np.ndarray]:
    """Beta profiles normalised exactly like treesliceR (not re-scaled to 1).

    Exposed for cross-validation and for reproducing R figures: ``values`` are
    ``series / total`` with ``total`` built from the whole tree, so the turnover
    cumulative curve ends below 1 (``sum_j min_slice(j) / min_full``) and the
    nestedness curve above 1.  :func:`pb_per_slice` returns the profile version,
    whose rows sum to exactly 1.  The returned dict also carries ``series``
    undivided, which is what ``treesliceR::r_phylo(index = "PB")`` returns.
    """
    return _pb_per_slice_core(
        tree,
        mat,
        adj,
        slices,
        component=component,
        approach=approach,
        weighted=weighted,
        multisite_domain=multisite_domain,
        normalisation="r_compatible",
    )


def _pb_per_slice_core(
    tree: TreeArray,
    mat,
    adj,
    slices: SliceStack,
    component: str,
    approach: str,
    weighted: bool,
    multisite_domain: str,
    normalisation: str,
) -> dict[str, np.ndarray]:
    """Shared engine of :func:`pb_per_slice` / :func:`r_compatible_profile`."""
    if approach not in ("multisite", "pairwise"):
        raise ValueError("approach must be 'multisite' or 'pairwise'")
    multisite_domain = _normalise_domain(multisite_domain)
    mat = _match_columns_to_tree(tree, mat, quantity="PB")
    _assert_slices_match_tree(tree, slices)
    adj_arr = np.asarray(adj)
    profiles = _neighbourhood_profiles(tree, mat, adj_arr, slices, weighted=weighted)

    n_focal = adj_arr.shape[0]
    k = slices.n_slices
    values = np.full((n_focal, k), np.nan)
    series = np.full((n_focal, k), np.nan)
    totals = np.full(n_focal, np.nan)
    status = np.array(["ok"] * n_focal, dtype=object)
    norm = np.array([""] * n_focal, dtype=object)

    for i, res in enumerate(profiles):
        if res["status"] != "ok":
            status[i] = res["status"]
            continue
        if res["all_absent"] or res["all_present"]:
            status[i] = "degenerate_neighbourhood"
            continue
        if res["n_species"] < 2:
            status[i] = "fewer_than_two_species"
            continue
        try:
            series_row, total = _pb_series(
                res, component=component, approach=approach, domain=multisite_domain
            )
        except ZeroDivisionError:  # pragma: no cover
            status[i] = "zero_denominator"
            continue
        row_sum = float(np.sum(series_row))
        if (
            not np.isfinite(series_row).all()
            or not np.isfinite(total)
            or total == 0
            or not np.isfinite(row_sum)
            or row_sum == 0
        ):
            status[i] = "degenerate_neighbourhood"
            continue
        if normalisation == "r_compatible":
            values[i] = series_row / total
            norm[i] = "full_tree"
        else:
            values[i] = series_row / row_sum
            # "full_tree" = the shares already partitioned the whole-tree index
            # (b + c is additive over slices); "profile" = the denominator had to
            # be rebuilt from the profile because min(b, c) is not additive.
            norm[i] = "full_tree" if np.isclose(row_sum, total, rtol=1e-12, atol=0.0) else "profile"
        series[i] = series_row
        totals[i] = total
    return {
        "values": values,
        "series": series,
        "total": totals,
        "status": status,
        "normalisation": norm,
    }


def pb_per_slice(
    tree: TreeArray,
    mat,
    adj,
    slices: SliceStack,
    component: str = "sorensen",
    approach: str = "multisite",
    weighted: bool = False,
    multisite_domain: str = "mixed",
) -> dict[str, np.ndarray]:
    """Per-slice relative phylogenetic beta diversity per focal neighbourhood.

    The tree must already be aligned with ``mat`` and ``slices`` must have been
    built from that aligned tree (align once, then slice - see
    ``phyloslicer.io.align_tree_matrix``).

    Returns a dict with:

    * ``values``: (n_focal, k) per-slice share of the neighbourhood's beta,
      normalised so that **each finite row sums to exactly 1**;
    * ``series``: (n_focal, k) the same shares **before** any row normalisation,
      i.e. the raw per-slice index in absolute (root -> tip) units.  This is the
      quantity ``treesliceR::r_phylo(index = "PB")`` returns, and
      ``values == series / row_sum``;
    * ``total``: (n_focal,) full-tree beta index value, directly comparable to
      :func:`pb` (it is *not* re-scaled by the profile);
    * ``status``: (n_focal,) strings explaining skipped neighbourhoods;
    * ``normalisation``: (n_focal,) ``"full_tree"`` when the shares already
      partition the whole-tree index (sorensen, whose ``b + c`` numerator is
      additive over slices, and any neighbourhood where the two denominators
      coincide numerically), or ``"profile"`` when the denominator had to be
      rebuilt from ``sum_j series_j`` because ``min(b, c)`` is not additive over
      slices.  Use :func:`r_compatible_profile` for the treesliceR curve
      (``series / total``), whose turnover rows sum to ~0.92-1.0 and whose
      nestedness rows were measured to overshoot by up to 54% on the committed
      ``bd20__random`` reference scenario at 30 slices.

    ``multisite_domain`` selects the summation domain of the shared term: the
    default ``"mixed"`` (accepted spelling: ``"reference"``) is the treesliceR
    arithmetic, whose value grows with the neighbourhood size - ``m / (m + 2)``
    on the construction of one shared skeleton branch plus one private terminal
    branch per site, i.e. 0.5 at m = 2 up to 0.909 at m = 20 - while
    ``"paired"`` sums every term over pairs, matching the pairwise aggregation
    for m > 2 and leaving m = 2 untouched.  See :data:`MULTISITE_DOMAINS`.  Use
    :func:`r_compatible_profile` for the exact treesliceR curve.
    """
    return _pb_per_slice_core(
        tree,
        mat,
        adj,
        slices,
        component=component,
        approach=approach,
        weighted=weighted,
        multisite_domain=multisite_domain,
        normalisation="profile",
    )


# ---------------------------------------------------------------------- #
# total beta (no slicing)
# ---------------------------------------------------------------------- #
def pb(
    tree: TreeArray,
    mat,
    adj,
    component: str = "sorensen",
    approach: str = "multisite",
    weighted: bool = False,
    multisite_domain: str = "mixed",
) -> _pd.DataFrame:
    """Full-tree phylogenetic beta diversity per focal neighbourhood.

    Returns ``site_id, PB, status``.  ``status`` is ``"ok"`` or the reason the
    neighbourhood was refused (``fewer_than_two_sites``,
    ``degenerate_neighbourhood``, ``fewer_than_two_species``,
    ``zero_denominator``) with ``PB`` set to NaN - the same refusals
    :func:`pb_per_slice` applies, which is what makes that function's ``total``
    column directly comparable to this one.

    ``multisite_domain`` must match the value passed to :func:`pb_per_slice`
    for the two functions to agree (``"reference"`` and ``"mixed"`` are the same
    domain); see :data:`MULTISITE_DOMAINS` for what the two domains mean and why
    the default value is not comparable across neighbourhoods of different sizes.

    The domain changes the whole-tree index reported here for every
    component.  It changes :func:`pb_per_slice`'s ``total`` the same way, but
    not that function's ``values`` profile - and so not the fitted ``CpB`` of
    ``cpb_rate`` - for ``sorensen`` and ``turnover``, because dividing the
    series by its own row sum cancels the scalar denominator the domain
    controls.  Only ``nestedness`` moves the rate, as its two terms carry
    different denominators.
    """
    multisite_domain = _normalise_domain(multisite_domain)
    if approach not in ("multisite", "pairwise"):
        raise ValueError("approach must be 'multisite' or 'pairwise'")
    mat_df = as_matrix(mat)
    tree2, mat2 = align_tree_matrix(tree, mat_df)
    if component not in ("sorensen", "turnover", "nestedness"):
        raise ValueError("component must be 'sorensen', 'turnover' or 'nestedness'")
    adj_arr = np.asarray(adj)
    A_all = site_edge_membership(tree2, mat2)
    if weighted:
        ranges = np.asarray(A_all.sum(axis=0)).ravel()
        ranges = np.where(ranges == 0, 1.0, ranges)
        vec = tree2.lengths / ranges
    else:
        vec = tree2.lengths

    rows = []
    for i in range(adj_arr.shape[0]):
        members = np.flatnonzero(adj_arr[i] == 1)
        if members.size < 2:
            rows.append({"site_id": i, "PB": np.nan, "status": "fewer_than_two_sites"})
            continue
        A_loc = A_all[members]
        occupancy = float(np.asarray(A_loc.sum()).ravel()[0])
        dense_size = A_loc.shape[0] * A_loc.shape[1]
        # the same three refusals pb_per_slice() applies, so that its `total`
        # really is "directly comparable to pb()": an all-empty or all-identical
        # neighbourhood, or a one-species matrix, is reported as NaN with the
        # reason instead of one function answering 0 and the other NaN
        if occupancy == 0.0:
            rows.append({"site_id": i, "PB": np.nan, "status": "degenerate_neighbourhood"})
            continue
        if occupancy == dense_size:
            rows.append({"site_id": i, "PB": np.nan, "status": "degenerate_neighbourhood"})
            continue
        if mat2.shape[1] < 2:
            rows.append({"site_id": i, "PB": np.nan, "status": "fewer_than_two_species"})
            continue
        pd_full = np.asarray(A_loc @ vec).ravel()
        ii, jj = np.triu_indices(members.size, k=1)
        indicator = _pair_indicator(ii, jj, members.size)
        A_pairs = ((indicator @ A_loc) > 0).astype(np.float64)
        pd_pairs = np.asarray(A_pairs @ vec).ravel()
        shared = pd_full[ii] + pd_full[jj] - pd_pairs
        a_ij = shared
        bc = pd_pairs - shared
        mn = pd_pairs - np.maximum(pd_full[ii], pd_full[jj])
        if approach == "multisite":
            if multisite_domain == "paired":
                a_tot = float(a_ij.sum())
            else:
                a_tot = float(
                    pd_full.sum()
                    - float((np.asarray(A_loc.sum(axis=0)).ravel() > 0).astype(float) @ vec)
                )
        if component == "sorensen":
            if approach == "multisite":
                bc_tot = float(bc.sum())
                val = bc_tot / (2.0 * a_tot + bc_tot)
            else:
                val = float(np.nanmean(bc / (2.0 * a_ij + bc)))
        elif component == "turnover":
            if approach == "multisite":
                turn = float(mn.sum())
                val = turn / (a_tot + turn)
            else:
                val = float(np.nanmean(mn / (a_ij + mn)))
        elif component == "nestedness":
            if approach == "multisite":
                bc_tot = float(bc.sum())
                turn = float(mn.sum())
                val = bc_tot / (2.0 * a_tot + bc_tot) - turn / (a_tot + turn)
            else:
                sor = np.nanmean(bc / (2.0 * a_ij + bc))
                tur = np.nanmean(mn / (a_ij + mn))
                val = float(sor - tur)
        else:
            raise ValueError("component must be 'sorensen', 'turnover' or 'nestedness'")
        if not np.isfinite(val):
            rows.append({"site_id": i, "PB": np.nan, "status": "zero_denominator"})
            continue
        rows.append({"site_id": i, "PB": float(val), "status": "ok"})
    return _pd.DataFrame(rows)
