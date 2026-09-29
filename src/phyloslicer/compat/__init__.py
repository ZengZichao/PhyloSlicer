"""treesliceR-compatible thin wrappers (eases migration from R).

Signatures *and* return conventions follow the R package, including its
quirks: ``squeeze_root``/``squeeze_tips`` accept ``criterion="my"/"pd"`` and
``dropNodes``; ``phylo_pieces`` uses the R ``n`` + ``method`` pair and can
return the ``timeSteps`` vector; ``prune_tips`` distinguishes ``method=1`` (keep
tips with long terminal branches) from ``method=2`` (drop them).

Where the main API deliberately differs from R, this layer re-establishes the R
behaviour rather than leaking the Python convention:

* :func:`phylo_pieces` with ``timeSteps=True`` returns R's increasing
  depth-from-root vector (``cumsum(n, ...)``), **not** the main API's
  decreasing "Ma before present" slice boundaries;
* :func:`squeeze_int` with ``invert=True`` returns R's ``(younger, older)``
  order, **not** the main API's ``(older, younger)``.

Both are asserted by ``tests/test_compat.py``.  Arguments that R validates are
validated loudly here instead of silently falling through to a different
branch: an unknown ``method`` or ``criterion`` raises :class:`ValueError`
(R returns ``NULL`` / silently skips the cut, which is worse for a scripted
migration than an error).
"""

from __future__ import annotations

from ..core.tree import TreeArray
from ..slicing.api import prune_tips as _prune_tips
from ..slicing.api import slice_interval as _slice_interval
from ..slicing.api import slice_pieces as _slice_pieces
from ..slicing.api import slice_rootward as _slice_rootward
from ..slicing.api import slice_tipward as _slice_tipward

__all__ = ["squeeze_root", "squeeze_tips", "squeeze_int", "phylo_pieces", "prune_tips"]

#: R's documented spellings plus the lowercase ones the R code actually tests.
_CRITERION_MAP = {
    "my": "time",
    "million years": "time",
    "time": "time",
    "pd": "pd",
}


def _criterion(criterion: str) -> str:
    """Map an R ``criterion`` string onto the Python one, failing loudly.

    R documents ``"PD"`` but its code only tests the lowercase ``"pd"``, so
    ``criterion = "PD"`` silently degenerates into "no cut at all" there.  Here
    the mapping is case-insensitive over the documented spellings and anything
    else is an error rather than a silent no-op.
    """
    if not isinstance(criterion, str):
        raise ValueError(f"criterion must be a string, got {type(criterion).__name__}")
    key = criterion.strip().lower()
    if key not in _CRITERION_MAP:
        raise ValueError(
            f"unknown criterion {criterion!r}; treesliceR accepts 'my' (million "
            "years) or 'PD' (accumulated phylogenetic diversity)"
        )
    return _CRITERION_MAP[key]


def squeeze_root(
    tree: TreeArray,
    time: float,
    criterion: str = "my",
    dropNodes: bool = False,
    allow_non_ultrametric: bool = False,
):
    """treesliceR ``squeeze_root``: keep the crown slice of width ``time``.

    ``time = 0`` returns the whole tree's crown slice of zero width (an empty
    slice) and ``time = total PD`` / ``time = tree depth`` the complete tree,
    exactly as R does.
    """
    return _slice_rootward(
        tree,
        time,
        criterion=_criterion(criterion),
        collapse=dropNodes,
        allow_non_ultrametric=allow_non_ultrametric,
    )


def squeeze_tips(
    tree: TreeArray,
    time: float,
    criterion: str = "my",
    dropNodes: bool = False,
    allow_non_ultrametric: bool = False,
):
    """treesliceR ``squeeze_tips``: drop the youngest ``time`` of the tree."""
    return _slice_tipward(
        tree,
        time,
        criterion=_criterion(criterion),
        collapse=dropNodes,
        allow_non_ultrametric=allow_non_ultrametric,
    )


def squeeze_int(
    tree: TreeArray,
    from_: float,
    to: float,
    invert: bool = False,
    criterion: str = "my",
    dropNodes: bool = False,
    allow_non_ultrametric: bool = False,
):
    """treesliceR ``squeeze_int`` (R's ``from`` is written ``from_`` in Python).

    With ``invert=True`` the return follows R's convention:
    ``[squeeze_root(to), squeeze_tips(from)]``, i.e. the **younger** complement
    first and the **older** complement second.  The main
    :func:`phyloslicer.slice_interval` returns the same two trees in the
    opposite ``(older, younger)`` order, so code ported from R should use this
    wrapper (or swap explicitly).
    """
    crit = _criterion(criterion)
    interior = _slice_interval(
        tree,
        from_,
        to,
        invert=invert,
        criterion=crit,
        collapse=dropNodes,
        allow_non_ultrametric=allow_non_ultrametric,
    )
    if invert:
        older, younger = interior
        return [younger, older]  # treesliceR: list(tree1 = root slice, tree2 = tip slice)
    return interior


def phylo_pieces(
    tree: TreeArray,
    n: int,
    criterion: str = "my",
    method: int = 1,
    timeSteps: bool = False,
    dropNodes: bool = False,
    returnTree: bool = False,
    timeSteps_as_depth: bool = True,
):
    """treesliceR ``phylo_pieces`` with R's ``n`` + ``method`` argument pair.

    Returns a list of slice trees, ordered root -> tips (as R's list is), or a
    ``[pieces, timeSteps(, tree)]`` list when ``timeSteps``/``returnTree`` are
    requested.

    ``timeSteps_as_depth`` keeps the R vector semantics: treesliceR builds its
    time steps as ``age <- cumsum(rep(n, nslices))``, i.e. the *increasing*
    cumulative depth measured from the root, and ``phylo_pieces`` returns that
    vector.  The default ``True`` returns the same thing
    (``stack.depth_windows[:, 1]`` - the tipward bound of each slice, so
    ``age[j]`` is the cumulative depth reached *by* piece ``j`` and the last
    element is the full tree depth).  Index-wise alignment with the returned
    pieces is therefore identical to R's: ``phylo_pieces(tr, n=4)`` on
    ``((A:4,B:4):4,(C:4,D:4):4);`` (``T=8``) yields ``[2, 4, 6, 8]``, exactly
    ``cumsum(rep(2, 4))``.  Under ``criterion="pd"`` R instead stores the
    time-mapped equal-PD cut points (``phylo_pieces.R:217-238``,
    ``AGE_t1 = age``), which are again the tipward bounds, so the same rule
    applies (``[3, 5, 6.5, 8]`` here).  Pass ``False`` to get the main API's
    :attr:`SliceStack.ages` instead (the older boundary of each slice in Ma
    before present, i.e. decreasing) - which is what
    :func:`phyloslicer.rate_profile` / the rate fitters regress against.
    """
    if method not in (1, 2):
        raise ValueError(
            f"method must be 1 (n is the number of pieces) or 2 (n is the piece "
            f"width); got {method!r}.  treesliceR returns NULL for any other value."
        )
    crit = _criterion(criterion)
    if method == 2:
        stack = _slice_pieces(tree, width=n, criterion=crit)
    else:
        stack = _slice_pieces(tree, n=n, criterion=crit)
    pieces = stack.to_sliced_trees(collapse=dropNodes)
    steps = stack.depth_windows[:, 1].copy() if timeSteps_as_depth else stack.ages.copy()
    if timeSteps and returnTree:
        return [pieces, steps, tree]
    if timeSteps:
        return [pieces, steps]
    if returnTree:
        return [pieces, tree]
    return pieces


def prune_tips(
    tree: TreeArray,
    time,
    qtl: bool = False,
    method: int = 1,
    allow_non_ultrametric: bool = False,
):
    """treesliceR ``prune_tips``: ``method=1`` keeps tips whose terminal branch
    is >= the threshold, ``method=2`` drops them.

    Any other ``method`` raises.  R's implementation neither validates nor
    falls through - it simply reaches the end of the function body and returns
    ``NULL``, which in a scripted migration turns into a confusing downstream
    error.
    """
    if method not in (1, 2):
        raise ValueError(
            f"method must be 1 or 2, got {method!r} (treesliceR returns NULL for any other value)"
        )
    side = "after" if method == 1 else "before"
    return _prune_tips(
        tree, time, quantiles=qtl, side=side, allow_non_ultrametric=allow_non_ultrametric
    )
