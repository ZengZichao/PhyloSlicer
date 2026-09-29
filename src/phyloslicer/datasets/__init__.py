"""Bundled example data: deterministic simulations only.

There is **no** download facility in this release.  The real treesliceR
tutorial data (Australian passerines, MIT licensed) can be exported from R
with a custom export script; :func:`fetch_tutorial_data` exists only to say so
loudly, because a placeholder that raises ``RuntimeError`` while the prose
promises "fetched on demand" is worse than an explicit
``NotImplementedError``.  Everything here runs offline and is reproducible from
the seed.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from ..simulate import make_dataset, yule_tree
from ..spatial import adjacency_from_bounds

__all__ = ["sim_passerines", "fetch_tutorial_data"]


def sim_passerines(
    n_sites: int = 208, n_species: int = 308, n_trees: int = 5, seed: int = 42
) -> dict:
    """Simulation mirroring the structure and scale of the treesliceR
    Australian-passerine case study (0.5-degree-style grid, 308 species,
    small posterior pool).

    It is *structurally* analogous, not numerically comparable: the R original
    ships 100 posterior trees (hence the small default ``n_trees = 5`` keeps
    examples fast rather than faithful), and the tree, matrix and coordinates
    are simulated.  Any figure produced from this bundle must not be presented
    as reproducing a published treesliceR result.
    """
    rng = np.random.default_rng(seed)
    tree = yule_tree(n_species, age=1.0, seed=seed)
    mat, coords = make_dataset(
        n_sites=n_sites, n_species=n_species, structure="clustered", seed=seed
    )
    # species columns must carry the tree's tip labels so tree and matrix can
    # be aligned directly (make_dataset numbers columns "1".."n_species")
    mat.columns = tree.tip_labels
    # bounds of a regular grid (0.5-degree cells)
    side = int(np.ceil(np.sqrt(n_sites)))
    gx, gy = np.meshgrid(np.arange(side, dtype=float), np.arange(side, dtype=float))
    bounds = pd.DataFrame(
        {
            "minx": 112.0 + gx.ravel()[:n_sites] * 0.5,
            "miny": -44.0 + gy.ravel()[:n_sites] * 0.5,
            "maxx": 112.5 + gx.ravel()[:n_sites] * 0.5,
            "maxy": -43.5 + gy.ravel()[:n_sites] * 0.5,
        }
    )
    adj = adjacency_from_bounds(bounds.to_numpy(), method="queen")
    trees = [tree] + [
        yule_tree(n_species, age=float(rng.uniform(0.9, 1.1)), seed=seed + i + 1)
        for i in range(n_trees - 1)
    ]
    return {"trees": trees, "mat": mat, "coords": coords, "bounds": bounds, "adj": adj}


def fetch_tutorial_data(dest: str = ".", confirm: bool = True) -> None:
    """Not implemented: PhyloSlicer bundles no remote tutorial data.

    Always raises :class:`NotImplementedError`.  No archive has been deposited
    and no DOI minted, so nothing could be fetched even if this function did
    something; :func:`sim_passerines` provides the offline stand-in.  This stub
    is kept (instead of silently deleted) so that any call tells the user the
    truth rather than returning the wrong data.
    """
    raise NotImplementedError(
        "phyloslicer.datasets ships no remote tutorial bundle and "
        "fetch_tutorial_data is not implemented; use "
        "sim_passerines() for a deterministic offline example, or export the "
        "treesliceR data from R yourself (data(pass_mat), data(pass_trees), "
        "data(AU_grid), data(AU_adj))"
    )
