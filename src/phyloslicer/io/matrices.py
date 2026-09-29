"""Presence/absence matrix construction and tree-matrix alignment."""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd
import scipy.sparse as sp

from ..core.tree import TreeArray

__all__ = ["presence_matrix", "align_tree_matrix", "as_matrix"]


def as_matrix(mat) -> pd.DataFrame:
    """Coerce a presence/absence input (DataFrame / ndarray / sparse) to a DataFrame.

    Rows are sites, columns are species.  Column labels are always coerced
    to ``str`` so downstream tree-matrix alignment (which compares column
    labels against string tip labels) works regardless of the input type:
    plain arrays get integer column labels, sparse matrices likewise, and
    DataFrames with numeric column names (e.g. species IDs read from CSV)
    are normalised here too.
    """
    if isinstance(mat, pd.DataFrame):
        out = mat.copy()
        out.columns = out.columns.astype(str)
        return out
    if sp.issparse(mat):
        dense = np.asarray(mat.todense())
    else:
        dense = np.asarray(mat)
    if dense.ndim == 1:
        dense = dense[None, :]
    out = pd.DataFrame(dense.astype(float))
    out.columns = out.columns.astype(str)
    return out


def presence_matrix(
    df: pd.DataFrame,
    site_col: str | None = None,
    species_col: str | None = None,
    value_col: str | None = None,
    site_index_col: bool = False,
) -> pd.DataFrame:
    """Build a sites x species presence/absence DataFrame from a (long) table.

    Parameters
    ----------
    df : DataFrame
        Either already a sites x species wide table (returned after light
        validation) or a long-format table with site/species columns.
    site_col, species_col, value_col : str, optional
        Long-format column names. If omitted, the first three columns are used
        (and ``value_col`` defaults to presence everywhere).
    site_index_col : bool
        Treat the DataFrame index as the site identifier in long format.

    Returns
    -------
    DataFrame with sites as rows, species as columns, values in {0, 1}.
    """
    if site_col is None and species_col is None:
        wide = df.copy()
        wide.index = wide.index.astype(str)
        wide.columns = wide.columns.astype(str)
        wide = (wide > 0).astype(float)
        return wide

    site_col = site_col or df.columns[0]
    species_col = species_col or df.columns[1]
    if value_col is None and len(df.columns) > 2:
        # an unrequested third column used to be read as an abundance by
        # position, so a latitude or elevation column (values > 0) turned every
        # record into a presence and silently fabricated the matrix
        warnings.warn(
            f"{df.columns[2]!r} and any later columns are ignored: pass "
            "value_col= explicitly to code presences from an abundance column",
            stacklevel=2,
        )

    if site_index_col:
        df = df.reset_index(names="__site__")
        site_col = "__site__"

    if value_col is None:
        values = pd.Series(1.0, index=df.index)
    else:
        values = (df[value_col] > 0).astype(float)

    wide = df.assign(__v__=values).pivot_table(
        index=site_col, columns=species_col, values="__v__", aggfunc="max", fill_value=0.0
    )
    wide.index = wide.index.astype(str)
    wide.columns = wide.columns.astype(str)
    return wide


def align_tree_matrix(tree: TreeArray, mat: pd.DataFrame) -> tuple[TreeArray, pd.DataFrame]:
    """Align a tree and a presence matrix, reporting every adjustment.

    Mirrors treesliceR's implicit ``keep.tip`` behaviour but *explicitly*:
    * species columns with zero occurrences are dropped (with a warning);
    * tree tips absent from the matrix are pruned (with a warning);
    * matrix species missing from the tree are dropped (with a warning);
    * the matrix rows/columns are reordered to (site order, tree tip order).
    """
    mat = as_matrix(mat)
    zero_cols = mat.columns[mat.sum(axis=0) == 0]
    if len(zero_cols) > 0:
        warnings.warn(
            f"removing {len(zero_cols)} species without any occurrence: "
            f"{', '.join(map(str, zero_cols[:5]))}{'...' if len(zero_cols) > 5 else ''}",
            stacklevel=2,
        )
        mat = mat.drop(columns=list(zero_cols))

    mat_species = set(mat.columns.astype(str))
    tree_species = set(tree.tip_labels)

    # fewer than two shared names -> nothing can be aligned; fail with the real
    # reason rather than letting the MRCA prune raise the misleading
    # "a subtree needs at least 2 tips"
    shared = mat_species & tree_species
    if len(shared) < 2:
        named = f": {next(iter(shared))!r}" if len(shared) == 1 else ""
        raise ValueError(
            f"only {len(shared)} of the {len(tree_species)} tree tip labels "
            f"occurs in the matrix's {len(mat_species)} columns{named}; a "
            "phylogenetic analysis needs at least two.  Check spelling, case and "
            "any tip-label prefix - align_tree_matrix prunes to the intersection."
        )

    extra_matrix = mat_species - tree_species
    if extra_matrix:
        warnings.warn(
            f"removing {len(extra_matrix)} matrix species absent from the tree: "
            f"{', '.join(sorted(extra_matrix)[:5])}{'...' if len(extra_matrix) > 5 else ''}",
            stacklevel=2,
        )
        mat = mat.drop(columns=list(extra_matrix))

    missing_from_matrix = tree_species - mat_species
    if missing_from_matrix:
        warnings.warn(
            f"pruning {len(missing_from_matrix)} tree tips absent from the matrix: "
            f"{', '.join(sorted(missing_from_matrix)[:5])}"
            f"{'...' if len(missing_from_matrix) > 5 else ''}",
            stacklevel=2,
        )
        keep_labels = [t for t in tree.tip_labels if t in mat_species]
        tree = tree.subarray(keep_labels)

    mat = mat[[t for t in tree.tip_labels if t in mat.columns]]
    # if some tips remain unrepresented (shouldn't happen), pad with zeros
    for t in tree.tip_labels:
        if t not in mat.columns:
            mat[t] = 0.0
    mat = mat[tree.tip_labels]
    return tree, mat
