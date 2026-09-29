#!/usr/bin/env python
"""Generate deterministic test data for PhyloSlicer testing suite.

All data is written to testing_suite/data/ and is fully reproducible from the
seeds used below.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Ensure the package is importable
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import phyloslicer as ps
from phyloslicer.simulate import birth_death_tree, to_newick_file, yule_tree

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
DATA_DIR.mkdir(parents=True, exist_ok=True)


def gen_small_tree():
    """4-tip symmetric tree with known analytic properties."""
    tree = ps.TreeArray.from_newick("((A:1,B:1):1,(C:1,D:1):1);")
    to_newick_file(tree, DATA_DIR / "small_tree.tre")
    return tree


def gen_medium_tree():
    """20-tip Yule tree, age=10, seed=42."""
    tree = yule_tree(20, age=10.0, seed=42)
    to_newick_file(tree, DATA_DIR / "medium_tree.tre")
    return tree


def gen_large_tree():
    """100-tip birth-death tree, age=10, seed=42."""
    tree = birth_death_tree(100, age=10.0, birth_rate=1.0, death_rate=0.2, seed=42)
    to_newick_file(tree, DATA_DIR / "large_tree.tre")
    return tree


def gen_matrix(tree, n_sites=30, seed=7, structure="clustered"):
    """Generate a presence/absence matrix aligned to the given tree."""
    rng = np.random.default_rng(seed)
    mat = (rng.random((n_sites, tree.n_tips)) < 0.3).astype(float)
    df = pd.DataFrame(mat, columns=tree.tip_labels)
    df.index = [f"site_{i + 1}" for i in range(n_sites)]
    # Ensure every species occurs at least once
    empty = df.sum(axis=0) == 0
    if empty.any():
        cols = df.columns[empty.to_numpy()]
        df.loc[df.index[: len(cols)], cols] = 1.0
    name = (
        f"{'medium' if tree.n_tips == 20 else 'large' if tree.n_tips == 100 else 'small'}_mat.csv"
    )
    df.to_csv(DATA_DIR / name)
    return df


def gen_posterior_trees():
    """5 posterior trees (50 tips each, slightly different ages)."""
    trees_dir = DATA_DIR / "posterior_trees"
    trees_dir.mkdir(exist_ok=True)
    rng = np.random.default_rng(99)
    for i in range(5):
        age = float(rng.uniform(9.5, 10.5))
        tree = yule_tree(50, age=age, seed=100 + i)
        to_newick_file(tree, trees_dir / f"post_{i + 1:04d}.tre")
    # Also a multi-tree Newick file
    trees = [yule_tree(50, age=10.0, seed=200 + i) for i in range(5)]
    nexus_text = ps.io.write_nexus(trees)
    (DATA_DIR / "posterior.nex").write_text(nexus_text)
    # Multi-tree Newick
    newick_text = "\n".join(t.to_newick() for t in trees)
    (DATA_DIR / "posterior.tre").write_text(newick_text + "\n")
    return trees


def gen_grid_adjacency(n_sites=30, grid_side=6):
    """Queen adjacency for a grid."""
    n = min(n_sites, grid_side * grid_side)
    xs = np.arange(grid_side, dtype=float)
    gx, gy = np.meshgrid(xs, xs)
    coords = pd.DataFrame({"x": gx.ravel()[:n], "y": gy.ravel()[:n]})
    bounds = pd.DataFrame(
        {
            "minx": coords["x"],
            "miny": coords["y"],
            "maxx": coords["x"] + 1.0,
            "maxy": coords["y"] + 1.0,
        }
    )
    bounds.to_csv(DATA_DIR / "grid_bounds.csv", index=False)
    adj = ps.spatial.adjacency_from_bounds(bounds.to_numpy(), method="queen")
    np.savetxt(DATA_DIR / "grid_adj.csv", adj, delimiter=",", fmt="%d")
    return adj, bounds


def gen_sim_passerines():
    """Use the bundled sim_passerines dataset (small scale)."""
    data = ps.datasets.sim_passerines(n_sites=50, n_species=50, n_trees=3, seed=42)
    to_newick_file(data["trees"][0], DATA_DIR / "passerines_tree.tre")
    data["mat"].to_csv(DATA_DIR / "passerines_mat.csv")
    data["bounds"].to_csv(DATA_DIR / "passerines_bounds.csv", index=False)
    np.savetxt(DATA_DIR / "passerines_adj.csv", data["adj"], delimiter=",", fmt="%d")
    # Save posterior trees
    pdir = DATA_DIR / "passerines_posterior"
    pdir.mkdir(exist_ok=True)
    for i, t in enumerate(data["trees"]):
        to_newick_file(t, pdir / f"tree_{i + 1:04d}.tre")
    return data


def main():
    print("=== Generating PhyloSlicer test data ===")
    print(f"Output directory: {DATA_DIR}")

    print("1. Small tree (4 tips)...")
    t1 = gen_small_tree()

    print("2. Medium tree (20 tips)...")
    t2 = gen_medium_tree()
    gen_matrix(t2, n_sites=30, seed=7)

    print("3. Large tree (100 tips)...")
    t3 = gen_large_tree()
    gen_matrix(t3, n_sites=50, seed=7)

    print("4. Posterior tree set...")
    gen_posterior_trees()

    print("5. Grid adjacency...")
    gen_grid_adjacency()

    print("6. Sim passerines dataset...")
    gen_sim_passerines()

    print("7. Small tree matrix...")
    gen_matrix(t1, n_sites=4, seed=3)

    print("\n=== Test data generation complete ===")
    print(f"Files in {DATA_DIR}:")
    for f in sorted(DATA_DIR.iterdir()):
        if f.is_file():
            print(f"  {f.name} ({f.stat().st_size} bytes)")
    for d in sorted(DATA_DIR.iterdir()):
        if d.is_dir():
            n = sum(1 for _ in d.iterdir())
            print(f"  {d.name}/ ({n} files)")


if __name__ == "__main__":
    main()
