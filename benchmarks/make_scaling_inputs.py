#!/usr/bin/env python
"""Generate the intermediate tree sizes for the fixed-k scaling benchmark.

The original cross-language benchmark (``benchmark_r.R`` /
``benchmark_phyloslicer.py``) covers 20- and 100-tip trees at k = 25/50 and
``make_large_inputs.py`` covers 1000 and 2000 tips at k = 100 only, so the
speed-up was never measured at a constant number of slices across tree
sizes.  This script supplies the missing 200- and 500-tip inputs, using the
same generator and seeding convention as ``make_large_inputs.py``
(``yule_tree(n, age=1, seed=100+n)``).

Usage (from the project root):
  python benchmarks/make_scaling_inputs.py
Writes benchmarks/scaling/large{n}__tree.tre (+ a 30-site matrix per size).
"""

from __future__ import annotations

from pathlib import Path

from phyloslicer.io.newick import write_newick
from phyloslicer.simulate import make_dataset, yule_tree

SIZES = (200, 500)


def main(out_dir: str = "benchmarks/scaling") -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for n_tips in SIZES:
        tree = yule_tree(n_tips, age=1.0, seed=100 + n_tips)
        (out / f"large{n_tips}__tree.tre").write_text(write_newick(tree) + "\n")
        mat, _ = make_dataset(
            n_sites=30, n_species=n_tips, structure="clustered", seed=200 + n_tips
        )
        mat.columns = tree.tip_labels
        mat.index = [f"s{i + 1}" for i in range(mat.shape[0])]
        mat.to_csv(out / f"large{n_tips}__mat.csv")
        print("wrote", n_tips, "tips,", mat.shape[0], "sites")


if __name__ == "__main__":
    main()
