#!/usr/bin/env python
"""Generate larger benchmark inputs (trees + presence matrices) for the
cross-language scaling benchmark.  Outputs benchmarks/large/{tree,mat} files.
"""

from __future__ import annotations

from pathlib import Path

from phyloslicer.io.newick import write_newick
from phyloslicer.simulate import make_dataset, yule_tree


def main(out_dir: str = "benchmarks/large") -> None:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    for n_tips in (1000, 2000):
        tree = yule_tree(n_tips, age=1.0, seed=100 + n_tips)
        (out / f"large{n_tips}__tree.tre").write_text(write_newick(tree) + "\n")
        mat, _ = make_dataset(
            n_sites=250, n_species=n_tips, structure="clustered", seed=200 + n_tips
        )
        mat.columns = tree.tip_labels
        mat.index = [f"s{i + 1}" for i in range(mat.shape[0])]
        mat.to_csv(out / f"large{n_tips}__mat.csv")
        print("wrote", n_tips, "tips,", mat.shape[0], "sites")


if __name__ == "__main__":
    main()
