"""Out-of-the-box PhyloSlicer run: no external data, no setup.

    python examples/quickstart.py            # 208 cells, 308 species, 12 trees
    python examples/quickstart.py --slices 25

Everything it reads ships inside the package (`phyloslicer.datasets`), so a
reader can run the worked example from the user guide with one command.
CI executes this file, which is what keeps the documented example honest.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from phyloslicer import cpd_rate, posterior_rate, slice_pieces  # noqa: E402
from phyloslicer.datasets import sim_passerines  # noqa: E402


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sites", type=int, default=208)
    ap.add_argument("--species", type=int, default=308)
    ap.add_argument("--trees", type=int, default=12)
    ap.add_argument("--seed", type=int, default=7)
    ap.add_argument("--slices", type=int, default=50)
    args = ap.parse_args(argv)

    d = sim_passerines(
        n_sites=args.sites, n_species=args.species, n_trees=args.trees, seed=args.seed
    )
    tree, mat = d["trees"][0], d["mat"]
    print(f"tree root age: {tree.root_age:.6f} (Ma before present)")

    sliced = slice_pieces(tree, n=args.slices, criterion="time")
    C = sliced.contribution_matrix()
    print(
        f"sliced into {sliced.n_slices} windows; contribution matrix {C.shape}"
        f" with {int((C > 0).sum())} non-zero cells"
    )

    rates = cpd_rate(tree, mat, n_slices=args.slices)
    cols = ["site_id", "CpD", "PD", "pDO", "converged", "r_squared", "origin_within_tree"]
    print(f"\nCpD fit: {int(rates['converged'].sum())} / {len(rates)} assemblages converged")
    print(rates[cols].head(3).to_string(index=False))
    extrapolated = int((~rates["origin_within_tree"]).sum())
    if extrapolated:
        print(
            f"note: {extrapolated} of {len(rates)} origin times fall before the root "
            f"(-ln(0.05)/rate > {tree.root_age:.3f}); they are model extrapolations, "
            "flagged by origin_within_tree=False"
        )

    post = posterior_rate(d["trees"], mat, rate="cpd", n_slices=args.slices)
    print("\nPosterior across the 12 pseudo-posterior trees:")
    print(
        post[["site_id", "rate_median", "hdi_low", "hdi_high", "n_trees"]]
        .head(3)
        .to_string(index=False)
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
