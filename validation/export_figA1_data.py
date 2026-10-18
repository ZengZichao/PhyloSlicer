#!/usr/bin/env python
"""Regenerate figA1_panel_data.csv from the archived validation inputs.

The counts behind panels (a) and (b) of Fig. A1 of the Supporting
information, plus the posterior medians behind panel (c):

  (a) CpD rates of the 24 assemblages that converged after the content of
      six assemblages of the 100-tip birth-and-death *clustered* validation
      scenario (bd100__clustered) had been set to zero (six rows drawn with
      ``numpy.random.default_rng(3).choice(mat.index, 6, replace=False)`` on the
      matrix's own 1-based site_id, matching the figure's draw exactly), at
      n_slices = 50;
  (b) CpD rates of six assemblages of the intact scenario (six rows drawn
      with numpy.random.default_rng(1).choice(30, 6, replace=False); the
      figure annotates them 1-based, i.e. 2, 5, 12, 14, 21 and 27) at
      n_slices = 10, 25, 50, 100 and 200;
  (b_all) the median CpD rate over all 30 assemblages of the intact scenario
      at the same five slice counts - the series behind the Appendix S1
      sentence "the median ... rises from 2.38 to 3.21" (distributed as
      rows since review P1-3; the percentage is computed from these
      unrounded values, and the rounded ones give 34.9 %);
  (c) posterior-median CpD rates of the first 20 assemblages of the bundled
      worked example (sim_passerines, 12 pseudo-posterior trees), plus the
      5-95% band over all 208 assemblages that Fig. A1c draws in grey.

Run it with no setup, from beside the file it writes:

    python export_figA1_data.py

It finds the package in this order: an installed ``phyloslicer``, then
``PHYLOSLICER_SRC``, then the sibling checkout ``PhyloSlicer-项目代码/src``.
Setting PHYLOSLICER_SRC remains supported for a checkout kept elsewhere.

Run it from the project root with ``python validation/export_figA1_data.py``;
it writes ``validation/figA1_panel_data.csv`` beside itself.  The rows are the
machine-readable form of what Fig. A1's caption asserts, and they are tabulated
in full as Table S5 of the Supporting information, so the figure's numbers can be
read on the page and re-derived from this script.

number that Fig. A1's caption asserts can be re-derived from this file:
the caption's 33.0-37.3 % rise across panel (b) is the per-assemblage
increase of ``value`` between slice_count 10 and 200, and the caption's
2.38 -> 3.21 median rise is the median over all 30 assemblages of the
intact scenario between the same two slice counts, now itself a row
series (panel b_all) of this file and recomputed by this script's
self-check.
"""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
os.chdir(HERE)
CHECKOUT_NAME = "PhyloSlicer-项目代码"


def _src_candidates():
    """Directories that may hold the ``phyloslicer`` package or its checkout."""
    cands = []
    env = os.environ.get("PHYLOSLICER_SRC")
    if env:
        pth = Path(env).expanduser()
        cands += [pth, pth / "src", pth.parent]
    try:
        spec = importlib.util.find_spec("phyloslicer")
    except (ImportError, ValueError):  # a broken entry on sys.path
        spec = None
    if spec is not None:
        locs = list(spec.submodule_search_locations or [])
        if spec.origin:
            locs.append(str(Path(spec.origin).parent))
        for loc in locs:
            pkg = Path(loc)
            cands += [pkg, pkg.parent, pkg.parent.parent]
    # this script ships inside the checkout, so its own parent is the root
    cands += [HERE.parent / "src", HERE.parent]
    cands += [HERE.parent / CHECKOUT_NAME / "src", HERE.parent / CHECKOUT_NAME]
    cands += sorted(q for pth in HERE.parent.glob(CHECKOUT_NAME + "*") for q in (pth / "src", pth))
    seen, out = set(), []
    for c in cands:
        if str(c) not in seen:
            seen.add(str(c))
            out.append(c)
    return out


def _ensure_phyloslicer() -> None:
    """Put ``phyloslicer`` on sys.path without hard-coding where it lives."""
    if importlib.util.find_spec("phyloslicer") is not None:
        return
    for cand in _src_candidates():
        if (cand / "phyloslicer").is_dir() and str(cand) not in sys.path:
            sys.path.insert(0, str(cand))
            if importlib.util.find_spec("phyloslicer") is not None:
                return
    raise RuntimeError(
        "cannot import phyloslicer.  Install it (pip install -e <checkout>) or "
        "point PHYLOSLICER_SRC at the checkout's src/ directory, e.g.  "
        "export PHYLOSLICER_SRC=/path/to/" + CHECKOUT_NAME + "/src"
    )


def _project_root(*required):
    """The checkout that ships the archived data files named in ``required``."""
    for cand in _src_candidates():
        root = cand.parent if cand.name == "src" else cand
        for probe in (cand, root):
            if all((probe / Path(*req)).is_file() for req in required):
                return probe
    raise RuntimeError(
        "cannot locate the PhyloSlicer source tree holding "
        + " and ".join("/".join(r) for r in required)
        + ".  Set PHYLOSLICER_SRC "
        "to a checkout's src/ directory, e.g.  export PHYLOSLICER_SRC=/path/to/"
        + CHECKOUT_NAME
        + "/src"
    )


_ensure_phyloslicer()

from phyloslicer import cpd_rate, posterior_rate  # noqa: E402
from phyloslicer.datasets import sim_passerines  # noqa: E402
from phyloslicer.io import presence_matrix, read_newick  # noqa: E402

PROJ = _project_root(
    ("validation", "reference", "bd100__tree.tre"),
    ("validation", "reference", "bd100__clustered__mat.csv"),
)
REFERENCE_DIR = Path(os.environ.get("PHYLOSLICER_REFERENCE_DIR", PROJ / "validation" / "reference"))


def main() -> int:
    tree = read_newick(str(REFERENCE_DIR / "bd100__tree.tre"))
    df = pd.read_csv(REFERENCE_DIR / "bd100__clustered__mat.csv", index_col=0)
    mat0 = presence_matrix(df, site_index_col=True)
    # Panel (a) must be the SAME draw as the figure, not merely a similar one:
    # PhyloSlicer-手稿附图/PhyloSlicer-附件图片-1.py seeds 3 and picks six of the
    # matrix's own 1-based site_id labels, so this draws on mat0.index too.  The
    # old code seeded 0 and picked 0-based positions, which emptied six different
    # rows and left Table S5 contradicting Fig. A1(a) row by row.
    zeroed_labels = list(
        np.random.default_rng(3).choice(np.asarray(mat0.index), size=6, replace=False)
    )
    # the site index is string-typed ('1'..'30'), so the drawn labels are kept
    # verbatim for .loc and only converted for display
    zeroed = sorted(int(x) for x in zeroed_labels)
    sampled = sorted(int(x) for x in np.random.default_rng(1).choice(30, 6, replace=False))

    rows = []
    mata = mat0.copy()
    mata.loc[zeroed_labels, :] = 0.0  # by label, matching the figure exactly
    ra = cpd_rate(tree, mata, n_slices=50)
    for _, r in ra.iterrows():
        rows.append(
            {
                "panel": "a",
                "series": "cpd_rate",
                "assemblage": int(r["site_id"]),
                "slice_count": 50,
                "value": round(float(r["CpD"]), 6),
                "converged": bool(r["converged"]),
                "note": "zeroed site_ids: " + " ".join(map(str, zeroed)),
            }
        )

    medians = {}
    for k in (10, 25, 50, 100, 200):
        rb = cpd_rate(tree, mat0, n_slices=k)
        medians[k] = float(rb["CpD"].median())
        rows.append(
            {
                "panel": "b_all",
                "series": "cpd_rate_median_all_30",
                "assemblage": "",
                "slice_count": k,
                "value": round(medians[k], 6),
                "converged": True,
                "note": "median over all 30 assemblages of the intact scenario",
            }
        )
        for s in sampled:
            rows.append(
                {
                    "panel": "b",
                    "series": "cpd_rate",
                    "assemblage": s,
                    "slice_count": k,
                    "value": round(float(rb["CpD"].iloc[s]), 6),
                    "converged": bool(rb["converged"].iloc[s]),
                    "note": "sampled with numpy default_rng(1); figure labels are 1-based",
                }
            )

    d = sim_passerines(n_sites=208, n_species=308, n_trees=12, seed=7)
    pc = posterior_rate(d["trees"], d["mat"], rate="cpd", n_slices=50)
    for i in range(20):
        rows.append(
            {
                "panel": "c",
                "series": "posterior_median_cpd",
                "assemblage": i,
                "slice_count": 50,
                "value": round(float(pc["rate_median"].iloc[i]), 6),
                "converged": True,
                "note": "12 pseudo-posterior Yule trees",
            }
        )
    rows.append(
        {
            "panel": "c",
            "series": "band_5_95_low_all_208",
            "assemblage": "",
            "slice_count": 50,
            "value": round(float(pc["rate_median"].quantile(0.05)), 6),
            "converged": True,
            "note": "grey band of Fig. A1c over all 208 assemblages",
        }
    )
    rows.append(
        {
            "panel": "c",
            "series": "band_5_95_high_all_208",
            "assemblage": "",
            "slice_count": 50,
            "value": round(float(pc["rate_median"].quantile(0.95)), 6),
            "converged": True,
            "note": "grey band of Fig. A1c over all 208 assemblages",
        }
    )

    out = HERE / "figA1_panel_data.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"wrote {len(rows)} rows -> {out}")

    # self-check against the caption's published numbers
    pb = pd.DataFrame([r for r in rows if r["panel"] == "b"])
    rises = [
        pb[pb.assemblage == s].sort_values("slice_count")["value"].pct_change().iloc[-1] + 1
        for s in sampled
    ]
    rises = [
        (
            pb[pb.assemblage == s].sort_values("slice_count")["value"].iloc[-1]
            / pb[pb.assemblage == s].sort_values("slice_count")["value"].iloc[0]
            - 1
        )
        * 100
        for s in sampled
    ]
    print(
        f"caption check: panel (b) rise range {min(rises):.1f}-{max(rises):.1f} % "
        f"(caption says 33.0-37.3 %)"
    )
    med_rise = (medians[200] / medians[10] - 1) * 100
    print(
        f"caption check: all-30 median {medians[10]:.2f} -> {medians[200]:.2f} "
        f"({med_rise:.1f} %; caption says 2.38 -> 3.21, 34.6 %, computed from the "
        f"unrounded medians - the rounded values give 34.9 %)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
