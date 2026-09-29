"""phyloslicer - high-performance temporal slicing of phylogenies.

PhyloSlicer implements the "slice phylogenies and infer phylogenetic
patterns over evolutionary time" paradigm introduced by treesliceR
(Araujo et al. 2025, Ecography) as an independent, vectorised Python
implementation with posterior-tree uncertainty quantification, robust
fitting diagnostics, and geospatial workflow integration.

Released under the BSD-3-Clause licence (see ``LICENSE``); the version of
record is :data:`__version__` and its release date is recorded in
``CHANGELOG.md``.
"""

from . import (
    compat,
    datasets,
    indices,
    io,
    rates,
    simulate,
    slicing,
    spatial,
    uncertainty,
    viz,
)
from .core import SlicedTree, TreeArray
from .rates import cpb_rate, cpd_rate, cpe_rate, fit_rate, origin_time, rate_profile, sensitivity
from .slicing import (
    SliceStack,
    prune_tips,
    slice_interval,
    slice_pieces,
    slice_rootward,
    slice_tipward,
)
from .uncertainty import bootstrap_sites, iter_trees, posterior_origin, posterior_rate

__version__ = "0.1.0"

__all__ = [
    "__version__",
    "TreeArray",
    "SlicedTree",
    "SliceStack",
    "slice_rootward",
    "slice_tipward",
    "slice_interval",
    "slice_pieces",
    "prune_tips",
    "fit_rate",
    "origin_time",
    "cpd_rate",
    "cpe_rate",
    "cpb_rate",
    "rate_profile",
    "sensitivity",
    "posterior_rate",
    "posterior_origin",
    "bootstrap_sites",
    "iter_trees",
    "core",
    "io",
    "slicing",
    "indices",
    "rates",
    "uncertainty",
    "spatial",
    "viz",
    "simulate",
    "compat",
    "datasets",
]
