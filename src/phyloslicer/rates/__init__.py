"""Rates layer."""

from .api import cpb_rate, cpd_rate, cpe_rate, rate_profile, sensitivity
from .fitting import (
    DEFAULT_RATE_BOUNDS,
    RateFit,
    adaptive_bounds,
    fit_rate,
    fit_rates_batch,
    origin_time,
)

__all__ = [
    "RateFit",
    "fit_rate",
    "fit_rates_batch",
    "origin_time",
    "DEFAULT_RATE_BOUNDS",
    "adaptive_bounds",
    "cpd_rate",
    "cpe_rate",
    "cpb_rate",
    "rate_profile",
    "sensitivity",
]
