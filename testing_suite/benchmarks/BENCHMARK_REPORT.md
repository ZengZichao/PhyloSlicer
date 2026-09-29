# PhyloSlicer Benchmark Report

**Version:** 0.1.0  
**Date:** 2026-09-26  
**Python:** 3.12.14  

## Methodology

Each configuration is run 3 times; the table reports the median.
Memory is peak traced memory (tracemalloc).

## Results

| n_tips | n_sites | n_slices | slice (s) | pd_per_slice (s) | cpd_rate (s) | cpd_MB | posterior (s) | n_converged |
|--------|---------|----------|-----------|------------------|--------------|--------|---------------|-------------|
| 20 | 10 | 50 | 0.0006 | 0.0042 | 0.0620 | 0.1 | 0.1648 | 10 |
| 50 | 20 | 50 | 0.0018 | 0.0065 | 0.0350 | 0.2 | 0.1803 | 20 |
| 100 | 30 | 50 | 0.0019 | 0.0081 | 0.0458 | 0.4 | 0.1651 | 30 |
| 100 | 30 | 100 | 0.0024 | 0.0075 | 0.0492 | 0.6 | 0.2092 | 30 |
| 100 | 30 | 200 | 0.0018 | 0.0091 | 0.0569 | 0.9 | 0.2080 | 30 |
| 200 | 50 | 100 | 0.0042 | 0.0127 | 0.0795 | 1.0 | 0.2710 | 50 |
| 500 | 100 | 100 | 0.0170 | 0.0551 | 0.0836 | 3.0 | 0.4725 | 100 |

## Scaling

The contribution matrix is `E x k` (edges x slices), so:
- Slicing cost is O(E * k) — linear in both tree size and slice count.
- PD per slice is O(S * E * k) via sparse matrix product.
- Rate fitting is vectorised: all sites solved simultaneously.
