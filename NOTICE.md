# Third-party notice

PhyloSlicer is an independent implementation. It contains no source code
derived from any other package, and it does not vendor, redistribute or link
against any third-party code. What it does carry are files that were produced
*with* third-party software, and those are recorded here.

## treesliceR

`validation/reference/*.csv` and `validation/reference/*.tre` are the outputs of
**treesliceR 1.1.0** (R), the reference implementation of the temporal-slicing
paradigm this package re-implements:

> treesliceR — MIT License, © 2023 treesliceR authors (see the package's own
> `LICENSE`: `YEAR: 2023`, `COPYRIGHT HOLDER: treesliceR authors`).
> Authors: Matheus Lima Araujo (aut, cre, cph; ORCID 0000-0002-9111-725X),
> Luiz Gabriel Souza e Souza Ferreira (ORCID 0009-0002-5881-9791), Gabriel
> Nakamura (ORCID 0000-0002-5144-5312), Marco Tulio Pacheco Coelho
> (ORCID 0000-0002-7831-3053), Thiago Fernando Rangel (ORCID
> 0000-0002-2001-7382).
> Package: <https://github.com/AraujoMat/treesliceR> ·
> Paper: <https://doi.org/10.1111/ecog.07364>

Those files are the *expected* side of a numerical cross-check, not a dependency:
PhyloSlicer runs and installs without R, and `validation/compare.py` only reads
the CSVs. `benchmarks/*_r.R` call the installed treesliceR to reproduce the
comparison tables; they are drivers, and they ship no treesliceR code.

Because the reference values were produced by an external program whose floating
point path differs from this one, the agreement tolerances in
[`validation/NOTES.md`](validation/NOTES.md) are part of the artifact: a rate
row that agrees to 1e-4 is a statement about two different optimisers solving
the same objective, not about identical bits.

## Runtime dependencies

numpy, scipy, pandas and matplotlib are required; dendropy, biopython, geopandas
and numba are optional extras. PhyloSlicer imports them at call time, never
copies them, and remains fully usable when the optional ones are absent.

## Documentation licence

The prose in `README*.md`, `docs/USER_GUIDE*.md`, `benchmarks/README*.md`,
`validation/NOTES*.md` and `testing_suite/*.md` is released under CC-BY 4.0; the
code in this repository is BSD-3-Clause (see [`LICENSE`](LICENSE)).
