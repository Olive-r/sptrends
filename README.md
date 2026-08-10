# sptrends <img src="man/figures/logo.png" align="right" height="139" alt="sptrends logo" />

**Statistical Inference for Spatiotemporal Trends in Gridded Data**

<!-- badges: start -->
[![R](https://img.shields.io/badge/R-%3E%3D4.1-blue?logo=r)](https://www.r-project.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21822842.svg)](https://doi.org/10.5281/zenodo.21822842)
[![pkgdown](https://img.shields.io/badge/docs-pkgdown-blue)](https://olive-r.github.io/sptrends/)
<!-- badges: end -->

sptrends provides a reproducible framework for statistical inference
of spatiotemporal trends in gridded environmental data. It offers modular
methods for serial-correlation treatment, trend testing, spatially explicit
inference, global and local spatial-autocorrelation diagnostics, robust slope
estimation,
multiple-testing correction, simulation, benchmarking, visualisation, mapping,
and reporting.

## Installation

Install the source package from a local archive or directory:

```r
remotes::install_local("path/to/sptrends")
```

## Quick start

A bundled annual NDVI dataset is available for reproducible examples:

```r
library(sptrends)

r <- read_ordered_stack(example_data("vhp_ndvi"))
result <- workflow_tst(r, report = FALSE, verbose = FALSE)
plot(result)
```

With your own data, replace the file path:

```r
# One raster file per time step, ordered from numbers in the file names.
r <- read_ordered_stack("path/to/rasters")

# Published True Significant Trends workflow
# (Gutierrez-Hernandez and Garcia, 2025, https://doi.org/10.1016/j.rsase.2024.101377).
result <- workflow_tst(r)

print(result)
summary(result)
plot(result)
```

## Workflows

- `workflow_tst()` implements True Significant Trends: selective
  prewhitening, Contextual Mann-Kendall, Theil-Sen and adaptive FDR
  ([Gutiérrez-Hernández and García, 2025](https://doi.org/10.1016/j.rsase.2024.101377)).
- `workflow_rta()` implements Robust Trend Analysis: Contextual
  Mann-Kendall, Theil-Sen and FDR-BH without prewhitening
  ([Gutiérrez-Hernández and García, 2024](https://doi.org/10.3390/rs16203886)).
- `workflow_trends()` combines user-selected preprocessing, trend, slope and
  FDR modules through the same validated interfaces.

Each analytical stage is also available independently:

```r
pw <- prewhiten(r, method = "TFPW_WS")
trend <- trend_test(pw$series, method = "CMK")
slope <- slope_estimator(pw$series, method = "TS")
fdr <- fdr_correction(trend$stats$p, method = "BKY")
```

`spatial_autocorrelation()` provides independent global and local permutation
inference for environmental variables, residuals, coefficients or inferential
fields. Local p-value rasters can be passed to `fdr_correction()` for BH, BKY
or BY control.

## Main modules

- `read_ordered_stack()` and `read_netcdf_stack()` — raster import with
  temporal-order checks.
- `compute_anomalies()` — deseasonalising, when the series has a
  known cycle.
- `prewhiten()` — serial-correlation treatment, trend-preserving.
- `trend_test()` — CMK with configurable odd neighbourhoods (3 by 3 by
  default), MK, OLS and modified MK.
- `slope_estimator()` — Theil-Sen, OLS and Siegel repeated median.
- `fdr_correction()` — BH, adaptive BKY and BY.
- `spatial_autocorrelation()` — permutation-based spatial diagnostics.
- `sim_trend_stack()`, `simulation_design()`, `compare_detections()`,
  `benchmark_methods()` and `benchmark_summary()` — known-truth
  simulation and reproducible method benchmarking, including S3
  performance plots and Monte Carlo summaries (empirical FDR and FWER
  where applicable) across scenario factors.

The current calculations are vectorised and selected stages can use parallel
processing, but analytical rasters are materialised in memory. Test
representative dimensions before processing very large datasets.

## Documentation

The function help pages describe assumptions, alternatives, computational
trade-offs, references and method-specific quality assurance. The vignettes
provide concise guides to preprocessing, trend tests, slope estimation, FDR
correction and complete workflows:

```r
browseVignettes("sptrends")
```

## Quality assurance

sptrends maintains 100% test coverage across every R source file, verified
with `covr::package_coverage()`, and passes `R CMD check` with 0 errors, 0
warnings and 0 notes on multiple platforms (local, and win-builder
R-devel/release/oldrelease) -- confirmed on every release. Style, spelling
and link integrity are checked periodically with `lintr`, `goodpractice`,
`spelling` and `urlchecker`.

Beyond the automated `tests/testthat/` suite, the package has been validated
through an external battery covering all 18 exported functions with
correctness checks against known-truth simulations, and a comprehensive
integral audit of code, citations and documentation.

`trend_test(method = "CMK")` has additionally been cross-checked against an
installed copy of [`ConMK`](https://github.com/antiphon/ConMK) (Antiphon,
GitHub, not on CRAN -- also available as a fork at
[`geoporttishare/ConMK`](https://github.com/geoporttishare/ConMK)), the
closest available external reference implementation of the contextual
Mann-Kendall test. The base statistic matched to floating-point precision,
and the optional `continuity = TRUE` argument reproduces `ConMK`'s own
p-values exactly at the specific edge case where the two implementations
would otherwise be expected to diverge. Full details in `?trend_test`'s
"External validation" section.

Every change is documented transparently in [NEWS.md](NEWS.md).

## Citation

To cite the package itself:

> Gutiérrez-Hernández, O., & García, L. V. (2026). sptrends: Statistical
> Inference for Spatiotemporal Trends in Gridded Data [R package]. Zenodo.
> <https://doi.org/10.5281/zenodo.21822842>

Or, from R:

```r
citation("sptrends")
```

Published workflows:

- Gutiérrez-Hernández, O. and García, L.V. (2025). Uncovering true
  significant trends in global greening. *Remote Sensing Applications:
  Society and Environment*, 37, 101377.
  <https://doi.org/10.1016/j.rsase.2024.101377>
- Gutiérrez-Hernández, O. and García, L.V. (2024). Robust Trend Analysis
  in Environmental Remote Sensing: A Case Study of Cork Oak Forest
  Decline. *Remote Sensing*, 16(20), 3886.
  <https://doi.org/10.3390/rs16203886>

## Authors

- Oliver Gutiérrez-Hernández
  ([ORCID](https://orcid.org/0000-0003-2580-5465)), Associate Professor,
  Department of Geography, University of Málaga.
- Luis V. García
  ([ORCID](https://orcid.org/0000-0002-5514-2941)), Tenured Scientist,
  Institute of Natural Resources and Agrobiology of Seville (IRNAS),
  Spanish National Research Council (CSIC).

## License

GPL (>= 3)
