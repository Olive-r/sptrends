# sptrends <img src="man/figures/logo.png" align="right" height="139" alt="sptrends logo" />

**Statistical Inference for Spatiotemporal Trends in Gridded Data**

<!-- badges: start -->
[![R](https://img.shields.io/badge/R-%3E%3D4.1-blue?logo=r)](https://www.r-project.org/)
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

```r
library(sptrends)

# One raster file per time step, ordered from numbers in the file names.
r <- read_ordered_stack("path/to/rasters")

# Published True Significant Trends workflow.
result <- workflow_tst(r)

print(result)
summary(result)
plot(result)
```

A bundled annual NDVI dataset is available for reproducible examples:

```r
r <- read_ordered_stack(example_data("vhp_ndvi"))
result <- workflow_tst(r, report = FALSE, verbose = FALSE)
plot(result)
```

## Workflows

- `workflow_tst()` implements True Significant Trends: selective
  prewhitening, Contextual Mann-Kendall, Theil-Sen and adaptive FDR.
- `workflow_rta()` implements Robust Trend Analysis: Contextual
  Mann-Kendall, Theil-Sen and FDR-BH without prewhitening.
- `workflow_trends()` combines user-selected preprocessing, trend, slope and
  FDR modules through the same validated interfaces.

Each analytical stage is also available independently:

```r
pw <- prewhiten(r, method = "TFPW_WS")
trend <- trend_test(pw$series, method = "CMK")
slope <- slope_estimator(pw$series, method = "TS")
fdr <- fdr_correction(trend$stats$p, method = c("BH", "BKY"))
```

`spatial_autocorrelation()` provides independent global and local permutation
inference for environmental variables, residuals, coefficients or inferential
fields. Local p-value rasters can be passed to `fdr_correction()` for BH, BKY
or BY control.

## Main modules

- `compute_anomalies()` and `prewhiten()` — temporal preprocessing.
- `trend_test()` — CMK with configurable odd neighbourhoods (3 by 3 by
  default), MK, OLS and modified MK.
- `slope_estimator()` — Theil-Sen, OLS and Siegel repeated median.
- `fdr_correction()` — BH, adaptive BKY and BY.
- `spatial_autocorrelation()` — permutation-based spatial diagnostics.
- `sim_trend_stack()`, `simulation_design()`, `compare_detections()`, and
  `benchmark_methods()` — known-truth simulation and reproducible method
  benchmarking, including S3 performance plots across scenario factors.
- `benchmark_summary()` — Monte Carlo summaries, including empirical FDR and
  FWER where applicable.
- `read_ordered_stack()` and `read_netcdf_stack()` — raster import with
  temporal-order checks.

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

## Citation

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
  doi:10.3390/rs16203886

## Authors

- Oliver Gutiérrez-Hernández
  ([ORCID](https://orcid.org/0000-0003-2580-5465)), University of Málaga.
- Luis V. García
  ([ORCID](https://orcid.org/0000-0002-5514-2941)), IRNAS-CSIC.

## License

GPL (>= 3)
