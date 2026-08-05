# Package index

## Published workflows

sptrends is not a collection of individual statistical functions – it is
a platform for **published workflows**. Each one below is a complete,
citable, peer-reviewed method with its own paper, not an ad hoc
combination this package invented: Robust Trend Analysis (RTA,
Gutiérrez-Hernández & García, 2024) and True Significant Trends (TST,
Gutiérrez-Hernández & García, 2025). See
[`?workflow_tst`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)
and
[`?workflow_rta`](https://olive-r.github.io/sptrends/reference/workflow_rta.md)
for how the two differ and why both are offered, rather than one
superseding the other.
[`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)/[`plot()`](https://rdrr.io/r/graphics/plot.default.html)
work the same way on every classed object this package returns – one
shared entry point per generic (the same convention `terra` itself
uses), not a differently-named function to remember for each result
type.

- [`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)
  : True Significant Trends (TST): the full pipeline in one call
- [`workflow_rta()`](https://olive-r.github.io/sptrends/reference/workflow_rta.md)
  : Robust Trend Analysis (RTA): the full pipeline in one call
- [`print(`*`<sptrends>`*`)`](https://olive-r.github.io/sptrends/reference/print.sptrends.md)
  : Print a sptrends result
- [`summary(`*`<sptrends>`*`)`](https://olive-r.github.io/sptrends/reference/summary.sptrends.md)
  : Summarise a sptrends result
- [`plot(`*`<sptrends>`*`)`](https://olive-r.github.io/sptrends/reference/plot.sptrends.md)
  : Plot a sptrends result

## Configurable trend workflow

[`workflow_trends()`](https://olive-r.github.io/sptrends/reference/workflow_trends.md)
is not itself a published method – it lets you assemble your own
combination of the same prewhitening, trend testing, slope estimation
and multiple-testing correction methods the two published workflows
above are built from, for the case where neither matches what a given
analysis needs. Returns the same kind of object
[`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)/[`plot()`](https://rdrr.io/r/graphics/plot.default.html)
recognise – see “Published workflows” above.

- [`workflow_trends()`](https://olive-r.github.io/sptrends/reference/workflow_trends.md)
  : Configure a monotonic or linear trend-analysis workflow

## Preprocessing

Steps applied to the raw raster time series *before* trend estimation or
significance testing – removing serial autocorrelation or seasonality
that would otherwise distort the steps that come after. Neither is one
of the core trend-analysis pillars itself; both prepare the input for
them. Both also return an object
[`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)/[`plot()`](https://rdrr.io/r/graphics/plot.default.html)
recognise – see “Published workflows” above.

- [`prewhiten()`](https://olive-r.github.io/sptrends/reference/prewhiten.md)
  : AR(1) prewhitening of raster time series
- [`compute_anomalies()`](https://olive-r.github.io/sptrends/reference/compute_anomalies.md)
  : Remove the seasonal cycle from raster time series

## Core methods

The classical, individually citable statistical techniques each
published workflow above is built from – every one of them usable
standalone, not only as a step inside
[`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)/[`workflow_rta()`](https://olive-r.github.io/sptrends/reference/workflow_rta.md).
Each function’s own help page traces it back to its original
publication. All four also return an object
[`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)/
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) recognise – see
“Published workflows” above.

- [`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md)
  : Trend tests for raster time series
- [`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md)
  : Slope estimators for raster time series
- [`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md)
  : Apply false discovery rate (FDR) correction to multiple p-values
- [`spatial_autocorrelation()`](https://olive-r.github.io/sptrends/reference/spatial_autocorrelation.md)
  : Permutation-based spatial autocorrelation tests

## Diagnostics

What remains standalone after the redesign: the interactive single-cell
inspector. Every other diagnostic previously listed here
(trend/Theil-Sen/prewhitening/FDR/Moran’s I summaries, histograms, maps,
category labels, and the FDR-masked direction map) is now reached via
[`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)/[`plot()`](https://rdrr.io/r/graphics/plot.default.html)
on the relevant object – see “Published workflows” and “Core methods”
above.

- [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md)
  : Inspect a single cell's (or area's) raw time series interactively

## Utilities

Everything else this platform needs around the science: reading gridded
time series in, simulating synthetic ground truth, benchmarking
detection methods against it, and the bundled real-world example
dataset.
[`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md)
also returns an object
[`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)/[`plot()`](https://rdrr.io/r/graphics/plot.default.html)
recognise – see “Published workflows” above – while remaining an
ordinary data frame underneath (indexing, `$`, and so on all still
work).

- [`read_ordered_stack()`](https://olive-r.github.io/sptrends/reference/read_ordered_stack.md)
  : Read and chronologically order a folder of raster files
- [`read_netcdf_stack()`](https://olive-r.github.io/sptrends/reference/read_netcdf_stack.md)
  : Read and chronologically order a single multi-temporal NetCDF file
- [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)
  : Generate a synthetic gridded time series with known true trends
- [`simulation_design()`](https://olive-r.github.io/sptrends/reference/simulation_design.md)
  : Build a factorial design of simulation scenarios
- [`example_data()`](https://olive-r.github.io/sptrends/reference/example_data.md)
  : Path to sptrends' bundled example dataset
- [`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md)
  : Compare detection methods against a known ground truth
- [`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md)
  : Benchmark statistical methods across known-truth simulation
  scenarios
- [`benchmark_summary()`](https://olive-r.github.io/sptrends/reference/benchmark_summary.md)
  : Summarise a method benchmark across Monte Carlo replicates
