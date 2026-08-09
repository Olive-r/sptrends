# Summarise a method benchmark across Monte Carlo replicates

Aggregates the replicate-level output of
[`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md)
while preserving scenarios and methods. For detection-like stages,
empirical FDR is the mean false-discovery proportion and empirical FWER
is the proportion of replicates containing at least one false positive.

## Usage

``` r
benchmark_summary(x, path = NULL, verbose = TRUE)
```

## Arguments

- x:

  Replicate-level result returned by
  [`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md).

- path:

  Character or `NULL`. If supplied, the summary is written as a CSV file
  at this path.

- verbose:

  Logical. If `TRUE`, reports progress, elapsed time and the estimated
  time remaining while scenario-method groups are summarised.

## Value

A data frame with one row per scenario and method, retained scenario
factors, the number of replicates, means and standard deviations of
numerical metrics, and `EmpiricalFDR`/`EmpiricalFWER` when detection
metrics are present.

## Details

**Function type:** **Benchmarking function** – summarises known-truth
experiments; it does not perform statistical inference on user data.

## Typical use

Run
[`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md)
and pass its result directly to this function.

## Methodological details

**Monte Carlo aggregation**

Replicates, rather than raster cells, are the independent Monte Carlo
units. Therefore FDR and FWER are aggregated across replicate-level
false discovery proportions and false-positive indicators, respectively.

**Computational considerations**

Aggregation operates on the retained benchmark table and does not rerun
simulations or methods.

**Limitations**

Summary precision depends on the number of independent replicates and
the range of scenarios evaluated. It does not generalise beyond those
designs.

**Quality assurance**

Tests verify grouping, scenario retention, means, standard deviations,
empirical FDR and FWER, CSV output and invalid input handling. The
retained external validation also verifies the summaries against 4,000
paired MK replicates and 33 independent simulation-cycle controls.

## See also

[`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md),
[`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md)

Other validation functions:
[`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md),
[`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md),
[`plot_detection_comparison()`](https://olive-r.github.io/sptrends/reference/plot_detection_comparison.md),
[`simulation_design()`](https://olive-r.github.io/sptrends/reference/simulation_design.md)

## Examples

``` r
x <- structure(
  data.frame(Scenario = rep("null", 2), Replicate = 1:2,
             Method = rep("method", 2), FP = c(1, 0),
             FalseDiscoveryProportion = c(1, 0),
             AnyFalsePositive = c(1, 0)),
  class = c("sptrends_benchmark", "data.frame"))
benchmark_summary(x, verbose = FALSE)
#>   Scenario Method n_replicates FP_mean     FP_sd FalseDiscoveryProportion_mean
#> 1     null method            2     0.5 0.7071068                           0.5
#>   FalseDiscoveryProportion_sd AnyFalsePositive_mean AnyFalsePositive_sd
#> 1                   0.7071068                   0.5           0.7071068
#>   EmpiricalFDR EmpiricalFWER
#> 1          0.5           0.5
```
