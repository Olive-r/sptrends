# Build a factorial design of simulation scenarios

Creates named argument lists for
[`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)
by crossing temporal, spatial, signal and noise conditions. It separates
experimental design from data generation so the complete scenario grid
can be inspected and retained.

## Usage

``` r
simulation_design(..., constants = list(), prefix = "scenario", verbose = TRUE)
```

## Arguments

- ...:

  Named vectors or lists of factor levels to cross.

- constants:

  Named list of arguments shared by every scenario.

- prefix:

  Character prefix used for generated scenario names.

- verbose:

  Logical. If `TRUE`, reports progress, elapsed time and the estimated
  time remaining while scenarios are assembled.

## Value

A named list of argument lists suitable for the `scenarios` argument of
[`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md),
with classes `"sptrends_simulation_design"` and `"sptrends"` for unified
printing, summaries, and plotting.

## Details

**Function type:** **Benchmarking function** – defines simulation
scenarios; it does not generate data or perform inference.

## Typical use

Define the factors that should vary, add shared settings through
`constants`, and pass the returned list to
[`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md).

## Methodological details

**Experimental design**

Every combination is retained. This makes comparisons across spatial and
temporal dependence explicit and prevents methods from being evaluated
on accidentally different scenario sets. Values that must remain
grouped, such as a two-number `signal_size`, should be wrapped in a
list.

**Computational considerations**

This function only constructs argument lists. Memory use grows with the
product of the numbers of supplied factor levels; data generation
remains deferred to
[`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md).

**Limitations**

A full factorial design can become unnecessarily large. Users should
vary scientifically relevant factors and keep fixed settings in
`constants`.

**Quality assurance**

Tests verify factorial completeness, deterministic ordering, grouped
values, validation failures, metadata retention and S3 presentation. The
complete external simulation-cycle validation passed all 33 prespecified
controls; see `inst/validation/` for the retained protocol and results.

## See also

[`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md),
[`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md)

Other validation functions:
[`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md),
[`benchmark_summary()`](https://olive-r.github.io/sptrends/reference/benchmark_summary.md),
[`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md),
[`plot_detection_comparison()`](https://olive-r.github.io/sptrends/reference/plot_detection_comparison.md)

## Examples

``` r
design <- simulation_design(
  spatial_model = c("independent", "exponential"),
  spatial_rho = c(0.3, 0.7), ar1 = c(0, 0.5),
  trend_strength = c(0, 0.05),
  constants = list(nrow = 20, ncol = 20, n_time = 20,
                   constant_block = FALSE), verbose = FALSE)
length(design)
#> [1] 16
```
