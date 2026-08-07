# Benchmark statistical methods across known-truth simulation scenarios

Coordinates reproducible Monte Carlo experiments without tying the
benchmark to sptrends implementations. Each method is an ordinary
function, so methods from other packages can be evaluated under exactly
the same simulated realisations and truth fields.

## Usage

``` r
benchmark_methods(
  scenarios,
  methods,
  n_replicates = 100L,
  stage = c("trend_test", "prewhitening", "slope", "fdr", "fwer", "custom", "detection"),
  simulator = sim_trend_stack,
  prepare = NULL,
  evaluator = NULL,
  seed = 1L,
  metrics = c("type_i", "type_ii", "type_iii", "field_power", "global_power",
    "within_image_power", "directional_power", "fdr"),
  evaluation_mask = NULL,
  verbose = TRUE
)
```

## Arguments

- scenarios:

  Named list of argument lists passed to `simulator`.

- methods:

  Named list of functions. Each receives `(input, simulation)` and
  returns a transformed raster series for `stage = "prewhitening"`, a
  detection object for trend/FDR/FWER stages, a slope raster/vector for
  `stage = "slope"`, or an object accepted by `evaluator`.

- n_replicates:

  Positive integer number of realisations per scenario.

- stage:

  One of `"prewhitening"`, `"trend_test"`, `"slope"`, `"fdr"`, `"fwer"`,
  or `"custom"`. `"detection"` is retained as an alias for
  `"trend_test"`.

- simulator:

  Function used to generate one known-truth realisation. It must accept
  the scenario arguments plus `seed`, and return at least `series`,
  `true_signal`, `true_direction`, and `true_slope` components when the
  corresponding built-in scorer is used.

- prepare:

  Optional function called once per replicate as
  `prepare(series, simulation)`. Its result becomes the identical
  `input` passed to every method. Use it to compute one common p-value
  or statistic map before comparing FDR or FWER methods.

- evaluator:

  Optional scoring function called as `evaluator(outputs, simulation)`.
  Required for `stage = "custom"`.

- seed:

  Integer seed that deterministically generates replicate seeds.

- metrics:

  Metrics passed to
  [`compare_detections()`](https://olivergh.github.io/sptrends/reference/compare_detections.md)
  for detection benchmarks.

- evaluation_mask:

  Optional common mask or a function of `simulation` returning one. The
  same mask is applied to every method so all methods are compared on
  exactly the same cells.

- verbose:

  Logical. If `TRUE`, report method-level progress, elapsed duration and
  estimated time remaining across the complete Monte Carlo experiment.
  Simulator-level verbosity is disabled automatically when the simulator
  exposes a `verbose` argument and the scenario does not override it.

## Value

A data frame with one row per scenario, replicate and method, including
explicit scenario factors, known-truth composition, elapsed time and
stage-appropriate accuracy metrics. It has classes
`"sptrends_benchmark"` and `"sptrends"`, providing unified
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
scenario-performance
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

## Details

**Function type:** **Benchmarking function** – coordinates simulation,
method execution, scoring and timing; it is not an inferential workflow.

## Typical use

Define named scenarios, define package or external method wrappers, and
run the same methods on every generated realisation. Built-in scoring is
available for prewhitening, trend tests, slopes, FDR and FWER. A custom
evaluator supports arbitrary result structures without changing the
Monte Carlo engine.

## Methodological details

**Paired method comparison**

A replicate seed is generated once and shared by every method within
that replicate. Consequently, method differences are paired within
identical data, known truth, prepared input and evaluation domain rather
than confounded with different random fields or cell subsets. The
returned table retains replicate-level results; uncertainty summaries
must be calculated across those independent replicates, not across
raster cells.

**Statistical assumptions**

Performance estimates are conditional on the simulated scenarios and the
common evaluation domain. Monte Carlo replicates, not raster cells, are
the independent units used to estimate repeated-sampling behaviour.

**Computational considerations**

Runtime grows with scenarios, replicates and methods. Scenario results
are retained at replicate level so expensive experiments can be
summarised or plotted without rerunning the methods.

**Limitations**

Built-in scorers require the documented truth components.
Method-specific objects or cluster-level targets require a custom
`evaluator`; incompatible method-specific evaluation domains are
deliberately rejected.

**Quality assurance**

Tests cover every supported stage, paired inputs, external simulators,
common masks, reproducible seeds, timing, failures, summaries and plots.
In an independent full run, 500 replicates were evaluated in each of
eight spatiotemporal scenarios. Cell-level MK decisions and directions
agreed exactly between sptrends and
[`Kendall::MannKendall()`](https://rdrr.io/pkg/Kendall/man/MannKendall.html)
for every retained performance metric, while paired seeds, summaries and
graphics passed all recorded controls. This validates the benchmark
orchestration and scoring; it does not imply that every method compared
by future users is equivalent.

## See also

[`simulation_design()`](https://olivergh.github.io/sptrends/reference/simulation_design.md),
[`sim_trend_stack()`](https://olivergh.github.io/sptrends/reference/sim_trend_stack.md),
[`compare_detections()`](https://olivergh.github.io/sptrends/reference/compare_detections.md),
[`benchmark_summary()`](https://olivergh.github.io/sptrends/reference/benchmark_summary.md)

Other validation functions:
[`benchmark_summary()`](https://olivergh.github.io/sptrends/reference/benchmark_summary.md),
[`compare_detections()`](https://olivergh.github.io/sptrends/reference/compare_detections.md),
[`plot_detection_comparison()`](https://olivergh.github.io/sptrends/reference/plot_detection_comparison.md),
[`simulation_design()`](https://olivergh.github.io/sptrends/reference/simulation_design.md)

## Examples

``` r
methods <- list(
  MK = function(series, simulation) {
    fit <- trend_test(series, method = "MK", report = FALSE,
                      verbose = FALSE)
    list(significant = fit$stats$p <= 0.05,
         direction = fit$stats$S)
  }
)
scenarios <- list(null = list(nrow = 6, ncol = 6, n_time = 8,
                              trend_fraction = 0))
result <- benchmark_methods(scenarios, methods, n_replicates = 2,
                            seed = 1, verbose = FALSE)
#> >> [compare_detections()] elapsed: 0.01 s
#> >> [compare_detections()] elapsed: 0.01 s
result
#> <sptrends benchmark>
#> Stage: trend_test | scenarios: 1 | methods: 1 | replicates: 2
#> Rows: 2 | elapsed method time: 0.083 s
```
