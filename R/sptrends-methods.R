# Adding a new "sptrends" subclass later (a future published workflow or
# core method): give its constructor `class(out) <- c("newclass",
# "sptrends")`, then write `.print_newclass()`, `.summary_newclass()`,
# and `.plot_newclass()` (see any existing *-methods.R file, or
# slope-methods.R, for the pattern) -- .sptrends_dispatch(), in
# utils-internal.R, finds them by that naming convention alone. The
# three generics below never need editing again; neither does the test
# file that exercises them (test-sptrends-methods.R), beyond adding the
# new class to its own object list. Only this file's own documentation
# (the @section blocks below, one per class) needs a new section by
# hand, since no amount of convention can auto-generate accurate prose.

#' Print a sptrends result
#'
#' A quick, one-line-per-detail overview of any classed object this
#' package returns -- `workflow_tst()`, `workflow_rta()`,
#' [workflow_trends()],
#' [trend_test()], [slope_estimator()],
#' [prewhiten()], and [fdr_correction()] all return an object
#' with `"sptrends"` as (one of) its classes, and `print()`,
#' [summary.sptrends()], and [plot.sptrends()] all work the same way
#' regardless of which one you have -- the specific one-line summary
#' shown depends on `x`'s own class, listed below. This mirrors the
#' convention used throughout `terra` itself (a single `print()`,
#' whether `x` is a `SpatRaster` or a `SpatVector`): one predictable
#' entry point per generic, not a different function name to remember
#' for each result type. The API is organised around the *object* a
#' function returns, not around remembering which reporting function
#' goes with which: `x <- workflow_tst(...)`, then `print(x)`,
#' `summary(x)`, `plot(x)`, regardless of what `x` actually is.
#'
#' | Generic | Purpose |
#' | --- | --- |
#' | `print()` | Quick overview |
#' | `summary()` | Detailed textual report |
#' | `plot()` | Visual exploration |
#'
#' `print()` itself is intended for a quick inspection of an object at
#' the console; use [summary.sptrends()] for more detailed textual
#' reporting and [plot.sptrends()] for graphical exploration.
#'
#' **Function type:** **Reporting/derived function** -- presents an existing
#' result and does not compute a new statistical estimate.
#'
#' @section Typical use:
#' `result <- workflow_tst(x); print(result)` for a concise console overview.
#'
#' @param x An object of class `"sptrends"` (also one of `"tst"`,
#'   `"rta"`, `"workflow_trends"`, `"trend_test"`, `"slope"`,
#'   `"prewhiten"`, `"fdr"`,
#'   `"spatial_autocorrelation"`, `"compare_detections"`,
#'   `"sptrends_simulation"`, `"sptrends_simulation_design"`, or
#'   `"sptrends_benchmark"`),
#'   from [workflow_tst()], [workflow_rta()], [workflow_trends()],
#'   [trend_test()],
#'   [slope_estimator()], [prewhiten()], [fdr_correction()],
#'   [spatial_autocorrelation()], [compare_detections()],
#'   [sim_trend_stack()], [simulation_design()], or [benchmark_methods()].
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @section Methodological details:
#' **Published workflow: `"tst"`.**
#' Whether prewhitening ran, how many cells were tested, and how many
#' are significant after FDR correction (if run).
#' **Published workflow: `"rta"`.**
#' The Theil-Sen slope range, the trend test's cell count, and how many
#' cells are significant after FDR-BH correction.
#' **Configurable workflow: `"workflow_trends"`.**
#' The selected preprocessing, trend-test, slope and FDR stages,
#' including skipped optional stages and the qualified Moran assessment
#' when requested.
#' **Trend estimation: `"trend_test"`.**
#' How many cells were tested, and how many are significant at the
#' conventional alpha=0.05 threshold, uncorrected.
#' **Trend estimation: `"slope"`.**
#' How many cells have a valid slope, and its range.
#' **Diagnostic: `"prewhiten"`.**
#' How many cells were prewhitened, out of how many valid cells.
#' **Diagnostic: `"fdr"`.**
#' How many cells are significant under each method that was requested
#' (raw, BH, BKY, and BY, if it was explicitly requested -- see
#' `?fdr_correction`'s own `method` argument for why `"BY"` is opt-in,
#' not part of its own default).
#' **Diagnostic: `"spatial_autocorrelation"`.**
#' Global results show the observed statistic (Moran's I or Getis-Ord
#' General G), its sign where applicable, and its permutation p-value.
#' Local results show the valid-cell count and the exploratory number of
#' cells below the raw `alpha`, followed by the route to
#' [fdr_correction()] for BH, BKY or BY.
#' **Validation: `"compare_detections"`.**
#' The comparison table itself, printed as a plain data frame (this
#' class also inherits from `"data.frame"`, so indexing, `$`, and so on
#' all work exactly as they would on any other one).
#' **Simulation and benchmarking.**
#' Simulation objects report their dimensions, dependence model and known
#' signal. Designs report scenario counts and varied factors. Benchmarks
#' report their stage, methods, scenarios, replicates and elapsed time.
#'
#' @seealso [summary.sptrends()] for detailed textual output and
#'   [plot.sptrends()] for graphical exploration.
#' @examples
#' \donttest{
#' # Annual mean NDVI from the bundled environmental dataset.
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' result <- workflow_tst(r, report = FALSE, verbose = FALSE)
#' print(result)  # dispatches to the "tst" case above
#' }
#'
#' @export
print.sptrends <- function(x, ...) {
  invisible(.sptrends_dispatch("print", x, ...))
}

#' Summarise a sptrends result
#'
#' The detailed textual report behind `print()`'s own one-line
#' overview -- see [print.sptrends()] for the full class list and the
#' rationale for one shared entry point per generic. Each class here
#' calls its own underlying reporting function directly, listed below,
#' rather than duplicating what `print()` already shows.
#'
#' **Function type:** **Reporting/derived function** -- summarises an existing
#' result and does not recompute its statistical analysis.
#'
#' @section Typical use:
#' `result <- workflow_tst(x); summary(result)` for a detailed textual report.
#'
#' @param object An object of class `"sptrends"` (also one of `"tst"`,
#'   `"rta"`, `"workflow_trends"`, `"trend_test"`, `"slope"`,
#'   `"prewhiten"`, `"fdr"`,
#'   `"spatial_autocorrelation"`, `"compare_detections"`,
#'   `"sptrends_simulation"`, `"sptrends_simulation_design"`, or
#'   `"sptrends_benchmark"`),
#'   from [workflow_tst()], [workflow_rta()], [workflow_trends()],
#'   [trend_test()],
#'   [slope_estimator()], [prewhiten()], [fdr_correction()],
#'   [spatial_autocorrelation()], [compare_detections()],
#'   [sim_trend_stack()], [simulation_design()], or [benchmark_methods()].
#' @param ... Passed on to the underlying reporting function (e.g.
#'   `path`; `alpha` for `"tst"`/`"rta"`/`"trend_test"` objects).
#' @return Invisibly, whatever the underlying reporting function itself
#'   returns -- see each section below.
#' @section Methodological details:
#' **Published workflow: `"tst"`.**
#' The full detail behind the `"tst"` case of [print.sptrends()]: the
#' uncorrected trend summary table and, if FDR correction was run, the
#' FDR summary. Returns a list with `trend` and `fdr` (or `NULL`).
#' **Published workflow: `"rta"`.**
#' The full detail behind the `"rta"` case of [print.sptrends()]: the
#' uncorrected trend summary table, the Theil-Sen slope summary, and the
#' FDR-BH summary. Returns a list with `trend` and `fdr`.
#' **Configurable workflow: `"workflow_trends"`.**
#' The uncorrected trend table, optional slope summary and selected FDR
#' summary.
#' **Trend estimation: `"trend_test"`.**
#' Cell counts and increase/decrease/no-change breakdown at multiple
#' alpha levels; calls `trend_summary()` internally.
#' **Trend estimation: `"slope"`.**
#' Valid cells, range, median, mean, and the increasing/decreasing/flat
#' breakdown; calls `slope_summary()` internally.
#' **Diagnostic: `"prewhiten"`.**
#' Valid cells, cells prewhitened, mean rho among them, and median
#' Durbin-Watson; calls `prewhiten_summary()` internally.
#' **Diagnostic: `"fdr"`.**
#' Significant/not-significant counts and percentages for every method
#' requested; calls `fdr_summary()` internally.
#' **Diagnostic: `"spatial_autocorrelation"`.**
#' Global results report the observed statistic, permutation distribution
#' and empirical summary statistics. Local results report the statistic
#' range, minimum permutation p-value and exploratory raw-significance
#' count.
#' **Validation: `"compare_detections"`.**
#' Which method scores best on each numeric metric in the table -- a
#' small table of its own, `metric`/`best_method`, not part of what
#' `compare_detections()` itself computes.
#' **Simulation and benchmarking.**
#' Simulation summaries quantify the true signal, true-null proportion and
#' slope range. Design summaries count levels per factor. Benchmark summaries
#' retain scenario factors and aggregate performance over independent Monte
#' Carlo replicates, including empirical FDR and FWER where available.
#'
#' @seealso [print.sptrends()] for a concise overview and
#'   [plot.sptrends()] for graphical exploration.
#' @examples
#' \donttest{
#' # Annual mean NDVI from the bundled environmental dataset.
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' result <- workflow_tst(r, report = FALSE, verbose = FALSE)
#' summary(result)  # dispatches to the "tst" case above
#' }
#'
#' @export
summary.sptrends <- function(object, ...) {
  .sptrends_dispatch("summary", object, ...)
}

#' Plot a sptrends result
#'
#' The graphical exploration counterpart to `print()`'s one-line
#' overview and `summary()`'s textual report -- see [print.sptrends()]
#' for the full class list and the rationale for one shared entry
#' point per generic. What `which` (and any other named
#' argument) accepts depends entirely on `x`'s own class -- see the
#' sections below, grouped by what kind of result each class
#' represents.
#'
#' **Function type:** **Reporting/derived function** -- visualises an existing
#' result and does not alter or recompute its statistical analysis.
#'
#' @section Typical use:
#' `result <- workflow_tst(x); plot(result)` draws the class-specific default;
#' use `which` for an alternative diagnostic view where supported.
#'
#' @param x An object of class `"sptrends"` (also one of `"tst"`,
#'   `"rta"`, `"workflow_trends"`, `"trend_test"`, `"slope"`,
#'   `"prewhiten"`, `"fdr"`,
#'   `"spatial_autocorrelation"`, `"compare_detections"`,
#'   `"sptrends_simulation"`, `"sptrends_simulation_design"`, or
#'   `"sptrends_benchmark"`),
#'   from [workflow_tst()], [workflow_rta()], [workflow_trends()],
#'   [trend_test()],
#'   [slope_estimator()], [prewhiten()], [fdr_correction()],
#'   [spatial_autocorrelation()], [compare_detections()],
#'   [sim_trend_stack()], [simulation_design()], or [benchmark_methods()].
#' @param ... Passed on to the underlying plotting logic -- see each
#'   section below for the arguments (typically `which`, and sometimes
#'   `method`, `smooth`, `alpha`, or `path`) `x`'s own class accepts.
#' @return `x`, invisibly.
#' @section Methodological details:
#' **Published workflow: `"tst"`.**
#' By default, draws a binarised trend map: increases and decreases from
#' the selected trend statistic are retained only where the selected
#' multiple-testing procedure rejects the null hypothesis (see
#' [direction_map()]).
#'
#' `which` -- the default and the three main views:
#'
#' | `which` | Draws |
#' | --- | --- |
#' | `"direction"` (default) | Binarised trend map after FDR correction (titled "TST direction map"; "RTA direction map" for the `"rta"` case below) |
#' | `"significance"` | FDR-BH/FDR-BKY significance maps side by side |
#' | `"trend"` | Uncorrected trend statistic/p-value/significance/direction |
#' | `"slope"` | Theil-Sen slope, masked to significant cells |
#'
#' `which = "trend"` draws **uncorrected** diagnostics only -- see the
#' "Warning" section of `?trend_test` before reporting significance
#' from this view. `which = "slope"` is the more reliable source for
#' per-cell direction/magnitude than the CMK statistic `Sm`, since
#' `Sm`'s neighbourhood averaging can occasionally disagree in sign
#' with a cell's own Theil-Sen estimate in a spatially heterogeneous
#' neighbourhood (see "Display smoothing" below for how this view
#' handles that visually).
#'
#' Eight further views, all **uncorrected** (no FDR, no significance
#' masking) -- diagnostics on the raw slope or the raw p-value, not a
#' final result to report:
#'
#' | `which` | Draws |
#' | --- | --- |
#' | `"slope_map"` | Continuous Theil-Sen slope, unmasked |
#' | `"slope_direction"` | Binary sign of the slope, unmasked |
#' | `"slope_hist"` | Histogram of slope values with a density curve |
#' | `"slope_bar"` | Bar chart of positive/negative/zero slope counts |
#' | `"pvalue_map"` | Continuous, uncorrected p-value |
#' | `"pvalue_significance"` | Binary significant/not at `alpha` |
#' | `"pvalue_hist"` | Classic binned histogram of p-values |
#' | `"pvalue_bar"` | Bar chart of significant vs. not |
#'
#' `slope_map`/`slope_direction`/`slope_hist`/`slope_bar` need
#' `x$theil_sen` (i.e. `workflow_tst()` must have been run with
#' `theil_sen = TRUE`); the four `pvalue_*` views always work, since
#' `x$trend` is never `NULL`. `slope_bar`/`pvalue_bar` accept
#' `probability = TRUE` for percentages instead of counts;
#' `pvalue_significance` uses `alpha` (default `0.05`), uncorrected
#' for multiple testing.
#'
#' `method`: `"BKY"` (default, matching `workflow_tst()`'s default
#' `fdr_method`), `"BH"` or `"BY"` -- which correction to use for
#' `which = "direction"` or `which = "slope"`. Ignored otherwise. If
#' `x` only has the other one (e.g. `workflow_tst()` was called with
#' `fdr_method = "BH"`), set `method` to match.
#' **Published workflow: `"rta"`.**
#' By default, draws the map you actually want to report: direction of
#' change masked by FDR-BH significance (see [direction_map()]).
#' Unlike the `"tst"` case above, there is no `method` argument --
#' `workflow_rta()` always uses FDR-BH (see `?workflow_rta`'s "How RTA
#' differs from TST,
#' and why"), so there is nothing to choose between.
#'
#' `which`: `"direction"` (default), `"significance"`, `"trend"`,
#' `"slope"`, and the eight further `slope_*`/`pvalue_*` diagnostic
#' views -- same meaning as the `"tst"` case above, but always
#' FDR-BH.
#' **Configurable workflow: `"workflow_trends"`.**
#' Provides the same trend, slope, p-value, significance and direction
#' views, using whichever optional slope and FDR stages were selected.
#' **Display smoothing (`"tst"`/`"rta"`/`"workflow_trends"`,
#' `which = "slope"` only).**
#' `smooth`: logical. If `TRUE`
#' (default) and `x$theil_sen` (`x$slope` for a `"workflow_trends"`
#' object -- same mechanism, different field name) was **not** already
#' smoothed at source
#' (via `workflow_tst()`'s own `theil_sen_args = list(smooth_neighbourhood =
#' TRUE)`, which is *not* `workflow_tst()`'s own default -- see
#' `?workflow_tst`'s
#' "Computational considerations" section for why `workflow_tst()`
#' and [workflow_rta()] deliberately do not differ here), the
#' significant-cells-only slope map is smoothed with a queen-3x3 median
#' filter **for this plot only**; if it was already smoothed at source,
#' this is not applied a second time. `x$theil_sen` itself is never
#' modified by plotting, and neither is any part of the significance
#' decision; this only changes what gets drawn. Because smoothing here
#' runs *after* masking to significant cells (a necessarily sparser set
#' of cells than the full raster), the result can look visually
#' blockier than smoothing a dense, unmasked raster would -- this is an
#' expected consequence of masking before smoothing, not a bug (see
#' `.mask_and_smooth_slope()`'s own `na.policy = "omit"` for the
#' related, and separate, fix ensuring a non-significant cell is never
#' itself painted with a colour). Setting `smooth = FALSE` when `workflow_tst()`
#' already smoothed at source cannot recover the unsmoothed values
#' (they were never kept) -- a message explains this if it happens. The
#' plot title always states which stage (if any) applied smoothing, so
#' the display is never ambiguous about what it shows. See
#' `?slope_estimator`'s "Optional queen-neighbourhood smoothing"
#' section for why this is a display convenience, not a validated
#' estimator, and is not applied to the value returned by [workflow_tst()].
#' **Trend estimation: `"trend_test"`.**
#' `which`: `"maps"` (default), all four uncorrected diagnostic maps
#' (trend statistic, p-value, significance, and direction), via
#' `trend_maps()`; `"histograms"`, histograms of the trend statistic and
#' p-value, via `trend_histograms()`. These are the **uncorrected**
#' result; see the "Warning" section of `?trend_test`
#' before reporting significance from them without a multiple-testing
#' correction (see [fdr_correction()] and the `"fdr"` case below).
#'
#' `alpha`: significance threshold, only used when `which = "maps"`.
#' **Trend estimation: `"slope"`.**
#' Draws the zero-centred diverging slope map; no `which` argument (only
#' one view exists for this class).
#' **Diagnostic: `"prewhiten"`.**
#' `which`: `"maps"` (default), the four diagnostic maps; `"histograms"`,
#' the two diagnostic histograms.
#' **Diagnostic: `"fdr"`.**
#' `which`: `"significance"` (default), significance maps (raster input
#' only); `"pvalue_histogram"`, histogram of the raw input p-values;
#' `"comparison"`, bar chart comparing significant counts across
#' whichever of raw/BH/BKY/BY were actually requested;
#' `"threshold"`, the BH vs. BKY step-up threshold plot (only if BKY was
#' requested).
#' **Diagnostic: `"spatial_autocorrelation"`.**
#' Global results draw the null distribution with the observed statistic
#' marked. Local results draw the statistic, empirical z, raw permutation
#' p-value and exploratory raw-significance map. Apply and plot
#' [fdr_correction()] separately for BH, BKY or BY inference.
#' **Validation: `"compare_detections"`.**
#' Intended mainly for simulation studies -- not typically the plot an
#' analyst runs on a real dataset, since it needs a known ground truth
#' to have been scored against in the first place (see
#' [compare_detections()]).
#'
#' A grouped bar chart of the comparison table -- one group of bars per
#' method, one bar per metric; no `which` argument (only one view
#' exists for this class). `metrics` chooses which columns to plot
#' (defaults to all of them).
#' **Simulation and benchmarking.**
#' Simulation plots expose the known signal, slope, direction, breaks, or
#' complete raster series. Design plots show the number of levels per factor.
#' Benchmark plots show performance against a varying scenario factor using
#' lines with uncertainty, grouped bars, replicate boxplots, two-factor
#' heatmaps, or multi-metric profiles. Use `metric`, `scenario`, `group`,
#' `facet`, `type`, `interval`, and `level` to configure these views.
#' Confidence intervals use the between-replicate standard error and a
#' Student-t critical value; intervals for rates and powers are clipped to
#' their admissible range from zero to one.
#'
#' @seealso [print.sptrends()] for a concise overview and
#'   [summary.sptrends()] for detailed textual output.
#' @examples
#' \donttest{
#' # Annual mean NDVI from the bundled environmental dataset.
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' result <- workflow_tst(r, report = FALSE, verbose = FALSE)
#'
#' # Default map: direction of change (greening/browning/no change),
#' # masked by FDR-BKY significance -- cells with grey are not
#' # significant, so they are left out of the coloured pattern.
#' plot(result)
#'
#' # The rate of change (not just direction) among significant cells,
#' # queen-3x3 median smoothed for display (the default for this view).
#' plot(result, which = "slope")
#' }
#'
#' @export
plot.sptrends <- function(x, ...) {
  .sptrends_dispatch("plot", x, ...)
}
