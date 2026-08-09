#' True Significant Trends (TST): the full pipeline in one call
#'
#' Implements the complete workflow for robust statistical inference of
#' monotonic trends in gridded raster time series. The main entry point
#' of sptrends: chains [prewhiten()], [trend_test()], [slope_estimator()],
#' and [fdr_correction()], in that order, into the True Significant
#' Trends (TST) workflow. Each step is optional and all underlying
#' parameters remain available directly on the individual functions --
#' `workflow_tst()` does not replace them, it saves wiring the calls
#' together for the common case, and returns a single `"tst"` object
#' with its own `print()`, `summary()`, and `plot()` methods.
#'
#' Unlike many trend-analysis workflows, TST separates preprocessing,
#' hypothesis testing, effect-size estimation, and multiple-testing
#' correction into independent but composable steps.
#' `workflow_tst()` simply orchestrates those steps into a reproducible
#' pipeline without hiding any of their parameters.
#'
#' Moran's I ([spatial_autocorrelation()]) is **not** part of this pipeline:
#' it is an independent, general spatial diagnostic. Its optional use on
#' inferential fields can reveal dependence relevant to FDR interpretation,
#' but cannot verify all assumptions of an FDR procedure -- see the package
#' vignette.
#'
#' **Function type:** **Core function** -- the complete published TST workflow.
#'
#' @section Typical use:
#' ```
#' raster time series
#'     |
#' workflow_tst()
#'     |
#' selective prewhitening -> CMK -> Theil-Sen -> adaptive FDR
#'     |
#' one `tst` result containing every stage
#' ```
#' For seasonal input, first use [compute_anomalies()] and pass its
#' `anomalies` raster. For long series, consider setting `max_pairs`
#' through `theil_sen_args`; see "Computational considerations" below.
#'
#' @section Methodological details:
#' **How it works.**
#' ```
#' Input raster
#'     |
#' (optional) prewhiten          -- removes serial autocorrelation
#'     |
#' trend_test                    -- is there a monotonic trend? (CMK)
#'     |
#' (optional) slope_estimator    -- how fast? (Theil-Sen or OLS)
#'     |
#' (optional) FDR correction     -- which cells survive multiple testing?
#'     |
#' "tst" object
#' ```
#' The four steps run in this fixed order because each one's own
#' assumptions depend on what came before it: prewhitening needs to
#' happen before the trend test, since the test assumes independent
#' observations; the slope is estimated on the same (optionally
#' prewhitened) series the test itself used, so the two describe the
#' same data; and FDR correction needs the test's own p-values to exist
#' first. Prewhitening, slope estimation, and FDR correction are each
#' individually optional -- not every analysis needs all three (a
#' quick exploratory look at significance alone might skip slope and
#' FDR entirely; data already known to have negligible serial
#' correlation might skip prewhitening) -- but the trend test itself is
#' not optional, since every other step either feeds into it or
#' consumes its output. The returned object keeps every intermediate
#' result it computed (see "Value" below), not only the final
#' one, so that any step's own output can be inspected or reused
#' without recomputing the whole pipeline.
#'
#' **Statistical assumptions: `alpha` and `q`.**
#' This is one of the most common mistakes in gridded trend analysis, so
#' it is worth spelling out: `alpha` (e.g. `0.05`) is defined for a
#' *single* hypothesis test. A raster is not one test -- it is one test
#' per cell, potentially thousands of them run simultaneously. Applying
#' `alpha` cell by cell, as if each cell were the only test being run, is
#' exactly the multiple-testing error this package's FDR step exists to
#' fix (see the "Warning" section of `?trend_test`); it is
#' not a harmless simplification.
#'
#' It is common practice to set `q` equal to `alpha` (e.g. both `0.05`) for
#' convenience and comparability between the uncorrected and corrected
#' results reported here -- this function does not enforce that, `alpha`
#' and `q` are independent arguments. But equal values does not mean equal
#' meaning: `alpha` bounds the error rate of each individual cell's test,
#' while `q` bounds the expected proportion of false positives among all
#' the cells called significant after correction (see [fdr_correction()] for the
#' full distinction). Do not read `trend_summary_table` (based on `alpha`,
#' uncorrected) and the FDR results (based on `q`) as answering the same
#' question just because the numbers match. Unlike `alpha`, which this
#' function reports at three conventional levels for context (see
#' `alpha` below), `q = 0.05` has little reason to change: it is the
#' standard target FDR level in the literature this package builds on,
#' and lowering it (e.g. to `0.01`) mainly costs statistical power rather
#' than offering a meaningfully different guarantee.
#'
#' **Statistical assumptions: monotonicity and seasonality.**
#' Every step here (prewhitening, the Contextual Mann-Kendall test,
#' Theil-Sen) is designed around a **monotonic** trend -- a consistent
#' tendency to increase or decrease -- not a periodic/seasonal cycle. If
#' `x` has a seasonal cycle (e.g. raw monthly data with an annual signal),
#' remove it first with [compute_anomalies()] and pass the anomalies to
#' `workflow_tst()`, not the raw seasonal series.
#'
#' **Computational considerations.**
#' Unlike the Mann-Kendall S statistic, the Theil-Sen slope needs every
#' pairwise slope in the series (`n*(n-1)/2` per cell) to take their
#' median -- it cannot be accumulated as a running sum, so it does not
#' scale the same way. For short series (a few dozen time steps) this is
#' unnoticeable; for long series (hundreds to thousands of steps, e.g.
#' multi-decade monthly data) it can become the slowest step in the whole
#' workflow. `theil_sen = TRUE` by default, to match the published TST
#' workflow, but for long series consider `theil_sen = FALSE`, or tune
#' `theil_sen_args` (`max_pairs` to subsample pairs, `n_cores` to
#' parallelise) -- see [slope_estimator()]. `smooth_neighbourhood` stays
#' at [slope_estimator()]'s own default (`FALSE`) unless you set it
#' yourself via `theil_sen_args = list(smooth_neighbourhood = TRUE)`:
#' that mechanism is not part of the published TST method (neither TST
#' nor RTA originally included any such smoothing -- it is this
#' package's own optional addition on top of both), so there is no
#' principled reason for `workflow_tst()` and [workflow_rta()] to
#' default to different
#' values for it -- see `?slope_estimator`'s "Optional queen-neighbourhood
#' smoothing" section for exactly what it does and why it is
#' off by default.
#'
#' **Limitations.** The workflow targets monotonic trends and does not
#' detect abrupt breaks, periodicity, or general nonlinear change. Optional
#' stages permit exploratory variants, but only the complete default
#' configuration reproduces the published TST workflow.
#'
#' **Quality assurance.** Each component is validated independently in its
#' own function. Workflow-level tests additionally verify stage order,
#' optional prewhitening, argument forwarding, shared-cluster behaviour,
#' sequential/parallel equivalence, timings, S3 structure, summaries,
#' plots, and propagation of raw and FDR-corrected results. See `?sptrends`
#' for the internal release protocol and external numerical controls.
#'
#' @param x A `terra::SpatRaster`; each layer is one time step, in
#'   increasing chronological order.
#' @param prewhiten **(Preprocessing)** Logical. If `TRUE` (default), run
#'   [prewhiten()] first. If `FALSE`, `x` is supplied to the trend
#'   test unmodified.
#' @param prewhiten_args A named list of extra arguments forwarded to
#'   [prewhiten()] (e.g. `list(dw_method = "test")`). Workflow-managed
#'   transport and reporting arguments cannot be overridden here.
#'
#'   **Known limitation with irregular time spacing and `TFPW_Y`:** this
#'   workflow has no top-level `t` argument of its own -- if your
#'   observation times are irregularly spaced (not one evenly-spaced
#'   step per layer) and you need a genuine time axis for
#'   [trend_test()]/[slope_estimator()] rather than the default
#'   `1:n` layer index, you must supply the correct `t` yourself in
#'   `prewhiten_args`, `cmk_args` and `theil_sen_args` separately, and
#'   keep them consistent by hand. `TFPW_Y` specifically drops the
#'   first observation (see [prewhiten()]'s own "Value" section), so
#'   the `t` you pass to `cmk_args`/`theil_sen_args` after using
#'   `TFPW_Y` must itself have its own first element dropped to stay
#'   aligned with the one-layer-shorter prewhitened series -- this
#'   workflow does not do that adjustment for you. Not an issue for
#'   `TFPW_WS` (default), which keeps every layer.
#' @param export_dw Logical. If `TRUE`, include the Durbin-Watson
#'   prewhitening diagnostics (`diagnostics` from [prewhiten()])
#'   in the return value, even though the trend test itself does not need
#'   them. Ignored if `prewhiten = FALSE`. Default `FALSE`.
#' @param cmk_args **(Trend detection)** A named list of extra arguments
#'   forwarded to
#'   [trend_test()] (e.g. `list(window_size = 5L)` for a broader CMK
#'   region, or `list(method = "MK", n_cores = 4)`). The default empty
#'   list preserves the 3 by 3 CMK region described by Neeti and Eastman
#'   (2011), as implemented in TerrSet's Kendall module;
#'   changing it creates a TST-inspired variant rather than an exact
#'   reproduction of the published workflow. See `prewhiten_args`
#'   above for a known limitation with irregular time spacing and
#'   `TFPW_Y` if you pass `t` here.
#' @param alpha Numeric vector of significance threshold(s) used for
#'   reporting the (uncorrected) trend result -- supplied to [trend_summary()]
#'   as `alpha`, and used as the single threshold for [trend_maps()]:
#'   `0.05` if it is one of the values in `alpha`, otherwise the
#'   strictest (smallest) value supplied. Defaults to a single `0.05`,
#'   matching `q`'s own default -- a normal analysis picks one alpha
#'   and one q, not several at once. Pass a vector explicitly (e.g.
#'   `alpha = c(0.1, 0.05, 0.01)`) to compare several thresholds side
#'   by side in `trend_summary_table` instead. See "Statistical
#'   assumptions: `alpha` and `q`" above.
#' @param theil_sen **(Slope estimation)** Logical. If `TRUE` (default),
#'   also compute the
#'   Theil-Sen slope (magnitude of change) via [slope_estimator()] --
#'   see "Computational considerations" above before leaving
#'   this on for long time series.
#' @param theil_sen_args A named list of extra arguments forwarded to
#'   [slope_estimator()] (e.g. `list(max_pairs = 20000, n_cores = 4)`).
#'   Ignored if `theil_sen = FALSE`. `smooth_neighbourhood` stays at
#'   [slope_estimator()]'s own default (`FALSE`) unless set here (e.g.
#'   `list(smooth_neighbourhood = TRUE)`) -- see "Theil-Sen has a real
#'   computational cost" above for why. See `prewhiten_args` above for
#'   a known limitation with irregular time spacing and `TFPW_Y` if
#'   you pass `t` here.
#' @param moran_check **(Multiple testing)** Logical. Forwarded to
#'   [fdr_correction()] as-is --
#'   if `TRUE`, runs [spatial_autocorrelation()] on the trend p-value
#'   raster as part of the FDR step. Default `FALSE` (Moran's I stays a
#'   separate, deliberate diagnostic -- see the package vignette).
#' @param fdr_method Character vector supplied as `method` to
#'   [fdr_correction()]. The Usage displays all supported values:
#'   `"BKY"`, `"BH"` and `"BY"`; when omitted, `"BKY"` remains the
#'   default originally used
#'   in the TST methodology (Gutiérrez-Hernández & García, 2025): it
#'   adapts to the estimated proportion of true nulls, gaining
#'   statistical power over BH while targeting false-discovery-rate
#'   control. Supply one method or a vector of methods; BY is the
#'   conservative option for arbitrary dependence. Set to `NULL` to
#'   skip FDR correction entirely.
#' @param q Numeric. Target FDR level, forwarded to [fdr_correction()] --
#'   not an adjusted or renamed `alpha`: it limits the expected false
#'   discovery proportion among rejected hypotheses. See [fdr_correction()]
#'   and "Statistical assumptions: `alpha` and `q`" above.
#' @param bky_implementation `"multtest"` (default, unchanged from previous
#'   versions) or `"original"`, forwarded to [fdr_bky()] via
#'   [fdr_correction()] -- see [fdr_bky()]'s own documentation for the
#'   difference between the two. Ignored if `fdr_method` does not
#'   include `"BKY"`.
#' @param report **(Reporting)** Logical. If `TRUE` (default), each step
#'   prints its own
#'   summary table and draws its diagnostic plots as it runs (same effect
#'   as calling each function directly with `report = TRUE`). Set to
#'   `FALSE` for silent, plot-free programmatic use.
#' @param verbose Logical. Print progress messages and elapsed time for
#'   the complete workflow. Per-stage times are also returned in `timing`.
#' @param n_cores Integer. `1` (default): every step runs sequentially.
#'   `> 1`: builds **one** `parallel::makeCluster()` PSOCK cluster,
#'   shared across every parallel step below (currently the CMK and
#'   Theil-Sen steps) instead of each step creating and tearing down
#'   its own -- avoiding the repeated process-spawn overhead of doing
#'   that two or three times in a row for what is, from the caller's
#'   own perspective, one parallel request. The cluster is closed
#'   automatically when `workflow_tst()` returns. Setting `n_cores`
#'   inside `cmk_args`/`theil_sen_args` directly still works exactly as
#'   before (each step falls back to building its own cluster from
#'   that value) -- but only when this top-level `n_cores` is left at
#'   its default of `1`; when both are set, this one wins and the
#'   per-step `n_cores` inside `cmk_args`/`theil_sen_args` is ignored,
#'   since a single shared cluster and per-step separate ones cannot
#'   both apply to the same call.
#'
#' @return An object of class `c("tst", "sptrends")` (the second, shared
#'   with [workflow_rta()]'s own return value, is for `print()`/`summary()`/
#'   `plot()` -- see [print.sptrends()]): a list
#'   with
#'   \item{prewhiten}{The full output of [prewhiten()], or
#'     `NULL` if `prewhiten = FALSE`.}
#'   \item{dw_diagnostics}{The prewhitening `diagnostics` raster, only if
#'     `export_dw = TRUE`.}
#'   \item{trend}{The `$stats` field of [trend_test()]'s own output (a
#'     `SpatRaster` with the trend statistics themselves) -- not the
#'     full [trend_test()] return value, which also includes
#'     `neighbourhood`/`window_size`.}
#'   \item{trend_summary_table}{The output of [trend_summary()] (invisible
#'     data frame).}
#'   \item{theil_sen}{The `$slope` field of [slope_estimator()]'s own
#'     output (a `SpatRaster`), or `NULL` if `theil_sen = FALSE`. Not
#'     the full [slope_estimator()] return value. Not smoothed by
#'     default (see the `theil_sen_args` note above) -- check
#'     `theil_sen_smoothed` if you need to know whether it was (e.g.
#'     because you set it yourself).}
#'   \item{theil_sen_smoothed}{Logical: whether `theil_sen` was computed
#'     with `smooth_neighbourhood = TRUE`. `FALSE` by default; `NULL` if
#'     `theil_sen = FALSE`.}
#'   \item{fdr}{The output of [fdr_correction()], or `NULL` if
#'     `fdr_method = NULL`.}
#'   \item{timing}{A named list of elapsed seconds per step actually run
#'     (`prewhiten`, `CMK`, `theil_sen`, `fdr` -- only the ones that ran
#'     are present), each measured with `Sys.time()` around that step's
#'     own function call. Coarse (whole-step, not line-by-line) and not a
#'     substitute for a proper profiler, but enough for noticing which
#'     step dominates on your own data.}
#'
#'   Use `print()` for a one-line summary, `summary()` for the full detail,
#'   and `plot()` for a map -- see [plot.sptrends()].
#'
#' @examples
#' \donttest{
#' # Annual mean NDVI from the bundled environmental dataset.
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#'
#' # Run the full workflow: prewhiten -> Contextual Mann-Kendall ->
#' # Theil-Sen -> FDR-BKY (only BKY, not BH -- see the "fdr_method"
#' # argument below for why).
#' result <- workflow_tst(r, report = FALSE, verbose = FALSE)
#'
#' # A "tst" object: printing it gives a one-line-per-step summary
#' # (cells modified by prewhitening, the trend test's cell count, the
#' # Theil-Sen slope range, and how many cells are significant after
#' # FDR-BKY correction).
#' result
#'
#' # Three reports worth seeing: how much of the map is significant at
#' # all, how fast the significant cells are changing, and which way.
#' plot(result, which = "significance")
#' plot(result, which = "slope")
#' plot(result)
#' }
#'
#' @references
#' Source of this workflow (primary reference -- this pipeline is a
#' direct implementation of it):
#' - Gutiérrez-Hernández, O. and García, L.V. (2025) Uncovering true
#'   significant trends in global greening. Remote Sensing Applications:
#'   Society and Environment, 37, 101377. \doi{10.1016/j.rsase.2024.101377}
#'
#' Step 1, selective AR(1) prewhitening (see [prewhiten()] for
#' the full reference list, including the Durbin-Watson gate):
#' - Wang, X.L. and Swail, V.R. (2001) Changes of Extreme Wave Heights in
#'   Northern Hemisphere Oceans and Related Atmospheric Circulation
#'   Regimes. Journal of Climate, 14(10), 2204-2221.
#'
#' Step 2, Contextual Mann-Kendall (see [trend_test()] for
#' the full reference list, including the foundational Mann-Kendall
#' statistic it builds on):
#' - Neeti, N. and Eastman, J.R. (2011) A Contextual Mann-Kendall Approach
#'   for the Assessment of Trend Significance in Image Time Series.
#'   Transactions in GIS, 15(5), 599-611. \doi{10.1111/j.1467-9671.2011.01280.x}
#'
#' Step 3, Theil-Sen slope (see [slope_estimator()] for the full
#' reference list):
#' - Theil, H. (1950) A rank-invariant method of linear and polynomial
#'   regression analysis. Indagationes Mathematicae, 12, 85-91 (Part I;
#'   published in three parts). No DOI available (pre-DOI-era publication).
#' - Sen, P.K. (1968) Estimates of the regression coefficient based on
#'   Kendall's tau. Journal of the American Statistical Association, 63,
#'   1379-1389. \doi{10.1080/01621459.1968.10480934}
#'
#' Step 4, FDR correction -- BKY (default) and BH (see [fdr_bky()] and
#' [fdr_correction()] for the full reference lists):
#' - Benjamini, Y., & Hochberg, Y. (1995) Controlling the False Discovery
#'   Rate: A Practical and Powerful Approach to Multiple Testing. Journal
#'   of the Royal Statistical Society: Series B, 57, 289-300.
#'   \doi{10.1111/j.2517-6161.1995.tb02031.x}
#' - Benjamini, Y., Krieger, A. M., & Yekutieli, D. (2006) Adaptive
#'   Linear Step-Up Procedures that Control the False Discovery Rate.
#'   Biometrika, 93(3), 491-507. \doi{10.1093/biomet/93.3.491}
#' @family pipeline functions
#' @export
workflow_tst <- function(x, prewhiten = TRUE, prewhiten_args = list(),
                 export_dw = FALSE, cmk_args = list(),
                 theil_sen = TRUE, theil_sen_args = list(),
                 alpha = 0.05, moran_check = FALSE,
                 fdr_method = c("BKY", "BH", "BY"), q = 0.05,
                 bky_implementation = c("multtest", "original"),
                 report = TRUE, verbose = TRUE, n_cores = 1) {
  finish_timer <- .sptrends_elapsed_timer("workflow_tst()", verbose)
  on.exit(finish_timer(), add = TRUE)
  if (missing(fdr_method)) fdr_method <- "BKY"
  if (!is.null(fdr_method)) {
    valid_fdr <- c("BH", "BKY", "BY")
    if (length(fdr_method) == 0L || !all(fdr_method %in% valid_fdr)) {
      stop("'fdr_method' must be NULL or one or more of: \"BH\", ",
           "\"BKY\", \"BY\".")
    }
  }
  bky_implementation <- match.arg(bky_implementation)
  if (!inherits(x, "SpatRaster")) stop("'x' must be a terra SpatRaster.")
  n_cores <- .validate_positive_integer(n_cores, "n_cores")
  alpha <- .validate_probability(alpha, "alpha", vector = TRUE)
  q <- .validate_probability(q, "q")
  .check_reserved_args(prewhiten_args, c("x", "report", "verbose"),
                       "prewhiten_args")
  .check_reserved_args(cmk_args,
                       c("x", "alpha", "report", "verbose",
                         "shared_cluster"),
                       "cmk_args")
  .check_reserved_args(theil_sen_args,
                       c("x", "report", "verbose", "shared_cluster"),
                       "theil_sen_args")

  # One PSOCK cluster shared across every parallel step below (CMK,
  # Theil-Sen), instead of each step creating and tearing down its own
  # -- see .sptrends_shared_cluster()'s own comment for why that
  # repeated overhead is worth avoiding. NULL (not a cluster) when
  # n_cores <= 1, so it can be passed straight through unconditionally.
  shared_cl <- .sptrends_shared_cluster(n_cores)
  if (!is.null(shared_cl)) {
    on.exit(parallel::stopCluster(shared_cl), add = TRUE)
  }

  n_steps <- 3 + isTRUE(theil_sen)
  step <- 0

  out <- list(prewhiten = NULL, dw_diagnostics = NULL, trend = NULL,
              trend_summary_table = NULL, theil_sen = NULL, fdr = NULL)
  timing <- list()

  series <- x
  if (isTRUE(prewhiten)) {
    step <- step + 1
    if (verbose) {
      message(sprintf("== Step %d/%d: prewhitening (Wang & Swail) ==",
                       step, n_steps))
    }
    t0 <- Sys.time()
    pw_args <- utils::modifyList(
      list(x = x, verbose = verbose, report = report), prewhiten_args)
    # "prewhiten" as a string, not a bare symbol: this function's own
    # 'prewhiten' parameter (TRUE/FALSE) shadows the function of the
    # same name in this scope. do.call() with a string name resolves it
    # via function lookup (like calling prewhiten(...) directly would),
    # skipping the logical local variable -- do.call(prewhiten, ...)
    # with the bare symbol would instead try to call TRUE/FALSE as a
    # function and error.
    pw_result <- do.call("prewhiten", pw_args)
    timing$prewhiten <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    out$prewhiten <- pw_result
    series <- pw_result$series
    # Wang-Swail keeps every layer of x (just modifies some in place);
    # Yue-Pilon's classic transform loses the first time step, so
    # series has one fewer layer than x -- names(x) would be the wrong
    # length for it. Use the tail of names(x) that actually lines up
    # with what series has, rather than assuming they always match.
    names(series) <- utils::tail(names(x), terra::nlyr(series))
    if (isTRUE(export_dw)) out$dw_diagnostics <- pw_result$diagnostics
  } else if (verbose) {
    step <- step + 1
    message(sprintf(
      "== Step %d/%d: prewhitening skipped (prewhiten = FALSE) ==",
      step, n_steps))
  }

  step <- step + 1
  if (verbose) {
    trend_method <- if (is.null(cmk_args$method)) "CMK" else cmk_args$method
    trend_label <- switch(trend_method,
      CMK = "Contextual Mann-Kendall",
      MK  = "classic Mann-Kendall",
      OLS = "OLS trend test",
      MMK = "modified Mann-Kendall (Hamed and Rao, 1998)",
      trend_method)
    message(sprintf("== Step %d/%d: %s ==", step, n_steps, trend_label))
  }
  t0 <- Sys.time()
  cmk_full_args <- utils::modifyList(
    list(x = series, alpha = alpha, report = report, verbose = verbose,
         shared_cluster = shared_cl),
    cmk_args)
  cmk_result <- do.call(trend_test, cmk_full_args)
  timing$CMK <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out$trend <- cmk_result$stats
  out$trend_summary_table <- trend_summary(cmk_result$stats, alpha = alpha,
                                            path = NULL, verbose = FALSE)

  if (isTRUE(theil_sen)) {
    step <- step + 1
    if (verbose) {
      message(sprintf(
        paste0("== Step %d/%d: Theil-Sen slope (this can be slow -- ",
               "see ?workflow_tst) =="),
        step, n_steps))
    }
    t0 <- Sys.time()
    ts_args <- utils::modifyList(
      list(x = series, report = report, verbose = verbose,
           shared_cluster = shared_cl),
      theil_sen_args)
    theil_result <- do.call(slope_estimator, ts_args)
    out$theil_sen <- theil_result$slope
    timing$theil_sen <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    out$theil_sen_smoothed <- theil_result$smoothed
  }

  if (!is.null(fdr_method)) {
    step <- step + 1
    if (verbose) {
      message(sprintf("== Step %d/%d: FDR correction ==", step, n_steps))
    }
    t0 <- Sys.time()
    out$fdr <- fdr_correction(
      out$trend$p, method = fdr_method, q = q,
      bky_implementation = bky_implementation,
      moran_check = moran_check, report = report, verbose = verbose)
    timing$fdr <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  } else if (verbose) {
    step <- step + 1
    message(sprintf(
      "== Step %d/%d: FDR correction skipped (fdr_method = NULL) ==",
      step, n_steps))
  }

  out$timing <- timing
  class(out) <- c("tst", "sptrends")
  out
}
