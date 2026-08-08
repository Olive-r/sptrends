#' Robust Trend Analysis (RTA): the full pipeline in one call
#'
#' Implements the complete Robust Trend Analysis workflow for monotonic
#' trends in gridded raster time series. One of sptrends' two entry
#' points: chains [slope_estimator()],
#' [trend_test()], and [fdr_correction()] (`method = "BH"`
#' only), in that order, into the **Robust Trend Analysis (RTA)**
#' workflow. Each step's underlying parameters remain available
#' directly on the individual functions -- `workflow_rta()` does not
#' replace them, it saves wiring the calls together for the common
#' case, and returns a single
#' `"rta"` object recognised by [print.sptrends()]. RTA is a
#' **different, shorter** published workflow from this package's other
#' integrated pipeline, [workflow_tst()] -- see "Methodological comparison
#' with TST" below before choosing between them. Both are implemented
#' here because both are genuinely published methods, not because one
#' supersedes the other -- this package aims to be a platform for
#' comparing such methods (see [compare_detections()]), not a vehicle
#' for only its own authors' preferred one.
#'
#' **When should I choose RTA?**: want to reproduce the published 2024
#' workflow exactly -> use `workflow_rta()`. Want the more recent
#' workflow, including selective prewhitening and the adaptive BKY
#' correction -> use [workflow_tst()]. See "Methodological comparison with
#' TST" below for the reasoning behind each design choice.
#'
#' Moran's I ([spatial_autocorrelation()]) is **not** part of this
#' pipeline, for the same reason it is not part of [workflow_tst()]: it is a
#' separate diagnostic for the spatial dependence assumption behind
#' FDR-BH, meant to be run independently (before or after) rather than
#' chained automatically, and pulling it in here would add a dependency
#' this function does not otherwise need -- see the package vignette.
#'
#' **Function type:** **Core function** -- one of the two integrated published
#' workflows in sptrends.
#'
#' @section Typical use:
#' ```
#' raster time series
#'     |
#' workflow_rta()
#'     |
#' Theil-Sen slope + CMK significance + FDR-BH
#'     |
#' one `rta` result containing every stage
#' ```
#' RTA analyses the supplied series without prewhitening. If the input
#' has a seasonal cycle, first use [compute_anomalies()] and pass its
#' `anomalies` raster.
#'
#' @section Methodological details:
#' **How it works.**
#' ```
#' Input raster
#'     |
#' slope_estimator    -- how fast? (Theil-Sen)
#'     |
#' trend_test         -- is there a monotonic trend? (CMK)
#'     |
#' FDR-BH correction  -- which cells survive multiple testing?
#'     |
#' "rta" object
#' ```
#' All three steps always run -- unlike [workflow_tst()], none of them
#' is optional here, matching the published RTA method exactly (see
#' "Comparison with TST" under "Methodological details" for what is
#' deliberately
#' different between the two workflows). Each step's own output is
#' kept in full on the returned object (see "Value" below) --
#' nothing is discarded once a later step begins.
#'
#' **Statistical assumptions: monotonicity and seasonality.**
#' Every step here (the Contextual Mann-Kendall test, Theil-Sen) is
#' designed around a **monotonic** trend -- a consistent tendency to
#' increase or decrease -- not a periodic/seasonal cycle. If `x` has a
#' seasonal cycle (e.g. raw monthly data with an annual signal), remove
#' it first with [compute_anomalies()] and pass the anomalies to
#' `workflow_rta()`,
#' not the raw seasonal series.
#'
#' **Computational considerations.**
#' Unlike the Mann-Kendall S statistic, the Theil-Sen slope needs every
#' pairwise slope in the series (`n*(n-1)/2` per cell) to take their
#' median -- it cannot be accumulated as a running sum, so it does not
#' scale the same way. For short series (a few dozen time steps) this is
#' unnoticeable; for long series (hundreds to thousands of steps, e.g.
#' multi-decade monthly data) it can become the slowest step in the whole
#' workflow. Tune `theil_sen_args` (`max_pairs` to subsample pairs,
#' `n_cores` to parallelise) -- see [slope_estimator()]. `smooth_neighbourhood`
#' is left at [slope_estimator()]'s own default (`FALSE`) unless you set
#' it yourself -- neither RTA nor TST originally included any such
#' smoothing (it is this package's own optional addition on top of both
#' published methods), so [workflow_tst()] does not default to it
#' either, for the
#' same reason.
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
#' question just because the numbers match. `q = 0.05` has little reason
#' to change: it is the standard target FDR level in the literature this
#' package builds on, and lowering it (e.g. to `0.01`) mainly costs
#' statistical power rather than offering a meaningfully different
#' guarantee.
#'
#' **Comparison with TST.** RTA (Gutiérrez-Hernández & García, 2024) and
#' TST (Gutiérrez-Hernández & García, 2025; see [workflow_tst()]) share two
#' pillars -- Theil-Sen and Contextual Mann-Kendall -- but differ in two
#' deliberate, independent ways that this function keeps faithful to the
#' published RTA method, rather than silently reusing TST's later choices.
#'
#' **Difference 1: prewhitening.** TST's first step, selective AR(1)
#' prewhitening, does not appear in RTA. Whether to prewhiten before a
#' Mann-Kendall-family test is a genuine, unresolved methodological debate,
#' not a settled question with one correct answer. Yue & Wang (2002) find
#' that prewhitening can substantially reduce power, particularly when a
#' real trend and real autocorrelation coexist: it can remove part of the
#' trend signal together with the autocorrelation. Bayazit & Önöz (2007)
#' argue the opposite case: skipping prewhitening when autocorrelation is
#' genuinely present can inflate the false-positive rate. RTA does not
#' prewhiten; TST prewhitens selectively, touching only cells whose own
#' Durbin-Watson statistic crosses a gating threshold (see [prewhiten()]),
#' specifically to limit unnecessary power loss. Neither position is
#' implemented here as universally correct.
#'
#' **Difference 2: FDR-BH only, not adaptive BKY.** [workflow_tst()]
#' defaults to the two-stage adaptive BKY correction, which estimates how
#' many tested hypotheses are likely genuinely non-null (\eqn{\hat\pi_0})
#' and relaxes its threshold when substantial real signal is detected; see
#' [fdr_bky()]. RTA instead uses the original, non-adaptive
#' Benjamini-Hochberg procedure, which does not estimate \eqn{\pi_0} and is
#' derived to control FDR under the global null: the least favourable case
#' in which every tested hypothesis could be truly null. It therefore keeps
#' a more conservative guarantee that does not relax as more signal is found.
#' See [fdr_correction()] for the procedure itself.
#'
#' These differences are independent: RTA's choice on one does not imply or
#' require its choice on the other. They happen to be the less adaptive
#' options in this package's two workflows, not because either dictates the
#' other.
#'
#' **Limitations.** RTA is intentionally faithful to its published method:
#' it does not prewhiten, it always estimates a Theil-Sen slope, and it uses
#' FDR-BH rather than exposing alternative corrections. Use
#' [workflow_trends()] when those stages or methods need to be configured.
#'
#' **Quality assurance.** CMK, slope estimation, and FDR correction are
#' validated independently in their module functions. Integration tests
#' verify the RTA stage sequence, argument forwarding, shared parallel
#' resources, timing fields, S3 return structure, reporting, and agreement
#' between direct module calls and workflow outputs. See `?sptrends` for the
#' complete internal and external quality-assurance strategy.
#'
#' @param x A `terra::SpatRaster`; each layer is one time step, in
#'   increasing chronological order. Used directly, unprewhitened -- see
#'   "Comparison with TST" under "Methodological details" below.
#' @param cmk_args A named list of extra arguments passed to
#'   [trend_test()] (e.g. `list(window_size = 5L)` for a broader CMK
#'   region, or `list(method = "MK", n_cores = 4)`). The default empty
#'   list preserves the 3 by 3 CMK region described by Neeti and Eastman
#'   (2011), as implemented in TerrSet's Kendall module;
#'   changing it creates an RTA-inspired variant rather than an exact
#'   reproduction of the published workflow.
#' @param theil_sen_args A named list of extra arguments passed to
#'   [slope_estimator()] (e.g. `list(max_pairs = 20000, n_cores = 4)`).
#'   `smooth_neighbourhood` is left at [slope_estimator()]'s own default
#'   (`FALSE`, same as [workflow_tst()]) unless you set it yourself -- see
#'   "Computational considerations" above.
#' @param alpha Numeric vector of significance thresholds used for
#'   reporting the (uncorrected) trend result -- passed to
#'   [trend_summary()] as `alpha`, and used as the single threshold for
#'   [trend_maps()]: `0.05` if it is one of the values in `alpha` (the
#'   default vector includes it), otherwise the strictest (smallest)
#'   value supplied. The three default values are not interchangeable:
#'   `0.05` is the conventional standard and the one actually used for
#'   the map; `0.1` is a more liberal threshold not unusual in
#'   exploratory trend studies; `0.01` is markedly more conservative.
#'   All three are shown side by side in `trend_summary_table` for
#'   context, not as equally valid choices.
#' @param q Numeric. Target FDR level, passed to [fdr_correction()] --
#'   not an adjusted or renamed `alpha`: it limits the expected false
#'   discovery proportion among rejected hypotheses. See [fdr_correction()]
#'   for the full distinction.
#' @param report Logical. If `TRUE` (default), each step prints its own
#'   summary table and draws its diagnostic plots as it runs (same
#'   effect as calling each function directly with `report = TRUE`).
#'   Set to `FALSE` for silent, plot-free programmatic use.
#' @param verbose Logical. Print progress messages and elapsed time for
#'   the complete workflow. Per-stage times are also returned in `timing`.
#' @param n_cores Integer. `1` (default): every step runs sequentially.
#'   `> 1`: builds **one** `parallel::makeCluster()` PSOCK cluster,
#'   shared across every parallel step below (currently the CMK and
#'   Theil-Sen steps) instead of each step creating and tearing down
#'   its own -- avoiding the repeated process-spawn overhead of doing
#'   that two or three times in a row for what is, from the caller's
#'   own perspective, one parallel request. The cluster is closed
#'   automatically when `workflow_rta()` returns. Setting `n_cores`
#'   inside `cmk_args`/`theil_sen_args` directly still works exactly as
#'   before (each step falls back to building its own cluster from
#'   that value) -- but only when this top-level `n_cores` is left at
#'   its default of `1`; when both are set, this one wins and the
#'   per-step `n_cores` inside `cmk_args`/`theil_sen_args` is ignored,
#'   since a single shared cluster and per-step separate ones cannot
#'   both apply to the same call. Mirrors [workflow_tst()]'s own
#'   identical `n_cores` argument exactly.
#'
#' @return An object of class `c("rta", "sptrends")` (the second, shared
#'   with [workflow_tst()]'s own return value, is for `print()`/`summary()`/
#'   `plot()` -- see [print.sptrends()]): a list
#'   with
#'   \item{theil_sen}{The output of [slope_estimator()]. Not smoothed
#'     by default -- check `theil_sen_smoothed` if you need to know
#'     whether it was (e.g. because you set it yourself).}
#'   \item{theil_sen_smoothed}{Logical: whether `theil_sen` was computed
#'     with `smooth_neighbourhood = TRUE`. `FALSE` by default.}
#'   \item{trend}{The output of [trend_test()].}
#'   \item{trend_summary_table}{The output of [trend_summary()]
#'     (invisible data frame).}
#'   \item{fdr}{The output of [fdr_correction()], `method = "BH"` only.}
#'   \item{timing}{A named list of elapsed seconds per step
#'     (`theil_sen`, `cmk`, `fdr`), each measured with `Sys.time()`
#'     around that step's own function call. Coarse (whole-step, not
#'     line-by-line) and not a substitute for a proper profiler, but
#'     enough for noticing which step dominates on your own data.}
#'
#'   Use [print.sptrends()] for a one-line summary.
#'
#' @examples
#' \donttest{
#' # Annual mean NDVI from the bundled environmental dataset.
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#'
#' # Run the full workflow: Theil-Sen -> Contextual Mann-Kendall ->
#' # FDR-BH (only BH -- see "Comparison with TST" above
#' # for why this function does not offer BKY as an option).
#' result <- workflow_rta(r, report = FALSE, verbose = FALSE)
#'
#' # An "rta" object: printing it gives a one-line-per-step summary
#' # (the Theil-Sen slope range, the trend test's cell count, and how
#' # many cells are significant after FDR-BH correction).
#' result
#'
#' # The same plot() methods used throughout this package, rather than
#' # reconstructing either view by hand: which cells are significant
#' # after FDR-BH, and how fast those cells are changing.
#' plot(result, which = "significance")
#' plot(result, which = "slope")
#' plot(result)
#' }
#'
#' @references
#' Source of this workflow (primary reference -- this pipeline is a
#' direct implementation of it):
#' - Gutiérrez-Hernández, O. and García, L.V. (2024) Robust Trend Analysis
#'   in Environmental Remote Sensing: A Case Study of Cork Oak Forest
#'   Decline. Remote Sensing, 16(20), 3886. \doi{10.3390/rs16203886}
#'
#' Step 1, Theil-Sen slope (see [slope_estimator()] for the full
#' reference list):
#' - Theil, H. (1950) A rank-invariant method of linear and polynomial
#'   regression analysis. Indagationes Mathematicae, 12, 85-91 (Part I;
#'   published in three parts). No DOI available (pre-DOI-era publication).
#' - Sen, P.K. (1968) Estimates of the regression coefficient based on
#'   Kendall's tau. Journal of the American Statistical Association, 63,
#'   1379-1389. \doi{10.1080/01621459.1968.10480934}
#'
#' Step 2, Contextual Mann-Kendall (see [trend_test()] for
#' the full reference list, including the foundational Mann-Kendall
#' statistic it builds on):
#' - Neeti, N. and Eastman, J.R. (2011) A Contextual Mann-Kendall Approach
#'   for the Assessment of Trend Significance in Image Time Series.
#'   Transactions in GIS, 15(5), 599-611. \doi{10.1111/j.1467-9671.2011.01280.x}
#'
#' Step 3, FDR-BH correction (see [fdr_correction()] for the complete
#' reference list):
#' - Benjamini, Y. and Hochberg, Y. (1995) Controlling the False
#'   Discovery Rate: A Practical and Powerful Approach to Multiple
#'   Testing. Journal of the Royal Statistical Society: Series B, 57,
#'   289-300. \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' On whether to prewhiten before a Mann-Kendall-family trend test (see
#' "Comparison with TST" under "Methodological details" above):
#' - Yue, S. and Wang, C.Y. (2002) Applicability of prewhitening to
#'   eliminate the influence of serial correlation on the Mann-Kendall
#'   test. Water Resources Research, 38(6), 4-1-4-6.
#'   \doi{10.1029/2001WR000861}
#' - Bayazit, M. and Önöz, B. (2007) To prewhiten or not to prewhiten in
#'   trend analysis? Hydrological Sciences Journal, 52(4), 611-624.
#'   \doi{10.1623/hysj.52.4.611}
#' @family pipeline functions
#' @export
workflow_rta <- function(x, cmk_args = list(), theil_sen_args = list(),
                 alpha = 0.05, q = 0.05, report = TRUE, verbose = TRUE,
                 n_cores = 1) {
  finish_timer <- .sptrends_elapsed_timer("workflow_rta()", verbose)
  on.exit(finish_timer(), add = TRUE)
  if (!inherits(x, "SpatRaster")) stop("'x' must be a terra SpatRaster.")
  n_cores <- .validate_positive_integer(n_cores, "n_cores")
  alpha <- .validate_probability(alpha, "alpha", vector = TRUE)
  q <- .validate_probability(q, "q")
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
  # repeated overhead is worth avoiding, and workflow_tst()'s own
  # identical mechanism, which this mirrors exactly.
  shared_cl <- .sptrends_shared_cluster(n_cores)
  if (!is.null(shared_cl)) {
    on.exit(parallel::stopCluster(shared_cl), add = TRUE)
  }

  out <- list(theil_sen = NULL, trend = NULL, trend_summary_table = NULL,
              fdr = NULL)
  n_steps <- 3
  step <- 0
  timing <- list()

  step <- step + 1
  if (verbose) {
    message(sprintf("== Step %d/%d: Theil-Sen slope (magnitude/direction) ==",
                     step, n_steps))
  }
  t0 <- Sys.time()
  ts_args <- utils::modifyList(
    list(x = x, report = report, verbose = verbose,
         shared_cluster = shared_cl),
    theil_sen_args)
  theil_result <- do.call(slope_estimator, ts_args)
  out$theil_sen <- theil_result$slope
  timing$theil_sen <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out$theil_sen_smoothed <- theil_result$smoothed

  step <- step + 1
  if (verbose) {
    trend_method <- if (is.null(cmk_args$method)) "CMK" else cmk_args$method
    trend_label <- switch(trend_method,
      CMK = "Contextual Mann-Kendall",
      MK  = "classic Mann-Kendall",
      OLS = "OLS trend test",
      MMK = "modified Mann-Kendall (Hamed and Rao, 1998)",
      trend_method)
    message(sprintf("== Step %d/%d: %s (significance) ==",
                     step, n_steps, trend_label))
  }
  t0 <- Sys.time()
  cmk_full_args <- utils::modifyList(
    list(x = x, alpha = alpha, report = report, verbose = verbose,
         shared_cluster = shared_cl),
    cmk_args)
  cmk_result <- do.call(trend_test, cmk_full_args)
  timing$CMK <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out$trend <- cmk_result$stats
  out$trend_summary_table <- trend_summary(cmk_result$stats, alpha = alpha,
                                            path = NULL, verbose = FALSE)

  step <- step + 1
  if (verbose) {
    message(sprintf("== Step %d/%d: FDR correction (BH) ==", step, n_steps))
  }
  t0 <- Sys.time()
  out$fdr <- fdr_correction(out$trend$p, method = "BH", q = q,
                             moran_check = FALSE, report = report,
                             verbose = verbose)
  timing$fdr <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  out$timing <- timing
  class(out) <- c("rta", "sptrends")
  out
}

#' @noRd
.print_rta <- function(x, ...) {
  cat("<Robust Trend Analysis (RTA) result>\n")

  slope_vals <- terra::values(x$theil_sen, mat = FALSE)
  n_valid_slope <- sum(!is.na(slope_vals))
  cat(sprintf("Theil-Sen slope: %d cells, range [%.4g, %.4g]\n",
              n_valid_slope, min(slope_vals, na.rm = TRUE),
              max(slope_vals, na.rm = TRUE)))

  s_name <- if ("Sm" %in% names(x$trend)) {
    "Sm"
  } else if ("S" %in% names(x$trend)) {
    "S"
  } else {
    "beta"
  }
  p_vals <- terra::values(x$trend$p, mat = FALSE)
  n_cells <- sum(!is.na(p_vals))
  cat(sprintf("Trend test: %d cells (%s statistic)\n", n_cells, s_name))

  n_sig <- sum(x$fdr$reject_BH, na.rm = TRUE)
  cat(sprintf("Significant after FDR-BH: %d (%.1f%%)\n",
              n_sig, 100 * n_sig / n_cells))

  invisible(x)
}

#' @noRd
.summary_rta <- function(object, ...) {
  cat("=== FDR correction (BH) -- the actual RTA result ===\n")
  fdr_tab <- fdr_summary(object$fdr)

  cat("\n=== Theil-Sen slope ===\n")
  print(summary(terra::values(object$theil_sen, mat = FALSE)))

  cat("\n=== Trend test, uncorrected (diagnostic only -- not the RTA result; see FDR correction above) ===\n")
  print(object$trend_summary_table)

  invisible(list(trend = object$trend_summary_table, fdr = fdr_tab))
}

#' @noRd
.plot_rta <- function(x, which = c("direction", "significance", "trend",
                                   "slope", "slope_map", "slope_direction",
                                   "slope_hist", "slope_bar", "pvalue_map",
                                   "pvalue_significance", "pvalue_hist",
                                   "pvalue_bar"),
                      smooth = TRUE, ...) {
  which <- match.arg(which)

  if (which == "direction") {
    already_smoothed_dir <- isTRUE(x$theil_sen_smoothed)
    slope_for_direction <- if (already_smoothed_dir) {
      x$theil_sen
    } else {
      suppressMessages(
        .smooth_slope_for_direction(x$theil_sen, isTRUE(smooth))
      )
    }
    direction <- direction_map(x$trend, x$fdr, slope = slope_for_direction,
                                    method = "BH",
                                    verbose = FALSE)
    fdr_direction_plot(direction, main = "RTA direction map")
  } else if (which == "significance") {
    fdr_significance_maps(x$fdr)
  } else if (which == "slope") {
    reject <- x$fdr$reject_BH
    already_smoothed <- isTRUE(x$theil_sen_smoothed)
    apply_smooth_now <- isTRUE(smooth) && !already_smoothed
    slope_sig <- .mask_and_smooth_slope(x$theil_sen, reject, apply_smooth_now)

    if (already_smoothed && !isTRUE(smooth)) {
      message(
        "Note: smooth = FALSE was requested, but this slope was already ",
        "computed with smoothing applied at source (theil_sen_args = ",
        "list(smooth_neighbourhood = TRUE) was requested) -- the unsmoothed ",
        "values were never kept to show instead. Rerun with theil_sen_args ",
        "= list(smooth_neighbourhood = FALSE) for a genuinely unsmoothed ",
        "result.")
    }
    subtitle <- "Theil-Sen slope, FDR-BH significant cells"

    range_lim <- .robust_diverging_range(slope_sig)
    # fill_range = TRUE: without it, cells beyond the robust range would
    # render as blank (NA) instead of saturating to the extreme colour --
    # see slope_map()'s own comment on this for the full reasoning.
    terra::plot(slope_sig, col = .sptrends_diverging_palette(50),
                range = range_lim, fill_range = TRUE, main = subtitle)
  } else if (which %in% c("slope_map", "slope_direction", "slope_hist",
                           "slope_bar")) {
    if (is.null(x$theil_sen)) {
      stop("No Theil-Sen slope in this rta object -- this should never ",
           "happen, since workflow_rta() always computes it; please ",
           "report this as a bug.")
    }
    switch(which,
      slope_map       = slope_map(x$theil_sen),
      slope_direction = slope_direction_map(x$theil_sen, ...),
      slope_hist      = slope_histogram(x$theil_sen, ...),
      slope_bar       = slope_direction_barplot(x$theil_sen, ...)
    )
  } else if (which %in% c("pvalue_map", "pvalue_significance", "pvalue_hist",
                           "pvalue_bar")) {
    switch(which,
      pvalue_map           = .plot_pvalue_map(x$trend, ...),
      pvalue_significance  = .plot_pvalue_significance(x$trend, ...),
      pvalue_hist          = .plot_pvalue_hist(x$trend, ...),
      pvalue_bar           = .plot_pvalue_bar(x$trend, ...)
    )
  } else {
    trend_maps(x$trend)
  }

  invisible(x)
}
