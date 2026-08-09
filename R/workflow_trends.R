#' Configure a monotonic or linear trend-analysis workflow
#'
#' This function analyses monotonic or linear temporal trends using a
#' user-selected trend test -- available methods include the
#' Mann-Kendall family (`"CMK"`, `"MK"`, `"MMK"`) and ordinary
#' least-squares regression (`"OLS"`); it does not detect abrupt
#' changes, periodicity, or general nonlinear temporal patterns. Like
#' [workflow_tst()] and [workflow_rta()], each of which implements one
#' specific published analytical workflow, this function analyses that
#' same class of trend. Unlike those two, `workflow_trends()` lets
#' you choose your own combination of prewhitening, trend testing,
#' slope estimation and multiple-testing correction methods -- useful
#' when no single published workflow matches what a given analysis
#' needs, or when comparing how sensitive a result is to that choice
#' (see the package's own validation framework, [sim_trend_stack()]/
#' [compare_detections()], for a way to evaluate that sensitivity
#' against data with a known answer, rather than only this function's
#' own runtime warnings below).
#'
#' **Function type:** **Core function** -- a configurable composition of the
#' same analytical stages used by the published workflows.
#'
#' @section Typical use:
#' ```
#' raster time series
#'     |
#' workflow_trends()
#'     |
#' optional prewhitening -> trend test -> optional slope -> optional FDR
#'     |
#' one `sptrends` result containing every computed stage
#' ```
#' Use this configurable workflow when the fixed published combinations
#' in [workflow_tst()] and [workflow_rta()] do not match the analysis.
#' For seasonal input, first pass the `$anomalies` component returned by
#' [compute_anomalies()].
#'
#' @section Methodological details:
#' **Relationship to TST.**
#' True Significant Trends (TST) is the methodological origin of
#' sptrends and supplies the architecture of the complete workflow:
#' temporal preprocessing, trend testing, slope estimation and
#' multiple-testing correction. TST includes all these stages, but not
#' every method now available here. `workflow_trends()` generalises the
#' architecture with additional methods and optional stages. A complete
#' configuration retains the analytical philosophy inherited from TST;
#' the name TST remains reserved for the published methods and sequence
#' reproduced by [workflow_tst()].
#'
#' Every method available on the individual functions this orchestrates
#' -- [prewhiten()], [trend_test()], [slope_estimator()],
#' [fdr_correction()] -- can be freely combined here, with two
#' exceptions this function warns about explicitly rather than
#' silently allowing (see "Combinations this function warns about"
#' below). The trend test and slope estimator answer different
#' questions and need not use the same statistical framework --
#' mixing a robust trend test with a parametric slope estimator (or
#' vice versa) is unconventional, but not warned about here.
#' Nevertheless, users remain responsible for ensuring the selected
#' combination is appropriate for their own data and scientific
#' objective; other combinations this function does not itself flag
#' can still be questionable for a given dataset (e.g. `"CMK"` on an
#' already spatially smoothed series -- see "Combinations this
#' function warns about" below for the one combination this function
#' does check for explicitly).
#'
#' **How it works.**
#' \enumerate{
#'   \item (optional) `prewhiten` -- removes serial autocorrelation.
#'   \item `trend_test` -- is there a monotonic trend?
#'   \item (optional) `slope_estimator` -- how fast? and
#'     (optional) FDR correction -- which cells survive multiple
#'     testing? These two are independent branches after the trend
#'     test, not a sequence: FDR correction depends only on the test's
#'     own p-values, not on whether a slope was estimated at all.
#' }
#' The same ordering logic as the two published workflows applies:
#' prewhitening (when requested) happens before the trend test, since
#' the test assumes independent observations; the slope (when
#' requested) is estimated on the same (optionally prewhitened) series
#' the test itself used; and FDR correction needs the test's own
#' p-values to exist first, but nothing else -- it does not depend on
#' whether a slope was estimated at all. Unlike prewhitening and FDR
#' correction, which were already optional here from the start, slope
#' estimation being skippable too (`slope_method = NULL`) was added
#' later, to bring this function to full parity with `workflow_tst()`'s
#' own `theil_sen = FALSE` -- see NEWS.md for when.
#'
#' **Statistical assumptions: `alpha` and `q`.**
#' As with `workflow_tst()`/`workflow_rta()`: *q* (inside `fdr_args`,
#' e.g. `fdr_args = list(q = 0.1)`) bounds the expected proportion of
#' false positives among cells called significant *after* correction --
#' it is not the same quantity as an uncorrected per-cell \eqn{\alpha}.
#' Applying a per-cell significance threshold across many cells
#' without multiplicity control can substantially increase the number
#' of false-positive detections, the multiple-testing problem
#' `fdr_method` exists to address. See `?workflow_tst`'s own "A note
#' on alpha and q" section for the fuller explanation, and
#' [fdr_correction()] for the distinction in full.
#'
#' **Statistical assumptions: monotonicity and seasonality.**
#' As with `workflow_tst()`/`workflow_rta()`: every method this function
#' can combine is designed around a **monotonic** trend, not a
#' periodic/seasonal cycle. If `x` has a seasonal cycle (e.g. raw
#' monthly data with an annual signal), remove it first with
#' [compute_anomalies()] and pass the anomalies here, not the raw
#' seasonal series.
#'
#' **Computational considerations.** Runtime depends primarily on the
#' selected methods and raster size. Theil-Sen and repeated-median slopes
#' are more expensive than OLS; CMK and parallel-capable stages can reuse
#' one workflow-level PSOCK cluster. Per-stage elapsed times are retained
#' in the returned object.
#'
#' **Limitations.** Configurability does not make every combination
#' scientifically interchangeable. The caller remains responsible for
#' matching the chosen assumptions to the data; the explicit safeguards
#' below cover important known conflicts, not every possible misuse.
#'
#' **Methodological safeguards.** Combining any prewhitening method with
#' `trend_method = "MMK"` is warned against because both address temporal
#' autocorrelation by different mechanisms. Applying both is redundant and
#' may compound the correction. The warning is issued before computation;
#' it does not forcibly stop an explicitly requested analysis.
#'
#' **Quality assurance.** Integration tests exercise every available method
#' family, optional stages, invalid and questionable combinations, protected
#' argument forwarding, BH/BKY/BY display paths, exact warning behaviour,
#' shared parallel resources, timing fields, and returned S3 structure.
#' Component results are compared with direct calls to [prewhiten()],
#' [trend_test()], [slope_estimator()], and [fdr_correction()]. See
#' `?sptrends` for the complete internal and external quality strategy.
#'
#' @param x A `terra::SpatRaster` with one layer per time step,
#'   ordered chronologically. All layers must share the same geometry.
#'   Temporal spacing is assumed regular unless explicit time
#'   coordinates are supplied through the relevant underlying
#'   function's own arguments (e.g. `trend_args = list(t = ...)`,
#'   `slope_args = list(t = ...)`).
#' @param prewhiten_method Which prewhitening method to apply before
#'   trend testing, or `"none"` to skip this step entirely. One of
#'   `"TFPW_WS"` (default; Wang and Swail, 2001), `"TFPW_Y"` (Yue,
#'   Pilon and Cavadias, 2002), `"TFPW_Z"` (Zhang, Vincent, Hogg and
#'   Niitsoo, 2000), `"VCTFPW"` (Wang, Chen, Becker and Liu, 2015), or
#'   `"none"`. See [prewhiten()]'s own `method` argument for the full
#'   comparison between these four, and "Combinations this function
#'   warns about" below for the one combination involving
#'   `trend_method` worth knowing about before combining.
#' @param trend_method Which trend test to run. One of `"CMK"`
#'   (default; Neeti and Eastman, 2011), `"MK"` (classic Mann-Kendall),
#'   `"OLS"` (tests whether an ordinary-least-squares regression slope
#'   differs from zero, rather than a rank-based statistic -- see
#'   [trend_test()]'s own `"OLS"` entry for the specific `p` this
#'   returns and the assumptions its own inference relies on), or
#'   `"MMK"` (Hamed and Rao, 1998). See [trend_test()]'s own `method`
#'   argument for the full comparison between these four.
#' @param slope_method Which slope estimator to run, or `NULL` to skip
#'   this step entirely (matching `workflow_tst()`'s own
#'   `theil_sen = FALSE`). One of `"TS"` (default; Theil-Sen), `"OLS"`,
#'   `"RM"` (Siegel's repeated median), or `NULL`. See
#'   [slope_estimator()]'s own `method` argument for the full
#'   comparison between the three.
#' @param fdr_method Which multiple-testing correction to apply, or
#'   `NULL` to skip this step. One of `"BKY"` (default; Benjamini,
#'   Krieger and Yekutieli, 2006), `"BH"` (Benjamini and Hochberg,
#'   1995), or `"BY"` (Benjamini and Yekutieli, 2001). `"BY"` is an
#'   opt-in safeguard when control under arbitrary dependence is
#'   scientifically required; it is generally more conservative than
#'   `"BH"` and adaptive `"BKY"`. See
#'   [fdr_correction()]'s own `method` argument, and `fdr_by()`'s own
#'   documentation, for the full reasoning.
#' @param prewhiten_args,trend_args,slope_args,fdr_args Named lists of
#'   extra arguments forwarded to the corresponding underlying
#'   function, beyond `method` itself (which this function's own
#'   arguments above already control) -- e.g.
#'   `trend_args = list(ties = TRUE, window_size = 5L)`,
#'   `fdr_args = list(q = 0.1)`. For CMK, omitting `window_size` preserves
#'   the 3 by 3 CMK region described by Neeti and Eastman (2011), as
#'   implemented in TerrSet's Kendall module.
#'   Arguments already controlled directly by `workflow_trends()`
#'   itself -- `x`, `method`, `report`, `verbose`, and cluster-
#'   management arguments (`n_cores`/`shared_cluster`) -- must not be
#'   repeated in these lists; a clear error is raised instead of silently
#'   overriding the workflow configuration.
#' @param n_cores Integer. Number of workers made available to stages
#'   that support parallel execution -- a single PSOCK cluster is
#'   created once and reused where possible, not one per step. Not
#'   every method actually uses it: `slope_method = "OLS"`'s own
#'   vectorised computation and the current `"RM"` implementation, for
#'   instance, do not. `1` (default): sequential. Requesting more cores
#'   than are actually available falls back to the number detected,
#'   with a message (see `?trend_test`'s own `n_cores` for the exact
#'   behaviour every step here shares).
#' @param report Logical. Forwarded to every step's own `report`
#'   argument -- if `TRUE`, each analytical stage may print its own
#'   component-level summary and draw its own diagnostic plots as it
#'   runs, not only a single consolidated summary of the whole
#'   workflow at the end. Set `FALSE` to suppress this intermediate
#'   output and inspect the returned object directly instead.
#' @param verbose Logical. Forwarded to every step's own `verbose`
#'   argument, and controls this function's own step-by-step progress
#'   messages and elapsed time. Per-stage times are returned in `timing`.
#'
#' @return A list of class
#'   `c("workflow_trends", "sptrends")`:
#'   \item{prewhiten}{The output of [prewhiten()], or `NULL` if
#'     `prewhiten_method = "none"`.}
#'   \item{trend}{The multilayer `SpatRaster` returned as `trend_test()`'s
#'     own `$stats` (the test statistic, raw `p`, and any
#'     method-specific layers -- see `?trend_test`'s own "Value"
#'     section for exactly which, since this differs by
#'     `trend_method`).}
#'   \item{trend_summary_table}{The output of [trend_summary()] on
#'     `trend`.}
#'   \item{slope}{The single-layer `SpatRaster` extracted from
#'     `slope_estimator()`'s own `$slope` (not the full `"slope"`
#'     object that function itself returns), or `NULL` if
#'     `slope_method = NULL`.}
#'   \item{fdr}{The output of [fdr_correction()], or `NULL` if
#'     `fdr_method = NULL`.}
#'   \item{timing}{A named list of per-step elapsed seconds, only for
#'     steps that actually ran.}
#'
#' @examples
#' \donttest{
#' # Annual mean NDVI from the bundled environmental dataset.
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#'
#' # A combination neither TST nor RTA offers on its own: Yue-Pilon
#' # prewhitening, classic Mann-Kendall, OLS slope, standard (BH) FDR.
#' result <- workflow_trends(r, prewhiten_method = "TFPW_Y",
#'                             trend_method = "MK", slope_method = "OLS",
#'                             fdr_method = "BH",
#'                             report = FALSE, verbose = FALSE)
#' result
#' plot(result, which = "significance")
#' plot(result, which = "slope")
#'
#' # Siegel's repeated median as the slope estimator instead of OLS --
#' # useful when a dataset is suspected to have enough outliers or
#' # leverage points that even Theil-Sen's own ~29% breakdown point
#' # might not fully resist them (see ?slope_estimator's own "RM" entry
#' # under `method`).
#' result_rm <- workflow_trends(r, prewhiten_method = "TFPW_Y",
#'                                trend_method = "MK", slope_method = "RM",
#'                                fdr_method = "BH",
#'                                report = FALSE, verbose = FALSE)
#' plot(result_rm, which = "slope")
#'
#' # Skipping optional stages entirely -- slope_method/fdr_method both
#' # NULL leaves $slope/$fdr NULL too, rather than some placeholder
#' # value, making the modularity concrete rather than just described.
#' result_test_only <- workflow_trends(r, prewhiten_method = "none",
#'                                       trend_method = "MMK",
#'                                       slope_method = NULL,
#'                                       fdr_method = NULL,
#'                                       report = FALSE, verbose = FALSE)
#' is.null(result_test_only$slope)
#' is.null(result_test_only$fdr)
#' }
#'
#' \dontrun{
#' # A combination this function warns about: MMK plus prewhitening.
#' # Not run automatically (unlike the examples above) specifically
#' # because its only purpose here is to demonstrate the warning
#' # itself, which would otherwise fire on every R CMD check and every
#' # example() call.
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' result_warns <- workflow_trends(r, trend_method = "MMK",
#'                                   report = FALSE, verbose = FALSE)
#' }
#'
#' @references
#' Primary method reference (`prewhiten_method = "TFPW_WS"`, the
#' default):
#' - Wang, X.L. and Swail, V.R. (2001) Changes of Extreme Wave Heights
#'   in Northern Hemisphere Oceans and Related Atmospheric Circulation
#'   Regimes. Journal of Climate, 14(10), 2204-2221.
#'
#' Primary method reference (`prewhiten_method = "TFPW_Y"`):
#' - Yue, S., Pilon, P., Phinney, B. and Cavadias, G. (2002) The
#'   influence of autocorrelation on the ability to detect trend in
#'   hydrological series. Hydrological Processes, 16(9), 1807-1829.
#'   \doi{10.1002/hyp.1095}
#'
#' Primary method reference (`prewhiten_method = "TFPW_Z"`):
#' - Zhang, X., Vincent, L.A., Hogg, W.D. and Niitsoo, A. (2000)
#'   Temperature and precipitation trends in Canada during the 20th
#'   century. Atmosphere-Ocean, 38(3), 395-429.
#'   \doi{10.1080/07055900.2000.9649654}
#'
#' Primary method reference (`prewhiten_method = "VCTFPW"`):
#' - Wang, W., Chen, Y., Becker, S. and Liu, B. (2015) Variance
#'   Correction Prewhitening Method for Trend Detection in
#'   Autocorrelated Data. Journal of Hydrologic Engineering, 20(12),
#'   04015033. \doi{10.1061/(ASCE)HE.1943-5584.0001234}
#'
#' Primary method reference (`trend_method = "CMK"`, the default):
#' - Neeti, N. and Eastman, J.R. (2011) A Contextual Mann-Kendall
#'   Approach for the Assessment of Trend Significance in Image Time
#'   Series. Transactions in GIS, 15(5), 599-611.
#'   \doi{10.1111/j.1467-9671.2011.01280.x}
#'
#' Primary method reference (`trend_method = "MMK"`):
#' - Hamed, K.H. and Rao, A.R. (1998) A modified Mann-Kendall trend
#'   test for autocorrelated data. Journal of Hydrology, 204(1-4),
#'   182-196. \doi{10.1016/S0022-1694(97)00125-X}
#'
#' Primary method reference (`fdr_method = "BH"`):
#' - Benjamini, Y. and Hochberg, Y. (1995) Controlling the False
#'   Discovery Rate: A Practical and Powerful Approach to Multiple
#'   Testing. Journal of the Royal Statistical Society: Series B, 57,
#'   289-300. \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Primary method reference (`fdr_method = "BKY"`, the default):
#' - Benjamini, Y., Krieger, A.M. and Yekutieli, D. (2006) Adaptive
#'   Linear Step-Up Procedures that Control the False Discovery Rate.
#'   Biometrika, 93(3), 491-507. \doi{10.1093/biomet/93.3.491}
#'
#' Primary method reference (`fdr_method = "BY"`, not recommended by
#' default -- see the `fdr_method` argument above for why):
#' - Benjamini, Y. and Yekutieli, D. (2001) The control of the false
#'   discovery rate in multiple testing under dependency. Annals of
#'   Statistics, 29(4), 1165-1188. \doi{10.1214/aos/1013699998}
#'
#' `trend_method = "MK"`, `trend_method = "OLS"`, and
#' `slope_method = "TS"`/`"OLS"`/`"RM"` do have their own foundational
#' references (Mann, 1945 and Kendall, 1975; Legendre, 1805 and Gauss,
#' 1809; Theil, 1950 and Sen, 1968; Siegel, 1982, respectively) --
#' cited in full under `?trend_test` and `?slope_estimator` rather than
#' repeated here.
#'
#' @family pipeline functions
#' @export
workflow_trends <- function(
    x,
    prewhiten_method = c("TFPW_WS", "TFPW_Y", "TFPW_Z", "VCTFPW",
                          "none"),
    trend_method = c("CMK", "MK", "OLS", "MMK"),
    slope_method = c("TS", "OLS", "RM"),
    fdr_method = c("BKY", "BH", "BY"),
    prewhiten_args = list(), trend_args = list(),
    slope_args = list(), fdr_args = list(),
    n_cores = 1, report = TRUE, verbose = TRUE) {
  finish_timer <- .sptrends_elapsed_timer("workflow_trends()", verbose)
  on.exit(finish_timer(), add = TRUE)
  prewhiten_method <- match.arg(prewhiten_method)
  trend_method <- match.arg(trend_method)
  if (missing(slope_method)) slope_method <- "TS"
  if (!is.null(slope_method)) {
    slope_method <- match.arg(slope_method, c("TS", "OLS", "RM"))
  }
  if (missing(fdr_method)) fdr_method <- "BKY"
  if (!is.null(fdr_method)) {
    valid_fdr <- c("BH", "BKY", "BY")
    if (length(fdr_method) == 0L || !all(fdr_method %in% valid_fdr)) {
      stop("'fdr_method' must be NULL or one or more of: \"BH\", ",
           "\"BKY\", \"BY\".")
    }
  }
  if (!inherits(x, "SpatRaster")) stop("'x' must be a terra SpatRaster.")
  n_cores <- .validate_positive_integer(n_cores, "n_cores")
  .check_reserved_args(
    prewhiten_args, c("x", "method", "report", "verbose"),
    "prewhiten_args")
  .check_reserved_args(
    trend_args, c("x", "method", "report", "verbose", "n_cores",
                  "shared_cluster"),
    "trend_args")
  .check_reserved_args(
    slope_args, c("x", "method", "report", "verbose", "n_cores",
                  "shared_cluster"),
    "slope_args")
  .check_reserved_args(
    fdr_args, c("p", "method", "report", "verbose"), "fdr_args")

  # One family of combinations this function warns about rather than silently
  # allowing -- see "Combinations this function warns about" above for
  # the reasoning behind it. Checked, and warned about, before any
  # step actually runs, ahead of any step progress messages.
  if (prewhiten_method != "none" && trend_method == "MMK") {
    warning(
      "prewhiten_method != \"none\" combined with trend_method = \"MMK\": ",
      "both correct temporal autocorrelation, by different means -- see ",
      "?trend_test's own \"MMK\" documentation. Proceeding anyway.",
      call. = FALSE)
  }

  # One PSOCK cluster shared across every parallel-capable step below,
  # instead of each step creating and tearing down its own -- see
  # .sptrends_shared_cluster()'s own comment for why that repeated
  # overhead is worth avoiding.
  shared_cl <- .sptrends_shared_cluster(n_cores)
  if (!is.null(shared_cl)) {
    on.exit(parallel::stopCluster(shared_cl), add = TRUE)
  }

  n_steps <- 4
  step <- 0
  timing <- list()
  out <- list(prewhiten = NULL, trend = NULL, trend_summary_table = NULL,
              slope = NULL, fdr = NULL)

  series <- x
  step <- step + 1
  if (prewhiten_method != "none") {
    if (verbose) {
      message(sprintf("== Step %d/%d: prewhitening (%s) ==",
                       step, n_steps, prewhiten_method))
    }
    t0 <- Sys.time()
    pw_full_args <- utils::modifyList(
      list(x = x, method = prewhiten_method, report = report,
           verbose = verbose),
      prewhiten_args)
    pw_result <- do.call(prewhiten, pw_full_args)
    timing$prewhiten <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    out$prewhiten <- pw_result
    series <- pw_result$series
    # As in workflow_tst(): most prewhitening methods keep every layer
    # of x, but yue_pilon-style transforms (TFPW_Y here) lose the
    # first time step, so series can have one fewer layer than x --
    # names(x) would be the wrong length for it in that case.
    names(series) <- utils::tail(names(x), terra::nlyr(series))
  } else if (verbose) {
    message(sprintf(
      "== Step %d/%d: prewhitening skipped (prewhiten_method = \"none\") ==",
      step, n_steps))
  }

  step <- step + 1
  if (verbose) {
    message(sprintf("== Step %d/%d: trend test (%s) ==",
                     step, n_steps, trend_method))
  }
  t0 <- Sys.time()
  trend_full_args <- utils::modifyList(
    list(x = series, method = trend_method, report = report,
         verbose = verbose, n_cores = n_cores, shared_cluster = shared_cl),
    trend_args)
  trend_result <- do.call(trend_test, trend_full_args)
  timing$trend <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  out$trend <- trend_result$stats
  summary_alpha <- if (is.null(fdr_args$q)) 0.05 else fdr_args$q
  out$trend_summary_table <- trend_summary(trend_result$stats, path = NULL,
                                            alpha = summary_alpha,
                                            verbose = FALSE)

  step <- step + 1
  if (!is.null(slope_method)) {
    if (verbose) {
      message(sprintf("== Step %d/%d: slope estimation (%s) ==",
                       step, n_steps, slope_method))
    }
    t0 <- Sys.time()
    slope_full_args <- utils::modifyList(
      list(x = series, method = slope_method, report = report,
           verbose = verbose, n_cores = n_cores, shared_cluster = shared_cl),
      slope_args)
    slope_result <- do.call(slope_estimator, slope_full_args)
    timing$slope <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    out$slope <- slope_result$slope
  } else if (verbose) {
    message(sprintf(
      "== Step %d/%d: slope estimation skipped (slope_method = NULL) ==",
      step, n_steps))
  }

  step <- step + 1
  if (!is.null(fdr_method)) {
    if (verbose) {
      message(sprintf("== Step %d/%d: FDR correction (%s) ==",
                       step, n_steps, paste(fdr_method, collapse = ", ")))
    }
    t0 <- Sys.time()
    fdr_full_args <- utils::modifyList(
      list(p = out$trend$p, method = fdr_method, report = report,
           verbose = verbose),
      fdr_args)
    out$fdr <- do.call(fdr_correction, fdr_full_args)
    timing$fdr <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  } else if (verbose) {
    message(sprintf(
      "== Step %d/%d: FDR correction skipped (fdr_method = NULL) ==",
      step, n_steps))
  }

  out$timing <- timing
  class(out) <- c("workflow_trends", "sptrends")
  out
}
