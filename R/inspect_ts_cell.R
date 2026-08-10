#' Inspect a single cell's (or area's) raw time series interactively
#'
#' Click a cell (or draw a polygon) on **whatever map is currently
#' displayed** -- a trend map, a significance map, a raw data map, it does
#' not matter what it shows, only that it shares the same spatial extent
#' as `x`.
#'
#' **Why inspect a single cell?**: raster trend maps summarise thousands
#' of time series into one image. Inspecting an individual location
#' helps determine whether an apparently unusual pixel reflects a
#' genuine temporal pattern, an isolated outlier, or a behaviour
#' representative of its surrounding neighbourhood -- a question no
#' summary map, on its own, can answer.
#'
#' **What it shows**: the raw time series behind the clicked location,
#' with a single fitted line overlaid (Theil-Sen or OLS -- see
#' `slope_method`), together with its confidence interval.
#'
#' **Prewhitening comparison**: optionally, a second panel shows the
#' same location's *prewhitened* series side by side (see `prewhitened`
#' below), so you can see the effect of that step at exactly the
#' location you are looking at.
#'
#' **When to use it**: to get oriented early on (click around right
#' after loading data) and to investigate a specific cell that looks
#' anomalous on a map you have already produced.
#'
#' **Built from the package's own pieces**: this function is not a
#' separate implementation of its own -- the fitted line calls
#' [slope_estimator()]'s own estimator, the neighbourhood aggregation
#' follows [trend_test()]'s own logic, and the prewhitened panel reads
#' [prewhiten()]'s own output directly. It is, in that sense, less a
#' standalone utility than a visual demonstration of how the rest of
#' sptrends' pieces fit together, applied to one location at a time.
#'
#' **Function type:** **Interactive exploration function** -- combines
#' existing sptrends estimators and neighbourhood logic at one selected
#' location; it introduces no new inferential method.
#'
#' @section Typical use:
#' ```
#' raster time series + a displayed map with matching geometry
#'     |
#' inspect_ts_cell()
#'     |
#' selected cell or area -> temporal plot + fitted slope
#' ```
#' Optionally supply the complete [prewhiten()] result to compare the
#' raw and transformed series at the same location.
#'
#' @param x The full time series stack (a `terra::SpatRaster`, one layer
#'   per time step) -- this is always where the plotted raw series comes
#'   from, regardless of what is currently displayed on screen.
#' @param prewhitened Optional. The full list returned by
#'   [prewhiten()] (not just its `$series`) -- if supplied, a
#'   second panel shows the same location's prewhitened series, its own
#'   Theil-Sen fit and confidence interval, and whether that location was
#'   actually modified by prewhitening (many cells are not, if their own
#'   Durbin-Watson statistic never crossed the gating threshold; see
#'   `?prewhiten`). Default `NULL`: only the raw-data panel is
#'   drawn.
#' @param selection_type `"point"` (default): click a single cell.
#'   `"polygon"`: draw a polygon (left-click to add vertices, press
#'   `Esc` or right-click to finish).
#' @param neighbourhood Logical, only used when `selection_type =
#'   "point"`. If
#'   `TRUE` (default), the clicked cell and its queen (or rook)
#'   neighbours are combined -- see the "How each mode aggregates its
#'   series" section below -- before estimating anything, borrowing
#'   spatial context the same way [trend_test()] does for
#'   significance. If `FALSE`, only the clicked cell's own series is
#'   used. Ignored when `selection_type = "polygon"` (a polygon is
#'   already an explicit choice of area). Intended for exploratory
#'   visualisation,
#'   at this one location, rather than a substitute for formal
#'   inference -- the slope shown here does not carry the same
#'   significance guarantee CMK's own variance-adjusted statistic does.
#' @param connectivity `"queen"` (default, 8 neighbours) or `"rook"` (4
#'   neighbours), passed to `terra::adjacent()`. Only relevant when
#'   `neighbourhood = TRUE`.
#' @param conf_level Confidence level for the fitted slope's confidence
#'   interval, reported in the legend as text (not drawn as a shaded
#'   band -- see "Confidence interval" below). Default `0.95`. Ignored
#'   when `slope_method = "RM"` or `compare_slopes = TRUE` (see both
#'   below).
#' @param t Finite numeric vector of unique, strictly increasing time
#'   points, one per layer. Defaults to `1:nlyr(x)`.
#' @param show_neighbours Logical, only used when `neighbourhood = TRUE`
#'   and `selection_type = "point"`. Answers, at a glance, a question any user
#'   looking at an unusual pixel asks: **is this cell's trend
#'   representative of its neighbourhood, or an outlier the aggregation
#'   is smoothing over?** If `TRUE`, draws a second figure after the
#'   main panel(s): a small-multiples grid, one mini-panel per cell
#'   actually aggregated (the clicked cell, highlighted, plus each of
#'   its individual queen/rook neighbours), each with its own raw
#'   series and fitted line (same `slope_method` as the main panel) --
#'   not the median-aggregated series
#'   the main panel shows, but each contributing cell on its own.
#'   Ignored (with a message) if the clicked cell has no
#'   neighbours with complete data (e.g. a corner cell with all-NA
#'   neighbours), since there would be nothing to compare against.
#'   Default `FALSE`.
#' @param slope_method Which estimator [slope_estimator()] uses for the
#'   fitted line(s) shown here. `"TS"` (default): robust to
#'   outliers, matches this package's own default trend workflow.
#'   `"OLS"`: ordinary least squares -- faster, but sensitive to
#'   outliers the way a single anomalous time step in the plot above
#'   can be. `"RM"`: Siegel's repeated median -- more robust than
#'   `"TS"`, but no confidence interval is shown (none implemented for
#'   it in this package). Ignored when `compare_slopes = TRUE`. This
#'   only affects this function's own quick single-cell
#'   fit; it has no bearing on which method a
#'   `workflow_tst()`/`workflow_rta()` run
#'   itself used.
#' @param compare_slopes Logical. If `TRUE`, ignores `slope_method` and
#'   draws all three estimators (`"TS"`, `"OLS"`, `"RM"`) as separate
#'   lines on the same panel(s), point estimates only -- no confidence
#'   intervals in this mode, since the three use genuinely different
#'   inferential frameworks (or none, for `"RM"`) and no single
#'   interval formula would honestly apply to all three. Default
#'   `FALSE`.
#' @param verbose Logical. If `TRUE` (default), print the click/draw
#'   instructions and the `show_neighbours` no-effect note; matches the
#'   `verbose` convention used throughout this package, though the
#'   interactive click/draw step itself still happens either way --
#'   this only silences the accompanying messages and elapsed time, not
#'   the interaction.
#' @param ... Ignored.
#'
#' @return Returns invisibly, a list with `cell` (the clicked/representative
#'   cell number) and `raw` (a list with `series`, `slope`, `ci_lower`,
#'   `ci_upper`, `conf_level` for the raw-data panel). If `prewhitened`
#'   was supplied, also `prewhitened` (the same fields for that panel)
#'   and `n_modified`/`n_total` (how many of the aggregated cells were
#'   actually modified by prewhitening). If `show_neighbours = TRUE` and
#'   at least one neighbour had a complete series, also `neighbours`: a
#'   list, one element per plotted cell (including the clicked one),
#'   each with `cell`, `slope`, and `is_centre`.
#'
#' @section Methodological details:
#' **Methods and method selection**
#'
#' `slope_method` selects Theil-Sen, OLS, or repeated median for the
#' exploratory fit. `compare_slopes = TRUE` displays all three point
#' estimates together; confidence intervals are then omitted because no
#' single inferential framework applies to all three estimators.
#'
#' **How each mode aggregates its series**
#'
#' All three modes -- point-only, point-with-neighbourhood, and polygon
#' -- follow the *same* two-step procedure, in the *same* order, so a
#' single confidence interval formula applies rigorously to all of them:
#' first take the per-time-step **median of raw values** across whichever
#' cells are included (just the one clicked cell, the clicked cell plus
#' its neighbours, or every cell inside the drawn polygon), producing one
#' aggregated series; *then* fit a single Theil-Sen slope to that series.
#' The same aggregation (same cells) is applied to the prewhitened series
#' too, when `prewhitened` is supplied, so the two panels are directly
#' comparable. Aggregating values first and estimating second is what
#' makes the Sen/Gilbert confidence interval below valid -- it is defined
#' for the pairwise slopes of a single series, not for a median of
#' several already-computed slope estimates (the reverse order), which
#' has no standard interval formula. This also means `neighbourhood =
#' TRUE`'s result is **not** the same number as `slope_estimator(x,
#' smooth_neighbourhood = TRUE)` at that same cell in a full-map run --
#' that other mechanism deliberately aggregates in the opposite order
#' (median of independently-computed per-cell slopes); see its own
#' documentation for why. The two exist for different purposes and are
#' not meant to match numerically.
#'
#' **Statistical assumptions and confidence intervals**
#'
#' **What it represents**: uncertainty in the *rate of change* itself,
#' not a prediction interval around the raw data points. Reported as
#' text in the legend (`ci_lower`/`ci_upper` in the returned value too),
#' not drawn as a shaded band on the plot -- only the single fitted
#' line and its slope/intercept are drawn. Not shown when
#' `slope_method = "RM"` or `compare_slopes = TRUE`.
#'
#' **How it is computed** (the standard rank-based method for the
#' Theil-Sen slope; Sen, 1968; Gilbert, 1987, pp. 217-219): given the
#' \eqn{N} pairwise slopes of the (already-aggregated) series, sorted,
#' the interval is \eqn{[\hat\beta_{(M1)}, \hat\beta_{(M2+1)}]}, where
#' \eqn{M1 = (N - C_\alpha)/2}, \eqn{M2 = (N + C_\alpha)/2}, and
#' \eqn{C_\alpha = z_{1-\alpha/2}\sqrt{Var(S)}} (the same Mann-Kendall
#' variance of `S` used throughout this package, tie-corrected). `M1`/
#' `M2` are rounded to the nearest integer rather than interpolated
#' between adjacent ranked slopes, unlike Gilbert's own recommendation --
#' a deliberate simplification, documented here rather than left silent.
#' This formula applies when `slope_method = "TS"` (the default);
#' for `slope_method = "OLS"`, the interval instead uses the standard
#' parametric simple-linear-regression interval (which assumes
#' normally-distributed, homoscedastic residuals) -- a genuinely
#' different formula, not the same computation reused, since it rests
#' on different assumptions matching its own point estimate.
#'
#' **Limitations**
#'
#' This is an exploratory, interactive view of selected locations, not a
#' raster-wide significance procedure. Neighbourhood or polygon medians
#' change the series being fitted, and repeated-median or multi-estimator
#' displays do not include confidence intervals.
#'
#' **Quality assurance**
#'
#' Tests cover point and polygon selection, neighbourhood aggregation,
#' raw/prewhitened comparisons, all slope choices, confidence intervals,
#' time metadata, incomplete series, interactive selection failures and
#' silent operation. Core numerical results are compared with direct
#' calculations and [slope_estimator()]. See `?sptrends` for the common
#' release-check protocol.
#'
#' @family reporting functions
#' @references
#' Confidence interval method:
#' - Sen, P.K. (1968) Estimates of the regression coefficient based on
#'   Kendall's tau. Journal of the American Statistical Association, 63,
#'   1379-1389. \doi{10.1080/01621459.1968.10480934}
#' - Gilbert, R.O. (1987) Statistical Methods for Environmental Pollution
#'   Monitoring. Van Nostrand Reinhold, New York, pp. 217-219.
#'
#' Origin of the Var(S) formula reused for the confidence interval (the
#' same one `trend_test()` builds on for CMK's adjusted
#' variance; see its own references for that extension):
#' - Mann, H.B. (1945) Nonparametric tests against trend. Econometrica,
#'   13(3), 245-259. \doi{10.2307/1907187}
#' - Kendall, M.G. (1975) Rank Correlation Methods (4th edn). Charles
#'   Griffin, London. No DOI available (pre-DOI-era publication).
#'
#' @seealso [prewhiten()] for transformed time series, [trend_test()] for
#'   trend significance, and [slope_estimator()] for raster-wide magnitude
#'   estimation.
#'
#' @examples
#' \dontrun{
#' # Interactive -- run this yourself, requires clicking on a plot.
#'
#' # Annual mean NDVI from the bundled environmental dataset.
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' terra::plot(r[[1]])         # any single-layer map with the same extent as r
#'
#' # That cell's own raw series, nothing borrowed from its surroundings.
#' inspect_ts_cell(r, neighbourhood = FALSE)
#'
#' # Compare against the function's own default -- the same cell combined
#' # with its queen neighbourhood -- to see how much borrowing spatial
#' # context changes the fit.
#' inspect_ts_cell(r)
#'
#' inspect_ts_cell(r, selection_type = "polygon")       # draw an area instead
#'
#' # Raw vs. prewhitened, side by side, at the clicked location
#' pw <- prewhiten(r, report = FALSE, verbose = FALSE)
#' inspect_ts_cell(r, prewhitened = pw)
#'
#' # Is the clicked cell's trend representative of its neighbourhood,
#' # or an outlier the median aggregation is smoothing over? Draws a
#' # second figure: one mini-panel per neighbour, plus the clicked cell
#' # itself (highlighted), each with its own Theil-Sen fit.
#' inspect_ts_cell(r, show_neighbours = TRUE)
#'
#' # Faster, but not robust to the effect a single anomalous time step
#' # can have on the fitted line -- see the "slope_method" argument.
#' inspect_ts_cell(r, slope_method = "OLS")
#'
#' # Siegel's repeated median -- more robust than Theil-Sen, but no
#' # confidence interval is shown for it (none is implemented here).
#' inspect_ts_cell(r, slope_method = "RM")
#'
#' # All three estimators at once, as three lines with no intervals.
#' inspect_ts_cell(r, compare_slopes = TRUE)
#' }
#'
#' @export
inspect_ts_cell <- function(x, prewhitened = NULL,
                             selection_type = c("point", "polygon"),
                             neighbourhood = TRUE,
                             connectivity = c("queen", "rook"),
                             conf_level = 0.95, t = NULL,
                             show_neighbours = FALSE,
                             slope_method = c("TS", "OLS", "RM"),
                             compare_slopes = FALSE,
                             verbose = TRUE, ...) {
  finish_timer <- .sptrends_elapsed_timer("inspect_ts_cell()", verbose)
  on.exit(finish_timer(), add = TRUE)
  selection_type <- match.arg(selection_type)
  connectivity <- match.arg(connectivity)
  slope_method <- match.arg(slope_method)
  if (!inherits(x, "SpatRaster")) stop("'x' must be a terra SpatRaster.")
  conf_level <- .validate_probability(conf_level, "conf_level")
  if (!is.null(prewhitened) &&
        (!is.list(prewhitened) || is.null(prewhitened$series) ||
           is.null(prewhitened$diagnostics))) {
    stop("'prewhitened' must be the full list returned by ",
         "prewhiten(), not just its $series.")
  }

  if (selection_type == "point") {
    if (verbose) message("Click a cell on the currently displayed map...")
    pt <- terra::click(x[[1]], n = 1, xy = TRUE, show = FALSE)
    if (is.null(pt) || nrow(pt) == 0) stop("No click registered.")
    cell <- terra::cellFromXY(x, as.matrix(pt[, c("x", "y")]))
    polygon <- NULL
  } else {
    if (verbose) {
      message("Draw a polygon on the currently displayed map (left-click ",
              "to add vertices, Esc or right-click to finish)...")
    }
    polygon <- terra::draw("polygon")
    cent <- terra::centroids(polygon)
    cell <- terra::cellFromXY(x, terra::crds(cent))
  }
  if (is.na(cell)) stop("Selected location falls outside 'x'.")

  .inspect_ts_cell_core(x, cell = cell, t = t, neighbourhood = neighbourhood,
                         connectivity = connectivity, conf_level = conf_level,
                         polygon = polygon, prewhitened = prewhitened,
                         show_neighbours = show_neighbours,
                         slope_method = slope_method,
                         compare_slopes = compare_slopes, verbose = verbose)
}

#' @noRd
.sen_confidence_interval <- function(series, t, conf_level = 0.95) {
  n <- length(series)
  pairs <- utils::combn(n, 2)
  slopes <- sort((series[pairs[2, ]] - series[pairs[1, ]]) /
                    (t[pairs[2, ]] - t[pairs[1, ]]))
  N <- length(slopes)

  tb <- table(series)
  tie_correction <- sum(tb * (tb - 1) * (2 * tb + 5))
  var_s <- (n * (n - 1) * (2 * n + 5) - tie_correction) / 18

  z <- stats::qnorm(1 - (1 - conf_level) / 2)
  c_alpha <- z * sqrt(var_s)

  m1 <- round((N - c_alpha) / 2)
  m2 <- round((N + c_alpha) / 2)
  lo_idx <- max(1, m1)
  hi_idx <- min(N, m2 + 1)

  list(lower = slopes[lo_idx], upper = slopes[hi_idx], conf_level = conf_level)
}

#' @noRd
.ols_confidence_interval <- function(series, t, slope, intercept,
                                      conf_level = 0.95) {
  # The standard textbook CI for a simple linear regression slope --
  # a genuinely different formula from .sen_confidence_interval()
  # above, not just the same computation under a different name: this
  # one is parametric (assumes normally-distributed, homoscedastic
  # residuals), Sen's is rank-based and distribution-free. Using Sen's
  # formula for an OLS slope, or vice versa, would silently mismatch
  # the point estimate to a CI built on different assumptions.
  n <- length(series)
  fitted <- intercept + slope * t
  resid <- series - fitted
  s2 <- sum(resid^2) / (n - 2)
  se_slope <- sqrt(s2 / sum((t - mean(t))^2))
  tcrit <- stats::qt(1 - (1 - conf_level) / 2, df = n - 2)
  list(lower = slope - tcrit * se_slope, upper = slope + tcrit * se_slope,
       conf_level = conf_level)
}

#' @noRd
.resolve_aggregation_cells <- function(x, cell, neighbourhood, connectivity,
                                        polygon) {
  if (!is.null(polygon)) return(NULL)  # polygon extraction handled separately
  if (isTRUE(neighbourhood)) {
    adj <- terra::adjacent(x, cells = cell, directions = connectivity,
                            pairs = FALSE)
    neighbour_cells <- as.vector(adj)
    neighbour_cells <- neighbour_cells[!is.na(neighbour_cells)]
    c(cell, neighbour_cells)
  } else {
    cell
  }
}

#' @noRd
.aggregate_series <- function(x, cells, polygon) {
  if (!is.null(polygon)) {
    ext <- terra::extract(x, polygon)[, -1, drop = FALSE]
    apply(ext, 2, stats::median, na.rm = TRUE)
  } else if (length(cells) > 1) {
    vals_mat <- terra::values(x)[cells, , drop = FALSE]
    apply(vals_mat, 2, stats::median, na.rm = TRUE)
  } else {
    as.numeric(x[cells])
  }
}

#' @noRd
.fit_slope_with_ci <- function(series, t, conf_level,
                                method = c("TS", "OLS", "RM")) {
  method <- match.arg(method)
  n <- length(series)
  mini <- terra::rast(nrows = 1, ncols = 1, xmin = 0, xmax = 1, ymin = 0,
                       ymax = 1, nlyrs = n)
  terra::values(mini) <- matrix(series, nrow = 1)
  slope_result <- slope_estimator(mini, method = method, t = t,
                                   smooth_neighbourhood = FALSE,
                                   report = FALSE, verbose = FALSE)
  slope <- terra::values(slope_result$slope)[1]

  if (method == "TS") {
    ci <- .sen_confidence_interval(series, t, conf_level = conf_level)
    # Sen's method: the intercept is the median of series - slope*t
    # (matches the rank-based slope above), not the OLS mean-based one.
    intercept <- stats::median(series - slope * t)
  } else if (method == "OLS") {
    intercept <- mean(series) - slope * mean(t)
    ci <- .ols_confidence_interval(series, t, slope, intercept,
                                    conf_level = conf_level)
  } else {
    # "RM" (Siegel's repeated median) has no confidence-interval
    # formula implemented in this package -- see ?slope_estimator's
    # own "Implementation notes" section. Its own intercept (already
    # computed by slope_estimator() itself, using the "direct"
    # convention verified against robslopes::RepeatedMedian()) is
    # reused directly rather than recomputed here, and no CI is
    # fabricated for it: ci_lower/ci_upper are NA, and
    # .plot_ts_panel() skips that legend line entirely when NA.
    intercept <- terra::values(slope_result$intercept)[1]
    ci <- list(lower = NA_real_, upper = NA_real_, conf_level = conf_level)
  }
  # t_center/y_center: the point on the fitted line at the middle of the
  # observed t range, not at t = 0 (the intercept's own anchor point,
  # which usually falls outside the data). Kept for the CI's own
  # mathematical property this represents -- the interval on the slope
  # pinches to a single point here and widens symmetrically away from
  # it -- even though it is no longer drawn as a shaded band (see
  # .plot_ts_panel(): the legend shows ci_lower/ci_upper as text only).
  t_center <- mean(t)
  y_center <- intercept + slope * t_center
  list(series = series, slope = slope, intercept = intercept,
       t_center = t_center, y_center = y_center,
       ci_lower = ci$lower, ci_upper = ci$upper, conf_level = conf_level,
       method = method)
}

#' @noRd
.fit_all_slopes <- function(series, t) {
  # No confidence intervals here, deliberately: TS/OLS/RM use
  # genuinely different inferential frameworks (rank-based,
  # parametric, and none implemented at all for RM respectively -- see
  # .fit_slope_with_ci()'s own comment on this), so a shared interval
  # across all three would either mismatch the point estimate to the
  # wrong assumptions or fabricate one that does not exist. This
  # function exists only to compare point estimates side by side.
  #
  # Each method's own intercept convention is used, not a single
  # shared formula: TS uses the median-based convention matching its
  # own rank-based slope (median(series - slope*t)); OLS uses the
  # mean-based convention matching its own least-squares slope; RM
  # reuses slope_estimator()'s own "direct" repeated-median intercept
  # directly, rather than recomputing it with either of the other two
  # formulas, which would not match the line slope_estimator() itself
  # reports for "RM".
  n <- length(series)
  mini <- terra::rast(nrows = 1, ncols = 1, xmin = 0, xmax = 1, ymin = 0,
                       ymax = 1, nlyrs = n)
  terra::values(mini) <- matrix(series, nrow = 1)

  ts_result <- slope_estimator(mini, method = "TS", t = t,
                                smooth_neighbourhood = FALSE,
                                report = FALSE, verbose = FALSE)
  ols_result <- slope_estimator(mini, method = "OLS", t = t,
                                 smooth_neighbourhood = FALSE,
                                 report = FALSE, verbose = FALSE)
  rm_result <- slope_estimator(mini, method = "RM", t = t,
                                smooth_neighbourhood = FALSE,
                                report = FALSE, verbose = FALSE)

  slopes <- c(TS = terra::values(ts_result$slope)[1],
              OLS = terra::values(ols_result$slope)[1],
              RM = terra::values(rm_result$slope)[1])

  intercepts <- c(
    TS = stats::median(series - slopes[["TS"]] * t),
    OLS = mean(series) - slopes[["OLS"]] * mean(t),
    RM = terra::values(rm_result$intercept)[1]
  )

  list(slope = slopes, intercept = intercepts)
}

#' @noRd
.plot_ts_panel_compare <- function(t, series, fits, label,
                                    xlab = "Time step") {
  graphics::plot(t, series, type = "b", pch = 19, col = "grey30",
                 xlab = xlab, ylab = "Value", main = label)
  cols <- c(TS = "steelblue", OLS = "firebrick", RM = "darkgreen")
  labels <- c(TS = "Theil-Sen", OLS = "OLS", RM = "Siegel (RM)")
  slopes <- fits$slope
  intercepts <- fits$intercept
  for (m in names(slopes)) {
    graphics::abline(a = intercepts[[m]], b = slopes[[m]],
                      col = cols[[m]], lwd = 2)
  }
  graphics::legend(
    "topleft",
    legend = sprintf("%s slope = %.4g", labels[names(slopes)], slopes),
    col = cols[names(slopes)], lwd = 2, bty = "n")
}

#' @noRd
.plot_ts_panel <- function(t, fit, label, xlab = "Time step") {
  graphics::plot(t, fit$series, type = "b", pch = 19, col = "grey30",
                 xlab = xlab, ylab = "Value", main = label)
  graphics::abline(a = fit$intercept, b = fit$slope, col = "steelblue", lwd = 2)
  method_label <- switch(fit$method, OLS = "OLS", RM = "Siegel (RM)",
                          "Theil-Sen")
  legend_text <- sprintf("%s slope = %.4g", method_label, fit$slope)
  legend_col <- "steelblue"
  legend_lwd <- 2
  if (!is.na(fit$ci_lower) && !is.na(fit$ci_upper)) {
    # "RM" has no confidence-interval formula implemented in this
    # package (see ?slope_estimator's own "Implementation notes") --
    # fit$ci_lower/fit$ci_upper are NA for that method specifically,
    # and this line is skipped rather than showing a fabricated or
    # missing-looking interval.
    legend_text <- c(legend_text,
                      sprintf("%.0f%% slope CI: [%.4g, %.4g]",
                              100 * fit$conf_level,
                              fit$ci_lower, fit$ci_upper))
    legend_col <- c(legend_col, NA)
    legend_lwd <- c(legend_lwd, NA)
  }
  graphics::legend(
    "topleft",
    legend = legend_text,
    col = legend_col, lwd = legend_lwd, bty = "n")
}

#' @noRd
.plot_ts_small_multiples <- function(x, t, cells, centre_cell,
                                      slope_method = "TS") {
  # One own-series-only fit per cell, no aggregation with its neighbours
  # -- the whole point is to see each contributing cell on its own,
  # unlike the median-aggregated panel(s) drawn before this one.
  fits <- vector("list", length(cells))
  for (i in seq_along(cells)) {
    series <- as.numeric(x[cells[i]])
    fits[[i]] <- if (anyNA(series)) {
      NULL
    } else {
      tryCatch(.fit_slope_with_ci(series, t, conf_level = 0.95,
                                   method = slope_method),
                error = function(e) NULL)
    }
  }
  names(fits) <- as.character(cells)

  n <- length(cells)
  ncol_grid <- ceiling(sqrt(n))
  nrow_grid <- ceiling(n / ncol_grid)
  old_par <- graphics::par(mfrow = c(nrow_grid, ncol_grid),
                            mar = c(2, 2, 2, 1))
  on.exit(graphics::par(old_par), add = TRUE)

  for (i in seq_along(cells)) {
    is_centre <- cells[i] == centre_cell
    title_col <- if (is_centre) "firebrick" else "grey20"
    if (is.null(fits[[i]])) {
      graphics::plot.new()
      graphics::title(main = sprintf("Cell %d (no data)", cells[i]),
                       col.main = title_col,
                       font.main = if (is_centre) 2 else 1)
      next
    }
    fit <- fits[[i]]
    graphics::plot(t, fit$series, type = "l", col = "grey50", xlab = "",
                   ylab = "",
                   main = sprintf("Cell %d%s", cells[i],
                                   if (is_centre) " (clicked)" else ""),
                   col.main = title_col,
                   font.main = if (is_centre) 2 else 1)
    graphics::abline(a = fit$intercept, b = fit$slope,
                      col = if (is_centre) "firebrick" else "steelblue",
                      lwd = 2)
    if (is_centre) graphics::box(col = "firebrick", lwd = 2)
  }

  invisible(lapply(seq_along(cells), function(i) {
    list(cell = cells[i],
         slope = if (is.null(fits[[i]])) NA_real_ else fits[[i]]$slope,
         is_centre = cells[i] == centre_cell)
  }))
}

#' @noRd
.inspect_ts_cell_core <- function(x, cell, t = NULL, neighbourhood = TRUE,
                                   connectivity = c("queen", "rook"),
                                   conf_level = 0.95, label = NULL,
                                   polygon = NULL, prewhitened = NULL,
                                   show_neighbours = FALSE,
                                   slope_method = c("TS", "OLS", "RM"),
                                   compare_slopes = FALSE,
                                   verbose = TRUE) {
  connectivity <- match.arg(connectivity)
  slope_method <- match.arg(slope_method)
  n <- terra::nlyr(x)
  # Prefer x's own real detected time (e.g. the years read_ordered_stack()
  # stores via terra::time(), rather than discarding them) over a
  # generic 1:n time-step index, when the caller hasn't supplied t
  # explicitly and x actually has valid, usable time metadata.
  used_real_time <- FALSE
  if (is.null(t)) {
    if (isTRUE(terra::has.time(x))) {
      real_t <- suppressWarnings(as.numeric(terra::time(x)))
      if (length(real_t) == n && !anyNA(real_t)) {
        t <- real_t
        used_real_time <- TRUE
      }
    }
    if (is.null(t)) t <- seq_len(n)
  }
  t <- .validate_time_axis(t, n)

  cells <- .resolve_aggregation_cells(x, cell, neighbourhood, connectivity,
                                       polygon)
  raw_series <- .aggregate_series(x, cells, polygon)
  if (anyNA(raw_series)) {
    stop("Selected location has no complete aggregated time series.")
  }

  if (is.null(label)) {
    label <- if (!is.null(polygon)) {
      "Polygon (median)"
    } else if (isTRUE(neighbourhood)) {
      sprintf("Cell %d + %d neighbours (median)", cell, length(cells) - 1)
    } else {
      sprintf("Cell %d (own series)", cell)
    }
  }

  raw_fit <- .fit_slope_with_ci(raw_series, t, conf_level,
                                 method = slope_method)

  draw_panel <- function(t_vals, series_vals, fit, panel_label, xlab) {
    if (isTRUE(compare_slopes)) {
      fits <- .fit_all_slopes(series_vals, t_vals)
      .plot_ts_panel_compare(t_vals, series_vals, fits, panel_label,
                              xlab = xlab)
      data.frame(method = names(fits$slope),
                 slope = unname(fits$slope),
                 intercept = unname(fits$intercept))
    } else {
      .plot_ts_panel(t_vals, fit, panel_label, xlab = xlab)
      NULL
    }
  }

  draw_neighbours <- isTRUE(show_neighbours) && isTRUE(neighbourhood) &&
    is.null(polygon) && length(cells) > 1
  if (!draw_neighbours && isTRUE(show_neighbours) && verbose) {
    message("show_neighbours = TRUE has no effect here: it needs ",
            "neighbourhood = TRUE, selection_type = \"point\", and at ",
            "least one neighbour cell.")
  }
  draw_neighbours_grid <- function() {
    if (!draw_neighbours) return(NULL)
    # A separate window, not a sequential draw on the same device --
    # otherwise this grid would simply overwrite the main panel(s)
    # drawn just before it, and the two are meant to be seen together.
    grDevices::dev.new()
    .plot_ts_small_multiples(x, t, cells, cell, slope_method = slope_method)
  }

  if (is.null(prewhitened)) {
    graphics::par(mfrow = c(1, 1))
    time_xlab <- if (used_real_time) "Year" else "Time step"
    slope_comparison <- draw_panel(t, raw_series, raw_fit,
                                    sprintf("Raw -- %s", label), time_xlab)
    neighbours_result <- draw_neighbours_grid()
    out <- list(
      cell = cell,
      raw = raw_fit[c("series", "slope", "ci_lower", "ci_upper", "conf_level")]
    )
    if (!is.null(slope_comparison)) out$slope_comparison <- slope_comparison
    if (!is.null(neighbours_result)) out$neighbours <- neighbours_result
    return(invisible(out))
  }

  pw_cells <- if (!is.null(polygon)) NULL else cells
  pw_series <- .aggregate_series(prewhitened$series, pw_cells, polygon)
  if (anyNA(pw_series)) {
    stop("Selected location has no complete prewhitened time series.")
  }
  # TFPW_WS keeps every time step (length(pw_series) == length(t));
  # TFPW_Y's classic transform loses the first one, so
  # length(pw_series) == length(t) - 1. Passing the full, longer t
  # alongside a shorter pw_series would silently misalign every time
  # value by one step (t[1:(n-1)] instead of the correct t[2:n]) rather
  # than error outright -- align explicitly here instead of assuming
  # they always match.
  pw_t <- utils::tail(t, length(pw_series))
  pw_fit <- .fit_slope_with_ci(pw_series, pw_t, conf_level,
                                method = slope_method)

  if (identical(prewhitened$method, "TFPW_Y")) {
    # No DW gate for this method -- every valid cell is prewhitened
    # unconditionally, so "modified" here just means "had a valid Rho
    # estimate", not a separate diagnostic field like Modified.
    rho_vals <- if (!is.null(polygon)) {
      terra::extract(prewhitened$diagnostics$Rho, polygon)[, -1]
    } else {
      terra::values(prewhitened$diagnostics$Rho)[cells]
    }
    n_modified <- sum(!is.na(rho_vals))
    n_total <- n_modified
  } else {
    modified_vals <- if (!is.null(polygon)) {
      terra::extract(prewhitened$diagnostics$Modified, polygon)[, -1]
    } else {
      terra::values(prewhitened$diagnostics$Modified)[cells]
    }
    n_modified <- sum(modified_vals == 1, na.rm = TRUE)
    n_total <- sum(!is.na(modified_vals))
  }
  pw_label <- sprintf("Prewhitened -- %s (%d/%d cells modified)",
                       label, n_modified, n_total)

  graphics::par(mfrow = c(1, 2))
  time_xlab <- if (used_real_time) "Year" else "Time step"
  raw_comparison <- draw_panel(t, raw_series, raw_fit,
                                sprintf("Raw -- %s", label), time_xlab)
  pw_comparison <- draw_panel(pw_t, pw_series, pw_fit, pw_label, time_xlab)
  neighbours_result <- draw_neighbours_grid()

  out <- list(
    cell = cell,
    raw = raw_fit[c("series", "slope", "ci_lower", "ci_upper", "conf_level")],
    prewhitened = pw_fit[c("series", "slope", "ci_lower", "ci_upper",
                            "conf_level")],
    n_modified = n_modified, n_total = n_total
  )
  if (!is.null(raw_comparison)) out$slope_comparison <- raw_comparison
  if (!is.null(pw_comparison)) out$pw_slope_comparison <- pw_comparison
  if (!is.null(neighbours_result)) out$neighbours <- neighbours_result
  invisible(out)
}
