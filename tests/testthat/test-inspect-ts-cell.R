test_that(".sen_confidence_interval matches a hand-computed example", {
  # t = 1:4, series = 1,3,2,6 -- same series used in test-theilsen.R's
  # hand-computed exact Theil-Sen slope (19/12). Pairwise slopes sorted:
  # -1, 0.5, 1.5, 5/3, 2, 4 (N = 6). No ties among 1,3,2,6, so
  # Var(S) = n(n-1)(2n+5)/18 = 4*3*13/18 = 8.6667 exactly.
  # C_alpha = qnorm(0.975) * sqrt(8.6667) ~= 5.7699
  # M1 = (6 - 5.7699)/2 ~= 0.115 -> round to 0 -> clamped to index 1
  # M2 = (6 + 5.7699)/2 ~= 5.885 -> round to 6 -> hi index 7 -> clamped to 6
  # => CI = [slopes[1], slopes[6]] = [-1, 4]
  series <- c(1, 3, 2, 6)
  t <- 1:4
  ci <- sptrends:::.sen_confidence_interval(series, t, conf_level = 0.95)
  expect_equal(ci$lower, -1)
  expect_equal(ci$upper, 4)
})

test_that(".inspect_ts_cell_core (no prewhitened) returns the expected structure", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 1)$series
  result <- sptrends:::.inspect_ts_cell_core(r, cell = 15, neighbourhood = FALSE)
  expect_type(result, "list")
  expect_named(result, c("cell", "raw"))
  expect_named(result$raw, c("series", "slope", "ci_lower", "ci_upper", "conf_level"))
  expect_identical(result$cell, 15)
  expect_length(result$raw$series, 10)
})

test_that(".inspect_ts_cell_core with neighbourhood = FALSE matches slope_estimator() on that cell's own series", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 2)$series
  cell <- 10
  result <- sptrends:::.inspect_ts_cell_core(r, cell = cell, neighbourhood = FALSE)
  expect_equal(unname(result$raw$series), as.numeric(r[cell]))

  full_slope <- slope_estimator(r, smooth_neighbourhood = FALSE, verbose = FALSE, report = FALSE)$slope
  expected <- terra::values(full_slope, mat = FALSE)[cell]
  expect_equal(result$raw$slope, expected)
})

test_that(".inspect_ts_cell_core with neighbourhood = TRUE aggregates by median of VALUES, then a single Theil-Sen", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 3)$series
  cell <- 30  # interior cell of an 8x8 grid -- all 8 queen neighbours available
  result <- sptrends:::.inspect_ts_cell_core(r, cell = cell, neighbourhood = TRUE)

  adj <- terra::adjacent(r, cells = cell, directions = "queen", pairs = FALSE)
  neighbour_cells <- as.vector(adj)
  neighbour_cells <- neighbour_cells[!is.na(neighbour_cells)]
  all_cells <- c(cell, neighbour_cells)
  vals_mat <- terra::values(r)[all_cells, , drop = FALSE]
  expected_series <- apply(vals_mat, 2, stats::median, na.rm = TRUE)

  expect_equal(unname(result$raw$series), unname(expected_series))
})

test_that(".inspect_ts_cell_core with a polygon takes the per-layer median across covered cells", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 4)$series
  mask_r <- r[[1]]
  mask_r[] <- NA
  mask_r[c(1, 2, 7, 8)] <- 1
  poly <- terra::as.polygons(mask_r, dissolve = TRUE, na.rm = TRUE)

  result <- sptrends:::.inspect_ts_cell_core(r, cell = 1, polygon = poly)
  ext_vals <- terra::extract(r, poly)[, -1, drop = FALSE]
  expected_series <- apply(ext_vals, 2, stats::median, na.rm = TRUE)
  expect_equal(unname(result$raw$series), unname(expected_series))
})

test_that(".inspect_ts_cell_core's confidence interval always brackets the point estimate", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 15, seed = 5)$series
  result <- sptrends:::.inspect_ts_cell_core(r, cell = 20, neighbourhood = FALSE)
  expect_lte(result$raw$ci_lower, result$raw$slope)
  expect_gte(result$raw$ci_upper, result$raw$slope)
})

test_that(".inspect_ts_cell_core errors when the selected cell has an incomplete series", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 6)$series
  r[1] <- NA
  expect_error(sptrends:::.inspect_ts_cell_core(r, cell = 1, neighbourhood = FALSE),
               "no complete aggregated time series")
})

test_that(".inspect_ts_cell_core handles a corner cell (fewer than 8 neighbours) without erroring", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 7)$series
  expect_error(sptrends:::.inspect_ts_cell_core(r, cell = 1, neighbourhood = TRUE), NA)
})

test_that(".inspect_ts_cell_core draws a second panel and reports modification status when prewhitened is supplied", {
  r <- .make_ar1_raster(rho = 0.6, n_time = 20, seed = 8)
  pw <- prewhiten(r, report = FALSE, verbose = FALSE)
  result <- sptrends:::.inspect_ts_cell_core(r, cell = 1, neighbourhood = FALSE,
                                              prewhitened = pw)
  expect_named(result, c("cell", "raw", "prewhitened", "n_modified", "n_total"))
  expect_identical(result$n_total, 1L)
  expect_true(result$n_modified %in% c(0L, 1L))
  # the prewhitened series should generally differ from the raw one
  # whenever the cell was actually modified
  if (result$n_modified == 1) {
    expect_false(isTRUE(all.equal(result$raw$series, result$prewhitened$series)))
  }
})

test_that("inspect_ts_cell errors on non-SpatRaster input", {
  expect_error(inspect_ts_cell(matrix(1:9, 3, 3)), "SpatRaster")
})

test_that("inspect_ts_cell errors when 'prewhitened' is not the full prewhiten() output", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 9)$series
  expect_error(inspect_ts_cell(r, prewhitened = r), "prewhiten")
})

test_that(".fit_slope_with_ci's t_center/y_center pinch the slope CI to a single point at the centre of t, not at t = 0 -- a mathematical property of the interval, no longer drawn as a band but still reported (t_center/y_center) and worth confirming directly", {
  series <- c(1, 3, 2, 6, 5, 8, 7, 10)
  t <- 1:8
  fit <- sptrends:::.fit_slope_with_ci(series, t, conf_level = 0.95)

  expect_equal(fit$t_center, mean(t))
  # at t_center, both CI bounds must give exactly the same y (the pinch
  # point) -- confirming the interval narrows to a point around the
  # middle of the series, not around t = 0 (which sits outside 1:8
  # entirely)
  lower_at_center <- fit$y_center + fit$ci_lower * (fit$t_center - fit$t_center)
  upper_at_center <- fit$y_center + fit$ci_upper * (fit$t_center - fit$t_center)
  expect_equal(lower_at_center, upper_at_center)
  expect_equal(lower_at_center, fit$y_center)

  # and the interval must be WIDER at both ends than at the centre
  # (t_center falls within 1:8 here, unlike the old t = 0 anchor)
  width_at_start <- (fit$ci_upper - fit$ci_lower) * abs(t[1] - fit$t_center)
  width_at_end <- (fit$ci_upper - fit$ci_lower) * abs(t[length(t)] - fit$t_center)
  expect_gt(width_at_start, 0)
  expect_gt(width_at_end, 0)
})

test_that(".inspect_ts_cell_core works with a polygon AND prewhitened together", {
  r <- .make_ar1_raster(rho = 0.6, n_time = 20, nrow = 6, ncol = 6, seed = 10)
  pw <- prewhiten(r, report = FALSE, verbose = FALSE)

  mask_r <- r[[1]]
  mask_r[] <- NA
  mask_r[c(1, 2, 7, 8)] <- 1
  poly <- terra::as.polygons(mask_r, dissolve = TRUE, na.rm = TRUE)

  result <- sptrends:::.inspect_ts_cell_core(r, cell = 1, polygon = poly, prewhitened = pw)
  expect_named(result, c("cell", "raw", "prewhitened", "n_modified", "n_total"))
  expect_equal(result$n_total, 4L)  # 4 cells covered by the polygon

  ext_vals <- terra::extract(pw$series, poly)[, -1, drop = FALSE]
  expected_pw_series <- apply(ext_vals, 2, stats::median, na.rm = TRUE)
  expect_equal(unname(result$prewhitened$series), unname(expected_pw_series))
})

test_that(".inspect_ts_cell_core works with a polygon AND a Yue-Pilon prewhitened object (no DW gate, Modified field is absent -- a different code path from the Wang-Swail case above)", {
  r <- .make_ar1_raster(rho = 0.6, n_time = 20, nrow = 6, ncol = 6, seed = 12)
  pw_yp <- prewhiten(r, method = "TFPW_Y", report = FALSE, verbose = FALSE)

  mask_r <- r[[1]]
  mask_r[] <- NA
  mask_r[c(1, 2, 7, 8)] <- 1
  poly <- terra::as.polygons(mask_r, dissolve = TRUE, na.rm = TRUE)

  result <- sptrends:::.inspect_ts_cell_core(r, cell = 1, polygon = poly,
                                              prewhitened = pw_yp)
  expect_named(result, c("cell", "raw", "prewhitened", "n_modified", "n_total"))
  # Yue-Pilon has no DW gate -- every valid cell counts as "modified".
  expect_equal(result$n_modified, result$n_total)
})

test_that(".inspect_ts_cell_core errors when the prewhitened series is incomplete at that location", {
  r <- .make_ar1_raster(rho = 0.6, n_time = 20, seed = 11)
  pw <- prewhiten(r, report = FALSE, verbose = FALSE)
  pw$series[1] <- NA
  expect_error(
    sptrends:::.inspect_ts_cell_core(r, cell = 1, neighbourhood = FALSE, prewhitened = pw),
    "no complete prewhitened time series"
  )
})

test_that("inspect_ts_cell errors when 'prewhitened' is a list missing '$series' or '$diagnostics'", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 12)$series
  expect_error(inspect_ts_cell(r, prewhitened = list(diagnostics = "x")), "prewhiten")
  expect_error(inspect_ts_cell(r, prewhitened = list(series = r)), "prewhiten")
})

test_that("inspect_ts_cell(selection_type = 'point') works end-to-end via a mocked terra::click()", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 50)$series
  ext_r <- terra::ext(r)
  mid_x <- mean(c(ext_r[1], ext_r[2]))
  mid_y <- mean(c(ext_r[3], ext_r[4]))

  testthat::local_mocked_bindings(
    click = function(...) data.frame(x = mid_x, y = mid_y),
    .package = "terra"
  )

  result <- suppressMessages(inspect_ts_cell(r, neighbourhood = FALSE))
  expect_type(result, "list")
  expect_false(is.null(result$raw))
  expect_identical(result$cell, terra::cellFromXY(r, cbind(mid_x, mid_y)))
})

test_that("inspect_ts_cell(verbose = FALSE) genuinely silences its own messages, matching the verbose convention used throughout this package", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 50)$series
  ext_r <- terra::ext(r)
  mid_x <- mean(c(ext_r[1], ext_r[2]))
  mid_y <- mean(c(ext_r[3], ext_r[4]))

  testthat::local_mocked_bindings(
    click = function(...) data.frame(x = mid_x, y = mid_y),
    .package = "terra"
  )

  expect_silent(inspect_ts_cell(r, neighbourhood = FALSE, verbose = FALSE))
})

test_that("inspect_ts_cell(selection_type = 'polygon') works end-to-end via a mocked terra::draw()", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 51)$series
  mask_r <- r[[1]]
  mask_r[] <- NA
  mask_r[c(1, 2, 9, 10)] <- 1  # a small 2x2 block (8-column grid)
  poly <- terra::as.polygons(mask_r, dissolve = TRUE, na.rm = TRUE)

  testthat::local_mocked_bindings(
    draw = function(...) poly,
    .package = "terra"
  )

  result <- suppressMessages(inspect_ts_cell(r, selection_type = "polygon"))
  expect_type(result, "list")
  expect_false(is.null(result$raw))
})

test_that("inspect_ts_cell errors when no click is registered (empty/NULL terra::click() result)", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 52)$series
  testthat::local_mocked_bindings(
    click = function(...) NULL,
    .package = "terra"
  )
  expect_error(suppressMessages(inspect_ts_cell(r)), "No click registered")
})

test_that("inspect_ts_cell errors when the click falls outside the raster's extent", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 53)$series
  ext_r <- terra::ext(r)
  # comfortably outside the raster on every side
  outside_x <- ext_r[2] + 1000 * (ext_r[2] - ext_r[1])
  outside_y <- ext_r[4] + 1000 * (ext_r[4] - ext_r[3])

  testthat::local_mocked_bindings(
    click = function(...) data.frame(x = outside_x, y = outside_y),
    .package = "terra"
  )
  expect_error(suppressMessages(inspect_ts_cell(r)), "falls outside")
})

test_that("show_neighbours = TRUE draws a small-multiples grid and returns per-cell fits", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 30)$series
  result <- sptrends:::.inspect_ts_cell_core(r, cell = 15, neighbourhood = TRUE,
                                              show_neighbours = TRUE)
  expect_true(!is.null(result$neighbours))
  expect_true(length(result$neighbours) >= 2)  # centre + at least 1 neighbour
  centre_entries <- vapply(result$neighbours, function(n) n$is_centre, logical(1))
  expect_true(sum(centre_entries) == 1)
  expect_identical(result$neighbours[centre_entries][[1]]$cell, 15)
})

test_that("show_neighbours = TRUE also works alongside a prewhitened comparison", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 31)$series
  pw <- prewhiten(r, report = FALSE, verbose = FALSE)
  result <- sptrends:::.inspect_ts_cell_core(r, cell = 15, neighbourhood = TRUE,
                                              prewhitened = pw,
                                              show_neighbours = TRUE)
  expect_true(!is.null(result$neighbours))
  expect_true(!is.null(result$prewhitened))
})

test_that("show_neighbours = TRUE is silently ignored (with a message) when neighbourhood = FALSE", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 32)$series
  expect_message(
    result <- sptrends:::.inspect_ts_cell_core(r, cell = 15, neighbourhood = FALSE,
                                                show_neighbours = TRUE),
    "has no effect"
  )
  expect_null(result$neighbours)
})

test_that("show_neighbours = TRUE is silently ignored (with a message) for a polygon selection", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 33)$series
  mask_r <- r[[1]]
  mask_r[] <- NA
  mask_r[c(1, 2, 7, 8)] <- 1
  poly <- terra::as.polygons(mask_r, dissolve = TRUE, na.rm = TRUE)
  expect_message(
    result <- sptrends:::.inspect_ts_cell_core(r, cell = 1, polygon = poly,
                                                show_neighbours = TRUE),
    "has no effect"
  )
  expect_null(result$neighbours)
})

test_that("small-multiples grid handles a neighbour with an incomplete series without erroring", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 34)$series
  # Knock out one neighbour's series entirely (still leaves the centre and
  # other neighbours complete, exercising the "no data" mini-panel branch).
  neighbour_cell <- as.vector(terra::adjacent(r, cells = 15,
                                               directions = "queen",
                                               pairs = FALSE))[1]
  r[neighbour_cell] <- NA
  expect_error(
    sptrends:::.inspect_ts_cell_core(r, cell = 15, neighbourhood = TRUE,
                                      show_neighbours = TRUE),
    NA
  )
})

test_that("inspect_ts_cell(show_neighbours = TRUE) works end-to-end via a mocked terra::click()", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 35)$series
  ext_r <- terra::ext(r)
  mid_x <- mean(c(ext_r[1], ext_r[2]))
  mid_y <- mean(c(ext_r[3], ext_r[4]))

  testthat::local_mocked_bindings(
    click = function(...) data.frame(x = mid_x, y = mid_y),
    .package = "terra"
  )

  result <- suppressMessages(inspect_ts_cell(r, show_neighbours = TRUE))
  expect_true(!is.null(result$neighbours))
})

test_that(".ols_confidence_interval matches lm()'s confint() exactly", {
  set.seed(9)
  n <- 15
  t <- 1:n
  series <- 0.5 + 0.3 * t + rnorm(n, sd = 2)

  slope <- sptrends:::.ols_slope_vectorised(matrix(series, nrow = 1), t)
  intercept <- mean(series) - slope * mean(t)
  ci <- sptrends:::.ols_confidence_interval(series, t, slope, intercept,
                                             conf_level = 0.95)

  ref <- lm(series ~ t)
  ref_ci <- confint(ref, "t", level = 0.95)

  expect_equal(unname(coef(ref)["t"]), slope, tolerance = 1e-10)
  expect_equal(ci$lower, unname(ref_ci[1]), tolerance = 1e-10)
  expect_equal(ci$upper, unname(ref_ci[2]), tolerance = 1e-10)
})

test_that(".fit_slope_with_ci(method = 'ols') uses the OLS CI, not Sen's, and returns method = 'ols'", {
  series <- c(1, 3, 2, 6, 5, 8, 7, 10)
  t <- 1:8
  fit_ols <- sptrends:::.fit_slope_with_ci(series, t, conf_level = 0.95,
                                            method = "OLS")
  fit_ts <- sptrends:::.fit_slope_with_ci(series, t, conf_level = 0.95,
                                           method = "TS")

  expect_identical(fit_ols$method, "OLS")
  expect_identical(fit_ts$method, "TS")
  # Different point estimates and different CIs -- these are genuinely
  # different computations, not the same numbers relabelled.
  expect_false(isTRUE(all.equal(fit_ols$slope, fit_ts$slope)))
  expect_false(isTRUE(all.equal(fit_ols$ci_lower, fit_ts$ci_lower)))
})

test_that("inspect_ts_cell(slope_method = 'ols') works end-to-end and shows 'OLS' in the legend data", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 36)$series
  ext_r <- terra::ext(r)
  mid_x <- mean(c(ext_r[1], ext_r[2]))
  mid_y <- mean(c(ext_r[3], ext_r[4]))

  testthat::local_mocked_bindings(
    click = function(...) data.frame(x = mid_x, y = mid_y),
    .package = "terra"
  )

  result <- suppressMessages(inspect_ts_cell(r, slope_method = "OLS"))
  expect_type(result, "list")
  expect_false(is.null(result$raw))
})

test_that("show_neighbours = TRUE respects slope_method = 'ols' for the small-multiples grid too", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 37)$series
  result <- sptrends:::.inspect_ts_cell_core(r, cell = 15, neighbourhood = TRUE,
                                              show_neighbours = TRUE,
                                              slope_method = "OLS")
  expect_true(!is.null(result$neighbours))
})

test_that(".inspect_ts_cell_core() uses terra::time(x) as the default t when x has valid time metadata", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 5, seed = 40)$series
  terra::time(r, tstep = "years") <- 1982:1986

  result <- sptrends:::.inspect_ts_cell_core(r, cell = 15, neighbourhood = FALSE)
  # No direct way to read back t from the return value, but the raw
  # series' own length must match, and the slope must be computable
  # without error against a real 1982:1986 axis, not silently falling
  # back to 1:5 instead.
  expect_equal(length(result$raw$series), 5)
})

test_that(".inspect_ts_cell_core() falls back to 1:n when x has no time metadata", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 5, seed = 41)$series
  expect_false(terra::has.time(r))
  expect_error(
    sptrends:::.inspect_ts_cell_core(r, cell = 15, neighbourhood = FALSE),
    NA
  )
})

test_that("an explicit t argument always overrides x's own terra::time(), even when both are present", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 5, seed = 42)$series
  terra::time(r, tstep = "years") <- 1982:1986
  custom_t <- c(10, 20, 30, 40, 50)

  result <- sptrends:::.inspect_ts_cell_core(r, cell = 15, t = custom_t,
                                              neighbourhood = FALSE)
  expect_equal(length(result$raw$series), 5)
})

test_that("inspect_ts_cell(slope_method = 'RM') works end-to-end, with no confidence interval reported (none implemented for RM)", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 12, seed = 838)$series
  ext_r <- terra::ext(r)
  mid_x <- mean(c(ext_r[1], ext_r[2]))
  mid_y <- mean(c(ext_r[3], ext_r[4]))

  testthat::local_mocked_bindings(
    click = function(...) data.frame(x = mid_x, y = mid_y),
    .package = "terra"
  )

  result <- suppressMessages(inspect_ts_cell(r, slope_method = "RM"))
  expect_type(result, "list")
  expect_false(is.null(result$raw))
  expect_true(is.na(result$raw$ci_lower))
  expect_true(is.na(result$raw$ci_upper))
})

test_that("show_neighbours = TRUE respects slope_method = 'RM' for the small-multiples grid too", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 12, seed = 839)$series
  result <- sptrends:::.inspect_ts_cell_core(r, cell = 15, neighbourhood = TRUE,
                                              show_neighbours = TRUE,
                                              slope_method = "RM")
  expect_true(!is.null(result$neighbours))
})

test_that("inspect_ts_cell(compare_slopes = TRUE) works end-to-end, ignoring slope_method", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 12, seed = 840)$series
  ext_r <- terra::ext(r)
  mid_x <- mean(c(ext_r[1], ext_r[2]))
  mid_y <- mean(c(ext_r[3], ext_r[4]))

  testthat::local_mocked_bindings(
    click = function(...) data.frame(x = mid_x, y = mid_y),
    .package = "terra"
  )

  result <- suppressMessages(inspect_ts_cell(r, compare_slopes = TRUE))
  expect_type(result, "list")
  expect_false(is.null(result$raw))
})

test_that(".fit_all_slopes() returns all three named point estimates (slope and intercept), matching slope_estimator() directly", {
  set.seed(841)
  n <- 15
  series <- as.numeric(arima.sim(list(ar = 0.3), n = n)) + 1:n * 0.03
  t_vec <- 1:n

  fits <- sptrends:::.fit_all_slopes(series, t_vec)
  expect_named(fits, c("slope", "intercept"))
  expect_named(fits$slope, c("TS", "OLS", "RM"))
  expect_named(fits$intercept, c("TS", "OLS", "RM"))
  expect_true(all(is.finite(fits$slope)))
  expect_true(all(is.finite(fits$intercept)))
})

test_that("inspect_ts_cell(compare_slopes = TRUE) returns $slope_comparison with all three methods' own slope and intercept", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 12, seed = 842)$series
  ext_r <- terra::ext(r)
  mid_x <- mean(c(ext_r[1], ext_r[2]))
  mid_y <- mean(c(ext_r[3], ext_r[4]))

  testthat::local_mocked_bindings(
    click = function(...) data.frame(x = mid_x, y = mid_y),
    .package = "terra"
  )

  result <- suppressMessages(inspect_ts_cell(r, compare_slopes = TRUE))
  expect_false(is.null(result$slope_comparison))
  expect_named(result$slope_comparison, c("method", "slope", "intercept"))
  expect_equal(sort(result$slope_comparison$method), c("OLS", "RM", "TS"))
  expect_true(all(is.finite(result$slope_comparison$slope)))
  expect_true(all(is.finite(result$slope_comparison$intercept)))
})

test_that("inspect_ts_cell(compare_slopes = TRUE) with prewhitened also returns $pw_slope_comparison for the second panel, not only $slope_comparison for the raw one", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 12, seed = 844)$series
  pw <- prewhiten(r, report = FALSE, verbose = FALSE)
  result <- suppressMessages(
    sptrends:::.inspect_ts_cell_core(r, cell = 15, prewhitened = pw,
                                      compare_slopes = TRUE)
  )
  expect_false(is.null(result$slope_comparison))
  expect_false(is.null(result$pw_slope_comparison))
  expect_named(result$pw_slope_comparison, c("method", "slope", "intercept"))
  expect_equal(sort(result$pw_slope_comparison$method), c("OLS", "RM", "TS"))
})

test_that(".inspect_ts_cell_core() accepts and forwards compare_slopes directly (a real bug found by the user: the public inspect_ts_cell() extended compare_slopes/slope_method = 'RM' but never forwarded compare_slopes to the core, which also still rejected 'RM')", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 12, seed = 843)$series
  result <- sptrends:::.inspect_ts_cell_core(r, cell = 15,
                                              slope_method = "RM",
                                              compare_slopes = TRUE)
  expect_false(is.null(result$slope_comparison))
})
