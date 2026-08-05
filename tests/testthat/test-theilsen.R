test_that("slope_estimator recovers a known linear trend", {
  n <- 15
  t <- seq_len(n)
  vals <- 2 + 0.5 * t
  r <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3)
  layers <- lapply(t, function(k) terra::setValues(r, rep(vals[k], 9)))
  x <- do.call(c, layers)

  slope <- slope_estimator(x, verbose = FALSE, report = FALSE)$slope
  expect_equal(unname(terra::values(slope, mat = FALSE)), rep(0.5, 9), tolerance = 1e-8)
})

test_that("slope_estimator subsamples pairs when max_pairs is exceeded", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 40, seed = 1)$series
  expect_message(slope_estimator(r, max_pairs = 50, seed = 1, verbose = TRUE, report = FALSE)$slope, "random sample")
})

test_that("slope_estimator gives the same result with n_cores = 1 and 2", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 10, seed = 2)$series
  slope_seq <- slope_estimator(r, n_cores = 1, verbose = FALSE, report = FALSE)$slope
  slope_par <- slope_estimator(r, n_cores = 2, verbose = FALSE, report = FALSE)$slope
  expect_equal(terra::values(slope_seq, mat = FALSE), terra::values(slope_par, mat = FALSE))
})

test_that("slope_estimator matches a hand-computed median of pairwise slopes", {
  # t = 1:4, values = 1,3,2,6. The 6 pairwise slopes are:
  # -1, 0.5, 1.5, 5/3, 2, 4 -- median (average of the 3rd/4th sorted) = 19/12.
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  layers <- lapply(c(1, 3, 2, 6), function(v) terra::setValues(r, rep(v, 4)))
  x <- do.call(c, layers)

  slope <- slope_estimator(x, max_pairs = Inf, verbose = FALSE, report = FALSE)$slope
  expect_equal(unname(terra::values(slope, mat = FALSE)), rep(19 / 12, 4), tolerance = 1e-10)
})

test_that("slope_estimator errors on fewer than 2 layers", {
  r <- terra::rast(nrows = 2, ncols = 2, xmin = 0, xmax = 2, ymin = 0, ymax = 2)
  x <- terra::setValues(r, 1:4)
  expect_error(slope_estimator(x, verbose = FALSE, report = FALSE)$slope, "at least 2 layers")
})

test_that("slope_estimator errors on non-SpatRaster input", {
  expect_error(slope_estimator(matrix(1:4, 2, 2), verbose = FALSE, report = FALSE)$slope, "SpatRaster")
})

test_that("slope_estimator validates the time axis for every method", {
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 6, seed = 602)$series
  for (method in c("TS", "OLS", "RM")) {
    expect_error(
      slope_estimator(r, method = method, t = 1:5,
                       report = FALSE, verbose = FALSE),
      "must have length"
    )
    expect_error(
      slope_estimator(r, method = method, t = c(1, 2, 3, 3, 5, 6),
                       report = FALSE, verbose = FALSE),
      "must not contain duplicate"
    )
    expect_error(
      slope_estimator(r, method = method, t = c(1, 2, 4, 3, 5, 6),
                       report = FALSE, verbose = FALSE),
      "strictly increasing"
    )
    expect_error(
      slope_estimator(r, method = method, t = c(1, 2, 3, 4, 5, Inf),
                       report = FALSE, verbose = FALSE),
      "finite"
    )
  }
})

test_that("slope_estimator respects a custom time vector", {
  # same values as the linear-trend test, but t spaced by 2 instead of 1 --
  # the slope per unit t should halve accordingly.
  n <- 15
  t_default <- seq_len(n)
  vals <- 2 + 0.5 * t_default
  r <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3)
  layers <- lapply(vals, function(v) terra::setValues(r, rep(v, 9)))
  x <- do.call(c, layers)

  slope_custom_t <- slope_estimator(x, t = seq(0, (n - 1) * 2, by = 2), verbose = FALSE, report = FALSE)$slope
  expect_equal(unname(terra::values(slope_custom_t, mat = FALSE)), rep(0.25, 9), tolerance = 1e-8)
})

test_that("slope_estimator messages 'Exact' when max_pairs is not exceeded", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 6, seed = 1)$series
  expect_message(slope_estimator(r, verbose = TRUE, report = FALSE)$slope, "Exact Theil-Sen")
})

test_that("slope_estimator subsamples and messages accordingly when max_pairs is exceeded", {
  # n=6 gives 15 possible pairs; forcing max_pairs=5 well below that.
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 6, seed = 1)$series
  expect_message(
    slope_estimator(r, max_pairs = 5, seed = 1, verbose = TRUE, report = FALSE)$slope,
    "random sample"
  )
})

test_that("slope_estimator's subsampled estimate is reproducible with the same seed", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 6, seed = 1)$series
  slope_a <- slope_estimator(r, max_pairs = 5, seed = 42, verbose = FALSE, report = FALSE)$slope
  slope_b <- slope_estimator(r, max_pairs = 5, seed = 42, verbose = FALSE, report = FALSE)$slope
  expect_equal(terra::values(slope_a, mat = FALSE), terra::values(slope_b, mat = FALSE))
})

test_that("slope_estimator leaves cells with missing values as NA", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 10, seed = 3)$series
  r[1] <- NA
  slope <- slope_estimator(r, verbose = FALSE, report = FALSE)$slope
  expect_true(is.na(unname(terra::values(slope, mat = FALSE))[1]))
})

test_that("slope_estimator messages the parallel-cores line when verbose = TRUE and n_cores > 1", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 8, seed = 1)$series
  expect_message(slope_estimator(r, n_cores = 2, verbose = TRUE, report = FALSE)$slope, "Parallel over")
})

test_that("slope_estimator handles n_cores greater than the number of valid cells", {
  # only 1 valid cell (rest NA), requesting 2 cores -- chunking must not error.
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 8, seed = 1)$series
  r[2:16] <- NA
  slope <- slope_estimator(r, n_cores = 2, verbose = FALSE, report = FALSE)$slope
  expect_equal(sum(!is.na(terra::values(slope, mat = FALSE))), 1)
})

test_that("slope_estimator's smooth_neighbourhood = FALSE (default) is unchanged from plain per-cell Theil-Sen", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 5)$series
  slope_default <- slope_estimator(r, verbose = FALSE, report = FALSE)$slope
  slope_explicit_false <- slope_estimator(r, smooth_neighbourhood = FALSE, verbose = FALSE, report = FALSE)$slope
  expect_equal(terra::values(slope_default), terra::values(slope_explicit_false))
})

test_that("slope_estimator's smooth_neighbourhood = TRUE changes the result and runs without error", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 6)$series
  slope_plain <- slope_estimator(r, verbose = FALSE, report = FALSE)$slope
  slope_smoothed <- slope_estimator(r, smooth_neighbourhood = TRUE, verbose = FALSE, report = FALSE)$slope
  expect_false(isTRUE(all.equal(terra::values(slope_plain), terra::values(slope_smoothed))))
  expect_named(slope_smoothed, "theilsen_slope")
})

test_that("slope_estimator's smooth_neighbourhood = TRUE messages a caveat when verbose = TRUE", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 7)$series
  expect_message(
    slope_estimator(r, smooth_neighbourhood = TRUE, verbose = TRUE, report = FALSE)$slope,
    "off by default"
  )
})

test_that("slope_estimator's smoothing matches a manual queen-3x3 median focal", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 8)$series
  slope_plain <- slope_estimator(r, verbose = FALSE, report = FALSE)$slope
  slope_smoothed <- slope_estimator(r, smooth_neighbourhood = TRUE, verbose = FALSE, report = FALSE)$slope

  w <- matrix(1, nrow = 3, ncol = 3)
  manual_smoothed <- terra::focal(slope_plain, w = w, fun = "median", na.rm = TRUE)
  expect_equal(unname(terra::values(slope_smoothed)), unname(terra::values(manual_smoothed)))
})

test_that("workflow_tst() forwards smooth_neighbourhood to slope_estimator() via theil_sen_args", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 9)$series
  # prewhiten = FALSE so workflow_tst() feeds slope_estimator() the same raw
  # series used directly below -- otherwise the comparison would be
  # between a prewhitened-then-smoothed and a raw-then-smoothed result.
  result <- workflow_tst(r, prewhiten = FALSE, theil_sen_args = list(smooth_neighbourhood = TRUE),
                report = FALSE, verbose = FALSE)
  slope_direct <- slope_estimator(r, smooth_neighbourhood = TRUE, verbose = FALSE, report = FALSE)$slope
  expect_equal(terra::values(result$theil_sen), terra::values(slope_direct))
})

test_that("slope_summary reports the expected metrics", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 40)$series
  slope <- slope_estimator(r, verbose = FALSE, report = FALSE)$slope
  tab <- suppressMessages(slope_summary(slope))
  expect_named(tab, c("metric", "value"))
  expect_true(all(c("valid_cells", "min_slope", "median_slope", "mean_slope",
                     "max_slope", "pct_increasing", "pct_decreasing",
                     "pct_flat") %in% tab$metric))
  vals <- terra::values(slope, mat = FALSE)
  expect_equal(as.numeric(tab$value[tab$metric == "valid_cells"]),
               sum(!is.na(vals)))
})

test_that("slope_summary prints informative messages", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 41)$series
  slope <- slope_estimator(r, verbose = FALSE, report = FALSE)$slope
  msgs <- testthat::capture_messages(slope_summary(slope))
  full <- paste(msgs, collapse = " ")
  expect_true(grepl("Valid cells", full, fixed = TRUE))
  expect_true(grepl("Slope range", full, fixed = TRUE))
  expect_true(grepl("Increasing", full, fixed = TRUE))
})

test_that("slope_map runs without error", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 42)$series
  slope <- slope_estimator(r, verbose = FALSE, report = FALSE)$slope
  expect_error(slope_map(slope), NA)
})

test_that("slope_estimator(report = TRUE) auto-prints the summary and draws the map", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 43)$series
  expect_message(
    slope_estimator(r, verbose = FALSE, report = TRUE)$slope,
    "Valid cells"
  )
})

test_that("slope_estimator(report = FALSE) stays silent regardless of verbose", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 44)$series
  expect_no_message(slope_estimator(r, verbose = FALSE, report = FALSE)$slope)
})

test_that("workflow_tst()'s and workflow_rta()'s own Theil-Sen step honours their report argument (no double-reporting)", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 45)$series
  # slope_summary()'s distinctive message text must not leak out when
  # the parent workflow itself was asked to stay quiet.
  expect_no_message(workflow_tst(r, report = FALSE, verbose = FALSE))
  expect_no_message(workflow_rta(r, report = FALSE, verbose = FALSE))
})

test_that("slope_summary can write a CSV", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 46)$series
  slope <- slope_estimator(r, verbose = FALSE, report = FALSE)$slope
  path <- tempfile(fileext = ".csv")
  suppressMessages(slope_summary(slope, path = path))
  expect_true(file.exists(path))
  tab_read <- utils::read.csv(path)
  expect_true("valid_cells" %in% tab_read$metric)
})

test_that("slope_map can write a PNG", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 47)$series
  slope <- slope_estimator(r, verbose = FALSE, report = FALSE)$slope
  path <- tempfile()
  slope_map(slope, path = path)
  expect_true(file.exists(paste0(path, "_theilsen_map.png")))
  unlink(paste0(path, "_theilsen_map.png"))
})

test_that("slope_map handles an all-NA slope without erroring (max_abs guard)", {
  r_na <- terra::rast(nrows = 4, ncols = 4, vals = NA_real_)
  expect_error(slope_map(r_na), NA)
})

test_that("slope_estimator(smooth_neighbourhood = TRUE) never fills a no-data cell from its neighbours (na.policy fix)", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 8, seed = 48)$series
  # Cell 13 (the centre of a 5x5 grid) has no data at all (every layer
  # NA there); every other cell has a complete, valid time series.
  # Before the na.policy = "omit" fix, terra::focal()'s own default
  # ("all") would assign cell 13 a smoothed value anyway, since its
  # queen neighbours all have valid slopes -- silently extending the
  # "has data" footprint into a cell that never had any observations.
  r[13] <- NA_real_

  slope <- slope_estimator(r, smooth_neighbourhood = TRUE, verbose = FALSE,
                            report = FALSE)$slope
  vals <- terra::values(slope, mat = FALSE)
  expect_true(is.na(vals[13]))
  expect_false(anyNA(vals[-13]))
})

test_that("slope_estimator() result has the 'theil'/'sptrends' classes", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 50)$series
  result <- slope_estimator(r, verbose = FALSE, report = FALSE)
  expect_identical(class(result), c("slope", "sptrends"))
  expect_s3_class(result, "sptrends")
  expect_true(inherits(result$slope, "SpatRaster"))
  expect_false(result$smoothed)
})

test_that("slope_estimator()'s $smoothed field reflects what actually happened", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 51)$series
  result_off <- slope_estimator(r, smooth_neighbourhood = FALSE, verbose = FALSE,
                                 report = FALSE)
  result_on <- slope_estimator(r, smooth_neighbourhood = TRUE, verbose = FALSE,
                                report = FALSE)
  expect_false(result_off$smoothed)
  expect_true(result_on$smoothed)
})

test_that("print.theil runs without error and returns invisibly", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 52)$series
  result <- slope_estimator(r, verbose = FALSE, report = FALSE)
  out <- capture.output(ret <- print(result))
  expect_true(any(grepl("Valid cells", out, fixed = TRUE)))
  expect_identical(ret, result)
})

test_that("print.theil notes when the slope was smoothed", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 53)$series
  result <- slope_estimator(r, smooth_neighbourhood = TRUE, verbose = FALSE,
                             report = FALSE)
  out <- capture.output(print(result))
  expect_true(any(grepl("smoothed", out, fixed = TRUE)))
})

test_that("summary.theil calls slope_summary() and returns its table", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 54)$series
  result <- slope_estimator(r, verbose = FALSE, report = FALSE)
  tab <- suppressMessages(summary(result))
  expect_true("valid_cells" %in% tab$metric)
})

test_that("plot.theil calls slope_map() without error", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 55)$series
  result <- slope_estimator(r, verbose = FALSE, report = FALSE)
  expect_error(plot(result), NA)
})

test_that("plot.theil(which = ...) dispatches to the other 3 reporting functions without error, not just the default map", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 55)$series
  result <- slope_estimator(r, verbose = FALSE, report = FALSE)
  expect_error(plot(result, which = "direction"), NA)
  expect_error(plot(result, which = "histogram"), NA)
  expect_error(plot(result, which = "bar"), NA)
})

test_that("workflow_tst()'s and workflow_rta()'s own $theil_sen field stays a bare SpatRaster (not the 'slope' wrapper)", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 56)$series
  result_tst <- workflow_tst(r, report = FALSE, verbose = FALSE)
  result_rta <- workflow_rta(r, report = FALSE, verbose = FALSE)
  expect_true(inherits(result_tst$theil_sen, "SpatRaster"))
  expect_false(inherits(result_tst$theil_sen, "slope"))
  expect_true(inherits(result_rta$theil_sen, "SpatRaster"))
  expect_false(inherits(result_rta$theil_sen, "slope"))
})

test_that("slope_estimator(method = 'ols') matches lm()'s slope exactly, cell by cell", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 60)$series
  result <- slope_estimator(r, method = "OLS", verbose = FALSE, report = FALSE)
  expect_identical(names(result$slope), "ols_slope")

  X <- terra::values(r, mat = TRUE)
  t <- 1:10
  expected <- apply(X, 1, function(row) {
    if (anyNA(row)) return(NA_real_)
    unname(stats::coef(stats::lm(row ~ t))[2])
  })
  actual <- terra::values(result$slope, mat = FALSE)
  expect_equal(actual, expected, tolerance = 1e-10)
})

test_that("slope_estimator(method = 'ols') handles NA cells the same way method = 'theilsen' does", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 8, seed = 61)$series
  r[1] <- NA  # knock out one cell entirely
  result_ols <- slope_estimator(r, method = "OLS", verbose = FALSE, report = FALSE)
  result_ts <- slope_estimator(r, method = "TS", verbose = FALSE, report = FALSE)
  na_ols <- is.na(terra::values(result_ols$slope, mat = FALSE))
  na_ts <- is.na(terra::values(result_ts$slope, mat = FALSE))
  expect_identical(na_ols, na_ts)
  expect_true(na_ols[1])
})

test_that("slope_estimator(method = 'ols') respects a custom t and ignores max_pairs/seed/n_cores", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 6, seed = 62)$series
  t_custom <- c(2000, 2001, 2003, 2004, 2005, 2009)  # uneven spacing
  expect_error(
    slope_estimator(r, method = "OLS", t = t_custom, max_pairs = 5,
                     seed = 1, n_cores = 2, verbose = FALSE, report = FALSE),
    NA
  )
})

test_that("slope_estimator(method = 'ols') supports smooth_neighbourhood, same class/layer-name convention as theilsen", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 63)$series
  result <- slope_estimator(r, method = "OLS", smooth_neighbourhood = TRUE,
                             verbose = FALSE, report = FALSE)
  expect_true(result$smoothed)
  expect_identical(names(result$slope), "ols_slope")
  expect_s3_class(result, "slope")
})

test_that("slope_estimator(method = 'ols') works with report = TRUE (reuses slope_summary()/slope_map())", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 64)$series
  expect_error(
    slope_estimator(r, method = "OLS", verbose = FALSE, report = TRUE),
    NA
  )
})

test_that("workflow_tst()/workflow_rta() forward method = 'ols' through theil_sen_args correctly", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 65)$series
  result_tst <- workflow_tst(r, theil_sen_args = list(method = "OLS"),
                     report = FALSE, verbose = FALSE)
  result_rta <- workflow_rta(r, theil_sen_args = list(method = "OLS"),
                     report = FALSE, verbose = FALSE)
  expect_true(inherits(result_tst$theil_sen, "SpatRaster"))
  expect_true(inherits(result_rta$theil_sen, "SpatRaster"))
})

test_that("slope_estimator errors clearly with fewer than 2 layers", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 1, seed = 310)$series
  expect_error(slope_estimator(r, verbose = FALSE), "at least 2 layers")
})

test_that("slope_estimator works at the exact minimum of 2 layers, without erroring", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 2, seed = 311)$series
  result <- slope_estimator(r, report = FALSE, verbose = FALSE)
  expect_true(inherits(result$slope, "SpatRaster"))
  expect_equal(terra::nlyr(result$slope), 1)
})

test_that("slope_estimator on a perfectly constant series gives exactly 0, not NA or an error", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 8, seed = 312)$series
  r[] <- 7  # every cell, every layer: exactly the same constant value
  result <- slope_estimator(r, report = FALSE, verbose = FALSE)
  slope_vals <- terra::values(result$slope, mat = FALSE)
  expect_true(all(slope_vals == 0, na.rm = TRUE))
  expect_false(anyNA(slope_vals))
})

test_that("slope_estimator does not error on a 1x1 raster", {
  r <- sim_trend_stack(nrow = 1, ncol = 1, n_time = 8, smooth_radius = 0,
                        seed = 313)$series
  expect_error(slope_estimator(r, report = FALSE, verbose = FALSE), NA)
})

test_that("slope_estimator errors clearly, not cryptically, when every cell is NA", {
  r_na <- terra::rast(nrows = 4, ncols = 4, nlyr = 5, vals = NA_real_)
  expect_error(slope_estimator(r_na, verbose = FALSE), "nothing to estimate")
})

test_that("slope_estimator handles a very long series (100 layers) without erroring", {
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 100, seed = 314)$series
  expect_error(slope_estimator(r, report = FALSE, verbose = FALSE), NA)
})

test_that("slope_estimator(method = 'ols', verbose = TRUE) prints its own progress messages, not just theilsen's (a genuinely separate code path)", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 8, seed = 315)$series
  expect_message(
    slope_estimator(r, method = "OLS", verbose = TRUE, report = FALSE),
    "Ordinary least squares"
  )
})

test_that("slope_estimator(method = 'ols', smooth_neighbourhood = TRUE, verbose = TRUE) prints the smoothing message on the ols path too, not just theilsen's", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 8, seed = 316)$series
  expect_message(
    slope_estimator(r, method = "OLS", smooth_neighbourhood = TRUE,
                     verbose = TRUE, report = FALSE),
    "queen-3x3 median smoothing"
  )
})
