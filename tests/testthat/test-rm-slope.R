test_that("slope_estimator(method = 'RM') matches a direct, unvectorised reimplementation of Siegel's (1982) repeated median formula (slope and intercept), using the exact upper-median order statistic robslopes::RepeatedMedian()'s own source confirms (floor((k + 2) / 2)), not R's own median()", {
  # Historical note (belongs here, in a test, not in the user-facing
  # help page): a real, direct comparison against
  # robslopes::RepeatedMedian() on live data confirmed the slope
  # matched to floating-point precision from the start, but the
  # intercept originally differed by 0.119 -- a real, meaningful
  # discrepancy, not numerical noise -- while this port used the
  # simpler "hierarchical" intercept convention. Switching to the
  # "direct" convention this reference implementation itself uses
  # brought the intercept to floating-point-level agreement too
  # (< 1e-8). See NEWS.md for the full account of that comparison.
  upper_median_ref <- function(v) {
    k <- length(v)
    sort(v)[floor((k + 2) / 2)]
  }
  rm_reference <- function(x, t) {
    n <- length(x)
    inner_slope <- numeric(n)
    inner_intercept <- numeric(n)
    for (i in seq_len(n)) {
      j <- setdiff(seq_len(n), i)
      slopes <- (x[j] - x[i]) / (t[j] - t[i])
      # "Direct"/"separate" intercept convention, verified empirically
      # to match robslopes::RepeatedMedian()'s own output -- the
      # simpler "hierarchical" convention (median(y - slope * t))
      # this reference implementation used in an earlier version of
      # this test did not match it.
      intercepts <- x[i] - slopes * t[i]
      inner_slope[i] <- upper_median_ref(slopes)
      inner_intercept[i] <- upper_median_ref(intercepts)
    }
    list(slope = upper_median_ref(inner_slope),
         intercept = upper_median_ref(inner_intercept))
  }

  set.seed(1201)
  n <- 11  # odd, so the upper-median vs ordinary-median distinction
           # for the OUTER median specifically does not mask a real
           # implementation difference in this particular test
  x <- as.numeric(arima.sim(list(ar = 0.3), n = n)) + 1:n * 0.05
  t_vec <- 1:n
  ref <- rm_reference(x, t_vec)

  r <- terra::rast(nrows = 1, ncols = 1, nlyrs = n)
  r <- terra::setValues(r, matrix(x, nrow = 1))
  result <- slope_estimator(r, method = "RM", t = t_vec, report = FALSE,
                             verbose = FALSE)
  slope_sptrends <- as.numeric(terra::extract(result$slope, 1))
  intercept_sptrends <- as.numeric(terra::extract(result$intercept, 1))

  expect_equal(slope_sptrends, ref$slope, tolerance = 1e-10)
  expect_equal(intercept_sptrends, ref$intercept, tolerance = 1e-10)
})

test_that("slope_estimator(method = 'RM') matches the same reference implementation for an even n as well, where the upper-median-vs-ordinary-median distinction genuinely matters for the outer median too", {
  upper_median_ref <- function(v) {
    k <- length(v)
    sort(v)[floor((k + 2) / 2)]
  }
  rm_reference <- function(x, t) {
    n <- length(x)
    inner_slope <- numeric(n)
    inner_intercept <- numeric(n)
    for (i in seq_len(n)) {
      j <- setdiff(seq_len(n), i)
      slopes <- (x[j] - x[i]) / (t[j] - t[i])
      intercepts <- x[i] - slopes * t[i]
      inner_slope[i] <- upper_median_ref(slopes)
      inner_intercept[i] <- upper_median_ref(intercepts)
    }
    list(slope = upper_median_ref(inner_slope),
         intercept = upper_median_ref(inner_intercept))
  }

  set.seed(1202)
  n <- 12  # even
  x <- as.numeric(arima.sim(list(ar = 0.4), n = n)) + 1:n * 0.03
  t_vec <- 1:n
  ref <- rm_reference(x, t_vec)

  r <- terra::rast(nrows = 1, ncols = 1, nlyrs = n)
  r <- terra::setValues(r, matrix(x, nrow = 1))
  result <- slope_estimator(r, method = "RM", t = t_vec, report = FALSE,
                             verbose = FALSE)
  slope_sptrends <- as.numeric(terra::extract(result$slope, 1))
  intercept_sptrends <- as.numeric(terra::extract(result$intercept, 1))

  # This is also, deliberately, a check that the port does NOT use
  # R's own median() (which would average the two middle inner-median
  # values instead of taking the specific upper one) -- if it did,
  # this test would fail even though the odd-n test above might still
  # pass by coincidence.
  expect_equal(slope_sptrends, ref$slope, tolerance = 1e-10)
  expect_equal(intercept_sptrends, ref$intercept, tolerance = 1e-10)
})

test_that("slope_estimator(method = 'RM') runs across multiple cells at once and returns the expected object structure", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 1203)$series
  result <- slope_estimator(r, method = "RM", report = FALSE, verbose = FALSE)

  expect_s3_class(result, c("slope", "sptrends"))
  expect_equal(result$method, "RM")
  expect_true("intercept" %in% names(result))
  expect_s4_class(result$intercept, "SpatRaster")
  expect_equal(names(result$slope), "rm_slope")
  expect_equal(names(result$intercept), "rm_intercept")

  slope_vals <- terra::values(result$slope, mat = FALSE)
  expect_true(all(is.finite(slope_vals[!is.na(slope_vals)])))
})

test_that("slope_estimator(method = 'RM') errors on fewer than 3 layers", {
  r <- sim_trend_stack(nrow = 2, ncol = 2, n_time = 2, seed = 1204)$series
  expect_error(
    slope_estimator(r, method = "RM", report = FALSE, verbose = FALSE),
    "at least 3 layers"
  )
})

test_that("slope_estimator(method = 'RM') rejects duplicate time values", {
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 10, seed = 1205)$series
  t_dup <- c(1, 2, 3, 3, 5, 6, 7, 8, 9, 10)  # duplicate at positions 3-4
  expect_error(
    slope_estimator(r, method = "RM", t = t_dup, report = FALSE,
                     verbose = FALSE),
    "must not contain duplicate"
  )
})

test_that("slope_estimator(method = 'RM') does not warn for the ordinary case of unique, monotonic t", {
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 10, seed = 1206)$series
  expect_warning(
    slope_estimator(r, method = "RM", report = FALSE, verbose = FALSE),
    NA
  )
})

test_that("slope_estimator(method = 'RM') handles a cell with missing values, returning NA without erroring on the others", {
  set.seed(1207)
  n <- 12
  vals <- rbind(
    as.numeric(arima.sim(list(ar = 0.3), n = n)),
    { v <- as.numeric(arima.sim(list(ar = 0.3), n = n)); v[3] <- NA; v }
  )
  r <- terra::rast(nrows = 1, ncols = 2, nlyrs = n)
  r <- terra::setValues(r, vals)

  result <- slope_estimator(r, method = "RM", report = FALSE, verbose = FALSE)
  slope_vals <- terra::values(result$slope, mat = FALSE)
  intercept_vals <- terra::values(result$intercept, mat = FALSE)

  expect_true(is.finite(slope_vals[1]))
  expect_true(is.na(slope_vals[2]))
  expect_true(is.finite(intercept_vals[1]))
  expect_true(is.na(intercept_vals[2]))
})

test_that("slope_estimator(method = 'RM', smooth_neighbourhood = TRUE) runs without error", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 10, seed = 1208)$series
  expect_error(
    slope_estimator(r, method = "RM", smooth_neighbourhood = TRUE,
                     report = FALSE, verbose = FALSE),
    NA
  )
})

test_that("slope_estimator(method = 'RM', report = TRUE) calls slope_summary()/slope_map() directly, the same generic functions TS/OLS already use -- not a 'no dedicated function yet' message, superseded once .summary_slope()/.plot_slope() were confirmed to already work generically with any method's own slope raster", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 10, seed = 1209)$series
  expect_message(
    slope_estimator(r, method = "RM", report = TRUE, verbose = FALSE),
    "Slope range"
  )
})

test_that("slope_estimator(method = 'RM', verbose = TRUE) prints its own progress message", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 10, seed = 1210)$series
  expect_message(
    slope_estimator(r, method = "RM", verbose = TRUE, report = FALSE),
    "Siegel \\(1982\\)"
  )
})

test_that("print() correctly labels an 'RM' result, and (regression) an 'OLS' result -- both were previously mislabelled 'Theil-Sen slope result' by a hardcoded print header", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 10, seed = 1211)$series

  res_rm <- slope_estimator(r, method = "RM", report = FALSE, verbose = FALSE)
  expect_output(print(res_rm), "Repeated median \\(Siegel\\)")

  res_ols <- slope_estimator(r, method = "OLS", report = FALSE,
                              verbose = FALSE)
  expect_output(print(res_ols), "Ordinary least squares")

  res_ts <- slope_estimator(r, method = "TS", report = FALSE, verbose = FALSE)
  expect_output(print(res_ts), "Theil-Sen")
})

test_that("summary()/plot() run without error on an 'RM' result, reusing the same generic functions TS/OLS already use", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 10, seed = 1212)$series
  res_rm <- slope_estimator(r, method = "RM", report = FALSE, verbose = FALSE)

  expect_error(suppressMessages(summary(res_rm)), NA)
  expect_error(plot(res_rm, which = "map"), NA)
  expect_error(plot(res_rm, which = "direction"), NA)
  expect_error(plot(res_rm, which = "histogram"), NA)
  expect_error(plot(res_rm, which = "bar"), NA)
})

test_that("slope_estimator(method = 'RM') ignores max_pairs/seed/n_cores without error (documented as not used by this method)", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 10, seed = 1213)$series
  expect_error(
    slope_estimator(r, method = "RM", max_pairs = 5, seed = 99,
                     n_cores = 2, report = FALSE, verbose = FALSE),
    NA
  )
})

test_that("slope_estimator(method = 'RM', smooth_neighbourhood = TRUE, verbose = TRUE) prints its own smoothing message -- missed by the existing smooth_neighbourhood test, which used verbose = FALSE", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 10, seed = 1215)$series
  expect_message(
    slope_estimator(r, method = "RM", smooth_neighbourhood = TRUE,
                     verbose = TRUE, report = FALSE),
    "Applying queen-3x3 median smoothing"
  )
})
