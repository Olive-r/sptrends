test_that("workflow_rta() errors on non-SpatRaster input", {
  expect_error(workflow_rta(matrix(1:9, 3, 3)), "SpatRaster")
})

test_that("workflow_rta() returns the expected structure with class 'rta'", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_rta(r, report = FALSE, verbose = FALSE)
  expect_s3_class(result, "rta")
  expect_named(result, c("theil_sen", "trend", "trend_summary_table", "fdr",
                          "theil_sen_smoothed", "timing"))
})

test_that("workflow_rta() does not prewhiten -- theil_sen and trend match calling the raw functions directly on x", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 2)$series
  result <- workflow_rta(r, report = FALSE, verbose = FALSE)

  direct_slope <- slope_estimator(r, verbose = FALSE, report = FALSE)$slope
  expect_equal(terra::values(result$theil_sen, mat = FALSE),
               terra::values(direct_slope, mat = FALSE))

  direct_trend <- trend_test(r, report = FALSE, verbose = FALSE)
  expect_equal(terra::values(result$trend$p, mat = FALSE),
               terra::values(direct_trend$stats$p, mat = FALSE))
})

test_that("workflow_rta()'s theil_sen step uses smooth_neighbourhood = FALSE by default (same as workflow_tst(), now aligned)", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 3)$series
  result <- workflow_rta(r, report = FALSE, verbose = FALSE)
  raw_slope <- slope_estimator(r, smooth_neighbourhood = FALSE, verbose = FALSE, report = FALSE)$slope
  expect_equal(terra::values(result$theil_sen, mat = FALSE),
               terra::values(raw_slope, mat = FALSE))
  expect_false(result$theil_sen_smoothed)
})

test_that("workflow_rta() only computes FDR-BH, not BKY", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 4)$series
  result <- workflow_rta(r, report = FALSE, verbose = FALSE)
  expect_null(result$fdr$reject_BKY)
  expect_false(is.null(result$fdr$reject_BH))
})

test_that("workflow_rta() never runs a Moran's I check", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 5)$series
  result <- workflow_rta(r, report = FALSE, verbose = FALSE)
  expect_null(result$fdr$moran)
  expect_null(result$fdr$moran_recommendation)
  expect_null(result$fdr$moran_assessment)
})

test_that("workflow_rta() forwards cmk_args and theil_sen_args correctly", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 20, seed = 6)$series
  result <- workflow_rta(r, cmk_args = list(method = "MK"),
                 theil_sen_args = list(max_pairs = 50),
                 report = FALSE, verbose = FALSE)
  # method = "MK" -> classic MK -> trend has "S"/"VarS", not "Sm"/"VarSm"
  expect_true("S" %in% names(result$trend))
  expect_false("Sm" %in% names(result$trend))
})

test_that("print.rta runs without error and reports the expected fields", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 7)$series
  result <- workflow_rta(r, report = FALSE, verbose = FALSE)
  expect_output(print(result), "Robust Trend Analysis")
  expect_output(print(result), "Theil-Sen slope")
  expect_output(print(result), "Significant after FDR-BH")
})

test_that("workflow_rta() messages every step label when verbose = TRUE", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 8)$series
  msgs <- character(0)
  withCallingHandlers(
    workflow_rta(r, report = FALSE, verbose = TRUE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  full <- paste(msgs, collapse = " ")
  expect_true(grepl("Theil-Sen slope", full, fixed = TRUE))
  expect_true(grepl("Contextual Mann-Kendall", full, fixed = TRUE))
  expect_true(grepl("FDR correction", full, fixed = TRUE))
})

test_that("summary.rta runs without error and returns trend/fdr tables", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 10)$series
  result <- workflow_rta(r, report = FALSE, verbose = FALSE)
  out <- capture.output(smry <- summary(result))
  expect_true(any(grepl("Trend test", out, fixed = TRUE)))
  expect_true(any(grepl("Theil-Sen slope", out, fixed = TRUE)))
  expect_true(any(grepl("FDR correction", out, fixed = TRUE)))
  expect_named(smry, c("trend", "fdr"))
})

test_that("plot.workflow_rta(which = 'direction') draws the default map", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 11)$series
  result <- workflow_rta(r, report = FALSE, verbose = FALSE)
  expect_error(plot(result), NA)
})

test_that("plot.rta supports which = 'significance' and which = 'trend'", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 12)$series
  result <- workflow_rta(r, report = FALSE, verbose = FALSE)
  expect_error(plot(result, which = "significance"), NA)
  expect_error(plot(result, which = "trend"), NA)
})

test_that("plot.workflow_rta(which = 'slope') runs without error, smoothed and unsmoothed", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 13)$series
  result <- workflow_rta(r, report = FALSE, verbose = FALSE)
  expect_error(plot(result, which = "slope"), NA)
  expect_error(plot(result, which = "slope", smooth = FALSE), NA)
})

test_that("plot.workflow_rta(which = 'slope', smooth = FALSE) messages when source smoothing was already applied", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 13)$series
  result <- workflow_rta(r, theil_sen_args = list(smooth_neighbourhood = TRUE),
                report = FALSE, verbose = FALSE)
  expect_true(result$theil_sen_smoothed)
  expect_message(plot(result, which = "slope", smooth = FALSE),
                  "already computed")
})

test_that("plot.workflow_rta(which = 'slope') never modifies the underlying rta object", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 14)$series
  result <- workflow_rta(r, report = FALSE, verbose = FALSE)
  before <- terra::values(result$theil_sen, mat = FALSE)
  plot(result, which = "slope")
  after <- terra::values(result$theil_sen, mat = FALSE)
  expect_identical(before, after)
})

test_that("plot.rta returns x invisibly", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 15)$series
  result <- workflow_rta(r, report = FALSE, verbose = FALSE)
  ret <- plot(result)
  expect_s3_class(ret, "rta")
})

test_that("workflow_rta() runs cleanly across every combination of cmk neighbourhood x theil_sen smoothing (a systematic sweep)", {
  grid <- expand.grid(
    cmk_neighbourhood = c(TRUE, FALSE),
    ts_smooth = c(TRUE, FALSE)
  )
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 101)$series

  for (i in seq_len(nrow(grid))) {
    cmk_i <- grid$cmk_neighbourhood[i]
    ts_i <- grid$ts_smooth[i]
    label <- sprintf("cmk_neighbourhood=%s, ts_smooth=%s", cmk_i, ts_i)

    result <- workflow_rta(r, cmk_args = list(method = if (cmk_i) "CMK" else "MK"),
                  theil_sen_args = list(smooth_neighbourhood = ts_i),
                  report = FALSE, verbose = FALSE)

    expect_s3_class(result, "rta")
    expect_error(print(result), NA, info = label)
    expect_error(capture.output(summary(result)), NA, info = label)
  }
})

test_that("plot.workflow_rta(which = 'slope') handles an all-NA/all-zero masked map without erroring (max_abs guard)", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 16)$series
  result <- workflow_rta(r, report = FALSE, verbose = FALSE)
  # Force zero significant cells, so the significance-masked slope map
  # is entirely NA -- the same guard already covered for trend_maps()
  # and slope_map(), exercised here for the "rta" case of
  # plot.sptrends()'s own version.
  result$fdr$reject_BH[] <- 0
  expect_error(plot(result, which = "slope"), NA)
})

test_that("plot.workflow_rta() supports the 8 new uncorrected slope_*/pvalue_* views", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 260)$series
  result <- workflow_rta(r, report = FALSE, verbose = FALSE)
  for (w in c("slope_map", "slope_direction", "slope_hist", "slope_bar",
              "pvalue_map", "pvalue_significance", "pvalue_hist", "pvalue_bar")) {
    expect_error(plot(result, which = w), NA, info = w)
  }
})

test_that("plot.workflow_rta(which = 'slope_map') errors defensively if x$theil_sen is somehow missing", {
  fake <- list(theil_sen = NULL)
  class(fake) <- c("rta", "sptrends")
  expect_error(plot(fake, which = "slope_map"), "should never happen")
})

test_that("plot.workflow_rta(which = 'direction') uses the already-smoothed slope as-is when theil_sen_smoothed = TRUE, without re-smoothing", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12,
                       trend_shape = "block", signal_size = c(8, 8),
                       trend_fraction = 1, trend_strength = 0.3,
                       noise_sd = 0.3, seed = 1)$series
  result <- workflow_rta(
    r, theil_sen_args = list(smooth_neighbourhood = TRUE),
    report = FALSE, verbose = FALSE
  )
  expect_true(isTRUE(result$theil_sen_smoothed))
  expect_error(plot(result, which = "direction"), NA)
})
