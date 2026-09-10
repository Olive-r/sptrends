test_that("workflow_trends() runs end to end with its own defaults and returns the expected structure", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 12, seed = 801)$series
  result <- workflow_trends(r, report = FALSE, verbose = FALSE)

  expect_s3_class(result, c("workflow_trends", "sptrends"))
  expect_identical(class(result), c("workflow_trends", "sptrends"))
  expect_true(all(c("prewhiten", "trend", "trend_summary_table", "slope",
                     "fdr", "timing") %in% names(result)))
  expect_false(is.null(result$prewhiten))
  expect_false(is.null(result$fdr))
  expect_true(all(c("prewhiten", "trend", "slope", "fdr") %in%
                    names(result$timing)))
})

test_that("workflow_trends() runs validly across a non-default combination of every method (TFPW_Y, MK, OLS, BH), with no warning", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 12, seed = 802)$series
  expect_warning(
    result <- workflow_trends(
      r, prewhiten_method = "TFPW_Y", trend_method = "MK",
      slope_method = "OLS", fdr_method = "BH",
      report = FALSE, verbose = FALSE),
    NA
  )
  expect_false(is.null(result$prewhiten))
  expect_identical(result$prewhiten$method, "TFPW_Y")
})

test_that("workflow_trends() warns when prewhiten_method != 'none' is combined with trend_method = 'MMK'", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 803)$series
  expect_warning(
    workflow_trends(r, trend_method = "MMK", report = FALSE,
                      verbose = FALSE),
    "both correct temporal autocorrelation"
  )
})

test_that("workflow_trends() does not warn about MMK when prewhiten_method = 'none'", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 804)$series
  expect_warning(
    workflow_trends(r, prewhiten_method = "none", trend_method = "MMK",
                      report = FALSE, verbose = FALSE),
    NA
  )
})

test_that("workflow_trends(prewhiten_method = 'none') skips prewhitening entirely", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 807)$series
  result <- workflow_trends(r, prewhiten_method = "none", report = FALSE,
                              verbose = FALSE)
  expect_true(is.null(result$prewhiten))
  expect_false("prewhiten" %in% names(result$timing))
})

test_that("workflow_trends(fdr_method = NULL) skips FDR correction entirely", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 808)$series
  result <- workflow_trends(r, fdr_method = NULL, report = FALSE,
                              verbose = FALSE)
  expect_true(is.null(result$fdr))
  expect_false("fdr" %in% names(result$timing))
})

test_that("workflow_trends() errors on non-SpatRaster input", {
  expect_error(
    workflow_trends(matrix(1:10, nrow = 2), report = FALSE,
                      verbose = FALSE),
    "must be a terra SpatRaster"
  )
})

test_that("workflow_trends(verbose = TRUE) prints its own step-by-step progress messages, naming the chosen method at each step", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 12, seed = 809)$series
  expect_message(
    workflow_trends(r, trend_method = "MK", verbose = TRUE, report = FALSE),
    "trend test \\(MK\\)"
  )
  expect_message(
    workflow_trends(r, prewhiten_method = "none", verbose = TRUE,
                      report = FALSE),
    "prewhitening skipped"
  )
  expect_message(
    workflow_trends(r, fdr_method = NULL, verbose = TRUE, report = FALSE),
    "FDR correction skipped"
  )
})

test_that("workflow_trends(n_cores = 2) gives the same trend result as the sequential path", {
  .skip_unless_parallel_tests()
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 12, seed = 810)$series
  res_seq <- workflow_trends(r, n_cores = 1, report = FALSE,
                               verbose = FALSE)
  res_par <- workflow_trends(r, n_cores = 2, report = FALSE,
                               verbose = FALSE)
  expect_equal(terra::values(res_seq$trend$Sm, mat = FALSE),
               terra::values(res_par$trend$Sm, mat = FALSE))
})

test_that("workflow_trends() forwards *_args correctly to the underlying functions", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 811)$series
  result <- workflow_trends(
    r, trend_method = "MK", trend_args = list(ties = TRUE),
    fdr_args = list(q = 0.1),
    report = FALSE, verbose = FALSE)
  expect_false(is.null(result$trend))
  expect_false(is.null(result$fdr))
})

test_that("print(), summary(), and plot() dispatch correctly for a workflow_trends object", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 812)$series
  result <- workflow_trends(r, report = FALSE, verbose = FALSE)

  expect_output(print(result), "<workflow_trends result>")
  expect_output(print(result), "Trend test:")

  smry <- summary(result)
  expect_true(is.list(smry))
  expect_true("trend" %in% names(smry))

  expect_error(plot(result, which = "trend"), NA)
  expect_error(plot(result, which = "slope"), NA)
  expect_error(plot(result, which = "significance"), NA)
})

test_that("plot(which = 'slope') errors clearly when there is no FDR correction to mask by", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 813)$series
  result <- workflow_trends(r, fdr_method = NULL, report = FALSE,
                              verbose = FALSE)
  expect_error(plot(result, which = "slope"), "No FDR correction")
})

test_that("plot(which = 'significance') errors clearly when there is no FDR correction", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 814)$series
  result <- workflow_trends(r, fdr_method = NULL, report = FALSE,
                              verbose = FALSE)
  expect_error(plot(result, which = "significance"), "No FDR correction")
})

test_that("print() covers every prewhitening branch: TFPW_Y, a Modified-based method, 'none', and FDR skipped/BH-only", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 815)$series

  res_tfpw_y <- workflow_trends(r, prewhiten_method = "TFPW_Y",
                                  report = FALSE, verbose = FALSE)
  expect_output(print(res_tfpw_y), "all \\d+ valid cells")

  res_ws <- workflow_trends(r, prewhiten_method = "TFPW_WS",
                              report = FALSE, verbose = FALSE)
  expect_output(print(res_ws), "cells modified")

  res_none <- workflow_trends(r, prewhiten_method = "none",
                                report = FALSE, verbose = FALSE)
  expect_output(print(res_none), "Prewhitening: skipped")

  res_no_fdr <- workflow_trends(r, fdr_method = NULL, report = FALSE,
                                  verbose = FALSE)
  expect_output(print(res_no_fdr), "FDR correction: skipped")

  res_bh <- workflow_trends(r, fdr_method = "BH", report = FALSE,
                              verbose = FALSE)
  expect_output(print(res_bh), "Significant after FDR-BH:")
})

test_that("print() shows the Moran's I assessment when moran_check = TRUE is forwarded via fdr_args", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 816)$series
  result <- workflow_trends(
    r, fdr_args = list(moran_check = TRUE), report = FALSE, verbose = FALSE)
  expect_output(print(result), "Moran's I assessment:")
})

test_that("summary() covers both the FDR-present and FDR-skipped paths", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 817)$series

  res_fdr <- workflow_trends(r, report = FALSE, verbose = FALSE)
  smry_fdr <- summary(res_fdr)
  expect_false(is.null(smry_fdr$fdr))

  res_no_fdr <- workflow_trends(r, fdr_method = NULL, report = FALSE,
                                  verbose = FALSE)
  smry_no_fdr <- summary(res_no_fdr)
  expect_true(is.null(smry_no_fdr$fdr))
})

test_that("plot(which = 'slope') falls back to reject_BH when only BH is available (no BKY)", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 818)$series
  result <- workflow_trends(r, fdr_method = "BH", report = FALSE,
                              verbose = FALSE)
  expect_error(plot(result, which = "slope"), NA)
})

test_that("plot(which = 'slope', smooth = TRUE) runs without error", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 819)$series
  result <- workflow_trends(r, report = FALSE, verbose = FALSE)
  expect_error(plot(result, which = "slope", smooth = TRUE), NA)
})

test_that("workflow_trends(slope_method = NULL) skips slope estimation entirely", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 820)$series
  result <- workflow_trends(r, slope_method = NULL, report = FALSE,
                              verbose = FALSE)
  expect_true(is.null(result$slope))
  expect_false("slope" %in% names(result$timing))
})

test_that("workflow_trends(slope_method = NULL, verbose = TRUE) prints its own skip message", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 12, seed = 821)$series
  expect_message(
    workflow_trends(r, slope_method = NULL, verbose = TRUE, report = FALSE),
    "slope estimation skipped"
  )
})

test_that("print() correctly reports the 'S' statistic for trend_method = 'MK', and 'beta' for 'OLS'", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 822)$series

  res_mk <- workflow_trends(r, trend_method = "MK", report = FALSE,
                              verbose = FALSE)
  expect_output(print(res_mk), "\\(S statistic\\)")

  res_ols <- workflow_trends(r, trend_method = "OLS", report = FALSE,
                               verbose = FALSE)
  expect_output(print(res_ols), "\\(beta statistic\\)")
})

test_that("plot(which = 'slope') errors clearly when slope_method = NULL", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 823)$series
  result <- workflow_trends(r, slope_method = NULL, report = FALSE,
                              verbose = FALSE)
  expect_error(plot(result, which = "slope"), "No slope")
})

test_that("plot() covers all twelve 'which' options, for full parity with workflow_tst()'s own plot()", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 12, seed = 824)$series
  result <- workflow_trends(r, report = FALSE, verbose = FALSE)

  every_which <- c("direction", "significance", "trend", "slope",
                    "slope_map", "slope_direction", "slope_hist",
                    "slope_bar", "pvalue_map", "pvalue_significance",
                    "pvalue_hist", "pvalue_bar")
  for (w in every_which) {
    expect_error(plot(result, which = w), NA, info = w)
  }
})

test_that("plot(which = 'direction') errors clearly with no FDR correction", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 825)$series
  result <- workflow_trends(r, fdr_method = NULL, report = FALSE,
                              verbose = FALSE)
  expect_error(plot(result, which = "direction"), "No FDR correction")
})

test_that("plot(which = 'slope_map'/'slope_direction'/'slope_hist'/'slope_bar') error clearly with slope_method = NULL", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 826)$series
  result <- workflow_trends(r, slope_method = NULL, report = FALSE,
                              verbose = FALSE)
  for (w in c("slope_map", "slope_direction", "slope_hist", "slope_bar")) {
    expect_error(plot(result, which = w), "No slope", info = w)
  }
})

test_that("plot(which = 'direction') falls back to the trend statistic's own sign when slope_method = NULL", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12,
                       trend_shape = "block", signal_size = c(8, 8),
                       trend_fraction = 1, trend_strength = 0.3,
                       noise_sd = 0.3, seed = 1)$series
  result <- workflow_trends(r, slope_method = NULL, report = FALSE,
                              verbose = FALSE)
  expect_true(is.null(result$slope))
  expect_error(plot(result, which = "direction"), NA)
})

test_that("plot(which = 'direction'/'slope') correctly picks reject_BH, not reject_BKY, when fdr_method = 'BH'", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 12, seed = 827)$series
  result <- workflow_trends(r, fdr_method = "BH", report = FALSE,
                              verbose = FALSE)
  expect_error(plot(result, which = "direction"), NA)
  expect_error(plot(result, which = "slope"), NA)
})

test_that("workflow_trends(slope_method = 'RM') works -- a real bug found during a documentation review: match.arg(slope_method, c('TS', 'OLS')) omitted 'RM' entirely, so this call errored before this fix despite slope_estimator() itself supporting 'RM' correctly", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 828)$series
  expect_error(
    result <- workflow_trends(r, slope_method = "RM", report = FALSE,
                                verbose = FALSE),
    NA
  )
  expect_false(is.null(result$slope))
})

test_that("workflow_trends(fdr_method = 'BY') works as an opt-in dependence safeguard", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 829)$series
  expect_error(
    result <- workflow_trends(r, fdr_method = "BY", report = FALSE,
                                verbose = FALSE),
    NA
  )
  expect_false(is.null(result$fdr))
  expect_false(is.null(result$fdr$reject_BY))
  expect_true(is.null(result$fdr$reject_BH))
  expect_true(is.null(result$fdr$reject_BKY))
  expect_output(print(result), "FDR-BY")
  expect_output(summary(result), "FDR correction")
  for (view in c("direction", "significance", "slope")) {
    expect_error(plot(result, which = view, smooth = FALSE), NA)
  }

  inconsistent <- result
  inconsistent$fdr$reject_BH <- inconsistent$fdr$reject_BY
  expect_error(
    plot(inconsistent, which = "direction"),
    "must contain exactly one"
  )
})

test_that("print() covers the VCTFPW prewhiten branch", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12, seed = 831)$series
  result <- workflow_trends(r, prewhiten_method = "VCTFPW",
                              report = FALSE, verbose = FALSE)
  expect_output(
    print(result),
    "Prewhitening \\(VCTFPW\\): [0-9]+ of [0-9]+ cells modified"
  )
})

test_that("workflow_trends forwards the configured CMK window", {
  r <- sim_trend_stack(nrow = 7, ncol = 7, n_time = 8, seed = 909)$series
  direct <- trend_test(r, window_size = 5L, report = FALSE,
                       verbose = FALSE)$stats
  result <- workflow_trends(
    r, prewhiten_method = "none", trend_method = "CMK",
    trend_args = list(window_size = 5L), slope_method = NULL,
    fdr_method = NULL, report = FALSE, verbose = FALSE
  )

  expect_equal(terra::values(result$trend, mat = TRUE),
               terra::values(direct, mat = TRUE), tolerance = 0)
})
