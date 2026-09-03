test_that("print.sptrends() correctly dispatches all six classes", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 200)$series
  pw <- .make_ar1_raster(rho = 0.6, n_time = 15, seed = 201)

  tst_result <- workflow_tst(r, report = FALSE, verbose = FALSE)
  rta_result <- workflow_rta(r, report = FALSE, verbose = FALSE)
  cmk_result <- trend_test(r, report = FALSE, verbose = FALSE)
  theil_result <- slope_estimator(r, report = FALSE, verbose = FALSE)
  prewhiten_result <- prewhiten(pw, report = FALSE, verbose = FALSE)
  fdr_result <- fdr_correction(cmk_result$stats$p, report = FALSE, verbose = FALSE)

  objects <- list(tst_result, rta_result, cmk_result, theil_result,
                   prewhiten_result, fdr_result)
  expected_labels <- c("True Significant Trends", "Robust Trend Analysis",
                        "Mann-Kendall", "Theil-Sen", "prewhitening",
                        "FDR correction")

  for (i in seq_along(objects)) {
    out <- capture.output(ret <- print(objects[[i]]))
    expect_true(any(grepl(expected_labels[i], out, fixed = TRUE)),
                info = expected_labels[i])
    expect_identical(ret, objects[[i]])
  }
})

test_that("summary.sptrends() correctly dispatches all six classes without error", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 202)$series
  pw <- .make_ar1_raster(rho = 0.6, n_time = 15, seed = 203)

  tst_result <- workflow_tst(r, report = FALSE, verbose = FALSE)
  rta_result <- workflow_rta(r, report = FALSE, verbose = FALSE)
  cmk_result <- trend_test(r, report = FALSE, verbose = FALSE)
  theil_result <- slope_estimator(r, report = FALSE, verbose = FALSE)
  prewhiten_result <- prewhiten(pw, report = FALSE, verbose = FALSE)
  fdr_result <- fdr_correction(cmk_result$stats$p, report = FALSE, verbose = FALSE)

  for (obj in list(tst_result, rta_result, cmk_result, theil_result,
                    prewhiten_result, fdr_result)) {
    expect_error(suppressMessages(summary(obj)), NA)
  }
})

test_that("plot.sptrends() correctly dispatches all six classes without error", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 204)$series
  pw <- .make_ar1_raster(rho = 0.6, n_time = 15, seed = 205)

  tst_result <- workflow_tst(r, report = FALSE, verbose = FALSE)
  rta_result <- workflow_rta(r, report = FALSE, verbose = FALSE)
  cmk_result <- trend_test(r, report = FALSE, verbose = FALSE)
  theil_result <- slope_estimator(r, report = FALSE, verbose = FALSE)
  prewhiten_result <- prewhiten(pw, report = FALSE, verbose = FALSE)
  fdr_result <- fdr_correction(cmk_result$stats$p, report = FALSE, verbose = FALSE)

  for (obj in list(tst_result, rta_result, cmk_result, theil_result,
                    prewhiten_result, fdr_result)) {
    expect_error(plot(obj), NA)
  }
})

test_that("print/summary/plot.sptrends() error clearly on an object with an unrecognised subclass", {
  fake <- list(a = 1)
  class(fake) <- c("not_a_real_class", "sptrends")
  expect_error(print(fake), "Unknown 'sptrends' object class")
  expect_error(summary(fake), "Unknown 'sptrends' object class")
  expect_error(plot(fake), "Unknown 'sptrends' object class")
})

test_that("the 18 previous per-class methods are gone from the public API; only the 3 unified ones remain", {
  removed <- c("print.tst", "print.rta", "print.cmk", "print.theil",
               "print.prewhiten", "print.fdr",
               "summary.tst", "summary.rta", "summary.cmk", "summary.theil",
               "summary.prewhiten", "summary.fdr",
               "plot.tst", "plot.rta", "plot.cmk", "plot.theil",
               "plot.prewhiten", "plot.fdr")
  exported <- getNamespaceExports("sptrends")
  expect_length(intersect(removed, exported), 0)

  # S3 methods registered only via S3method(), with no matching export(),
  # are intentionally absent from getNamespaceExports()/ls("package:...")
  # -- that's what makes them reachable only via generic dispatch
  # (print(x), never print.sptrends(x) directly), not a bug. The correct
  # way to confirm they're properly registered is getS3method().
  expect_false(is.null(getS3method("print", "sptrends", optional = TRUE)))
  expect_false(is.null(getS3method("summary", "sptrends", optional = TRUE)))
  expect_false(is.null(getS3method("plot", "sptrends", optional = TRUE)))
})

test_that("plot.sptrends() forwards ... correctly to each class's own which= options", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 206)$series
  cmk_result <- trend_test(r, report = FALSE, verbose = FALSE)
  expect_error(plot(cmk_result, which = "histograms"), NA)

  fdr_result <- fdr_correction(cmk_result$stats$p, report = FALSE, verbose = FALSE)
  expect_error(plot(fdr_result, which = "comparison"), NA)

  pw <- prewhiten(.make_ar1_raster(rho = 0.6, n_time = 15, seed = 207),
                              report = FALSE, verbose = FALSE)
  expect_error(plot(pw, which = "histograms"), NA)
})


test_that("print()/summary()/plot() dispatch correctly for 'spatial_autocorrelation' and 'compare_detections' too", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 1, seed = 210)$series[[1]]
  moran_result <- spatial_autocorrelation(r, nperm = 9, seed = 210,
                                           verbose = FALSE, report = FALSE)
  expect_s3_class(moran_result, "sptrends")
  expect_error(print(moran_result), NA)
  expect_error(suppressMessages(summary(moran_result)), NA)
  expect_error(plot(moran_result), NA)

  sim <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 10, seed = 211)
  trend_mk  <- trend_test(sim$series, method = "MK",
                                        report = FALSE, verbose = FALSE)
  trend_cmk <- trend_test(sim$series, method = "CMK",
                                        report = FALSE, verbose = FALSE)
  comparison <- compare_detections(
    detections = list(MK = trend_mk$stats$p <= 0.05,
                       CMK = trend_cmk$stats$p <= 0.05),
    ground_truth = sim$true_slope
  )
  expect_s3_class(comparison, "sptrends")
  expect_s3_class(comparison, "data.frame")
  expect_error(print(comparison), NA)
  expect_error(plot(comparison), NA)
  tab <- summary(comparison)
  expect_true(all(c("metric", "best_method") %in% names(tab)))
})

test_that("compare_detections() results keep working like an ordinary data frame", {
  sim <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 8, seed = 212)
  trend <- trend_test(sim$series, report = FALSE, verbose = FALSE)
  comparison <- compare_detections(
    detections = list(CMK = trend$stats$p <= 0.05),
    ground_truth = sim$true_slope
  )
  expect_identical(comparison$Method, "CMK")
  expect_equal(nrow(comparison), 1)
  expect_true("Sensitivity" %in% names(comparison))
})

test_that("spatial_autocorrelation_summary()/spatial_autocorrelation_null_plot()/classify_moran()/workflow_summary()/direction_map() are no longer part of the public API", {
  exported <- getNamespaceExports("sptrends")
  expect_false("spatial_autocorrelation_summary" %in% exported)
  expect_false("spatial_autocorrelation_null_plot" %in% exported)
  expect_false("classify_moran" %in% exported)
  expect_false("workflow_summary" %in% exported)
  # direction_map() was reconsidered and internalised: the binarised
  # trend direction is the same underlying computation as the
  # uncorrected one, just with a significance filter on top -- not two
  # genuinely different results needing two entry points. Reachable
  # via plot(x, which = "direction") for workflow_tst()/workflow_rta() results.
  expect_false("direction_map" %in% exported)
})

test_that("plot_detection_comparison() is part of the public API (exported for callers who want to visualise compare_detections() output directly, per CRAN feedback on ':::' usage in examples)", {
  exported <- getNamespaceExports("sptrends")
  expect_true("plot_detection_comparison" %in% exported)
})

test_that(".moran_category() covers all three descriptive categories (low/moderate/strong)", {
  expect_identical(sptrends:::.moran_category(0.05), "low")
  expect_identical(sptrends:::.moran_category(-0.05), "low")
  expect_identical(sptrends:::.moran_category(0.2), "moderate")
  expect_identical(sptrends:::.moran_category(-0.2), "moderate")
  expect_identical(sptrends:::.moran_category(0.5), "strong")
  expect_identical(sptrends:::.moran_category(-0.5), "strong")
})

test_that(".summary_compare_detections() handles a metric that is NA for every method (e.g. no method ever predicted positive)", {
  n <- 20
  ground_truth <- c(rep(TRUE, 5), rep(FALSE, n - 5))
  comparison <- compare_detections(
    detections = list(A = rep(FALSE, n), B = rep(FALSE, n)),
    ground_truth = ground_truth
  )
  # Precision = TP / (TP + FP); both methods never predict positive, so
  # TP + FP = 0 for both -- Precision is NA for every row.
  expect_true(all(is.na(comparison$Precision)))
  tab <- summary(comparison)
  precision_row <- tab[tab$metric == "Precision", ]
  expect_true(is.na(precision_row$best_method))
})
