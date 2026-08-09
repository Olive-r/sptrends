test_that("the 15 functions folded into S3 methods are no longer part of the public API", {
  internalised <- c(
    "trend_maps", "trend_histograms", "trend_summary",
    "slope_map", "slope_summary",
    "prewhiten_maps", "prewhiten_histograms", "prewhiten_summary",
    "fdr_summary", "fdr_pvalue_histogram", "fdr_significance_maps",
    "fdr_comparison_barplot", "fdr_threshold_plot",
    "fdr_direction_plot", "fdr_direction_summary"
  )
  exported <- getNamespaceExports("sptrends")
  still_exported <- intersect(internalised, exported)
  expect_length(still_exported, 0)
})

test_that("the 15 internalised functions still exist and work when called from within the package", {
  # testthat's own environment inherits from the package namespace, so
  # these remain directly callable by bare name in this file -- see
  # https://testthat.r-lib.org/reference/test_package.html. This test
  # exists to make that assumption explicit and pinned down, not to
  # exercise each function's own logic again (already covered by their
  # existing, dedicated test files).
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 90)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  expect_error(suppressMessages(trend_summary(trend$stats)), NA)
  expect_error(trend_histograms(trend$stats), NA)
  expect_error(trend_maps(trend$stats), NA)

  slope <- slope_estimator(r, report = FALSE, verbose = FALSE)
  expect_error(suppressMessages(slope_summary(slope$slope)), NA)
  expect_error(slope_map(slope$slope), NA)

  pw <- .make_ar1_raster(rho = 0.6, n_time = 15, seed = 91)
  pw_result <- prewhiten(pw, report = FALSE, verbose = FALSE)
  expect_error(suppressMessages(prewhiten_summary(pw_result$diagnostics)), NA)
  expect_error(prewhiten_histograms(pw_result$diagnostics), NA)
  expect_error(prewhiten_maps(pw_result$diagnostics), NA)

  fdr_result <- fdr_correction(trend$stats$p, report = FALSE, verbose = FALSE)
  expect_error(suppressMessages(fdr_summary(fdr_result)), NA)
  expect_error(fdr_pvalue_histogram(fdr_result$p), NA)
  expect_error(fdr_significance_maps(fdr_result), NA)
  expect_error(fdr_comparison_barplot(fdr_result), NA)
  expect_error(fdr_threshold_plot(fdr_result), NA)

  direction <- direction_map(trend$stats, fdr_result, method = "BH")
  expect_error(fdr_direction_plot(direction), NA)
  expect_error(suppressMessages(fdr_direction_summary(trend$stats, fdr_result)), NA)
})

test_that("fdr_bh() and fdr_bky() are no longer part of the public API either", {
  exported <- getNamespaceExports("sptrends")
  expect_false("fdr_bh" %in% exported)
  expect_false("fdr_bky" %in% exported)
})

test_that("fdr_bh()/fdr_bky() still exist and work when called from within the package, and fdr_correction() still calls them correctly", {
  p <- c(0.001, 0.01, 0.02, 0.5, 0.8)
  expect_error(fdr_bh(p), NA)
  expect_error(fdr_bky(p), NA)

  res <- fdr_correction(p, report = FALSE, verbose = FALSE)
  expect_true(!is.null(res$reject_BH))
  expect_true(!is.null(res$reject_BKY))
})

test_that("workflow_summary() and prepare_cmk_neighbourhood() are no longer part of the public API (method_citation(), classify_moran(), and workflow_summary() itself were removed entirely, not just hidden)", {
  exported <- getNamespaceExports("sptrends")
  expect_false("workflow_summary" %in% exported)
  expect_false("prepare_cmk_neighbourhood" %in% exported)
  expect_false("method_citation" %in% exported)
  expect_false("classify_moran" %in% exported)
  expect_false(exists("method_citation", where = asNamespace("sptrends"),
                       inherits = FALSE))
  expect_false(exists("workflow_summary", where = asNamespace("sptrends"),
                       inherits = FALSE))
})

test_that("prepare_cmk_neighbourhood() still exists and trend_test() still uses it internally when called from within the package", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 6, seed = 221)$series
  X <- terra::values(r, mat = TRUE)
  ok <- stats::complete.cases(X)
  expect_error(prepare_cmk_neighbourhood(r, ok), NA)

  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  expect_s3_class(trend, "trend_test")
})

test_that("fdr_direction_map() was renamed to direction_map(), and direction_map() is now internal (folded into plot(x, which = 'direction'))", {
  exported <- getNamespaceExports("sptrends")
  expect_false("fdr_direction_map" %in% exported)
  expect_false("direction_map" %in% exported)

  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 222)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  fdr_result <- fdr_correction(trend$stats$p, report = FALSE, verbose = FALSE)
  expect_error(sptrends:::direction_map(trend$stats, fdr_result, method = "BH"), NA)
})
