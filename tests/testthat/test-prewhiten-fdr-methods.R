test_that("prewhiten() and fdr_correction() results have the 'sptrends' superclass", {
  r <- .make_ar1_raster(rho = 0.6, n_time = 15, seed = 70)
  pw <- prewhiten(r, report = FALSE, verbose = FALSE)
  expect_identical(class(pw), c("prewhiten", "sptrends"))
  expect_s3_class(pw, "sptrends")

  r2 <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 71)$series
  trend <- trend_test(r2, report = FALSE, verbose = FALSE)
  fdr_res <- fdr_correction(trend$stats$p, report = FALSE, verbose = FALSE)
  expect_identical(class(fdr_res), c("fdr", "sptrends"))
  expect_s3_class(fdr_res, "sptrends")
})

test_that("print.prewhiten runs without error and returns invisibly", {
  r <- .make_ar1_raster(rho = 0.6, n_time = 15, seed = 72)
  pw <- prewhiten(r, report = FALSE, verbose = FALSE)
  out <- capture.output(ret <- print(pw))
  expect_true(any(grepl("Prewhitened", out, fixed = TRUE)))
  expect_identical(ret, pw)
})

test_that("summary.prewhiten calls prewhiten_summary() and returns its table", {
  r <- .make_ar1_raster(rho = 0.6, n_time = 15, seed = 73)
  pw <- prewhiten(r, report = FALSE, verbose = FALSE)
  tab <- suppressMessages(summary(pw))
  expect_true("valid_cells" %in% tab$metric)
})

test_that("plot.prewhiten supports both 'maps' and 'histograms'", {
  r <- .make_ar1_raster(rho = 0.6, n_time = 15, seed = 74)
  pw <- prewhiten(r, report = FALSE, verbose = FALSE)
  expect_error(plot(pw), NA)
  expect_error(plot(pw, which = "histograms"), NA)
})

test_that("print.fdr runs without error and returns invisibly", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 75)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  fdr_res <- fdr_correction(trend$stats$p, report = FALSE, verbose = FALSE)
  out <- capture.output(ret <- print(fdr_res))
  expect_true(any(grepl("significant", out, fixed = TRUE)))
  expect_identical(ret, fdr_res)
})

test_that("summary.fdr calls fdr_summary() and returns its table", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 76)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  fdr_res <- fdr_correction(trend$stats$p, report = FALSE, verbose = FALSE)
  tab <- suppressMessages(summary(fdr_res))
  expect_true("method" %in% names(tab))
})

test_that("plot.fdr supports all four 'which' options", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 77)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  fdr_res <- fdr_correction(trend$stats$p, report = FALSE, verbose = FALSE)
  expect_error(plot(fdr_res), NA)
  expect_error(plot(fdr_res, which = "pvalue_histogram"), NA)
  expect_error(plot(fdr_res, which = "comparison"), NA)
  expect_error(plot(fdr_res, which = "threshold"), NA)
})

test_that("plot.fdr(which = 'significance') errors clearly on a plain numeric vector input", {
  p <- c(0.01, 0.2, 0.5, 0.03)
  fdr_res <- fdr_correction(p, report = FALSE, verbose = FALSE)
  expect_error(plot(fdr_res, which = "significance"), "no rasters")
})

test_that("plot.fdr(which = 'threshold') errors clearly when BKY was not requested", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 78)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  fdr_res <- fdr_correction(trend$stats$p, method = "BH", report = FALSE, verbose = FALSE)
  expect_error(plot(fdr_res, which = "threshold"), "no BKY results")
})

test_that("fdr_correction() stores the raw p-values in $p, for both vector and raster input", {
  p_vec <- c(0.01, 0.2, 0.5, 0.03)
  fdr_vec <- fdr_correction(p_vec, report = FALSE, verbose = FALSE)
  expect_equal(fdr_vec$p, p_vec)

  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 79)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  fdr_r <- fdr_correction(trend$stats$p, report = FALSE, verbose = FALSE)
  expect_equal(fdr_r$p, terra::values(trend$stats$p, mat = FALSE))
})
