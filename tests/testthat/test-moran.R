test_that("spatial_autocorrelation returns a valid result structure", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 1, seed = 1)$series[[1]]
  res <- spatial_autocorrelation(
    r, nperm = 19, seed = 1, verbose = FALSE, report = FALSE
  )

  expect_true(is.numeric(res$statistic))
  expect_true(res$p >= 0)
  expect_true(res$p <= 1)
  expect_length(res$null_dist, 19)
  expect_true(res$sign %in% c("positive", "negative"))
  expect_identical(res$scope, "global")
  expect_identical(
    class(res),
    c("spatial_autocorrelation_global", "spatial_autocorrelation",
      "sptrends")
  )
})

test_that("the explicit global scope preserves the previous default result", {
  r <- sim_trend_stack(nrow = 7, ncol = 7, n_time = 1, seed = 2)$series[[1]]
  implicit <- spatial_autocorrelation(
    r, nperm = 19, seed = 1, verbose = FALSE, report = FALSE
  )
  explicit <- spatial_autocorrelation(
    r, scope = "global", nperm = 19, seed = 1,
    verbose = FALSE, report = FALSE
  )

  expect_equal(implicit$statistic, explicit$statistic)
  expect_equal(implicit$p, explicit$p)
  expect_equal(implicit$null_dist, explicit$null_dist)
})

test_that("local scope returns raw permutation inference rasters", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 1, seed = 3)$series[[1]]

  local <- spatial_autocorrelation(
    r, scope = "local", nperm = 19, seed = 1,
    verbose = FALSE, report = FALSE
  )
  expect_s3_class(local, "spatial_autocorrelation_local")
  expect_s4_class(local$statistic, "SpatRaster")
  expect_s4_class(local$z, "SpatRaster")
  expect_s4_class(local$p, "SpatRaster")
  expect_s4_class(local$significant_raw, "SpatRaster")
  expect_false("p_adjusted" %in% names(local))
  expect_false(any(c("adjustment", "q") %in%
                     names(formals(spatial_autocorrelation))))
})

test_that("scope and alpha are validated", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 1, seed = 4)$series[[1]]

  expect_error(
    spatial_autocorrelation(
      r, scope = "regional", verbose = FALSE, report = FALSE
    ),
    "should be one of"
  )
  for (bad in list(0, 1, -0.1, 1.1, NA_real_, Inf, c(0.01, 0.05), "0.05")) {
    expect_error(
      spatial_autocorrelation(
        r, alpha = bad, verbose = FALSE, report = FALSE
      ),
      "'alpha'"
    )
  }
})

test_that("spatial_autocorrelation is stable across core counts", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 1, seed = 5)$series[[1]]

  res_seq <- spatial_autocorrelation(
    r, nperm = 9, seed = 1, n_cores = 1,
    verbose = FALSE, report = FALSE
  )
  res_par <- spatial_autocorrelation(
    r, nperm = 9, seed = 1, n_cores = 2,
    verbose = FALSE, report = FALSE
  )

  expect_equal(res_seq$statistic, res_par$statistic)
})

test_that("spatial_autocorrelation errors on multi-layer input", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 3, seed = 1)$series
  expect_error(spatial_autocorrelation(r, verbose = FALSE), "single layer")
})

test_that("non-finite raster values are not treated as valid cells", {
  r <- terra::rast(nrows = 3, ncols = 3)
  terra::values(r) <- c(1:8, Inf)
  expect_error(
    spatial_autocorrelation(r, verbose = FALSE, report = FALSE),
    "finite values or NA"
  )
})

test_that("classify_moran returns a label without erroring", {
  expect_silent(cat_label <- suppressMessages(classify_moran(0.5)))
  expect_true(cat_label %in% c("low", "moderate", "strong"))
})

test_that("report = TRUE prints the summary and draws the null plot", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 1, seed = 20)$series[[1]]
  expect_message(
    spatial_autocorrelation(
      r, nperm = 9, seed = 1, verbose = FALSE, report = TRUE
    ),
    "Observed I"
  )
})

test_that("report = FALSE stays silent", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 1, seed = 21)$series[[1]]
  expect_no_message(
    spatial_autocorrelation(
      r, nperm = 9, seed = 1, verbose = FALSE, report = FALSE
    )
  )
})

test_that("connectivity is validated via match.arg", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 1, seed = 22)$series[[1]]
  expect_error(
    spatial_autocorrelation(r, connectivity = "diagonal", nperm = 9,
                            verbose = FALSE, report = FALSE),
    "should be one of"
  )
})

test_that("print and summary of Moran results show the category", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 1, seed = 220)$series[[1]]
  moran_result <- spatial_autocorrelation(r, nperm = 9, seed = 220,
                                          verbose = FALSE, report = FALSE)
  out <- capture.output(print(moran_result))
  expect_true(any(grepl("category:", out, fixed = TRUE)))

  tab <- suppressMessages(summary(moran_result))
  expect_true("category" %in% tab$metric)
})
