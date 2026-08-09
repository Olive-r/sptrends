.make_pattern_raster <- function(values, nrow = 2, ncol = 2) {
  ncell_ <- nrow * ncol
  r <- terra::rast(nrows = nrow, ncols = ncol, xmin = 0, xmax = ncol, ymin = 0, ymax = nrow)
  layers <- lapply(values, function(v) terra::setValues(r, rep(v, ncell_)))
  do.call(c, layers)
}

test_that("compute_anomalies gives exactly zero anomalies for a purely repeating cycle", {
  # Same 3-value pattern repeated twice, identical in every cell -- the
  # climatology at each position exactly matches every observed value at
  # that position, so the anomaly must be exactly zero everywhere.
  r <- .make_pattern_raster(c(10, 20, 30, 10, 20, 30))
  res <- compute_anomalies(r, cycle = 3, verbose = FALSE)

  expect_equal(terra::nlyr(res$anomalies), 6)
  expect_equal(terra::nlyr(res$climatology), 3)
  expect_true(all(abs(terra::values(res$anomalies, mat = TRUE)) < 1e-10))
  expect_equal(unname(terra::values(res$climatology, mat = TRUE)[1, ]), c(10, 20, 30))
})

test_that("compute_anomalies computes the exact expected anomaly for a hand-worked case", {
  # position climatology: (10+12)/2=11, (20+22)/2=21, (30+32)/2=31
  # expected anomalies:    -1, -1, -1,   1, 1, 1
  r <- .make_pattern_raster(c(10, 20, 30, 12, 22, 32))
  res <- compute_anomalies(r, cycle = 3, verbose = FALSE)

  vals <- unname(terra::values(res$anomalies, mat = TRUE)[1, ])
  expect_equal(vals, c(-1, -1, -1, 1, 1, 1))
})

test_that("compute_anomalies standardise = TRUE gives z-scores with the expected sign and magnitude", {
  # position 1 values across cycles: 10, 12 -> mean 11, sd = sqrt(2)
  # anomaly at layer 1 = (10 - 11) / sqrt(2) = -1/sqrt(2)
  r <- .make_pattern_raster(c(10, 20, 30, 12, 22, 32))
  res <- compute_anomalies(r, cycle = 3, standardise = TRUE, verbose = FALSE)

  expect_equal(terra::nlyr(res$climatology_sd), 3)
  val1 <- unname(terra::values(res$anomalies, mat = TRUE)[1, 1])
  expect_equal(val1, -1 / sqrt(2), tolerance = 1e-8)
})

test_that("compute_anomalies gives NA (not Inf/NaN) where a cycle position has zero variance", {
  # position 1 is constant (10, 10) across both cycles -> sd = 0 there.
  r <- .make_pattern_raster(c(10, 20, 30, 10, 22, 32))
  res <- compute_anomalies(r, cycle = 3, standardise = TRUE, verbose = FALSE)

  val_const_position <- unname(terra::values(res$anomalies, mat = TRUE)[1, 1])
  expect_true(is.na(val_const_position))
  expect_false(is.nan(val_const_position))
})

test_that("compute_anomalies errors on cycle < 2", {
  r <- .make_pattern_raster(c(1, 2, 3, 4))
  expect_error(compute_anomalies(r, cycle = 1, verbose = FALSE), "cycle.*>= 2")
})

test_that("compute_anomalies errors when there are fewer layers than the cycle length", {
  r <- .make_pattern_raster(c(1, 2, 3))
  expect_error(compute_anomalies(r, cycle = 12, verbose = FALSE), "Fewer layers")
})

test_that("compute_anomalies errors on non-SpatRaster input", {
  expect_error(compute_anomalies(matrix(1:4, 2, 2), verbose = FALSE), "SpatRaster")
})

test_that("compute_anomalies messages about a partial final cycle when verbose = TRUE", {
  r <- .make_pattern_raster(c(10, 20, 30, 10))
  expect_message(compute_anomalies(r, cycle = 3, verbose = TRUE), "not a multiple of cycle")
})

test_that("compute_anomalies stays silent about partial cycles when verbose = FALSE", {
  r <- .make_pattern_raster(c(10, 20, 30, 10))
  expect_silent(compute_anomalies(r, cycle = 3, verbose = FALSE))
})

test_that("compute_anomalies messages the standardising step when verbose = TRUE", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 12, seed = 1)$series
  expect_message(
    compute_anomalies(r, cycle = 12, standardise = TRUE, verbose = TRUE),
    "Standardising"
  )
})
