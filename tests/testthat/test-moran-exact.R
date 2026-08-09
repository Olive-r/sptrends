.make_grid <- function(values, nrow, ncol) {
  r <- terra::rast(nrows = nrow, ncols = ncol, xmin = 0, xmax = ncol, ymin = 0, ymax = nrow)
  terra::setValues(r, values)
}

test_that("spatial_autocorrelation computes the exact expected I under queen connectivity", {
  # 2x2 "checkerboard": A=1 (top-left), B=10 (top-right), C=10 (bottom-left), D=1 (bottom-right).
  # Hand-computed: I_queen = -1/3 (complete graph, all 4 cells mutually adjacent).
  r <- .make_grid(c(1, 10, 10, 1), nrow = 2, ncol = 2)
  res <- spatial_autocorrelation(r, connectivity = "queen", nperm = 9, seed = 1, verbose = FALSE, report = FALSE)
  expect_equal(res$statistic, -1 / 3, tolerance = 1e-8)
  expect_identical(res$sign, "negative")
})

test_that("spatial_autocorrelation computes the exact expected I under rook connectivity", {
  # Same checkerboard: I_rook = -1 exactly (only edge-adjacent pairs count,
  # every rook-neighbour pair is maximally dissimilar here).
  r <- .make_grid(c(1, 10, 10, 1), nrow = 2, ncol = 2)
  res <- spatial_autocorrelation(r, connectivity = "rook", nperm = 9, seed = 1, verbose = FALSE, report = FALSE)
  expect_equal(res$statistic, -1, tolerance = 1e-8)
})

test_that("queen and rook give genuinely different results on the same data", {
  r <- .make_grid(c(1, 10, 10, 1), nrow = 2, ncol = 2)
  res_queen <- spatial_autocorrelation(r, connectivity = "queen", nperm = 9, seed = 1, verbose = FALSE, report = FALSE)
  res_rook  <- spatial_autocorrelation(r, connectivity = "rook", nperm = 9, seed = 1, verbose = FALSE, report = FALSE)
  expect_false(isTRUE(all.equal(res_queen$statistic, res_rook$statistic)))
})

test_that("spatial_autocorrelation computes the exact expected I for a simple positive-autocorrelation chain", {
  # 1x4 row, values 1,1,10,10 -- similar values adjacent, hand-computed I = 1/3.
  r <- .make_grid(c(1, 1, 10, 10), nrow = 1, ncol = 4)
  res <- spatial_autocorrelation(r, nperm = 9, seed = 1, verbose = FALSE, report = FALSE)
  expect_equal(res$statistic, 1 / 3, tolerance = 1e-8)
  expect_identical(res$sign, "positive")
})

test_that("spatial_autocorrelation's alternative argument behaves sensibly for a strongly positive I", {
  r <- .make_grid(c(1, 1, 1, 1, 1, 1, 10, 10, 10), nrow = 3, ncol = 3)
  res_greater <- spatial_autocorrelation(r, alternative = "greater", nperm = 99, seed = 1, verbose = FALSE, report = FALSE)
  res_less    <- spatial_autocorrelation(r, alternative = "less", nperm = 99, seed = 1, verbose = FALSE, report = FALSE)
  res_two     <- spatial_autocorrelation(r, alternative = "two.sided", nperm = 99, seed = 1, verbose = FALSE, report = FALSE)

  # same observed I regardless of which tail is tested
  expect_equal(res_greater$statistic, res_less$statistic)
  expect_equal(res_greater$statistic, res_two$statistic)
  # a strongly positive I should look significant against "greater" but not "less"
  expect_lt(res_greater$p, res_less$p)
  expect_true(all(c(res_greater$p, res_less$p, res_two$p) >= 0 & c(res_greater$p, res_less$p, res_two$p) <= 1))
})

test_that("spatial_autocorrelation errors with fewer than 3 valid cells", {
  r <- .make_grid(c(1, NA, NA, NA), nrow = 2, ncol = 2)
  expect_error(spatial_autocorrelation(r, verbose = FALSE), "Too few valid cells")
})

test_that("spatial_autocorrelation validates permutation controls", {
  r <- .make_grid(1:9, nrow = 3, ncol = 3)
  for (bad in list(0, -1, 1.5, NA_real_, Inf, c(9, 19))) {
    expect_error(
      spatial_autocorrelation(r, nperm = bad, verbose = FALSE,
                               report = FALSE),
      "'nperm'"
    )
  }
  for (bad in list(0, -1, 1.5, NA_real_, Inf, c(1, 2))) {
    expect_error(
      spatial_autocorrelation(r, n_cores = bad, verbose = FALSE,
                               report = FALSE),
      "'n_cores'"
    )
  }
  expect_error(
    spatial_autocorrelation(r, seed = c(1, 2), verbose = FALSE,
                             report = FALSE),
    "'seed'"
  )
  expect_error(
    spatial_autocorrelation(
      r, precomputed_neighbourhood = list(W = matrix(1, 2, 2)),
      verbose = FALSE, report = FALSE
    ),
    "must be an object returned by"
  )
})

test_that("spatial_autocorrelation errors when no valid cell has a valid neighbour", {
  # A 1x5 row with only cells 1, 3, 5 valid -- none of them are adjacent to
  # each other (cells 2 and 4, their only possible neighbours, are NA).
  r <- .make_grid(c(1, NA, 5, NA, 9), nrow = 1, ncol = 5)
  expect_error(spatial_autocorrelation(r, verbose = FALSE), "No valid cell has valid neighbours")
})

test_that("spatial_autocorrelation excludes NA cells from N correctly", {
  r <- .make_grid(c(1, 10, 10, NA), nrow = 2, ncol = 2)
  res <- spatial_autocorrelation(r, nperm = 9, seed = 1, verbose = FALSE, report = FALSE)
  expect_identical(res$N, 3L)
})

test_that("classify_moran labels the three categories correctly, including boundaries", {
  expect_identical(suppressMessages(classify_moran(0.05)), "low")
  expect_identical(suppressMessages(classify_moran(0.1)), "moderate")   # boundary: not < 0.1
  expect_identical(suppressMessages(classify_moran(0.29)), "moderate")
  expect_identical(suppressMessages(classify_moran(0.3)), "strong")     # boundary: not < 0.3
  expect_identical(suppressMessages(classify_moran(-0.5)), "strong")    # uses abs(I)
})

test_that("spatial_autocorrelation_summary returns the expected fields and can write a CSV", {
  r <- .make_grid(c(1, 10, 10, 1), nrow = 2, ncol = 2)
  res <- spatial_autocorrelation(r, nperm = 19, seed = 1, verbose = FALSE, report = FALSE)
  tab <- suppressMessages(spatial_autocorrelation_summary(res))

  expect_s3_class(tab, "data.frame")
  expect_true(all(c("N", "I", "sign", "p_value", "null_mean", "null_sd") %in% tab$metric))

  path <- tempfile(fileext = ".csv")
  suppressMessages(spatial_autocorrelation_summary(res, path = path))
  expect_true(file.exists(path))
  unlink(path)
})

test_that("spatial_autocorrelation_null_plot runs without error, including the out-of-range annotation branch", {
  r <- .make_grid(c(1, 10, 10, 1), nrow = 2, ncol = 2)
  res <- spatial_autocorrelation(r, nperm = 19, seed = 1, verbose = FALSE, report = FALSE)
  before <- graphics::par("mar")
  expect_silent(spatial_autocorrelation_null_plot(res))
  expect_equal(graphics::par("mar"), before)

  # force the "observed I outside the null histogram range" branch
  res_fake <- res
  res_fake$null_dist <- rep(0, 19)
  res_fake$statistic <- 5
  expect_error(suppressWarnings(spatial_autocorrelation_null_plot(res_fake)), NA)
})

test_that("spatial_autocorrelation errors on non-SpatRaster input", {
  expect_error(spatial_autocorrelation(matrix(1:9, 3, 3)), "SpatRaster")
})

test_that("spatial_autocorrelation messages the parallel line when verbose = TRUE and n_cores > 1", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 1, seed = 1)$series[[1]]
  expect_message(
    spatial_autocorrelation(r, nperm = 9, n_cores = 2, verbose = TRUE, report = FALSE),
    "Parallel permutations"
  )
})

test_that("spatial_autocorrelation_summary notes when p is at the mathematical floor", {
  # a strongly clustered raster with very few permutations makes the
  # observed I more extreme than every single permutation very likely.
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 1, smooth_radius = 4, seed = 1)$series[[1]]
  res <- spatial_autocorrelation(r, nperm = 9, seed = 1, verbose = FALSE, report = FALSE)
  skip_if(res$p != 1 / (res$nperm + 1), "observed I was not at the floor with this seed -- try another")
  expect_message(spatial_autocorrelation_summary(res), "mathematical floor")
})

test_that("spatial_autocorrelation_null_plot can write a PNG via path", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 1, seed = 1)$series[[1]]
  res <- spatial_autocorrelation(r, nperm = 19, seed = 1, verbose = FALSE, report = FALSE)
  path <- tempfile(fileext = ".png")
  spatial_autocorrelation_null_plot(res, path = path)
  expect_true(file.exists(path))
  unlink(path)
})

test_that("spatial_autocorrelation(method = 'getis_ord') computes the exact expected G, matching a hand-computed reference", {
  # Same 2x2 queen-connectivity checkerboard as the Moran exact tests
  # above -- independently verified in Python (direct double-sum
  # definition of General G, Getis & Ord 1992) before writing this
  # test: G = 1.0 for this exact configuration.
  r <- .make_grid(c(1, 10, 10, 1), nrow = 2, ncol = 2)
  res <- spatial_autocorrelation(r, method = "getis_ord", connectivity = "queen",
                                  nperm = 9, seed = 1, verbose = FALSE,
                                  report = FALSE)
  expect_equal(res$statistic, 1.0, tolerance = 1e-8)
  expect_identical(res$method, "getis_ord")
  expect_true(is.na(res$sign))  # direction is not a General G concept
})

test_that("spatial_autocorrelation(method = 'getis_ord') rejects negative values with an informative error", {
  r <- .make_grid(c(1, -5, 10, 1), nrow = 2, ncol = 2)
  expect_error(
    spatial_autocorrelation(r, method = "getis_ord", verbose = FALSE),
    "non-negative"
  )
})

test_that("moran and getis_ord give genuinely different statistics on the same data (not the same number relabelled)", {
  r <- .make_grid(c(1, 10, 10, 1), nrow = 2, ncol = 2)
  res_moran <- spatial_autocorrelation(r, method = "moran", nperm = 9, seed = 1,
                                        verbose = FALSE, report = FALSE)
  res_getis <- spatial_autocorrelation(r, method = "getis_ord", nperm = 9,
                                        seed = 1, verbose = FALSE, report = FALSE)
  expect_false(isTRUE(all.equal(res_moran$statistic, res_getis$statistic)))
})

test_that("print()/summary()/plot() work on a getis_ord result the same way they do for moran", {
  r <- .make_grid(c(1, 10, 10, 1), nrow = 2, ncol = 2)
  res <- spatial_autocorrelation(r, method = "getis_ord", nperm = 9, seed = 1,
                                  verbose = FALSE, report = FALSE)
  expect_error(print(res), NA)
  expect_error(suppressMessages(summary(res)), NA)
  expect_error(plot(res), NA)
})

test_that("spatial_autocorrelation(method = 'moran') errors clearly on a perfectly constant raster (zero variance, undefined denominator)", {
  r <- terra::rast(nrows = 5, ncols = 5, vals = 3)  # every cell: exactly 3
  expect_error(
    spatial_autocorrelation(r, method = "moran", verbose = FALSE),
    "zero variance"
  )
})

test_that("spatial_autocorrelation(method = 'getis_ord') errors clearly when every cell is exactly zero (undefined denominator, a different degenerate case from negative values)", {
  r <- terra::rast(nrows = 5, ncols = 5, vals = 0)  # every cell: exactly 0
  expect_error(
    spatial_autocorrelation(r, method = "getis_ord", verbose = FALSE),
    "exactly zero everywhere"
  )
})

test_that("spatial_autocorrelation(method = 'getis_ord') does NOT error on a constant *positive* raster (only all-zero is degenerate, not any constant)", {
  r <- terra::rast(nrows = 5, ncols = 5, vals = 3)  # every cell: exactly 3
  expect_error(
    spatial_autocorrelation(r, method = "getis_ord", nperm = 9, seed = 1,
                             verbose = FALSE, report = FALSE),
    NA
  )
})

test_that("spatial_autocorrelation(method = 'getis_ord', verbose = TRUE) prints its own 'Observed G' message", {
  r <- abs(sim_trend_stack(nrow = 6, ncol = 6, n_time = 1, seed = 1)$series[[1]])
  expect_message(
    spatial_autocorrelation(r, method = "getis_ord", nperm = 9, seed = 1,
                             verbose = TRUE, report = FALSE),
    "Observed G"
  )
})

test_that("the shared .prepare_spatial_neighbourhood() helper gives results identical to the two independent implementations it replaced", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 6, seed = 60)$series
  X <- terra::values(r, mat = TRUE)
  ok <- stats::complete.cases(X)

  nb_shared <- sptrends:::.prepare_spatial_neighbourhood(r, ok)
  nb_cmk <- prepare_cmk_neighbourhood(r, ok)
  nb_moran <- sptrends:::.prepare_moran_neighbourhood(r, ok)

  expect_equal(as.matrix(nb_shared$W), as.matrix(nb_cmk$W))
  expect_equal(as.matrix(nb_shared$W), as.matrix(nb_moran$W))
  expect_equal(nb_shared$nb_count, nb_cmk$nb_count)
  expect_equal(nb_shared$S0, nb_moran$S0)
  expect_equal(nb_shared$S0, sum(nb_shared$nb_count))
})

test_that("spatial_autocorrelation(precomputed_neighbourhood = ...) gives identical results whether it is skipped, or supplied from either prepare_cmk_neighbourhood() or the shared internal helper", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 1, seed = 61)$series[[1]]
  r_full <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 6,
                             seed = 61)$series
  X <- terra::values(r_full, mat = TRUE)
  ok <- stats::complete.cases(X)

  result_default <- spatial_autocorrelation(r, nperm = 49, seed = 1,
                                             verbose = FALSE, report = FALSE)

  nb_cmk <- prepare_cmk_neighbourhood(r_full, ok)
  result_from_cmk <- spatial_autocorrelation(
    r, precomputed_neighbourhood = nb_cmk, nperm = 49, seed = 1,
    verbose = FALSE, report = FALSE)

  nb_shared <- sptrends:::.prepare_spatial_neighbourhood(r_full, ok)
  result_from_shared <- spatial_autocorrelation(
    r, precomputed_neighbourhood = nb_shared, nperm = 49, seed = 1,
    verbose = FALSE, report = FALSE)

  expect_equal(result_default$statistic, result_from_cmk$statistic)
  expect_equal(result_default$statistic, result_from_shared$statistic)
  expect_equal(result_default$p, result_from_cmk$p)
  expect_equal(result_default$p, result_from_shared$p)
})
