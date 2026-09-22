test_that("formal simulation helpers cover defensive covariance paths", {
  independent <- sptrends:::.simulation_correlation(
    c(0, 1), "independent", rho = 0, range = 1, smoothness = 0.5
  )
  expect_identical(independent, c(1, 0))

  expect_error(
    sptrends:::.simulate_spatial_innovation(
      2, 2, "independent", rho = 0, range = 1,
      smoothness = 0.5, sd = 1, dist = "t", df = 2
    ),
    "must be > 2"
  )
  innovation <- sptrends:::.simulate_spatial_innovation(
    2, 2, "independent", rho = 0, range = 1,
    smoothness = 0.5, sd = 2, dist = "t", df = 4
  )
  expect_identical(dim(innovation), c(2L, 2L))
  expect_true(all(is.finite(innovation)))

  expect_error(sptrends:::.sim_noise(3, 1, "t", 2), "must be > 2")
  heavy_tailed <- sptrends:::.sim_noise(20, 1, "t", 4)
  expect_length(heavy_tailed, 20L)
  expect_true(all(is.finite(heavy_tailed)))
})

test_that("exact covariance simulation rejects a non-positive matrix", {
  testthat::local_mocked_bindings(
    .simulation_correlation = function(...) {
      matrix(c(1, 2, 2, 1), 2, 2)
    },
    .package = "sptrends"
  )
  expect_error(
    sptrends:::.simulate_gaussian_field_exact(
      1, 2, "matern", rho = 0, range = 1,
      smoothness = 1, sd = 1
    ),
    "not positive semidefinite"
  )
})

test_that("large invalid embeddings fail before exact decomposition", {
  testthat::local_mocked_bindings(
    .simulation_correlation = function(distance, ...) {
      result <- matrix(0, nrow(distance), ncol(distance))
      result[1, 2] <- 1
      result
    },
    .package = "sptrends"
  )
  expect_error(
    sptrends:::.simulate_gaussian_field(
      51, 50, "matern", rho = 0, range = 1,
      smoothness = 1, sd = 1, max_expansions = 0
    ),
    "too large for the exact covariance fallback"
  )
})

test_that("small invalid embeddings use the exact covariance fallback", {
  calls <- 0L
  testthat::local_mocked_bindings(
    .simulation_correlation = function(distance, ...) {
      calls <<- calls + 1L
      if (calls == 1L) {
        result <- matrix(0, nrow(distance), ncol(distance))
        result[1, 2] <- 1
        return(result)
      }
      diag(nrow(distance))
    },
    .package = "sptrends"
  )
  result <- sptrends:::.simulate_gaussian_field(
    2, 2, "matern", rho = 0, range = 1,
    smoothness = 1, sd = 1, max_expansions = 0
  )
  expect_identical(dim(result), c(2L, 2L))
  expect_equal(calls, 2L)
})

test_that("custom and rotated signal masks cover every input form", {
  raster_mask <- terra::rast(nrows = 3, ncols = 4)
  terra::values(raster_mask) <- c(rep(1, 4), rep(0, 8))
  raster_result <- sptrends:::.simulation_signal_mask(
    3, 4, "custom", 1, "centre", 0, 1, raster_mask
  )
  expect_equal(sum(raster_result), 4)

  expect_error(
    sptrends:::.simulation_signal_mask(
      3, 4, "custom", 1, "centre", 0, 1, NULL
    ),
    "custom_mask.*required"
  )
  expect_error(
    sptrends:::.simulation_signal_mask(
      3, 4, "custom", 1, "centre", 0, 1, matrix(1, 2, 2)
    ),
    "dimensions"
  )
  expect_error(
    sptrends:::.simulation_signal_mask(
      3, 4, "custom", 1, "centre", 0, 1, 1:5
    ),
    "exactly"
  )

  rotated <- sptrends:::.simulation_signal_mask(
    9, 9, "rectangle", c(3, 7), "centre", 45, 1
  )
  expect_true(any(rotated))
  expect_true(any(!rotated))
  expect_true(all(sptrends:::.simulation_signal_mask(
    2, 3, "radial", 1, "centre", 0, 1
  )))
})

test_that("formal simulator validation covers malformed controls", {
  expect_error(sim_trend_stack(smooth_radius = -1), "non-negative integer")
  expect_error(sim_trend_stack(trend_fraction = -0.1), "trend_fraction")
  expect_error(sim_trend_stack(break_fraction = 1.1), "break_fraction")
  expect_error(sim_trend_stack(spatial_rho = NA_real_), "spatial_rho")
  expect_error(sim_trend_stack(signal_angle = Inf), "signal_angle")
  expect_error(sim_trend_stack(signal_axis_ratio = NA_real_),
               "signal_axis_ratio")
  expect_error(sim_trend_stack(signal_size = "large"), "signal_size")
  expect_error(
    sim_trend_stack(
      n_time = 4, break_type = "mean", break_time = 4,
      constant_block = FALSE
    ),
    "break_time"
  )
})

test_that("simulator reports progress, remaining time and duration", {
  output <- capture.output(
    sim_trend_stack(
      nrow = 2, ncol = 2, n_time = 3,
      trend_fraction = 0, constant_block = FALSE,
      spatial_model = "independent", seed = 1, verbose = TRUE
    ),
    type = "message"
  )
  text <- paste(output, collapse = " ")
  expect_match(text, "progress:", fixed = TRUE)
  expect_match(text, "remaining:", fixed = TRUE)
  expect_match(text, "elapsed:", fixed = TRUE)
})
