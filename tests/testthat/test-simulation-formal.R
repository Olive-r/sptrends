test_that("formal covariance functions satisfy their analytical contracts", {
  distance <- c(0, 1, 2)
  expect_equal(
    sptrends:::.simulation_correlation(
      distance, "exponential", rho = 0.6, range = 1, smoothness = 0.5),
    c(1, 0.6, 0.6^2))
  expect_equal(
    sptrends:::.simulation_correlation(
      distance, "gaussian", rho = 0.6, range = 1, smoothness = 0.5),
    c(1, 0.6, 0.6^4))
  expect_equal(
    sptrends:::.simulation_correlation(
      distance, "matern", rho = 0.6, range = 2, smoothness = 0.5),
    exp(-distance / 2), tolerance = 1e-12)

  distance_matrix <- matrix(distance, nrow = 1L)
  matrix_result <- sptrends:::.simulation_correlation(
    distance_matrix, "matern", rho = 0.6, range = 2,
    smoothness = 0.5)
  expect_identical(dim(matrix_result), dim(distance_matrix))
})

test_that("Matérn correlation agrees with the independent fields package", {
  skip_if_not_installed("fields")
  distance <- seq(0, 4, length.out = 21)
  for (smoothness in c(0.5, 1, 1.5, 2.5)) {
    observed <- sptrends:::.simulation_correlation(
      distance, "matern", rho = 0.5, range = 1.7,
      smoothness = smoothness)
    reference <- fields::Matern(
      distance, aRange = 1.7, smoothness = smoothness)
    expect_equal(observed, reference, tolerance = 1e-10)
  }
})

test_that("formal fields reproduce target variance and lag correlation", {
  set.seed(23)
  fields <- replicate(
    300,
    sptrends:::.simulate_gaussian_field(
      12, 12, "exponential", rho = 0.7, range = 1,
      smoothness = 0.5, sd = 2),
    simplify = FALSE)
  centre <- vapply(fields, function(x) x[6, 6], numeric(1))
  neighbour <- vapply(fields, function(x) x[6, 7], numeric(1))
  expect_equal(stats::var(centre), 4, tolerance = 0.7)
  expect_equal(stats::cor(centre, neighbour), 0.7, tolerance = 0.1)
})

test_that("exact Gaussian fields honour dimensions and finite covariance", {
  set.seed(71)
  field <- sptrends:::.simulate_gaussian_field_exact(
    4, 5, "matern", rho = 0, range = 1.5,
    smoothness = 1, sd = 2)
  expect_identical(dim(field), c(4L, 5L))
  expect_true(all(is.finite(field)))
})

test_that("exact signal geometries and truth rasters are internally coherent", {
  sim <- sim_trend_stack(
    nrow = 20, ncol = 20, n_time = 5,
    trend_shape = "ellipse", signal_size = c(8, 12),
    signal_axis_ratio = 0.5, trend_fraction = 1,
    constant_block = FALSE, spatial_model = "independent", seed = 9)
  slope <- terra::values(sim$true_slope, mat = FALSE)
  signal <- terra::values(sim$true_signal, mat = FALSE)
  direction <- terra::values(sim$true_direction, mat = FALSE)
  expect_identical(signal, as.numeric(slope != 0))
  expect_identical(direction, sign(slope))
  expect_true(sum(signal) > 0)
  expect_true(sum(signal) < length(signal))
  expect_identical(sim$parameters$spatial_model, "independent")
})

test_that("centred square signal_size is the exact requested cell count", {
  sim <- sim_trend_stack(
    nrow = 100, ncol = 100, n_time = 2, trend_shape = "square",
    signal_size = 5, trend_fraction = 1, constant_block = FALSE,
    spatial_model = "independent", seed = 8)
  expect_equal(sum(terra::values(sim$true_signal, mat = FALSE)), 25)
  nonzero <- terra::values(sim$true_slope, mat = FALSE)
  expect_true(all(nonzero[nonzero != 0] > 0))
})

test_that("custom signal masks are preserved exactly", {
  mask <- matrix(FALSE, 6, 7)
  mask[2:4, 3:6] <- TRUE
  sim <- sim_trend_stack(
    nrow = 6, ncol = 7, n_time = 4, trend_shape = "custom",
    custom_mask = mask, trend_fraction = 1, constant_block = FALSE,
    spatial_model = "independent", seed = 10)
  expect_identical(
    terra::values(sim$true_signal, mat = FALSE),
    as.numeric(as.vector(t(mask))))
})

test_that("random signal placement and orientation are reproducible", {
  arguments <- list(
    nrow = 12, ncol = 14, n_time = 3, trend_shape = "ellipse",
    signal_size = c(5, 8), signal_location = "random",
    signal_angle = "random", trend_fraction = 1,
    constant_block = FALSE, spatial_model = "independent", seed = 91)
  first <- do.call(sim_trend_stack, arguments)
  second <- do.call(sim_trend_stack, arguments)
  expect_identical(
    terra::values(first$true_signal, mat = FALSE),
    terra::values(second$true_signal, mat = FALSE))
  expect_identical(
    terra::values(first$series, mat = TRUE),
    terra::values(second$series, mat = TRUE))
  expect_identical(first$parameters, second$parameters)
})

test_that("simulation objects retain the complete generating specification", {
  mask <- matrix(FALSE, 5, 6)
  mask[2:4, 2:5] <- TRUE
  simulation <- sim_trend_stack(
    nrow = 5, ncol = 6, n_time = 8,
    trend_shape = "custom", custom_mask = mask,
    trend_fraction = 1, ar1 = 0.4, noise_sd = 2,
    noise_dist = "t", t_df = 5,
    spatial_model = "matern", spatial_rho = 0,
    spatial_range = 2, spatial_smoothness = 1.5,
    signal_location = "centre", signal_angle = 20,
    signal_axis_ratio = 0.7, constant_block = FALSE,
    break_type = "mean", break_time = 4,
    break_fraction = 1, break_magnitude = 1.2, seed = 19)

  expected <- c(
    "nrow", "ncol", "n_time", "trend_strength", "trend_shape",
    "trend_fraction", "ar1", "noise_sd", "noise_dist", "t_df",
    "smooth_radius", "spatial_model", "spatial_rho", "spatial_range",
    "spatial_smoothness", "signal_size", "signal_location",
    "signal_angle", "signal_axis_ratio", "custom_mask",
    "constant_block", "break_type", "break_time", "break_fraction",
    "break_magnitude", "seed")
  expect_identical(names(simulation$parameters), expected)
  expect_identical(simulation$parameters$custom_mask, mask)
  expect_equal(simulation$parameters$break_time, 4)
})

test_that("formal spatial parameters reject impossible values", {
  expect_error(
    sim_trend_stack(spatial_model = "exponential", spatial_rho = 0),
    "must be in")
  expect_error(sim_trend_stack(spatial_range = 0), "spatial_range")
  expect_error(sim_trend_stack(spatial_smoothness = 0),
               "spatial_smoothness")
  expect_error(sim_trend_stack(signal_size = c(2, 3, 4)), "signal_size")
  expect_error(sim_trend_stack(signal_axis_ratio = 0), "signal_axis_ratio")
  expect_error(sim_trend_stack(signal_angle = "fixed"), "signal_angle")
})

test_that("Matérn covariance is controlled by range and smoothness", {
  expect_error(
    sim_trend_stack(
      nrow = 4, ncol = 4, n_time = 2,
      spatial_model = "matern", spatial_rho = 0,
      spatial_range = 1.5, spatial_smoothness = 1,
      constant_block = FALSE, seed = 8),
    NA)
})
