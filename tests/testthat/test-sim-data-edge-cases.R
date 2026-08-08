test_that("sim_trend_stack works with n_time = 1 (regression test)", {
  sim <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 1, seed = 1)
  expect_equal(terra::nlyr(sim$series), 1)
  expect_false(anyNA(terra::values(sim$series, mat = FALSE)))
})

test_that("sim_trend_stack works with n_time = 2", {
  sim <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 2, seed = 1)
  expect_equal(terra::nlyr(sim$series), 2)
})

test_that("trend_test errors clearly on a single-layer raster", {
  sim <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 1, seed = 1)
  expect_error(trend_test(sim$series, verbose = FALSE), "at least 2 layers")
})

test_that("prewhiten errors clearly on too few layers", {
  sim <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 2, seed = 1)
  expect_error(prewhiten(sim$series, verbose = FALSE), "at least 3 layers")
})

test_that("slope_estimator errors clearly on a single-layer raster", {
  sim <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 1, seed = 1)
  expect_error(slope_estimator(sim$series, verbose = FALSE, report = FALSE)$slope, "at least 2 layers")
})

test_that("sim_trend_stack returns series and true_slope with matching geometry", {
  sim <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 5, seed = 1)
  expect_s4_class(sim$series, "SpatRaster")
  expect_s4_class(sim$true_slope, "SpatRaster")
  expect_equal(terra::nlyr(sim$true_slope), 1)
  expect_equal(terra::ncell(sim$series), terra::ncell(sim$true_slope))
})

test_that("trend_fraction = 0 gives a complete null field (true_slope all zero)", {
  sim <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 5, trend_fraction = 0, seed = 1)
  expect_true(all(terra::values(sim$true_slope, mat = FALSE) == 0))
})

test_that("trend_fraction = 1 gives a non-zero true slope everywhere trend_shape assigns one", {
  sim <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 5, trend_fraction = 1,
                          trend_shape = "block", constant_block = FALSE, seed = 1)
  vals <- terra::values(sim$true_slope, mat = FALSE)
  expect_true(any(vals != 0))
})

test_that("trend_shape = 'gradient' produces a left-right varying true slope", {
  sim <- sim_trend_stack(nrow = 6, ncol = 12, n_time = 3, trend_shape = "gradient",
                          trend_fraction = 1, constant_block = FALSE, seed = 1)
  m <- matrix(terra::values(sim$true_slope, mat = FALSE), nrow = 6, byrow = TRUE)
  expect_lt(mean(m[, 1]), mean(m[, ncol(m)]))
})

test_that("noise_sd scales the noise magnitude", {
  sim_low <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 20, noise_sd = 0.1,
                              smooth_radius = 0, seed = 1)
  sim_high <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 20, noise_sd = 5,
                               smooth_radius = 0, seed = 1)
  sd_low <- stats::sd(terra::values(sim_low$series, mat = FALSE))
  sd_high <- stats::sd(terra::values(sim_high$series, mat = FALSE))
  expect_gt(sd_high, sd_low)
})

test_that("sim_trend_stack rejects an out-of-range trend_fraction", {
  expect_error(sim_trend_stack(trend_fraction = 1.5), "trend_fraction")
  expect_error(sim_trend_stack(trend_fraction = -0.1), "trend_fraction")
})

test_that("noise_dist = 't' gives heavier tails than 'gaussian' at the same noise_sd", {
  sim_g <- sim_trend_stack(nrow = 20, ncol = 20, n_time = 30, noise_dist = "gaussian",
                            noise_sd = 1, smooth_radius = 0, trend_fraction = 0, seed = 1)
  sim_t <- sim_trend_stack(nrow = 20, ncol = 20, n_time = 30, noise_dist = "t", t_df = 3,
                            noise_sd = 1, smooth_radius = 0, trend_fraction = 0, seed = 1)

  vals_g <- terra::values(sim_g$series, mat = FALSE) - 10
  vals_t <- terra::values(sim_t$series, mat = FALSE) - 10

  # similar standard deviation (both rescaled to noise_sd = 1) ...
  expect_equal(stats::sd(vals_g), stats::sd(vals_t), tolerance = 0.15)
  # ... but the t-distributed noise has a heavier tail (higher kurtosis proxy).
  kurtosis_proxy <- function(x) mean((x - mean(x))^4) / stats::sd(x)^4
  expect_gt(kurtosis_proxy(vals_t), kurtosis_proxy(vals_g))
})

test_that("sim_trend_stack rejects t_df <= 2", {
  expect_error(sim_trend_stack(noise_dist = "t", t_df = 2), "t_df")
  expect_error(sim_trend_stack(noise_dist = "t", t_df = 1), "t_df")
})

test_that("sim_trend_stack(break_type = 'none') is unchanged: true_break is all 0, break_time is NULL", {
  sim <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 8, seed = 300)
  vals <- terra::values(sim$true_break, mat = FALSE)
  expect_true(all(vals == 0))
  expect_null(sim$break_time)
})

test_that("sim_trend_stack rejects an out-of-range break_fraction", {
  expect_error(
    sim_trend_stack(nrow = 8, ncol = 8, n_time = 6, break_type = "mean",
                     break_fraction = 1.5),
    "break_fraction"
  )
})

test_that("sim_trend_stack rejects a break_time outside [1, n_time - 1]", {
  expect_error(
    sim_trend_stack(nrow = 8, ncol = 8, n_time = 6, break_type = "mean",
                     break_time = 6, seed = 1),
    "break_time"
  )
  expect_error(
    sim_trend_stack(nrow = 8, ncol = 8, n_time = 6, break_type = "mean",
                     break_time = 0, seed = 1),
    "break_time"
  )
})

test_that("sim_trend_stack(break_type = 'mean') gives an exact step of break_magnitude at break_time, matching a hand-computed reference", {
  sim <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 10,
                          trend_fraction = 0, break_type = "mean",
                          break_time = 5, break_fraction = 1,
                          break_magnitude = 3, noise_sd = 0, ar1 = 0,
                          smooth_radius = 0, constant_block = FALSE,
                          seed = 301)
  X <- terra::values(sim$series, mat = TRUE)
  # Every cell has a break (break_fraction = 1, no constant_block to
  # exclude any); no noise, no trend -- the only thing that can differ
  # between step 5 and step 6 is the break itself.
  step_size <- X[, 6] - X[, 5]
  expect_true(all(abs(step_size - 3) < 1e-8))
  expect_equal(sim$break_time, 5)
})

test_that("sim_trend_stack(break_type = 'slope') gives the correct pre/post slope, matching a hand-computed reference", {
  sim <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 10,
                          trend_fraction = 0, break_type = "slope",
                          break_time = 5, break_fraction = 1,
                          break_magnitude = 2, noise_sd = 0, ar1 = 0,
                          smooth_radius = 0, constant_block = FALSE,
                          seed = 302)
  X <- terra::values(sim$series, mat = TRUE)
  t <- 1:10
  # unname(): terra::values(..., mat = TRUE) retains the raster's own
  # layer names as column names on the returned matrix -- a scalar
  # extraction like X[1, 5] carries that name through into any
  # arithmetic done on it, which expect_equal() then treats as a
  # mismatch against an unnamed reference value even when the numbers
  # themselves agree.
  slope_before <- unname((X[1, 5] - X[1, 1]) / (5 - 1))
  slope_after <- unname((X[1, 10] - X[1, 6]) / (10 - 6))
  # trend_fraction = 0 -> base slope is 0 for every cell; the only
  # slope present is the break's own break_magnitude, kicking in after
  # break_time.
  expect_equal(slope_before, 0, tolerance = 1e-8)
  expect_equal(slope_after, 2, tolerance = 1e-8)
})

test_that("a cell can have both a trend and a break at the same time (composition, not mutual exclusion)", {
  sim <- sim_trend_stack(nrow = 12, ncol = 12, n_time = 10,
                          trend_fraction = 1, break_type = "mean",
                          break_fraction = 1, break_time = 5,
                          break_magnitude = 3, seed = 303)
  slope_vals <- terra::values(sim$true_slope, mat = FALSE)
  break_vals <- terra::values(sim$true_break, mat = FALSE)
  # With both fractions at 1 (minus the small constant_block corner),
  # the overwhelming majority of cells should show up as having BOTH a
  # nonzero true slope and a true break -- if they were mutually
  # exclusive, this overlap would be zero or near it instead.
  both <- sum(slope_vals != 0 & break_vals == 1, na.rm = TRUE)
  expect_true(both > 0.5 * length(slope_vals))
})

test_that("break_fraction rounding down to zero active blocks gives true_break all 0 and break_time NULL, not an error", {
  sim <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, break_type = "mean",
                          break_fraction = 0.001, seed = 304)
  vals <- terra::values(sim$true_break, mat = FALSE)
  expect_true(all(vals == 0))
  expect_null(sim$break_time)
})

test_that("break_time = NULL defaults to round(n_time / 2), reflected in the returned value", {
  sim <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 11, break_type = "mean",
                          break_fraction = 1, seed = 305)
  expect_equal(sim$break_time, round(11 / 2))
})

test_that("cells with no break are completely unaffected by break_magnitude, regardless of break_type", {
  sim <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 10,
                          trend_fraction = 0, break_type = "mean",
                          break_fraction = 0.3, break_magnitude = 100,
                          noise_sd = 0, ar1 = 0, smooth_radius = 0,
                          constant_block = FALSE, seed = 306)
  X <- terra::values(sim$series, mat = TRUE)
  break_vals <- terra::values(sim$true_break, mat = FALSE)
  no_break_rows <- which(break_vals == 0)
  expect_true(length(no_break_rows) > 0)  # sanity: some cells have no break
  # With no trend, no noise, and no break, an unaffected cell's series
  # must be constant across all 10 time steps -- not shifted by even a
  # fraction of break_magnitude = 100.
  for (i in no_break_rows[seq_len(min(5, length(no_break_rows)))]) {
    expect_equal(diff(range(X[i, ])), 0, tolerance = 1e-8)
  }
})

test_that("constant_block cells never show a break, even when randomly selected for one", {
  # Try several seeds -- constant_block always overrides whatever
  # break_fraction's random block draw happened to pick for that corner.
  for (s in 1:5) {
    sim <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 8, break_type = "mean",
                            break_fraction = 1, constant_block = TRUE, seed = s)
    corner_cells <- which(
      (rep(seq_len(8), each = 8) <= 2) & (rep(seq_len(8), times = 8) <= 2)
    )
    break_vals <- terra::values(sim$true_break, mat = FALSE)
    expect_true(all(break_vals[corner_cells] == 0))
  }
})

test_that("true_break's total count roughly matches break_fraction of the spatial blocks (not exact due to rounding, but in the right ballpark)", {
  sim <- sim_trend_stack(nrow = 20, ncol = 20, n_time = 8, break_type = "mean",
                          break_fraction = 0.5, constant_block = FALSE, seed = 307)
  break_vals <- terra::values(sim$true_break, mat = FALSE)
  pct_break <- mean(break_vals)
  # Block-based selection, not cell-based, so this won't be exactly
  # 50% -- but should be reasonably close for a 20x20 grid with several
  # blocks.
  expect_true(pct_break > 0.25 && pct_break < 0.75)
})

test_that("break selection is reproducible with the same seed", {
  sim1 <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 8, break_type = "mean",
                           break_fraction = 0.3, seed = 308)
  sim2 <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 8, break_type = "mean",
                           break_fraction = 0.3, seed = 308)
  expect_identical(terra::values(sim1$true_break, mat = FALSE),
                    terra::values(sim2$true_break, mat = FALSE))
  expect_identical(terra::values(sim1$series, mat = FALSE),
                    terra::values(sim2$series, mat = FALSE))
})

test_that("break_fraction = 0 (the default) leaves true_break all 0 even when break_type != 'none'", {
  sim <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 8, break_type = "slope",
                          break_fraction = 0, seed = 309)
  vals <- terra::values(sim$true_break, mat = FALSE)
  expect_true(all(vals == 0))
  expect_null(sim$break_time)
})

test_that("true_break can be used directly as compare_detections() ground truth, matching how true_slope already is", {
  sim <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 8, trend_fraction = 0,
                          break_type = "mean", break_fraction = 0.4,
                          break_magnitude = 5, seed = 310)
  ground_truth <- as.logical(terra::values(sim$true_break, mat = FALSE))
  # A trivial "detector" that just guesses TRUE everywhere, purely to
  # confirm compare_detections() accepts this ground truth vector
  # without needing any real change-point method to exist yet.
  guess_all_true <- rep(TRUE, length(ground_truth))
  comparison <- compare_detections(
    detections = list(naive = guess_all_true),
    ground_truth = ground_truth
  )
  expect_s3_class(comparison, "compare_detections")
  expect_equal(comparison$Sensitivity, 1)
})

test_that("sim_trend_stack works on a 1x1 grid for every trend_shape, not just the default (regression test for a real division-by-zero bug: max(dist_centre) and half_width were both exactly 0 for a single cell)", {
  for (shape in c("radial", "gradient", "block")) {
    sim <- sim_trend_stack(nrow = 1, ncol = 1, n_time = 8,
                            trend_shape = shape, smooth_radius = 0, seed = 1)
    expect_true(inherits(sim$series, "SpatRaster"))
    expect_equal(terra::ncell(sim$series), 1)
    expect_false(anyNA(terra::values(sim$true_slope, mat = FALSE)))
  }
})

test_that("sim_trend_stack(smooth_radius > 0) on a raster too small for the window warns and falls back to unsmoothed noise, instead of letting terra::focal()'s own cryptic error propagate", {
  expect_warning(
    sim <- sim_trend_stack(nrow = 1, ncol = 1, n_time = 8,
                            smooth_radius = 1, seed = 1),
    "skipped"
  )
  expect_false(anyNA(terra::values(sim$series, mat = FALSE)))
})

test_that("sim_trend_stack(smooth_radius > 0) works correctly after being batched into a single multi-layer terra::focal() call (regression test for the loop-to-batch optimisation)", {
  sim <- sim_trend_stack(nrow = 20, ncol = 20, n_time = 6, noise_sd = 2,
                          smooth_radius = 2, trend_fraction = 0, seed = 400)
  X <- terra::values(sim$series, mat = TRUE)
  expect_false(anyNA(X))
  expect_equal(dim(X), c(400, 6))

  # The rescaling step is designed to preserve each layer's own
  # cell-to-cell standard deviation after smoothing (smoothing alone
  # would shrink it, since averaging neighbours flattens variance) --
  # this is the one property the refactor could plausibly have broken
  # if the multi-layer batching handled columns differently from the
  # original per-layer loop.
  per_layer_sd <- apply(X, 2, sd)
  expect_true(all(abs(per_layer_sd - 2) < 0.5))

  # Every layer should differ from every other -- a coarse sanity
  # check that the multi-layer batching didn't accidentally smooth
  # all layers together into one shared result instead of each
  # independently (terra::focal() on a multi-layer input processes
  # each layer's own neighbourhood, not a cross-layer one).
  all_different <- all(vapply(2:ncol(X), function(k) any(X[, k] != X[, 1]),
                               logical(1)))
  expect_true(all_different)
})

test_that("sim_trend_stack separates spatial scale from intensity while preserving its established default", {
  args <- list(nrow = 20, ncol = 20, n_time = 6,
               trend_fraction = 0, constant_block = FALSE,
               smooth_radius = 2, seed = 401)
  legacy <- do.call(sim_trend_stack, args)
  explicit_one <- do.call(sim_trend_stack,
                          c(args, list(spatial_rho = 1)))
  expect_equal(terra::values(legacy$series),
               terra::values(explicit_one$series))

  independent_by_rho <- do.call(sim_trend_stack,
                                c(args, list(spatial_rho = 0)))
  independent_by_radius <- do.call(
    sim_trend_stack,
    utils::modifyList(args, list(smooth_radius = 0)))
  expect_equal(terra::values(independent_by_rho$series),
               terra::values(independent_by_radius$series))

  intermediate_a <- do.call(sim_trend_stack,
                            c(args, list(spatial_rho = 0.5)))
  intermediate_b <- do.call(sim_trend_stack,
                            c(args, list(spatial_rho = 0.5)))
  expect_equal(terra::values(intermediate_a$series),
               terra::values(intermediate_b$series))
  expect_false(isTRUE(all.equal(terra::values(intermediate_a$series),
                                terra::values(explicit_one$series))))
})

test_that("sim_trend_stack validates spatial_rho as a finite [0, 1] scalar", {
  for (bad in list(-0.1, 1.1, NA_real_, Inf, c(0, 1), "high")) {
    expect_error(sim_trend_stack(spatial_rho = bad), "spatial_rho")
  }
  expect_no_error(sim_trend_stack(nrow = 5, ncol = 5, n_time = 3,
                                   spatial_rho = 0, seed = 1))
  expect_no_error(sim_trend_stack(nrow = 5, ncol = 5, n_time = 3,
                                   spatial_rho = 1, seed = 1))
})

test_that("sim_trend_stack enforces the documented smooth_radius contract", {
  for (bad in list(-1, 0.5, NA_real_, Inf, c(0, 1), "wide")) {
    expect_error(sim_trend_stack(smooth_radius = bad), "smooth_radius")
  }
  expect_no_error(sim_trend_stack(nrow = 5, ncol = 5, n_time = 3,
                                   smooth_radius = 0, seed = 1))
})
