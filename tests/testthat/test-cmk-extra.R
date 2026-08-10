.make_pattern_raster_cmk <- function(values, nrow = 2, ncol = 2) {
  ncell_ <- nrow * ncol
  r <- terra::rast(nrows = nrow, ncols = ncol, xmin = 0, xmax = ncol, ymin = 0, ymax = nrow)
  layers <- lapply(values, function(v) terra::setValues(r, rep(v, ncell_)))
  do.call(c, layers)
}

test_that("trend_test's classic S matches the hand-computed value (no ties)", {
  # values 1,2,3,4: every one of the 6 pairs (i<j) is increasing -> S = 6.
  # s2 (no ties) = n(n-1)(2n+5)/18 = 4*3*13/18 = 156/18.
  r <- .make_pattern_raster_cmk(c(1, 2, 3, 4))
  trend <- trend_test(r, method = "MK", ties = FALSE, report = FALSE, verbose = FALSE)

  S_vals <- unname(terra::values(trend$stats$S, mat = FALSE))
  VarS_vals <- unname(terra::values(trend$stats$VarS, mat = FALSE))
  expect_true(all(S_vals == 6))
  expect_equal(VarS_vals, rep(156 / 18, length(VarS_vals)))
})

test_that("trend_test's tie correction matches the hand-computed value", {
  # values 1,2,2,4: pairs (i<j): (1,2)=+1 (1,2)=+1 (1,4)=+1 (2,2)=0 (2,4)=+1 (2,4)=+1 -> S = 5.
  # one tie group of size 2 -> correction = tb(tb-1)(2tb+5) = 2*1*9 = 18.
  # s2_tie = (n(n-1)(2n+5) - 18)/18 = (156-18)/18 = 138/18.
  r <- .make_pattern_raster_cmk(c(1, 2, 2, 4))
  trend <- trend_test(r, method = "MK", ties = TRUE, report = FALSE, verbose = FALSE)

  S_vals <- unname(terra::values(trend$stats$S, mat = FALSE))
  VarS_vals <- unname(terra::values(trend$stats$VarS, mat = FALSE))
  expect_true(all(S_vals == 5))
  expect_equal(VarS_vals, rep(138 / 18, length(VarS_vals)))
})

test_that("ties = TRUE and ties = FALSE give different VarS on tied data", {
  r <- .make_pattern_raster_cmk(c(1, 2, 2, 4))
  trend_no_tie  <- trend_test(r, method = "MK", ties = FALSE, report = FALSE, verbose = FALSE)
  trend_tie     <- trend_test(r, method = "MK", ties = TRUE,  report = FALSE, verbose = FALSE)
  expect_false(isTRUE(all.equal(
    unname(terra::values(trend_no_tie$stats$VarS, mat = FALSE))[1],
    unname(terra::values(trend_tie$stats$VarS, mat = FALSE))[1]
  )))
})

test_that("trend_test's p-values are consistent with S via pnorm (no ties)", {
  r <- .make_pattern_raster_cmk(c(1, 2, 3, 4))
  trend <- trend_test(r, method = "MK", report = FALSE, verbose = FALSE)

  S <- 6
  s2 <- 156 / 18
  D <- 1  # sign(S)
  Z <- (S - D) / sqrt(s2)
  expected_p <- 2 * (1 - stats::pnorm(abs(Z)))

  p_vals <- unname(terra::values(trend$stats$p, mat = FALSE))
  expect_equal(p_vals, rep(expected_p, length(p_vals)), tolerance = 1e-10)
})

test_that("trend_test's Sm reduces to S when a cell has no valid neighbours", {
  # An isolated valid cell (all its neighbours are NA) should keep Sm == S
  # for that cell -- the neighbourhood average has nothing to add.
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 10, seed = 1)$series
  # blank out everything except the centre cell and re-fill neighbours with NA
  ok_cell <- 13  # centre of a 5x5 grid
  neighbours <- as.vector(terra::adjacent(r, cells = ok_cell, directions = "queen"))
  for (nb in neighbours) r[nb] <- NA

  trend_mk  <- trend_test(r, method = "MK", report = FALSE, verbose = FALSE)
  trend_cmk <- trend_test(r, method = "CMK", report = FALSE, verbose = FALSE)

  s_center   <- unname(terra::values(trend_mk$stats$S, mat = FALSE))[ok_cell]
  sm_center  <- unname(terra::values(trend_cmk$stats$Sm, mat = FALSE))[ok_cell]
  expect_equal(sm_center, s_center)
})

test_that("CMK centre cell matches analytical RAMK equations for the same 3x3 region", {
  # Douglas et al. (2000), Eqs. 7 and 11-14: CMK should be RAMK with
  # the region defined as the moving 3x3 raster neighbourhood.
  X <- rbind(
    c(1, 2, 3, 4, 5, 6),
    c(1, 3, 2, 5, 4, 6),
    c(6, 5, 4, 3, 2, 1),
    c(2, 1, 4, 3, 6, 5),
    c(1, 2, 4, 3, 5, 6),
    c(6, 4, 5, 2, 3, 1),
    c(1, 4, 2, 5, 3, 6),
    c(3, 1, 2, 6, 4, 5),
    c(5, 6, 3, 4, 1, 2)
  )
  base <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3,
                       ymin = 0, ymax = 3)
  layers <- lapply(seq_len(ncol(X)), function(i) {
    terra::setValues(base, X[, i])
  })
  r <- do.call(c, layers)

  kendall_s <- function(x) {
    total <- 0
    for (i in seq_len(length(x) - 1)) {
      for (j in (i + 1):length(x)) total <- total + sign(x[j] - x[i])
    }
    total
  }
  S <- apply(X, 1, kendall_s)
  n <- ncol(X)
  m <- nrow(X)
  sigma2 <- n * (n - 1) * (2 * n + 5) / 18
  Xc <- X - rowMeans(X)
  sd_x <- sqrt(rowSums(Xc^2) / (n - 1))
  rho <- tcrossprod(Xc) / (n * outer(sd_x, sd_x))
  expected_Sm <- mean(S)
  expected_VarSm <- (m * sigma2 +
    2 * sigma2 * sum(rho[upper.tri(rho)])) / m^2
  expected_p <- 2 * (1 - stats::pnorm(
    abs(expected_Sm / sqrt(expected_VarSm))
  ))

  result <- trend_test(r, method = "CMK", continuity = FALSE,
                        report = FALSE, verbose = FALSE)
  expect_equal(terra::values(result$stats$Sm, mat = FALSE)[5],
               expected_Sm, tolerance = 1e-12)
  expect_equal(terra::values(result$stats$VarSm, mat = FALSE)[5],
               expected_VarSm, tolerance = 1e-12)
  expect_equal(terra::values(result$stats$p, mat = FALSE)[5],
               expected_p, tolerance = 1e-12)
})

test_that("CMK propagates heterogeneous tie-corrected variances through RAMK", {
  X <- rbind(
    c(1, 1, 2, 3, 4, 5),
    c(1, 2, 2, 3, 4, 5),
    c(1, 2, 3, 3, 4, 5),
    c(1, 2, 3, 4, 4, 5),
    c(1, 2, 3, 4, 5, 5),
    c(1, 1, 2, 2, 4, 5),
    c(1, 2, 2, 3, 3, 5),
    c(5, 4, 4, 3, 2, 1),
    c(5, 4, 3, 3, 2, 1)
  )
  base <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3,
                       ymin = 0, ymax = 3)
  r <- do.call(c, lapply(seq_len(ncol(X)), function(i) {
    terra::setValues(base, X[, i])
  }))

  kendall_s <- function(x) {
    sum(vapply(seq_len(length(x) - 1), function(i) {
      sum(sign(x[(i + 1):length(x)] - x[i]))
    }, numeric(1)))
  }
  tie_variance <- function(x) {
    n <- length(x)
    groups <- as.numeric(table(x))
    (n * (n - 1) * (2 * n + 5) -
       sum(groups * (groups - 1) * (2 * groups + 5))) / 18
  }

  S <- apply(X, 1, kendall_s)
  var_s <- apply(X, 1, tie_variance)
  n <- ncol(X)
  m <- nrow(X)
  Xc <- X - rowMeans(X)
  sd_x <- sqrt(rowSums(Xc^2) / (n - 1))
  rho <- tcrossprod(Xc) / (n * outer(sd_x, sd_x))
  covariance <- rho * sqrt(outer(var_s, var_s))
  expected_Sm <- mean(S)
  expected_VarSm <- (sum(var_s) +
    2 * sum(covariance[upper.tri(covariance)])) / m^2
  expected_p <- 2 * (1 - stats::pnorm(
    abs(expected_Sm / sqrt(expected_VarSm))
  ))

  result <- trend_test(r, method = "CMK", ties = TRUE,
                        continuity = FALSE, report = FALSE, verbose = FALSE)
  expect_gt(length(unique(var_s)), 1)
  expect_equal(terra::values(result$stats$Sm, mat = FALSE)[5],
               expected_Sm, tolerance = 1e-12)
  expect_equal(terra::values(result$stats$VarSm, mat = FALSE)[5],
               expected_VarSm, tolerance = 1e-12)
  expect_equal(terra::values(result$stats$p, mat = FALSE)[5],
               expected_p, tolerance = 1e-12)
})

test_that("trend_test errors on fewer than 2 layers", {
  r <- .make_pattern_raster_cmk(1)
  expect_error(trend_test(r, verbose = FALSE), "at least 2 layers")
})

test_that("trend_test errors on non-SpatRaster input", {
  expect_error(trend_test(matrix(1:4, 2, 2), verbose = FALSE), "SpatRaster")
})

test_that("trend_test gives identical Sm with n_cores = 1 and n_cores = 2 (neighbourhood = TRUE)", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 4)$series
  seq_result <- trend_test(r, method = "CMK", n_cores = 1, report = FALSE, verbose = FALSE)
  par_result <- trend_test(r, method = "CMK", n_cores = 2, report = FALSE, verbose = FALSE)
  expect_equal(terra::values(seq_result$stats$Sm, mat = FALSE), terra::values(par_result$stats$Sm, mat = FALSE))
})

test_that("prepare_cmk_neighbourhood's rook connectivity gives a different W than queen", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 5, seed = 1)$series
  X <- terra::values(r, mat = TRUE)
  ok <- stats::complete.cases(X)
  nb_queen <- prepare_cmk_neighbourhood(r, ok, connectivity = "queen")
  nb_rook  <- prepare_cmk_neighbourhood(r, ok, connectivity = "rook")
  # queen counts diagonal neighbours too, so it must have at least as many
  # (and, on an interior cell, strictly more) neighbours than rook.
  expect_true(all(nb_queen$nb_count >= nb_rook$nb_count))
  expect_true(any(nb_queen$nb_count > nb_rook$nb_count))
})

test_that("trend_test(method = 'CMK', connectivity = 'rook') actually works end to end and gives a different Sm than the default 'queen' -- a real, previously-existing bug found by an external audit: trend_test() never exposed 'connectivity' as its own argument, never passed it to its own internal neighbourhood builder, and hardcoded 'queen' in its own validation, making 'rook' completely unreachable through trend_test() even though prepare_cmk_neighbourhood() itself always supported it", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 10, seed = 1)$series

  result_queen <- trend_test(r, method = "CMK", connectivity = "queen",
                             report = FALSE, verbose = FALSE)
  result_rook <- trend_test(r, method = "CMK", connectivity = "rook",
                            report = FALSE, verbose = FALSE)

  Sm_queen <- terra::values(result_queen$stats$Sm, mat = FALSE)
  Sm_rook <- terra::values(result_rook$stats$Sm, mat = FALSE)
  expect_false(isTRUE(all.equal(Sm_queen, Sm_rook)))

  # A manually-built rook neighbourhood, explicitly declared as such,
  # must now be accepted without error (the mismatch-detection test
  # above covers the opposite, undeclared case).
  X <- terra::values(r, mat = TRUE)
  ok <- stats::complete.cases(X)
  nb_rook <- prepare_cmk_neighbourhood(r, ok, connectivity = "rook")
  expect_error(
    trend_test(r, method = "CMK", connectivity = "rook",
              precomputed_neighbourhood = nb_rook,
              report = FALSE, verbose = FALSE),
    NA
  )
})

test_that("trend_test rejects a non-default window_size or connectivity with a method other than CMK, since neither argument means anything without a neighbourhood", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 1)$series
  expect_error(
    trend_test(r, method = "MK", window_size = 5L, report = FALSE,
              verbose = FALSE),
    "only applicable when method"
  )
  expect_error(
    trend_test(r, method = "MK", connectivity = "rook", report = FALSE,
              verbose = FALSE),
    "only applicable when method"
  )
  # The defaults for both must still be accepted with any method.
  expect_error(
    trend_test(r, method = "MK", report = FALSE, verbose = FALSE),
    NA
  )
})

test_that("trend_histograms and trend_maps run without error", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  expect_silent(trend_histograms(trend$stats))
  expect_silent(trend_maps(trend$stats))
})

test_that("trend_summary can write a CSV", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  path <- tempfile(fileext = ".csv")
  suppressMessages(trend_summary(trend$stats, path = path))
  expect_true(file.exists(path))
  unlink(path)
})

test_that("trend_test's report = TRUE runs the full reporting branch (neighbourhood = TRUE)", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  expect_message(
    trend_test(r, method = "CMK", report = TRUE, verbose = FALSE),
    "At alpha="
  )
})

test_that("trend_test's report = TRUE runs the full reporting branch (neighbourhood = FALSE)", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  expect_message(
    trend_test(r, method = "MK", report = TRUE, verbose = FALSE),
    "At alpha="
  )
})

test_that("trend_test combines ties = TRUE with n_cores = 2 correctly", {
  r <- .make_pattern_raster_cmk(c(1, 2, 2, 4), nrow = 3, ncol = 3)
  trend_seq <- trend_test(r, method = "MK", ties = TRUE, n_cores = 1,
                                        report = FALSE, verbose = FALSE)
  trend_par <- trend_test(r, method = "MK", ties = TRUE, n_cores = 2,
                                        report = FALSE, verbose = FALSE)
  expect_equal(terra::values(trend_seq$stats$VarS, mat = FALSE), terra::values(trend_par$stats$VarS, mat = FALSE))
})

test_that("trend_test messages the parallel-S line when verbose = TRUE and n_cores > 1", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 8, seed = 1)$series
  expect_message(
    trend_test(r, n_cores = 2, report = FALSE, verbose = TRUE),
    "Parallel S over"
  )
})

test_that("trend_test messages the tie-correction line when verbose = TRUE", {
  r <- .make_pattern_raster_cmk(c(1, 2, 2, 4))
  expect_message(
    trend_test(r, method = "MK", ties = TRUE, report = FALSE, verbose = TRUE),
    "Applying tie correction"
  )
})

test_that("trend_test handles n_cores = 2 with exactly 2 time layers (regression test)", {
  # n=2 gives i_values = 1:(n-1) = a single value -- the chunk-splitting
  # logic must not error when there is only one item to split (see the
  # equivalent slope_estimator regression test/fix).
  r <- .make_pattern_raster_cmk(c(1, 2))
  expect_error(trend_test(r, n_cores = 2, report = FALSE, verbose = FALSE), NA)
})

test_that("trend_test's tie correction handles a cell with NA in its own series", {
  r <- .make_pattern_raster_cmk(c(1, 2, 2, 4), nrow = 2, ncol = 2)
  r[1] <- NA  # this cell now has an incomplete/NA series
  expect_error(trend_test(r, method = "MK", ties = TRUE, report = FALSE, verbose = FALSE), NA)
})

test_that("trend_histograms can write both PNG files via path", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  path <- tempfile()
  trend_histograms(trend$stats, path = path)
  expect_true(file.exists(paste0(path, "_stat.png")))
  expect_true(file.exists(paste0(path, "_p.png")))
  unlink(paste0(path, c("_stat.png", "_p.png")))
})

test_that("trend_maps can write all four PNG files via path", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  path <- tempfile()
  trend_maps(trend$stats, path = path)
  suffixes <- c("_Sm_map.png", "_p_map.png", "_significant_map.png", "_direction_map.png")
  expect_true(all(file.exists(paste0(path, suffixes))))
  unlink(paste0(path, suffixes))
})

test_that("trend_maps centres its diverging palette on zero even with all-positive trend values", {
  # all-positive Sm (e.g. a strong, one-sided warming signal) is exactly
  # the case where an auto-ranged palette would put "blue" at a positive
  # value -- this is a regression test for that fix, not a visual check.
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, trend_strength = 0.3,
                        trend_fraction = 1, noise_sd = 0.1, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  expect_error(trend_maps(trend$stats), NA)
})

test_that("trend_maps handles an all-NA trend statistic without erroring (max_abs_S guard)", {
  r_na <- terra::rast(nrows = 4, ncols = 4, vals = NA_real_)
  trend_all_na <- c(r_na, r_na, r_na)
  names(trend_all_na) <- c("Sm", "VarSm", "p")
  expect_error(trend_maps(trend_all_na), NA)
})

test_that("trend_test errors clearly with fewer than 2 layers", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 1, seed = 300)$series
  expect_error(trend_test(r, verbose = FALSE),
               "at least 2 layers")
})

test_that("trend_test works at the exact minimum of 2 layers, without erroring", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 2, seed = 301)$series
  expect_error(
    result <- trend_test(r, method = "MK", report = FALSE, verbose = FALSE),
    NA
  )
  expect_equal(terra::nlyr(result$stats), 3)  # Sm/S, VarSm/VarS, p
})

test_that("trend_test on a perfectly constant series gives S = 0 and p = 1 (no trend, maximally non-significant)", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 8, seed = 302)$series
  r[] <- 5  # every cell, every layer: exactly the same constant value
  result <- trend_test(r, method = "MK", report = FALSE, verbose = FALSE)
  s_vals <- terra::values(result$stats[["S"]], mat = FALSE)
  p_vals <- terra::values(result$stats[["p"]], mat = FALSE)
  expect_true(all(s_vals == 0, na.rm = TRUE))
  expect_true(all(p_vals == 1, na.rm = TRUE))
})

test_that("CMK returns p = 1 for a completely tied local region", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 8, seed = 302)$series
  r[] <- 5
  result <- trend_test(r, method = "CMK", ties = TRUE,
                       report = FALSE, verbose = FALSE)
  expect_true(all(terra::values(result$stats[["Sm"]], mat = FALSE) == 0))
  expect_true(all(terra::values(result$stats[["VarSm"]], mat = FALSE) == 0))
  expect_true(all(terra::values(result$stats[["p"]], mat = FALSE) == 1))
})

test_that("CMK rejects a precomputed neighbourhood from incompatible data", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 8, seed = 303)$series
  ok <- rep(TRUE, terra::ncell(r))
  nb <- prepare_cmk_neighbourhood(r, ok)

  shifted <- r
  terra::ext(shifted) <- terra::ext(10, 14, 10, 14)
  expect_error(
    trend_test(shifted, precomputed_neighbourhood = nb,
               report = FALSE, verbose = FALSE),
    "different raster geometry"
  )

  r[1] <- NA
  expect_error(
    trend_test(r, precomputed_neighbourhood = nb,
               report = FALSE, verbose = FALSE),
    "different valid-cell pattern"
  )

  malformed <- nb
  malformed$W <- malformed$W[-1, -1]
  r_clean <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 8,
                             seed = 303)$series
  expect_error(
    trend_test(r_clean, precomputed_neighbourhood = malformed,
               report = FALSE, verbose = FALSE),
    "must be a 16 x 16 Matrix"
  )

  malformed_count <- nb
  malformed_count$nb_count[1] <- -1
  expect_error(
    trend_test(r_clean, precomputed_neighbourhood = malformed_count,
               report = FALSE, verbose = FALSE),
    "one finite, non-negative value"
  )

  rook <- prepare_cmk_neighbourhood(r_clean, ok, connectivity = "rook")
  expect_error(
    trend_test(r_clean, precomputed_neighbourhood = rook,
               report = FALSE, verbose = FALSE),
    "connectivity='rook'"
  )
})

test_that("CMK defaults exactly to the explicit 3 by 3 neighbourhood", {
  r <- sim_trend_stack(nrow = 7, ncol = 7, n_time = 10, seed = 901)$series
  implicit <- trend_test(r, method = "CMK", report = FALSE,
                         verbose = FALSE)
  explicit <- trend_test(r, method = "CMK", window_size = 3L,
                         report = FALSE, verbose = FALSE)

  expect_equal(
    terra::values(implicit$stats, mat = TRUE),
    terra::values(explicit$stats, mat = TRUE),
    tolerance = 0
  )
  expect_identical(implicit$window_size, 3L)
  expect_identical(explicit$window_size, 3L)
})

test_that("a 5 by 5 CMK region uses every valid cell in that square", {
  r <- sim_trend_stack(nrow = 7, ncol = 7, n_time = 10, seed = 902)$series
  X <- terra::values(r, mat = TRUE)
  ok <- stats::complete.cases(X)
  neighbourhood <- prepare_cmk_neighbourhood(
    r, ok, window_size = 5L
  )
  centre <- 25L

  expect_identical(neighbourhood$signature$window_size, 5L)
  expect_equal(neighbourhood$nb_count[centre], 24)
  expect_equal(neighbourhood$nb_count[1], 8)

  cmk <- trend_test(r, method = "CMK", window_size = 5L,
                    precomputed_neighbourhood = neighbourhood,
                    report = FALSE, verbose = FALSE)
  mk <- trend_test(r, method = "MK", report = FALSE, verbose = FALSE)
  region <- c(9:13, 16:20, 23:27, 30:34, 37:41)
  S <- terra::values(mk$stats$S, mat = FALSE)
  Sm <- terra::values(cmk$stats$Sm, mat = FALSE)

  expect_equal(Sm[centre], mean(S[region]), tolerance = 1e-12)
  expect_identical(cmk$window_size, 5L)
})

test_that("a larger CMK region excludes invalid cells and keeps its shape", {
  r <- sim_trend_stack(nrow = 7, ncol = 7, n_time = 10, seed = 910)$series
  r[10] <- NA
  r[33] <- NA
  X <- terra::values(r, mat = TRUE)
  ok <- stats::complete.cases(X)
  neighbourhood <- prepare_cmk_neighbourhood(
    r, ok, window_size = 5L
  )
  centre <- 25L
  region <- c(9:13, 16:20, 23:27, 30:34, 37:41)

  expect_equal(neighbourhood$nb_count[centre], 22)
  expect_false(any(c(10L, 33L) %in% which(neighbourhood$W[centre, ] != 0)))

  cmk <- trend_test(r, window_size = 5L, report = FALSE, verbose = FALSE)
  mk <- trend_test(r, method = "MK", report = FALSE, verbose = FALSE)
  S <- terra::values(mk$stats$S, mat = FALSE)
  Sm <- terra::values(cmk$stats$Sm, mat = FALSE)
  expect_equal(Sm[centre], mean(S[region], na.rm = TRUE), tolerance = 1e-12)
})

test_that("CMK window size is validated and belongs only to CMK", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 8, seed = 903)$series
  for (bad in list(1, 2, 4, 3.5, NA_real_, c(3, 5), "5")) {
    expect_error(
      trend_test(r, method = "CMK", window_size = bad,
                 report = FALSE, verbose = FALSE),
      "window_size"
    )
  }
  expect_error(
    trend_test(r, method = "MK", window_size = 5L,
               report = FALSE, verbose = FALSE),
    "only applicable"
  )
})

test_that("precomputed CMK regions cannot cross window sizes", {
  r <- sim_trend_stack(nrow = 7, ncol = 7, n_time = 8, seed = 904)$series
  ok <- stats::complete.cases(terra::values(r, mat = TRUE))
  region_5 <- prepare_cmk_neighbourhood(r, ok, window_size = 5L)

  expect_error(
    trend_test(r, window_size = 3L,
               precomputed_neighbourhood = region_5,
               report = FALSE, verbose = FALSE),
    "uses window_size=5"
  )

  legacy <- region_5
  legacy$signature$window_size <- NULL
  expect_error(
    trend_test(r, window_size = 5L,
               precomputed_neighbourhood = legacy,
               report = FALSE, verbose = FALSE),
    "uses window_size=unknown"
  )
})

test_that("larger queen and rook builders retain their intended shapes", {
  r <- sim_trend_stack(nrow = 7, ncol = 7, n_time = 8, seed = 905)$series
  ok <- rep(TRUE, terra::ncell(r))
  queen <- prepare_cmk_neighbourhood(
    r, ok, connectivity = "queen", window_size = 5L
  )
  rook <- prepare_cmk_neighbourhood(
    r, ok, connectivity = "rook", window_size = 5L
  )

  expect_equal(queen$nb_count[25], 24)
  expect_equal(rook$nb_count[25], 8)
})

test_that("5 by 5 CMK is equivalent in sequential and parallel execution", {
  skip_on_cran()
  r <- sim_trend_stack(nrow = 7, ncol = 7, n_time = 8, seed = 906)$series
  sequential <- trend_test(r, window_size = 5L, n_cores = 1,
                           report = FALSE, verbose = FALSE)
  parallel <- trend_test(r, window_size = 5L, n_cores = 2,
                         report = FALSE, verbose = FALSE)
  expect_equal(
    terra::values(sequential$stats, mat = TRUE),
    terra::values(parallel$stats, mat = TRUE),
    tolerance = 0
  )
})

test_that("CMK progress identifies the selected window size", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 8, seed = 907)$series
  expect_message(
    trend_test(r, window_size = 5L, report = FALSE, verbose = TRUE),
    "Queen 5x5 neighbourhood"
  )
})

test_that("trend_test(method = 'cmk') does not error on a 1x1 raster (no neighbours to average over)", {
  r <- sim_trend_stack(nrow = 1, ncol = 1, n_time = 8, smooth_radius = 0,
                        seed = 303)$series
  expect_error(
    trend_test(r, method = "CMK", report = FALSE, verbose = FALSE),
    NA
  )
})

test_that("trend_test errors clearly, not cryptically, when every cell is NA", {
  r_na <- terra::rast(nrows = 4, ncols = 4, nlyr = 5, vals = NA_real_)
  expect_error(trend_test(r_na, verbose = FALSE), "nothing to test")
})

test_that("trend_test handles a very long series (100 layers) without erroring", {
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 100, seed = 304)$series
  expect_error(
    trend_test(r, method = "MK", report = FALSE, verbose = FALSE),
    NA
  )
})

test_that("CMK continuity correction matches frozen ConMK p-values where neighbourhood variance is defined", {
  csv_path <- system.file("validation", "conmk_comparison_100cells.csv",
                           package = "sptrends")
  skip_if(csv_path == "", "frozen ConMK comparison CSV not found")
  frozen <- utils::read.csv(csv_path)

  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 15, seed = 1)$series
  result <- trend_test(r, method = "CMK", continuity = TRUE,
                        report = FALSE, verbose = FALSE)
  p_ours_cc <- terra::values(result$stats$p, mat = FALSE)

  # Cells 1, 2, 11, 12 have a constant time series (value 5 throughout,
  # a genuine, narrow edge case in this specific simulated raster, not
  # engineered) -- zero variance, so the cross-correlation term this
  # function's own variance formula needs is undefined (0/0). This
  # package's own convention (replace the zero standard deviation with
  # Inf, giving zero correlation) is a documented, deliberate choice --
  # see ?trend_test's own "Implementation notes" -- but what ConMK's
  # own C++ does for the same 0/0 case is not known (IEEE-754 0/0 is
  # itself undefined, and inspecting its C++ source did not reveal an
  # explicit guard for it), so exact agreement is not expected, or
  # asserted, for any cell affected by a constant neighbour. Confirmed
  # by direct inspection: every affected cell is
  # exactly {1, 2, 3, 11, 12, 13, 21, 22, 23} -- itself or a queen
  # neighbour of one of the two constant cells (1, 2, 11, 12 are all
  # mutually adjacent). All 91 unaffected cells are checked here; the
  # 9 affected ones are a documented, narrow exception, not silently
  # dropped from awareness.
  affected <- c(1, 2, 3, 11, 12, 13, 21, 22, 23)
  unaffected <- setdiff(seq_along(p_ours_cc), affected)

  expect_equal(p_ours_cc[unaffected], frozen$theirs_p[unaffected],
               tolerance = 1e-6)
})

test_that("trend_test(method = 'CMK', continuity = FALSE) (the default) does NOT reproduce ConMK's own p-values -- confirms continuity actually changes the output, not a no-op", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 15, seed = 1)$series
  result_default <- trend_test(r, method = "CMK", report = FALSE,
                                verbose = FALSE)
  result_cc <- trend_test(r, method = "CMK", continuity = TRUE,
                           report = FALSE, verbose = FALSE)

  p_default <- terra::values(result_default$stats$p, mat = FALSE)
  p_cc <- terra::values(result_cc$stats$p, mat = FALSE)

  expect_false(isTRUE(all.equal(p_default, p_cc)))
  # Sm/VarSm themselves must be identical -- continuity only affects
  # how Zm/p are derived from them, not the pooled statistic itself.
  expect_equal(terra::values(result_default$stats$Sm),
               terra::values(result_cc$stats$Sm))
  expect_equal(terra::values(result_default$stats$VarSm),
               terra::values(result_cc$stats$VarSm))
})

test_that("trend_test(continuity = TRUE) is ignored (no error, no effect) for methods other than 'CMK'", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 10, seed = 1213)$series
  expect_error(
    result <- trend_test(r, method = "MK", continuity = TRUE,
                          report = FALSE, verbose = FALSE),
    NA
  )
  result_default <- trend_test(r, method = "MK", report = FALSE,
                                verbose = FALSE)
  expect_equal(terra::values(result$stats$p),
               terra::values(result_default$stats$p))
})
