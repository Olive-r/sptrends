.local_moran_null_reference <- function(r, result, connectivity = "queen") {
  values <- terra::values(r, mat = FALSE)
  ok <- !is.na(values)
  valid <- which(ok)
  W <- sptrends:::.prepare_moran_neighbourhood(r, ok, connectivity)$W
  vapply(result$permutation_seeds, function(seed) {
    set.seed(seed)
    permuted <- values
    permuted[valid] <- values[valid][sample.int(length(valid))]
    centred <- permuted
    centred[ok] <- centred[ok] - mean(centred[ok])
    centred[!ok] <- 0
    centred[valid] * as.numeric(W %*% centred)[valid] /
      mean(centred[valid]^2)
  }, numeric(length(valid)))
}

test_that("local Moran statistics satisfy the global-local identity", {
  r <- terra::rast(nrows = 3, ncols = 3)
  terra::values(r) <- c(1, 2, 4, 3, 6, 8, 5, 7, 9)
  local <- spatial_autocorrelation(
    r, scope = "local", connectivity = "rook", nperm = 19,
    seed = 1, verbose = FALSE, report = FALSE
  )
  global <- spatial_autocorrelation(
    r, scope = "global", connectivity = "rook", nperm = 19,
    seed = 1, verbose = FALSE, report = FALSE
  )
  ok <- rep(TRUE, terra::ncell(r))
  nb <- sptrends:::.prepare_moran_neighbourhood(r, ok, "rook")
  observed <- terra::values(local$statistic, mat = FALSE)
  expect_equal(sum(observed), nb$S0 * global$statistic, tolerance = 1e-12)
})

test_that("local p-values feed BH, BKY and BY through fdr_correction", {
  r <- sim_trend_stack(
    nrow = 6, ncol = 6, n_time = 1, seed = 11
  )$series[[1]]
  local <- spatial_autocorrelation(
    r, scope = "local", nperm = 39, seed = 2,
    verbose = FALSE, report = FALSE
  )
  fdr <- fdr_correction(
    local$p, method = c("BH", "BKY", "BY"), q = 0.05,
    verbose = FALSE, report = FALSE
  )
  raw <- terra::values(local$p, mat = FALSE)
  expect_equal(
    terra::values(fdr$rasters$q_BH, mat = FALSE),
    stats::p.adjust(raw, method = "BH")
  )
  expect_equal(
    terra::values(fdr$rasters$q_BY, mat = FALSE),
    stats::p.adjust(raw, method = "BY")
  )
  expect_equal(
    terra::values(fdr$rasters$q_BKY, mat = FALSE),
    sptrends:::fdr_bky(raw)$q_value
  )
  expect_equal(
    terra::values(local$significant_raw, mat = FALSE),
    as.numeric(raw <= local$alpha)
  )
})

test_that("local Moran matches a direct matrix calculation", {
  r <- terra::rast(nrows = 2, ncols = 3)
  values <- c(1, 5, 2, 8, 3, 7)
  terra::values(r) <- values
  result <- spatial_autocorrelation(
    r, scope = "local", connectivity = "queen", nperm = 19,
    seed = 21, verbose = FALSE, report = FALSE
  )
  ok <- rep(TRUE, length(values))
  W <- sptrends:::.prepare_moran_neighbourhood(r, ok, "queen")$W
  centred <- values - mean(values)
  expected <- centred * as.numeric(W %*% centred) /
    mean(centred^2)
  expect_equal(
    terra::values(result$statistic, mat = FALSE), expected,
    tolerance = 1e-12
  )
  observed <- terra::values(result$statistic, mat = FALSE)
  null_mean <- terra::values(result$null_mean, mat = FALSE)
  null_sd <- terra::values(result$null_sd, mat = FALSE)
  expect_equal(
    terra::values(result$z, mat = FALSE),
    (observed - null_mean) / null_sd
  )
})

test_that("local Getis-Ord Gi* includes the focal cell", {
  r <- terra::rast(
    nrows = 1, ncols = 3,
    xmin = 0, xmax = 3, ymin = 0, ymax = 1
  )
  terra::values(r) <- c(1, 2, 7)
  result <- spatial_autocorrelation(
    r, method = "getis_ord", scope = "local", connectivity = "rook",
    nperm = 19, seed = 3, verbose = FALSE, report = FALSE
  )
  observed <- terra::values(result$statistic, mat = FALSE)
  expect_equal(observed, c(3, 10, 9) / 10, tolerance = 1e-12)
  expect_output(print(result), "Getis-Ord Gi\\*")
  expect_message(
    spatial_autocorrelation(
      r, method = "getis_ord", scope = "local", connectivity = "rook",
      nperm = 9, seed = 3, verbose = TRUE, report = FALSE
    ),
    "Computing local Gi\\* permutation inference"
  )
})

test_that("local results are reproducible and parallel invariant", {
  .skip_unless_parallel_tests()
  r <- sim_trend_stack(
    nrow = 6, ncol = 6, n_time = 1, seed = 12
  )$series[[1]]
  sequential <- spatial_autocorrelation(
    r, scope = "local", nperm = 29,
    seed = 4, n_cores = 1, verbose = FALSE, report = FALSE
  )
  parallel <- spatial_autocorrelation(
    r, scope = "local", nperm = 29,
    seed = 4, n_cores = 2, verbose = FALSE, report = FALSE
  )
  expect_equal(sequential$permutation_seeds, parallel$permutation_seeds)
  expect_equal(
    terra::values(sequential$p, mat = FALSE),
    terra::values(parallel$p, mat = FALSE)
  )
})

test_that("local alternatives use their stated permutation tails", {
  r <- terra::rast(nrows = 3, ncols = 3)
  terra::values(r) <- c(1, 1, 2, 1, 9, 8, 2, 8, 9)
  alternatives <- c("greater", "less", "two.sided")
  results <- lapply(alternatives, function(alternative) {
    spatial_autocorrelation(
      r, scope = "local", alternative = alternative, nperm = 39,
      seed = 22, verbose = FALSE, report = FALSE
    )
  })
  observed <- terra::values(results[[1]]$statistic, mat = FALSE)
  null <- .local_moran_null_reference(r, results[[1]])
  expected_greater <- (rowSums(null >= observed) + 1) / 40
  expected_less <- (rowSums(null <= observed) + 1) / 40
  null_mean <- rowMeans(null)
  expected_two <- (
    rowSums(abs(null - null_mean) >= abs(observed - null_mean)) + 1
  ) / 40
  expect_equal(
    terra::values(results[[1]]$p, mat = FALSE), expected_greater
  )
  expect_equal(terra::values(results[[2]]$p, mat = FALSE), expected_less)
  expect_equal(terra::values(results[[3]]$p, mat = FALSE), expected_two)
  expect_true(all(expected_two >= 1 / 40))
})

test_that("isolated valid cells are not treated as local tests", {
  r <- terra::rast(
    nrows = 3, ncols = 3,
    xmin = 0, xmax = 3, ymin = 0, ymax = 3
  )
  terra::values(r) <- c(1, NA, 2, NA, NA, NA, NA, NA, 3)
  expect_error(
    spatial_autocorrelation(
      r, scope = "local", connectivity = "rook", nperm = 19,
      verbose = FALSE,
      report = FALSE
    ),
    "No valid cell has valid neighbours"
  )

  terra::values(r) <- c(1, 2, NA, NA, NA, NA, NA, NA, 3)
  result <- spatial_autocorrelation(
    r, scope = "local", connectivity = "rook", nperm = 19,
    seed = 24, verbose = FALSE, report = FALSE
  )
  statistic <- terra::values(result$statistic, mat = FALSE)
  p <- terra::values(result$p, mat = FALSE)
  expect_true(is.na(statistic[9]))
  expect_true(is.na(p[9]))
  expect_true(all(is.finite(statistic[c(1, 2)])))
  expect_equal(result$N, 3)
  expect_equal(result$N_tested, 2)
  fdr <- fdr_correction(
    result$p, method = c("BH", "BKY", "BY"),
    report = FALSE, verbose = FALSE
  )
  expect_equal(fdr$summary_bky$m, 2)
  expect_true(is.na(terra::values(fdr$rasters$q_BH, mat = FALSE)[9]))
  expect_true(is.na(terra::values(fdr$rasters$q_BKY, mat = FALSE)[9]))
  expect_true(is.na(terra::values(fdr$rasters$q_BY, mat = FALSE)[9]))
})

test_that("local precomputed neighbourhood preserves the complete result", {
  r <- sim_trend_stack(
    nrow = 5, ncol = 5, n_time = 1, seed = 25
  )$series[[1]]
  ok <- !is.na(terra::values(r, mat = FALSE))
  neighbourhood <- sptrends:::.prepare_moran_neighbourhood(r, ok, "queen")
  direct <- spatial_autocorrelation(
    r, scope = "local", nperm = 29,
    seed = 26, verbose = FALSE, report = FALSE
  )
  reused <- spatial_autocorrelation(
    r, scope = "local", nperm = 29,
    seed = 26, precomputed_neighbourhood = neighbourhood,
    verbose = FALSE, report = FALSE
  )
  expect_equal(direct$permutation_seeds, reused$permutation_seeds)
  expect_equal(
    terra::values(direct$p, mat = FALSE),
    terra::values(reused$p, mat = FALSE)
  )
})

test_that("local NA geometry and S3 reporting are preserved", {
  r <- terra::rast(nrows = 4, ncols = 4)
  terra::values(r) <- seq_len(16)
  r[6] <- NA
  result <- spatial_autocorrelation(
    r, scope = "local", nperm = 19,
    seed = 5, verbose = FALSE, report = FALSE
  )
  expect_true(is.na(terra::values(result$statistic, mat = FALSE)[6]))
  expect_output(print(result), "fdr_correction")
  expect_s3_class(summary(result), "data.frame")
  expect_error(plot(result), NA)
})

test_that("local edge branches and reporting paths remain usable", {
  r <- terra::rast(nrows = 3, ncols = 3)
  terra::values(r) <- seq_len(9)
  one <- spatial_autocorrelation(
    r, scope = "local", nperm = 1, seed = 33,
    verbose = FALSE, report = FALSE
  )
  expect_true(all(is.na(terra::values(one$null_sd, mat = FALSE))))

  ok <- rep(TRUE, 9)
  W <- sptrends:::.prepare_moran_neighbourhood(r, ok, "queen")$W
  gi <- sptrends:::.spatial_local_statistic(
    terra::values(r, mat = FALSE), W, ok, "getis_ord"
  )
  expect_length(gi, 9)

  result <- spatial_autocorrelation(
    r, scope = "local", nperm = 9,
    seed = 34, verbose = FALSE, report = FALSE
  )
  csv <- tempfile(fileext = ".csv")
  png <- tempfile(fileext = ".png")
  expect_output(summary(result, path = csv), "significant_raw")
  expect_true(file.exists(csv))
  expect_error(plot(result, path = png), NA)
  expect_true(file.exists(png))
  expect_output(
    spatial_autocorrelation(
      r, scope = "local", nperm = 9, seed = 35,
      verbose = FALSE, report = TRUE
    ),
    "significant_raw"
  )
})
