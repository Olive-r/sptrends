test_that(".validate_time_axis enforces one finite increasing value per layer", {
  expect_equal(sptrends:::.validate_time_axis(2001:2004, 4), 2001:2004)
  expect_error(sptrends:::.validate_time_axis(letters[1:4], 4), "numeric")
  expect_error(sptrends:::.validate_time_axis(1:3, 4), "length")
  expect_error(
    sptrends:::.validate_time_axis(c(1, 2, NA, 4), 4),
    "finite"
  )
  expect_error(
    sptrends:::.validate_time_axis(c(1, 2, 2, 4), 4),
    "duplicate"
  )
  expect_error(
    sptrends:::.validate_time_axis(c(1, 3, 2, 4), 4),
    "strictly increasing"
  )
})

test_that(".with_timer passes through the result and reports timing when verbose", {
  result <- sptrends:::.with_timer("test step", {
    Sys.sleep(0)
    42
  }, verbose = FALSE)
  expect_identical(result, 42)

  expect_message(
    sptrends:::.with_timer("test step", 1 + 1, verbose = TRUE),
    "completed in"
  )
})

test_that(".sptrends_progress supports one step and respects quiet mode", {
  output <- capture.output({
    pb <- sptrends:::.sptrends_progress(1, "label", verbose = TRUE)
    sptrends:::.sptrends_progress_step(pb, 1)
    sptrends:::.sptrends_progress_close(pb)
  }, type = "message")
  expect_s3_class(pb, "sptrends_progress")
  expect_match(paste(output, collapse = " "), "100.0%", fixed = TRUE)
  expect_null(sptrends:::.sptrends_progress(10, "label", verbose = FALSE))
  expect_null(sptrends:::.sptrends_progress(0, "label", verbose = TRUE))
})

test_that(".sptrends_progress returns a usable progress bar object for multiple steps", {
  output <- capture.output({
    pb <- sptrends:::.sptrends_progress(5, "label", verbose = TRUE)
    sptrends:::.sptrends_progress_step(pb, 3)
    sptrends:::.sptrends_progress_close(pb)
  }, type = "message")
  expect_false(is.null(pb))
  expect_s3_class(pb, "sptrends_progress")
  text <- paste(output, collapse = " ")
  expect_match(text, "progress:", fixed = TRUE)
  expect_match(text, "elapsed:", fixed = TRUE)
  expect_match(text, "remaining:", fixed = TRUE)
})

test_that(".sptrends_progress_step/.sptrends_progress_close do nothing on a NULL bar", {
  expect_null(sptrends:::.sptrends_progress_step(NULL, 1))
  expect_null(sptrends:::.sptrends_progress_close(NULL))
  expect_null(sptrends:::.sptrends_progress_step(new.env(), 1))
  expect_null(sptrends:::.sptrends_progress_close(new.env()))

  output <- capture.output({
    pb <- sptrends:::.sptrends_progress(3, "closed", verbose = TRUE)
    sptrends:::.sptrends_progress_step(pb, 0)
    pb$last_step <- 1
    pb$last_update <- proc.time()[["elapsed"]]
    expect_null(sptrends:::.sptrends_progress_step(pb, 2))
    sptrends:::.sptrends_progress_step(pb, 4)
    sptrends:::.sptrends_progress_close(pb)
    expect_null(sptrends:::.sptrends_progress_step(pb, 3))
    expect_null(sptrends:::.sptrends_progress_close(pb))
  }, type = "message")
  expect_match(paste(output, collapse = " "), "estimating", fixed = TRUE)
})

test_that(".sptrends_format_duration formats elapsed and remaining time", {
  expect_identical(sptrends:::.sptrends_format_duration(0), "00:00")
  expect_identical(sptrends:::.sptrends_format_duration(65), "01:05")
  expect_identical(sptrends:::.sptrends_format_duration(3661), "01:01:01")
  expect_identical(sptrends:::.sptrends_format_duration(NA_real_), "unknown")
  expect_identical(sptrends:::.sptrends_format_duration(Inf), "unknown")
  expect_identical(sptrends:::.sptrends_format_duration(-1), "unknown")
  expect_identical(sptrends:::.sptrends_format_duration("1"), "unknown")
  expect_identical(sptrends:::.sptrends_format_duration(c(1, 2)), "unknown")
})

test_that(".safe_categorical_plot skips the map and messages when no values are present", {
  r <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3)
  r <- terra::setValues(r, rep(NA_real_, 9))
  expect_message(
    sptrends:::.safe_categorical_plot(r, values = c(-1, 0, 1), colours = c("red", "grey", "blue"),
                                       labels = c("dec", "none", "inc"), main = "test map"),
    "no valid cells"
  )
})

test_that(".safe_categorical_plot draws normally when at least one category is present", {
  r <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3)
  r <- terra::setValues(r, c(rep(1, 5), rep(0, 4)))
  expect_silent(
    sptrends:::.safe_categorical_plot(r, values = c(-1, 0, 1), colours = c("red", "grey", "blue"),
                                       labels = c("dec", "none", "inc"), main = "test map")
  )
})

test_that(".save_current_plot writes a PNG file at the given path", {
  r <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3)
  r <- terra::setValues(r, 1:9)
  terra::plot(r)
  path <- tempfile(fileext = ".png")
  sptrends:::.save_current_plot(path)
  expect_true(file.exists(path))
  unlink(path)
})

test_that(".sptrends_parallel_lapply with n_cores <= 1 matches plain lapply", {
  result <- sptrends:::.sptrends_parallel_lapply(1:5, function(x) x^2, n_cores = 1)
  expect_identical(result, lapply(1:5, function(x) x^2))
})

test_that(".sptrends_parallel_lapply with n_cores > 1 gives the same results as sequential", {
  seq_result <- sptrends:::.sptrends_parallel_lapply(1:6, function(x) x * 2, n_cores = 1)
  par_result <- sptrends:::.sptrends_parallel_lapply(1:6, function(x) x * 2, n_cores = 2)
  expect_identical(seq_result, par_result)
})

test_that(".sptrends_parallel_lapply exports variables from the caller's environment to workers", {
  multiplier <- 10
  result <- sptrends:::.sptrends_parallel_lapply(
    1:3, function(x) x * multiplier, n_cores = 2,
    export_vars = "multiplier", export_env = environment()
  )
  expect_identical(result, list(10, 20, 30))
})

test_that(".sptrends_parallel_lapply's seed gives reproducible parallel results", {
  draw_random <- function(x) stats::runif(1)
  result_a <- sptrends:::.sptrends_parallel_lapply(1:4, draw_random, n_cores = 2, seed = 123)
  result_b <- sptrends:::.sptrends_parallel_lapply(1:4, draw_random, n_cores = 2, seed = 123)
  expect_identical(result_a, result_b)
})

test_that(".sptrends_parallel_lapply's packages argument makes a package available on workers", {
  # median() is unambiguous, but this exercises the same clusterCall() path
  # used to load a namespaced package on each PSOCK worker before FUN runs.
  result <- sptrends:::.sptrends_parallel_lapply(
    list(1:5, 6:10), function(x) stats::median(x), n_cores = 2, packages = "stats"
  )
  expect_identical(result, list(3L, 8L))
})

test_that(".sptrends_parallel_lapply's n_cores-capping message fires safely via a mocked detectCores()", {
  skip_if_not_installed("testthat", minimum_version = "3.2.0")
  # Never actually request more cores than the machine has (that trips
  # R CMD check's own core-limit protection, independent of this
  # package's own capping logic -- see the NEWS entry for the earlier
  # mistake this exact approach avoids). Instead, mock detectCores() to
  # report fewer cores than we ask for, while the real request (2) stays
  # safely within CRAN's own limit regardless of what the mock reports.
  testthat::local_mocked_bindings(
    detectCores = function(...) 1L,
    .package = "parallel"
  )
  expect_message(
    result <- sptrends:::.sptrends_parallel_lapply(1:4, function(x) x * 2, n_cores = 2),
    "using"
  )
  expect_identical(result, as.list(c(2, 4, 6, 8)))
})

test_that(".robust_diverging_range() caps the range at 2 SD rather than the single most extreme cell, when that cell is a genuine outlier", {
  # 24 cells clustered near 0-1, one wild outlier at 100 -- the naive
  # max_abs range would stretch from -100 to 100, compressing all the
  # other 24 cells' colour into a sliver near the palette's midpoint.
  vals <- c(runif(24, -1, 1), 100)
  r <- terra::rast(nrows = 5, ncols = 5, vals = vals)

  range_lim <- sptrends:::.robust_diverging_range(r)
  sd_val <- terra::global(r, "sd", na.rm = TRUE)$sd

  expect_equal(range_lim, c(-2 * sd_val, 2 * sd_val), tolerance = 1e-8)
  expect_lt(range_lim[2], 100)  # genuinely capped, not just echoing max_abs
})

test_that(".robust_diverging_range() falls back to max_abs when standard deviation would not actually clip anything (no real outlier)", {
  # A raster with no genuine outlier: 2 SD comfortably exceeds max_abs
  # here, so min(max_abs, 2*sd) should just be max_abs itself -- the
  # cap should never make the range WIDER than the real data.
  vals <- seq(-1, 1, length.out = 25)
  r <- terra::rast(nrows = 5, ncols = 5, vals = vals)

  range_lim <- sptrends:::.robust_diverging_range(r)
  expect_equal(range_lim, c(-1, 1), tolerance = 1e-8)
})

test_that(".robust_diverging_range() does not error on an all-NA or all-zero raster (same degenerate cases as the max_abs guards it replaces)", {
  r_na <- terra::rast(nrows = 4, ncols = 4, vals = NA_real_)
  expect_equal(sptrends:::.robust_diverging_range(r_na), c(-1, 1))

  r_zero <- terra::rast(nrows = 4, ncols = 4, vals = 0)
  expect_equal(sptrends:::.robust_diverging_range(r_zero), c(-1, 1))
})

test_that(".robust_diverging_range() falls back to max_abs on a constant, non-zero raster (sd = 0 but max_abs != 0 -- a different degenerate case from the all-zero one above, which returns earlier at the max_abs == 0 check without ever reaching the sd_val guard)", {
  r_const <- terra::rast(nrows = 4, ncols = 4, vals = 5)
  expect_equal(sptrends:::.robust_diverging_range(r_const), c(-5, 5))
})

test_that(".sptrends_brand is a fixed, named 5-colour palette with valid hex colours", {
  brand <- sptrends:::.sptrends_brand
  expect_named(brand, c("navy", "royal_blue", "azure", "cyan", "white"))
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", unlist(brand))))
  expect_identical(brand$cyan, "#22CFF5")
})

test_that(".sptrends_shared_cluster() returns NULL for n_cores <= 1, and a real, usable PSOCK cluster for n_cores > 1", {
  expect_null(sptrends:::.sptrends_shared_cluster(1))
  expect_null(sptrends:::.sptrends_shared_cluster(0))

  cl <- sptrends:::.sptrends_shared_cluster(2)
  on.exit(parallel::stopCluster(cl))
  expect_s3_class(cl, "SOCKcluster")
  expect_equal(length(cl), 2)
  # A real, minimal round-trip through the cluster it built, not just
  # checking the class of the object.
  expect_equal(parallel::parSapply(cl, 1:3, function(x) x * 2), c(2, 4, 6))
})

test_that(".sptrends_shared_cluster() warns and caps at the number of detected cores when n_cores exceeds it", {
  testthat::local_mocked_bindings(
    detectCores = function(...) 2L,
    .package = "parallel"
  )
  expect_message(
    cl <- sptrends:::.sptrends_shared_cluster(99),
    "only 2 logical cores"
  )
  on.exit(parallel::stopCluster(cl))
  expect_equal(length(cl), 2)
})

test_that(".sptrends_parallel_lapply()'s shared_cluster argument reuses an already-running cluster instead of building its own, and gives identical results to building one from n_cores directly", {
  cl <- parallel::makeCluster(2, type = "PSOCK")
  on.exit(parallel::stopCluster(cl))

  x <- 10
  result_shared <- sptrends:::.sptrends_parallel_lapply(
    1:5, function(i) i + x, export_vars = "x", export_env = environment(),
    shared_cluster = cl)
  result_own <- sptrends:::.sptrends_parallel_lapply(
    1:5, function(i) i + x, n_cores = 2,
    export_vars = "x", export_env = environment())

  expect_equal(unlist(result_shared), unlist(result_own))
  expect_equal(unlist(result_shared), 11:15)

  # The cluster passed in must still be alive and usable afterwards --
  # .sptrends_parallel_lapply() must not have stopped it, unlike the
  # cluster it builds and tears down itself when shared_cluster is NULL.
  expect_equal(parallel::parSapply(cl, 1:2, function(x) x), c(1, 2))
})

test_that(".sptrends_parallel_lapply()'s shared_cluster argument respects its own seed argument, giving reproducible results across two separate calls sharing the same cluster", {
  cl <- parallel::makeCluster(2, type = "PSOCK")
  on.exit(parallel::stopCluster(cl))

  draw_random <- function(i) stats::runif(1)
  result_a <- sptrends:::.sptrends_parallel_lapply(
    1:4, draw_random, shared_cluster = cl, seed = 42)
  result_b <- sptrends:::.sptrends_parallel_lapply(
    1:4, draw_random, shared_cluster = cl, seed = 42)

  expect_equal(unlist(result_a), unlist(result_b))
})

test_that(".sptrends_parallel_lapply()'s shared_cluster argument respects its own packages argument, loading namespaces on the shared cluster's workers", {
  cl <- parallel::makeCluster(2, type = "PSOCK")
  on.exit(parallel::stopCluster(cl))

  result <- sptrends:::.sptrends_parallel_lapply(
    1:2, function(i) requireNamespace("stats", quietly = TRUE),
    shared_cluster = cl, packages = "stats")
  expect_true(all(unlist(result)))
})

test_that(".sptrends_elapsed_timer reports elapsed time only when verbose", {
  expect_message({
    finish <- sptrends:::.sptrends_elapsed_timer("core step", TRUE)
    expect_null(finish())
  }, "\\[core step\\] elapsed: [0-9]+\\.[0-9]{2} s")

  expect_silent({
    finish <- sptrends:::.sptrends_elapsed_timer("quiet step", FALSE)
    expect_null(finish())
  })
})
