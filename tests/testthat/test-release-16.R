# Small numerical regressions; no real parallel workers are started here.
.release16_raster <- function(values) {
  if (is.null(dim(values))) values <- matrix(values, nrow = 1)
  r <- terra::rast(nrows = 1, ncols = nrow(values), nlyrs = ncol(values),
                   xmin = 0, xmax = nrow(values), ymin = 0, ymax = 1)
  terra::values(r) <- values
  r
}

test_that("Yue-Pilon retains constant residuals without NaN", {
  values <- rbind(rep(4, 6), 2 * (1:6) + 3)
  result <- prewhiten(.release16_raster(values), method = "TFPW_Y",
                      report = FALSE, verbose = FALSE)
  expect_equal(unname(terra::values(result$series, mat = TRUE)),
                unname(values[, -1, drop = FALSE]))
  expect_equal(as.numeric(terra::values(result$diagnostics$Rho)), c(0, 0))
})

test_that("prewhitening histograms tolerate absent DW and unmodified NA cells", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  diagnostics <- .release16_raster(rbind(c(NA, 0.2, 1, 0), c(NA, NA, NA, NA)))
  names(diagnostics) <- c("DW_initial", "Rho", "Modified", "Clamped")
  expect_no_error(prewhiten_histograms(diagnostics))
  diagnostics$Modified <- terra::setValues(diagnostics$Modified, c(0, NA))
  expect_no_error(suppressMessages(prewhiten_histograms(diagnostics)))
})

test_that("published workflows forward n_cores without starting clusters in this test", {
  original_trend <- trend_test
  original_slope <- slope_estimator
  seen <- integer()
  testthat::local_mocked_bindings(
    .sptrends_shared_cluster = function(n_cores) NULL,
    trend_test = function(...) {
      args <- list(...); seen <<- c(seen, args$n_cores)
      args$n_cores <- 1L
      do.call(original_trend, args)
    },
    slope_estimator = function(...) {
      args <- list(...); seen <<- c(seen, args$n_cores)
      args$n_cores <- 1L
      do.call(original_slope, args)
    }
  )
  r <- .release16_raster(rbind(c(1, 3, 2, 4), c(4, 2, 3, 1)))
  workflow_tst(r, prewhiten = FALSE, n_cores = 2L,
               cmk_args = list(method = "MK"), report = FALSE, verbose = FALSE)
  workflow_rta(r, n_cores = 2L, cmk_args = list(method = "MK"),
               report = FALSE, verbose = FALSE)
  expect_equal(seen, rep(2L, 4))
})

test_that("original BKY plotting thresholds reproduce the actual decisions", {
  for (p in list(c(0.5, 0.7), c(0.001, 0.002), c(0.01, 0.049, 0.5))) {
    result <- fdr_bky(p, implementation = "original")
    crossing <- which(result$p_sorted <= result$thresh_bky)
    count <- if (length(crossing)) max(crossing) else 0L
    expect_equal(sum(result$reject), count)
  }
})

test_that("F1 is zero when errors occur without true positives", {
  result <- compare_detections(list(m = c(FALSE, TRUE)), c(TRUE, FALSE),
                               metrics = "f1", verbose = FALSE)
  expect_equal(result$F1, 0)
})

test_that("detection comparisons require the same evaluated cells", {
  detections <- list(a = c(NA, TRUE, FALSE), b = c(TRUE, TRUE, FALSE))
  truth <- c(TRUE, FALSE, FALSE)
  expect_error(compare_detections(detections, truth, verbose = FALSE),
                "different valid evaluation domains")
  result <- compare_detections(detections, truth,
                               evaluation_mask = c(FALSE, TRUE, TRUE),
                               verbose = FALSE)
  expect_equal(result$FP, c(1, 1))
  expect_equal(result$TN, c(1, 1))
})

test_that("slope and prewhitening benchmarks respect the common mask", {
  truth <- .release16_raster(matrix(c(1, 1), ncol = 1))
  simulation <- list(true_slope = truth, true_signal = c(TRUE, TRUE))
  slopes <- list(a = c(1, 100), b = c(1, -100))
  scored <- .benchmark_score(slopes, simulation, "slope", NULL,
                             character(), c(TRUE, FALSE))
  expect_equal(scored$RMSE, c(0, 0))
  outputs <- list(a = .release16_raster(rbind(1:4, 100 * (1:4))),
                  b = .release16_raster(rbind(1:4, -100 * (1:4))))
  scored <- .benchmark_prewhitening(outputs, simulation, c(TRUE, FALSE))
  expect_equal(scored$SlopeRMSE, c(0, 0))
  info <- .benchmark_scenario_data(list(), simulation, c(TRUE, FALSE))
  expect_equal(info$DomainCells, 1)
  expect_error(.benchmark_common_domain(list(c(TRUE, FALSE), c(TRUE, TRUE)),
                                         c(1, 1), NULL), "different valid")
})

test_that("global spatial diagnostics handle minimal permutation counts", {
  r <- .release16_raster(matrix(c(1, 2, 4), ncol = 1))
  result <- spatial_autocorrelation(r, nperm = 1, seed = 1,
                                    report = FALSE, verbose = FALSE)
  expect_no_error(suppressMessages(summary(result)))
  invalid <- .release16_raster(matrix(c(1, 0, 0), ncol = 1))
  expect_error(spatial_autocorrelation(invalid, method = "getis_ord",
                                        report = FALSE, verbose = FALSE),
                "at least two positive")
})

test_that("slope histograms accept a single valid cell", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_no_error(slope_histogram(.release16_raster(matrix(1))))
})

test_that("missing neighbour panels retain list positions", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  r <- .release16_raster(rbind(rep(NA, 4), c(1, 2, 4, 3)))
  result <- .plot_ts_small_multiples(r, 1:4, 1:2, 2)
  expect_length(result, 2)
  expect_true(is.na(result[[1]]$slope))
  expect_true(is.finite(result[[2]]$slope))
})

test_that("TST direction supports BY without a separately estimated slope", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  r <- .release16_raster(rbind(c(1, 3, 2, 4), c(4, 2, 3, 1)))
  result <- workflow_tst(r, prewhiten = FALSE, theil_sen = FALSE,
                         cmk_args = list(method = "MK"), fdr_method = "BY",
                         report = FALSE, verbose = FALSE)
  expect_no_error(plot(result, which = "direction", method = "BY"))
})
