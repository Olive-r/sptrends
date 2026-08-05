test_that("direction_map classifies into -1/0/1 only", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  fdr_result <- fdr_correction(trend$stats$p, report = FALSE, verbose = FALSE)

  direction <- direction_map(trend$stats, fdr_result, method = "BH", verbose = FALSE)
  expect_equal(names(direction), "binarised_trend_map")
  vals <- unique(terra::values(direction, mat = FALSE))
  expect_true(all(vals %in% c(-1, 0, 1, NA)))
})

test_that("direction_map errors when the requested method is absent", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  fdr_result <- fdr_correction(trend$stats$p, method = "BH", report = FALSE, verbose = FALSE)

  expect_error(direction_map(trend$stats, fdr_result, method = "BKY"), "rerun fdr_correction")
})

test_that("fdr_direction_summary returns one row per available method", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  fdr_result <- fdr_correction(trend$stats$p, method = c("BH", "BKY"), report = FALSE, verbose = FALSE)

  tab <- fdr_direction_summary(trend$stats, fdr_result)
  expect_setequal(tab$method, c("raw", "BH", "BKY"))
  expect_true(all(c("n_increase", "n_decrease", "n_not_significant") %in%
                    names(tab)))
})

test_that("direction_map does not error on rasters with NA cells (regression test)", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 1)$series
  # introduce a real NA cell, as a masked/incomplete-series cell would produce
  r[1] <- NA
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  fdr_result <- fdr_correction(trend$stats$p, method = "BH", report = FALSE, verbose = FALSE)

  direction <- direction_map(trend$stats, fdr_result, method = "BH", verbose = FALSE)
  vals <- terra::values(direction, mat = FALSE)
  expect_true(anyNA(vals))
  expect_true(all(vals %in% c(-1, 0, 1, NA)))
})

test_that("fdr_direction_summary can write a CSV", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  fdr_result <- fdr_correction(trend$stats$p, method = c("BH", "BKY"), report = FALSE, verbose = FALSE)

  path <- tempfile(fileext = ".csv")
  suppressMessages(fdr_direction_summary(trend$stats, fdr_result, path = path))
  expect_true(file.exists(path))
  unlink(path)
})
