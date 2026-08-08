.write_year_files <- function(dir, years, nrow = 4, ncol = 4, prefix = "year", ext = "tif") {
  r <- sim_trend_stack(nrow = nrow, ncol = ncol, n_time = length(years), seed = 1)$series
  for (i in seq_along(years)) {
    terra::writeRaster(r[[i]], file.path(dir, sprintf("%s%d.%s", prefix, years[i], ext)), overwrite = TRUE)
  }
  invisible(dir)
}

test_that("read_ordered_stack auto-detects a simple year pattern and orders correctly", {
  dir <- tempfile("sptrends_ros_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  .write_year_files(dir, years = c(3, 1, 2))  # written out of order on purpose

  s <- read_ordered_stack(dir, report = FALSE, verbose = FALSE)
  expect_equal(terra::nlyr(s), 3)
  expect_named(s, c("year1", "year2", "year3"))
})

test_that("read_ordered_stack reads a non-GeoTIFF terra raster format by default", {
  dir <- tempfile("sptrends_ros_grd_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  .write_year_files(dir, years = c(3, 1, 2), ext = "grd")

  s <- read_ordered_stack(dir, report = FALSE, verbose = FALSE)
  expect_equal(terra::nlyr(s), 3)
  expect_named(s, c("year1", "year2", "year3"))
})

test_that("read_ordered_stack messages the auto-detected pattern when verbose = TRUE", {
  dir <- tempfile("sptrends_ros_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  .write_year_files(dir, years = c(1, 2, 3))

  expect_message(read_ordered_stack(dir, report = FALSE, verbose = TRUE), "auto-detected")
})

test_that("read_ordered_stack accepts an explicit order_regex", {
  dir <- tempfile("sptrends_ros_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  .write_year_files(dir, years = c(2001, 2002, 2003), prefix = "obs_")

  s <- read_ordered_stack(dir, order_regex = "obs_([0-9]+)", report = FALSE, verbose = FALSE)
  expect_equal(terra::nlyr(s), 3)
})

test_that("read_ordered_stack errors when the supplied order_regex fails to extract a unique number", {
  dir <- tempfile("sptrends_ros_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  .write_year_files(dir, years = c(1, 2, 3))

  expect_error(
    read_ordered_stack(dir, order_regex = "not_present_([0-9]+)", report = FALSE, verbose = FALSE),
    "NOT VERIFIABLE"
  )
})

test_that("read_ordered_stack errors when no automatic pattern extracts a unique number", {
  dir <- tempfile("sptrends_ros_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  # file names with no digits at all -- every candidate pattern requires
  # extracting at least one digit, so none of them can succeed.
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 3, seed = 1)$series
  for (nm in c("alpha", "beta", "gamma")) {
    terra::writeRaster(r[[which(nm == c("alpha", "beta", "gamma"))]],
                        file.path(dir, paste0(nm, ".tif")), overwrite = TRUE)
  }

  expect_error(read_ordered_stack(dir, report = FALSE, verbose = FALSE), "NOT VERIFIABLE")
})

test_that("read_ordered_stack errors with no matching files in the folder", {
  dir <- tempfile("sptrends_ros_empty_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  expect_error(read_ordered_stack(dir, report = FALSE, verbose = FALSE), "No files matching")
})

test_that("read_ordered_stack notes gaps in a non-consecutive sequence when verbose = TRUE", {
  dir <- tempfile("sptrends_ros_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  .write_year_files(dir, years = c(1, 2, 5))  # gap between 2 and 5

  expect_message(read_ordered_stack(dir, report = FALSE, verbose = TRUE), "gaps in the numbering")
})

test_that("read_ordered_stack builds the stack silently when verbose = FALSE", {
  dir <- tempfile("sptrends_ros_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  .write_year_files(dir, years = c(1, 2, 3))

  expect_silent(read_ordered_stack(dir, report = FALSE, verbose = FALSE))
})

test_that("read_ordered_stack de-duplicates layer names that genuinely collide after tools::file_path_sans_ext(), even though order_regex extracted unique numbers", {
  # A real, constructible collision: tools::file_path_sans_ext() only
  # strips the LAST "." onward, so two files named "sample.v1" and
  # "sample.v2" both strip to the identical stem "sample" -- while a
  # custom order_regex matching digits specifically inside that final
  # ".vN" segment extracts genuinely unique numbers (1, 2) from each,
  # passing try_regex()'s own uniqueness check well before layer names
  # are ever built. filetype is passed explicitly to terra::writeRaster()
  # since ".v1"/".v2" are not real GeoTIFF extensions terra could
  # otherwise auto-detect from.
  dir <- tempfile("sptrends_ros_collision_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 2, seed = 1)$series
  terra::writeRaster(r[[1]], file.path(dir, "sample.v1"), filetype = "GTiff",
                      overwrite = TRUE)
  terra::writeRaster(r[[2]], file.path(dir, "sample.v2"), filetype = "GTiff",
                      overwrite = TRUE)

  s <- read_ordered_stack(dir, pattern = "\\.v[0-9]$",
                           order_regex = "v([0-9]+)$",
                           report = FALSE, verbose = FALSE)
  expect_equal(terra::nlyr(s), 2)
  # Both layers would otherwise be named "sample" -- verify they were
  # actually made unique, not silently overwriting one another.
  expect_equal(length(unique(names(s))), 2)
  expect_true(all(grepl("^sample", names(s))))
})

test_that("read_ordered_stack messages about the same name collision when verbose = TRUE", {
  dir <- tempfile("sptrends_ros_collision_msg_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 2, seed = 2)$series
  terra::writeRaster(r[[1]], file.path(dir, "sample.v1"), filetype = "GTiff",
                      overwrite = TRUE)
  terra::writeRaster(r[[2]], file.path(dir, "sample.v2"), filetype = "GTiff",
                      overwrite = TRUE)

  expect_message(
    read_ordered_stack(dir, pattern = "\\.v[0-9]$", order_regex = "v([0-9]+)$",
                        report = FALSE, verbose = TRUE),
    "duplicate layer names"
  )
})

test_that("read_netcdf_stack errors when the file does not exist", {
  expect_error(read_netcdf_stack(tempfile(fileext = ".nc"), verbose = FALSE), "File not found")
})

test_that("read_netcdf_stack errors on a multi-variable file with var = NULL", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  r1 <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 3, seed = 1)$series
  r2 <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 3, seed = 2)$series
  terra::time(r1) <- as.Date("2000-01-01") + (0:2) * 365
  terra::time(r2) <- as.Date("2000-01-01") + (0:2) * 365
  names(r1) <- rep("var_a", 3)
  names(r2) <- rep("var_b", 3)
  ds <- terra::sds(list(var_a = r1, var_b = r2))
  terra::writeCDF(ds, path, overwrite = TRUE)

  expect_error(read_netcdf_stack(path, verbose = FALSE), "has 2 variables")
  unlink(path)
})

test_that("read_netcdf_stack reads the requested variable correctly", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  r1 <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 3, seed = 1)$series
  r2 <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 3, seed = 2)$series
  terra::time(r1) <- as.Date("2000-01-01") + (0:2) * 365
  terra::time(r2) <- as.Date("2000-01-01") + (0:2) * 365
  names(r1) <- rep("var_a", 3)
  names(r2) <- rep("var_b", 3)
  ds <- terra::sds(list(var_a = r1, var_b = r2))
  terra::writeCDF(ds, path, overwrite = TRUE)

  s <- read_netcdf_stack(path, var = "var_b", report = FALSE, verbose = FALSE)
  expect_equal(terra::nlyr(s), 3)
  unlink(path)
})

test_that("read_netcdf_stack errors on a non-increasing time sequence", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 3, seed = 1)$series
  terra::time(r) <- as.Date(c("2000-01-01", "2000-01-01", "2002-01-01"))  # duplicate date
  terra::writeCDF(r, path, varname = "trend_var", overwrite = TRUE)

  expect_error(read_netcdf_stack(path, report = FALSE, verbose = FALSE), "INVALID TEMPORAL ORDER")
  unlink(path)
})

test_that("read_ordered_stack reads a specific NetCDF variable via 'var' across multiple files", {
  skip_if_not_installed("ncdf4")
  dir <- tempfile("sptrends_ros_var_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  for (yr in 1:3) {
    r1 <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 1, seed = yr)$series
    r2 <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 1, seed = yr + 10)$series
    names(r1) <- "var_a"
    names(r2) <- "var_b"
    ds <- terra::sds(list(var_a = r1, var_b = r2))
    terra::writeCDF(ds, file.path(dir, sprintf("year%d.nc", yr)), overwrite = TRUE)
  }

  s <- read_ordered_stack(dir, var = "var_b", report = FALSE, verbose = FALSE)
  expect_equal(terra::nlyr(s), 3)
})

test_that("read_netcdf_stack errors when the file has no usable time dimension", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 3, seed = 1)$series
  # deliberately leave terra::time() unset (all NA)
  terra::writeCDF(r, path, varname = "trend_var", overwrite = TRUE)

  expect_error(read_netcdf_stack(path, report = FALSE, verbose = FALSE), "no usable time dimension")
  unlink(path)
})

test_that("read_ordered_stack draws the diagnostic plot when report = TRUE", {
  dir <- tempfile("sptrends_ros_plot_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  .write_year_files(dir, years = c(1, 2, 3))

  expect_error(read_ordered_stack(dir, report = TRUE, verbose = FALSE), NA)
})

test_that("read_ordered_stack accepts a custom candidate_regex list", {
  dir <- tempfile("sptrends_ros_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  .write_year_files(dir, years = c(1, 2, 3), prefix = "custom_pattern_")

  s <- read_ordered_stack(dir, candidate_regex = "custom_pattern_([0-9]+)", report = FALSE, verbose = FALSE)
  expect_equal(terra::nlyr(s), 3)
})

test_that("read_ordered_stack auto-detects the MODIS-style '_A' date pattern", {
  dir <- tempfile("sptrends_ros_modis_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 3, seed = 1)$series
  for (i in 1:3) {
    terra::writeRaster(r[[i]], file.path(dir, sprintf("MOD13Q1_A200%d065.tif", i)), overwrite = TRUE)
  }

  s <- read_ordered_stack(dir, report = FALSE, verbose = FALSE)
  expect_equal(terra::nlyr(s), 3)
})

test_that("read_ordered_stack rejects a candidate pattern that extracts duplicate numbers", {
  dir <- tempfile("sptrends_ros_dup_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  # different file names, but both contain only the single digit "1" --
  # the generic "any run of digits" candidate extracts "1" from both,
  # a duplicate, so that candidate is rejected too (not just the ones
  # with no digits at all) before the overall function gives up.
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 2, seed = 1)$series
  terra::writeRaster(r[[1]], file.path(dir, "sceneA1.tif"), overwrite = TRUE)
  terra::writeRaster(r[[2]], file.path(dir, "sceneB1.tif"), overwrite = TRUE)

  expect_error(read_ordered_stack(dir, report = FALSE, verbose = FALSE), "NOT VERIFIABLE")
})

test_that("read_netcdf_stack prints the order-verification table and final message when verbose = TRUE", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 4, seed = 1)$series
  terra::time(r) <- as.Date("2000-01-01") + (0:3) * 365
  terra::writeCDF(r, path, varname = "trend_var", overwrite = TRUE)

  expect_message(
    read_netcdf_stack(path, report = FALSE, verbose = TRUE),
    "Stack built from NetCDF"
  )
  unlink(path)
})

test_that("read_netcdf_stack draws the diagnostic plot when report = TRUE", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 4, seed = 1)$series
  terra::time(r) <- as.Date("2000-01-01") + (0:3) * 365
  terra::writeCDF(r, path, varname = "trend_var", overwrite = TRUE)

  expect_error(read_netcdf_stack(path, report = TRUE, verbose = FALSE), NA)
  unlink(path)
})

test_that("read_netcdf_stack handles a numeric (non-Date) time dimension", {
  skip_if_not_installed("ncdf4")
  path <- tempfile(fileext = ".nc")
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 4, seed = 1)$series
  terra::time(r, tstep = "years") <- 2000:2003  # plain numeric years, not a Date/POSIXct
  terra::writeCDF(r, path, varname = "trend_var", overwrite = TRUE)

  expect_error(read_netcdf_stack(path, report = FALSE, verbose = FALSE), NA)
  unlink(path)
})

test_that("read_ordered_stack() stores the detected order numbers as terra::time() metadata, not just layer names", {
  dir <- tempfile("sptrends_ros_time_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  .write_year_files(dir, years = c(1982, 1983, 1984))

  s <- read_ordered_stack(dir, report = FALSE, verbose = FALSE)
  expect_true(terra::has.time(s))
  expect_equal(as.numeric(terra::time(s)), c(1982, 1983, 1984))
})

test_that("read_netcdf_stack()'s sub-annual warning falls back to the generic 'sub-annual' wording when gap_days cannot be computed (single time step, mocked monthly step classification)", {
  skip_if_not_installed("ncdf4")
  # A single time step gives diff(ordered_time) length 0, so gap_days
  # (its median) is genuinely NA regardless of the actual date used --
  # this isolates the fallback branch (no "~N days apart" text possible)
  # without depending on terra::timeInfo()'s own real-world date-parsing
  # heuristics for what actually counts as "monthly" data, which this
  # package does not control. terra::timeInfo() is mocked to report
  # step = "months" specifically, the same established pattern already
  # used for mocking parallel::detectCores() elsewhere in this suite.
  path <- tempfile(fileext = ".nc")
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 1, seed = 1)$series
  terra::time(r) <- as.Date("2000-01-01")
  terra::writeCDF(r, path, varname = "trend_var", overwrite = TRUE)

  testthat::local_mocked_bindings(
    timeInfo = function(...) list(step = "months"),
    .package = "terra"
  )
  expect_warning(
    read_netcdf_stack(path, report = FALSE, verbose = FALSE),
    "sub-annual"
  )
  unlink(path)
})

test_that("read_ordered_stack()'s INVALID TEMPORAL ORDER safety net triggers correctly when the internal permutation step misbehaves (mocking this package's own .order_for_stacking(), not base::order() itself)", {
  # order_numbers is already guaranteed unique and non-NA by try_regex()
  # before .order_for_stacking() is ever called -- under a correct
  # permutation this branch cannot trigger, which is exactly why it is
  # a safety net rather than normal control flow. Exercised here by
  # mocking this package's own internal wrapper to return a
  # deliberately wrong (identity) permutation instead of the real sort
  # order, rather than mocking base::order() itself (used throughout R
  # internally, and not safely mockable without risking the whole test
  # run) or trying to construct a "naturally" invalid sequence, which
  # the uniqueness guarantee upstream makes impossible by construction.
  #
  # Years deliberately span 1-11 (not just single digits): list.files()
  # returns files in *alphabetical* order regardless of write order, and
  # for single-digit years (e.g. 1, 2, 3) alphabetical and numerical
  # order happen to coincide -- order_numbers would already arrive
  # correctly sorted, making the identity mock below a no-op that
  # accidentally still "worked" rather than a genuine wrong permutation.
  # With "year10"/"year11" in the mix, list.files()'s own alphabetical
  # order ("year1", "year10", "year11", "year2", ...) is genuinely NOT
  # numerically sorted, so identity really is the wrong permutation here.
  dir <- tempfile("sptrends_ros_badorder_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  .write_year_files(dir, years = 1:11)

  testthat::local_mocked_bindings(
    .order_for_stacking = function(x) seq_along(x),  # wrong: identity, not sorted
    .package = "sptrends"
  )
  expect_error(
    read_ordered_stack(dir, report = FALSE, verbose = FALSE),
    "INVALID TEMPORAL ORDER"
  )
})
