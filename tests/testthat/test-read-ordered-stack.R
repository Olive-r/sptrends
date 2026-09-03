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

# ================================================================
# Modo explicito: files / time / cycle_type -- diseñado en sesion
# con revision externa (7 rondas), verificado con calculo directo
# antes de implementar cada pieza delicada (semimonthly/10day y
# los limites de mes/bisiesto, la formula de centro con duracion
# par, la deteccion de series incompletas con 'end').
# ================================================================

.write_plain_files <- function(dir, n, nrow = 4, ncol = 4, prefix = "f") {
  r <- sim_trend_stack(nrow = nrow, ncol = ncol, n_time = n, seed = 1)$series
  paths <- file.path(dir, sprintf("%s%02d.tif", prefix, seq_len(n)))
  for (i in seq_len(n)) {
    terra::writeRaster(r[[i]], paths[i], overwrite = TRUE)
  }
  paths
}

test_that("read_ordered_stack requires exactly one of 'dir' or 'files'", {
  expect_error(read_ordered_stack(), "exactly one of 'dir' or 'files'")
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- .write_plain_files(dir, 2)
  expect_error(
    read_ordered_stack(dir = dir, files = files, report = FALSE,
                       verbose = FALSE),
    "exactly one of 'dir' or 'files'"
  )
})

test_that("read_ordered_stack rejects incompatible argument combinations", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- .write_plain_files(dir, 2)

  expect_error(
    read_ordered_stack(files = files, time = Sys.Date(),
                       cycle_type = "monthly", start = Sys.Date(),
                       report = FALSE, verbose = FALSE),
    "mutually exclusive"
  )
  expect_error(
    read_ordered_stack(files = files, cycle_type = "monthly",
                       report = FALSE, verbose = FALSE),
    "requires 'start'"
  )
  expect_error(
    read_ordered_stack(files = files, time = Sys.Date() + 0:1,
                       report = FALSE, verbose = FALSE),
    NA
  )
  expect_error(
    read_ordered_stack(files = files,
                       start = as.Date("2001-01-01"),
                       report = FALSE, verbose = FALSE),
    "only applicable with 'cycle_type'"
  )
  expect_error(
    read_ordered_stack(files = files, cycle_type = "monthly",
                       start = as.Date("2001-01-01"),
                       time_anchor = "middle",
                       report = FALSE, verbose = FALSE)
  )
})

test_that("read_ordered_stack rejects 'cycle_type'/'time' given with 'dir' instead of 'files', and rejects a 'time_anchor' given without 'cycle_type' in the automatic-mode path", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  .write_year_files(dir, years = c(1, 2, 3))

  expect_error(
    read_ordered_stack(dir, cycle_type = "monthly",
                       start = as.Date("2001-01-01"),
                       report = FALSE, verbose = FALSE),
    "requires 'files'"
  )
  expect_error(
    read_ordered_stack(dir, time = Sys.Date(), report = FALSE,
                       verbose = FALSE),
    "requires 'files'"
  )
  expect_error(
    read_ordered_stack(dir, time_anchor = "centre", report = FALSE,
                       verbose = FALSE),
    "only applicable with 'cycle_type'"
  )
})

test_that("read_ordered_stack rejects 'time_anchor' explicitly set to a non-default value with a cycle_type that has no start/centre/end ambiguity ('annual'/'daily')", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- .write_plain_files(dir, 1)

  expect_error(
    read_ordered_stack(files = files, cycle_type = "annual",
                       start = as.Date("2001-01-01"),
                       time_anchor = "end",
                       report = FALSE, verbose = FALSE),
    "only applies to"
  )
  expect_error(
    read_ordered_stack(files = files, cycle_type = "daily",
                       start = as.Date("2001-01-01"),
                       time_anchor = "end",
                       report = FALSE, verbose = FALSE),
    "only applies to"
  )
})

test_that("read_ordered_stack validates 'start' and 'end' as single non-missing Date values, and rejects 'end' before 'start'", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- .write_plain_files(dir, 1)

  for (bad in list("2001-01-01", 1L, NA, as.Date(c("2001-01-01", "2001-02-01")),
                   as.Date(NA))) {
    expect_error(
      read_ordered_stack(files = files, cycle_type = "monthly",
                         start = bad, report = FALSE, verbose = FALSE),
      "'start' must be one non-missing Date"
    )
  }
  for (bad in list("2001-01-01", 1L, as.Date(c("2001-01-01", "2001-02-01")),
                   as.Date(NA))) {
    expect_error(
      read_ordered_stack(files = files, cycle_type = "monthly",
                         start = as.Date("2001-01-01"), end = bad,
                         report = FALSE, verbose = FALSE),
      "'end' must be one non-missing Date"
    )
  }
  expect_error(
    read_ordered_stack(files = files, cycle_type = "monthly",
                       start = as.Date("2001-02-01"),
                       end = as.Date("2001-01-01"),
                       report = FALSE, verbose = FALSE),
    "must not be before"
  )
})

test_that("read_ordered_stack(files) alone gives ordinal time (1, 2, 3...)", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- .write_plain_files(dir, 3)

  s <- read_ordered_stack(files = files, report = FALSE, verbose = FALSE)
  expect_equal(terra::nlyr(s), 3)
  expect_equal(as.numeric(terra::time(s)), c(1, 2, 3))
})

test_that("read_ordered_stack(files, time) validates time fully -- length, class, NA, duplicates, strictly increasing", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- .write_plain_files(dir, 3)
  fechas_ok <- as.Date(c("2001-01-01", "2001-02-01", "2001-03-01"))

  expect_error(
    read_ordered_stack(files = files, time = fechas_ok[1:2],
                       report = FALSE, verbose = FALSE),
    "layers in"
  )
  expect_error(
    read_ordered_stack(files = files,
                       time = c(fechas_ok[1], NA, fechas_ok[3]),
                       report = FALSE, verbose = FALSE),
    "must not contain NA"
  )
  expect_error(
    read_ordered_stack(files = files,
                       time = c(fechas_ok[1], fechas_ok[1], fechas_ok[3]),
                       report = FALSE, verbose = FALSE),
    "duplicates"
  )
  expect_error(
    read_ordered_stack(files = files,
                       time = rev(fechas_ok),
                       report = FALSE, verbose = FALSE),
    "strictly increasing"
  )

  s <- read_ordered_stack(files = files, time = fechas_ok,
                          report = FALSE, verbose = FALSE)
  expect_equal(as.Date(terra::time(s)), fechas_ok)
})

test_that("read_ordered_stack(files, cycle_type, start) generates correct calendar dates -- monthly, semimonthly, 10day (with real leap-year arithmetic)", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  files_12 <- .write_plain_files(dir, 12, prefix = "m")
  s_m <- read_ordered_stack(files = files_12, cycle_type = "monthly",
                            start = as.Date("1982-01-01"),
                            time_anchor = "start",
                            report = FALSE, verbose = FALSE)
  expect_equal(as.Date(terra::time(s_m)),
              seq(as.Date("1982-01-01"), by = "month", length.out = 12))

  files_4 <- .write_plain_files(dir, 4, prefix = "sm")
  s_sm <- read_ordered_stack(files = files_4, cycle_type = "semimonthly",
                             start = as.Date("1982-01-01"),
                             time_anchor = "start",
                             report = FALSE, verbose = FALSE)
  expect_equal(as.Date(terra::time(s_sm)),
              as.Date(c("1982-01-01", "1982-01-16",
                        "1982-02-01", "1982-02-16")))

  # bisiesto real: 2020 es bisiesto -- el segundo periodo
  # semimonthly de febrero debe terminar el 29, no el 28, verificado
  # via time_anchor="end"
  files_2 <- .write_plain_files(dir, 2, prefix = "leap")
  s_leap <- read_ordered_stack(files = files_2, cycle_type = "semimonthly",
                               start = as.Date("2020-02-01"),
                               time_anchor = "end",
                               report = FALSE, verbose = FALSE)
  expect_equal(as.Date(terra::time(s_leap))[2], as.Date("2020-02-29"))
})

test_that("time_anchor='centre' picks the earlier of two central days for an even-duration period -- deterministic tie-break, verified numerically", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- .write_plain_files(dir, 1, prefix = "c")

  # periodo "10-day" dia1-dia10 (10 dias, duracion par): candidatos
  # centrales son dia5 y dia6 -- la regla elige el primero (dia5)
  s <- read_ordered_stack(files = files, cycle_type = "10-day",
                          start = as.Date("1982-01-01"),
                          time_anchor = "centre",
                          report = FALSE, verbose = FALSE)
  expect_equal(as.Date(terra::time(s)), as.Date("1982-01-05"))
})

test_that("read_ordered_stack detects an incomplete series when 'end' is supplied -- the exact bug found and fixed during design (v4 -> v5)", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files_11 <- .write_plain_files(dir, 11, prefix = "incomplete")

  # 11 archivos, pero el rango declarado (año completo, mensual)
  # exige 12 -- debe fallar, no aceptar silenciosamente la serie
  # incompleta (el fallo real que la comprobacion anterior con
  # 'end' no detectaba: solo comprobaba si la ULTIMA fecha superaba
  # 'end', no si el CONTEO coincidia)
  expect_error(
    read_ordered_stack(files = files_11, cycle_type = "monthly",
                       start = as.Date("1982-01-01"),
                       end = as.Date("1982-12-31"),
                       report = FALSE, verbose = FALSE),
    "expects 12"
  )

  # con los 12 completos, si debe pasar
  files_12 <- .write_plain_files(dir, 12, prefix = "complete")
  expect_error(
    read_ordered_stack(files = files_12, cycle_type = "monthly",
                       start = as.Date("1982-01-01"),
                       end = as.Date("1982-12-31"),
                       report = FALSE, verbose = FALSE),
    NA
  )
})

test_that("read_ordered_stack rejects an 'end' that does not fall exactly on the end of a real period -- a real bug found by an external review: 12 monthly layers with end='1982-12-15' (mid-December, not a real month end) used to be silently accepted, since the old check only compared LAYER COUNT, not whether the last generated period actually ended on 'end'", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files_12b <- .write_plain_files(dir, 12, prefix = "midmonth")

  expect_error(
    read_ordered_stack(files = files_12b, cycle_type = "monthly",
                       start = as.Date("1982-01-01"),
                       end = as.Date("1982-12-15"),
                       report = FALSE, verbose = FALSE),
    "does not fall exactly on the end of a period"
  )
})

test_that("read_ordered_stack's gap detection does not crash with 'object pattern_used not found' when order_regex is supplied manually -- a real bug found by an external review: pattern_used was only initialised inside the automatic-detection branch, never in the manual order_regex branch", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 2, seed = 1)$series
  terra::writeRaster(r[[1]], file.path(dir, "pos100.tif"), overwrite = TRUE)
  terra::writeRaster(r[[2]], file.path(dir, "pos205.tif"), overwrite = TRUE)

  expect_error(
    read_ordered_stack(dir, order_regex = "pos([0-9]+)",
                       report = FALSE, verbose = TRUE),
    NA
  )
})

test_that("read_ordered_stack rejects a 'start' that is not a valid boundary for the declared cycle_type -- internally incoherent declaration, not a guess about the real product", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- .write_plain_files(dir, 1)

  expect_error(
    read_ordered_stack(files = files, cycle_type = "monthly",
                       start = as.Date("1982-01-17"),
                       report = FALSE, verbose = FALSE),
    "not a valid period boundary"
  )
  expect_error(
    read_ordered_stack(files = files, cycle_type = "semimonthly",
                       start = as.Date("1982-01-05"),
                       report = FALSE, verbose = FALSE),
    "not a valid period boundary"
  )
  expect_error(
    read_ordered_stack(files = files, cycle_type = "10-day",
                       start = as.Date("1982-01-15"),
                       report = FALSE, verbose = FALSE),
    "not a valid period boundary"
  )
  # "daily" acepta cualquier fecha
  expect_error(
    read_ordered_stack(files = files, cycle_type = "daily",
                       start = as.Date("1982-01-17"),
                       report = FALSE, verbose = FALSE),
    NA
  )
})

test_that("read_ordered_stack(files) validates files itself -- existence, non-empty, no duplicates", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- .write_plain_files(dir, 2)

  expect_error(
    read_ordered_stack(files = character(0), report = FALSE,
                       verbose = FALSE),
    "non-empty"
  )
  expect_error(
    read_ordered_stack(files = c(files, files[1]), report = FALSE,
                       verbose = FALSE),
    "duplicate"
  )
  expect_error(
    read_ordered_stack(files = c(files, "no_existe_este_archivo.tif"),
                       report = FALSE, verbose = FALSE),
    "do not exist"
  )
})

test_that("time_anchor defaults to 'centre', matching Earth Trends Modeler's own convention ('the Julian day is the center of the period') -- not 'start'", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- .write_plain_files(dir, 1)

  s_default <- read_ordered_stack(files = files, cycle_type = "monthly",
                                  start = as.Date("1982-01-01"),
                                  report = FALSE, verbose = FALSE)
  s_centre <- read_ordered_stack(files = files, cycle_type = "monthly",
                                 start = as.Date("1982-01-01"),
                                 time_anchor = "centre",
                                 report = FALSE, verbose = FALSE)
  s_start <- read_ordered_stack(files = files, cycle_type = "monthly",
                                start = as.Date("1982-01-01"),
                                time_anchor = "start",
                                report = FALSE, verbose = FALSE)

  expect_equal(as.Date(terra::time(s_default)), as.Date(terra::time(s_centre)))
  expect_false(isTRUE(all.equal(
    as.Date(terra::time(s_default)), as.Date(terra::time(s_start))
  )))
})

test_that("changing time_anchor's own default to 'centre' does not break calls that never touch it -- the exact bug that broke a vignette: automatic mode (no cycle_type) with time_anchor left untouched must not error, and neither should cycle_type='annual'/'daily' with time_anchor left untouched", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  .write_year_files(dir, years = c(1, 2, 3))

  # el caso exacto que rompio a-getting-started.Rmd
  expect_error(
    read_ordered_stack(dir, report = FALSE, verbose = FALSE),
    NA
  )

  files <- .write_plain_files(dir, 1, prefix = "annual_check")
  expect_error(
    read_ordered_stack(files = files, cycle_type = "annual",
                       start = as.Date("1982-01-01"),
                       report = FALSE, verbose = FALSE),
    NA
  )
  expect_error(
    read_ordered_stack(files = files, cycle_type = "daily",
                       start = as.Date("1982-01-01"),
                       report = FALSE, verbose = FALSE),
    NA
  )
})

test_that("read_ordered_stack(cycle_type = '16-day') generates fixed annual-reset periods, matching real MODIS convention (verified against real MOD13Q1 file names earlier: year-day 17 is the second 16-day period)", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  files_23 <- .write_plain_files(dir, 23, prefix = "d16")
  s16 <- read_ordered_stack(files = files_23, cycle_type = "16-day",
                            start = as.Date("2001-01-01"),
                            time_anchor = "start",
                            report = FALSE, verbose = FALSE)
  dates16 <- as.Date(terra::time(s16))
  expect_equal(dates16[1], as.Date("2001-01-01"))
  expect_equal(dates16[2], as.Date("2001-01-17"))
  expect_equal(length(dates16), 23L)
})

test_that("read_ordered_stack(cycle_type = '8-day'/'weekly') generates the exact counts and boundaries matching Earth Trends Modeler's own table, verified with direct calculation before implementing", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  files_46 <- .write_plain_files(dir, 46, prefix = "d8")
  s8 <- read_ordered_stack(files = files_46, cycle_type = "8-day",
                           start = as.Date("2001-01-01"),
                           time_anchor = "start",
                           report = FALSE, verbose = FALSE)
  dates8 <- as.Date(terra::time(s8))
  expect_equal(dates8[1], as.Date("2001-01-01"))
  expect_equal(dates8[2], as.Date("2001-01-09"))
  expect_equal(length(dates8), 46L)

  files_52 <- .write_plain_files(dir, 52, prefix = "wk")
  swk <- read_ordered_stack(files = files_52, cycle_type = "weekly",
                            start = as.Date("2001-01-01"),
                            time_anchor = "start",
                            report = FALSE, verbose = FALSE)
  dateswk <- as.Date(terra::time(swk))
  expect_equal(dateswk[1], as.Date("2001-01-01"))
  expect_equal(dateswk[2], as.Date("2001-01-08"))
  expect_equal(length(dateswk), 52L)
  # ETM's own number is 52, not 53 -- confirms the last day of the
  # year (31 December) deliberately belongs to no period here.
  swk_end <- read_ordered_stack(files = files_52, cycle_type = "weekly",
                                start = as.Date("2001-01-01"),
                                time_anchor = "end",
                                report = FALSE, verbose = FALSE)
  expect_equal(as.Date(terra::time(swk_end))[52], as.Date("2001-12-30"))
})

test_that("16day's last period of the year respects leap years in its own end date, even though the yearly count stays fixed at 23", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- .write_plain_files(dir, 1, prefix = "leap16")

  s_normal <- read_ordered_stack(
    files = files, cycle_type = "16-day",
    start = as.Date("2001-12-19"), time_anchor = "end",
    report = FALSE, verbose = FALSE
  )
  expect_equal(as.Date(terra::time(s_normal)), as.Date("2001-12-31"))

  s_leap <- read_ordered_stack(
    files = files, cycle_type = "16-day",
    start = as.Date("2020-12-18"), time_anchor = "end",
    report = FALSE, verbose = FALSE
  )
  expect_equal(as.Date(terra::time(s_leap)), as.Date("2020-12-31"))
})

test_that("read_ordered_stack rejects a 'start' that is not one of the 16-day fixed year-day boundaries", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- .write_plain_files(dir, 1)

  expect_error(
    read_ordered_stack(files = files, cycle_type = "16-day",
                       start = as.Date("2001-01-10"),
                       report = FALSE, verbose = FALSE),
    "not a valid period boundary"
  )
  expect_error(
    read_ordered_stack(files = files, cycle_type = "16-day",
                       start = as.Date("2001-01-17"),
                       report = FALSE, verbose = FALSE),
    NA
  )
})

test_that("automatic mode's gap detection is scoped to simple annual sequences, not combined codes -- the exact false-positive risk found during design", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  # patron generico (no un año simple de 4 digitos): la brecha entre
  # valores no deberia dispararse como "gap" -- solo un aviso
  # informativo, nunca un error duro
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 2, seed = 1)$series
  terra::writeRaster(r[[1]], file.path(dir, "code_100.tif"), overwrite = TRUE)
  terra::writeRaster(r[[2]], file.path(dir, "code_205.tif"), overwrite = TRUE)

  expect_message(
    read_ordered_stack(dir, report = FALSE, verbose = TRUE),
    "not evenly spaced"
  )
})

test_that("read_ordered_stack(files, var) reads a specific NetCDF variable in explicit mode too", {
  skip_if_not_installed("ncdf4")
  dir <- tempfile("sptrends_ros_var_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))

  files <- character(3)
  for (yr in 1:3) {
    r1 <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 1, seed = yr)$series
    r2 <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 1, seed = yr + 10)$series
    names(r1) <- "var_a"; names(r2) <- "var_b"
    ds <- terra::sds(list(var_a = r1, var_b = r2))
    files[yr] <- file.path(dir, sprintf("year%d.nc", yr))
    terra::writeCDF(ds, files[yr], overwrite = TRUE)
  }

  s <- read_ordered_stack(files = files, var = "var_b",
                          report = FALSE, verbose = FALSE)
  expect_equal(terra::nlyr(s), 3)
})

test_that("read_ordered_stack(files, time) validates time's class and finiteness", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- .write_plain_files(dir, 2)

  expect_error(
    read_ordered_stack(files = files, time = c("2001-01-01", "2001-02-01"),
                       report = FALSE, verbose = FALSE),
    "must be Date, POSIXct or numeric"
  )
  expect_error(
    read_ordered_stack(files = files, time = c(1, Inf),
                       report = FALSE, verbose = FALSE),
    "must be finite"
  )
})

test_that("read_ordered_stack's explicit mode prints the expected verbose messages, for both files+time and files+cycle_type", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- .write_plain_files(dir, 2)

  expect_message(
    read_ordered_stack(files = files, time = as.Date(c("2001-01-01",
                                                        "2001-02-01")),
                       report = FALSE, verbose = TRUE),
    "treated as authoritative"
  )
  expect_message(
    read_ordered_stack(files = files, cycle_type = "monthly",
                       start = as.Date("2001-01-01"),
                       report = FALSE, verbose = TRUE),
    "Declared series"
  )
  expect_message(
    read_ordered_stack(files = files, report = FALSE, verbose = TRUE),
    "Temporal order verification"
  )
  expect_message(
    read_ordered_stack(files = files, report = FALSE, verbose = TRUE),
    "Stack built"
  )
})

test_that("read_ordered_stack draws the diagnostic plot in explicit mode when report = TRUE", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- .write_plain_files(dir, 2)

  expect_error(
    read_ordered_stack(files = files, report = TRUE, verbose = FALSE),
    NA
  )
})

test_that("read_ordered_stack de-duplicates layer names in explicit mode when two files share the same stripped name", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  subdir_a <- file.path(dir, "a"); dir.create(subdir_a)
  subdir_b <- file.path(dir, "b"); dir.create(subdir_b)
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 2, seed = 1)$series
  terra::writeRaster(r[[1]], file.path(subdir_a, "same.tif"), overwrite = TRUE)
  terra::writeRaster(r[[2]], file.path(subdir_b, "same.tif"), overwrite = TRUE)
  files <- c(file.path(subdir_a, "same.tif"), file.path(subdir_b, "same.tif"))

  expect_message(
    s <- read_ordered_stack(files = files, report = FALSE, verbose = TRUE),
    "duplicate layer names"
  )
  expect_equal(names(s), c("same", "same_1"))
})

test_that("read_ordered_stack(cycle_type = '10-day') correctly rolls through all three within-month periods (day 1, day 11, day 21) across a full month", {
  dir <- tempfile("sptrends_ros_"); dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  files <- .write_plain_files(dir, 6, prefix = "d10roll")

  s <- read_ordered_stack(files = files, cycle_type = "10-day",
                          start = as.Date("2001-01-01"),
                          time_anchor = "start",
                          report = FALSE, verbose = FALSE)
  dates <- as.Date(terra::time(s))
  expect_equal(dates, as.Date(c("2001-01-01", "2001-01-11", "2001-01-21",
                                "2001-02-01", "2001-02-11", "2001-02-21")))
})
