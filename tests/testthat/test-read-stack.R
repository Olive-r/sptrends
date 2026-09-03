test_that("read_netcdf_stack warns on monthly (sub-annual) data", {
  skip_if_not_installed("ncdf4")

  path <- tempfile(fileext = ".nc")
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 24, seed = 1)$series
  terra::time(r) <- as.Date("2000-01-15") + round(seq(0, 23) * 30.44)
  terra::writeCDF(r, path, varname = "trend_var", overwrite = TRUE)

  expect_warning(
    read_netcdf_stack(path, report = FALSE, verbose = FALSE),
    "sub-annual"
  )

  unlink(path)
})

test_that("read_netcdf_stack does not warn on yearly data", {
  skip_if_not_installed("ncdf4")

  path <- tempfile(fileext = ".nc")
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 10, seed = 1)$series
  terra::time(r) <- as.Date("2000-06-15") + round(seq(0, 9) * 365.25)
  terra::writeCDF(r, path, varname = "trend_var", overwrite = TRUE)

  warnings_caught <- character(0)
  withCallingHandlers(
    read_netcdf_stack(path, report = FALSE, verbose = FALSE),
    warning = function(w) {
      warnings_caught <<- c(warnings_caught, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warnings_caught, 0)

  unlink(path)
})

test_that("read_netcdf_stack errors clearly when the file does not exist", {
  expect_error(
    read_netcdf_stack(tempfile(fileext = ".nc"), report = FALSE, verbose = FALSE),
    "File not found"
  )
})

test_that("read_netcdf_stack errors clearly when the file has no usable time dimension", {
  skip_if_not_installed("ncdf4")

  path <- tempfile(fileext = ".nc")
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 5, seed = 1)$series
  # Deliberately no terra::time(r) <- ... here -- writeCDF() then has
  # nothing to write into the file's own time dimension, reproducing
  # the "all NA" case this check exists for.
  terra::writeCDF(r, path, varname = "trend_var", overwrite = TRUE)

  expect_error(
    read_netcdf_stack(path, report = FALSE, verbose = FALSE),
    "TEMPORAL ORDER NOT VERIFIABLE"
  )

  unlink(path)
})

test_that("read_netcdf_stack errors clearly when the file's own time values have duplicates (a real, reachable case -- unlike read_ordered_stack()'s file-name-derived numbers, a NetCDF's time dimension is never pre-validated for uniqueness before this check)", {
  skip_if_not_installed("ncdf4")

  path <- tempfile(fileext = ".nc")
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 5, seed = 1)$series
  # Two layers sharing the same date -- a genuinely possible real-world
  # NetCDF authoring mistake (e.g. a processing script re-running the
  # same time step twice), not a contrived setup.
  terra::time(r) <- as.Date(c("2000-01-01", "2000-01-01", "2001-01-01",
                               "2002-01-01", "2003-01-01"))
  terra::writeCDF(r, path, varname = "trend_var", overwrite = TRUE)

  expect_error(
    read_netcdf_stack(path, report = FALSE, verbose = FALSE),
    "INVALID TEMPORAL ORDER"
  )

  unlink(path)
})

test_that("read_netcdf_stack reads back exactly what was written (round-trip fidelity)", {
  skip_if_not_installed("ncdf4")

  path <- tempfile(fileext = ".nc")
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 8, seed = 1)$series
  fechas <- as.Date("2000-06-15") + round(seq(0, 7) * 365.25)
  terra::time(r) <- fechas
  terra::writeCDF(r, path, varname = "trend_var", overwrite = TRUE)

  leido <- read_netcdf_stack(path, report = FALSE, verbose = FALSE)

  expect_equal(terra::nlyr(leido), terra::nlyr(r))
  expect_equal(as.character(terra::time(leido)), as.character(fechas))
  expect_equal(
    terra::values(leido, mat = FALSE),
    terra::values(r, mat = FALSE),
    tolerance = 1e-6
  )

  unlink(path)
})
