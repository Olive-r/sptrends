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

test_that("fdr_direction_summary supports 'BY' when present in fdr_result -- a real gap found by an external audit, not previously tested", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  fdr_result <- fdr_correction(trend$stats$p, method = c("BH", "BKY", "BY"),
                               report = FALSE, verbose = FALSE)

  tab <- suppressMessages(fdr_direction_summary(trend$stats, fdr_result))
  expect_true("BY" %in% tab$method)
  expect_setequal(tab$method, c("raw", "BH", "BKY", "BY"))
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

test_that("direction_map(slope = ...) uses the slope's own sign, not the trend statistic's", {
  # trend_shape="block" con signal_size = todo el raster: senal completa y
  # uniforme, sin el desvanecimiento hacia los bordes de la forma "radial"
  # por defecto -- en un raster tan pequeno (5x5), casi todas las celdas
  # estan cerca del borde, y ese desvanecimiento diluia la senal efectiva
  # incluso con trend_strength alto (fallo real visto en la ronda anterior).
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 12,
                       trend_shape = "block", signal_size = c(5, 5),
                       trend_fraction = 1, trend_strength = 0.3,
                       noise_sd = 0.05, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  fdr_result <- fdr_correction(trend$stats$p, method = "BH", report = FALSE,
                               verbose = FALSE)

  # Direccion por defecto (sin slope): viene del estadistico del test
  direction_stat <- direction_map(trend$stats, fdr_result, method = "raw",
                                  verbose = FALSE)

  # Pendiente artificial con el signo INVERTIDO al del estadistico, en
  # todas las celdas significativas -- si el parametro 'slope' de verdad
  # cambia la fuente del signo, la clasificacion resultante debe ser la
  # opuesta exactamente donde antes habia +1/-1 (0 se queda en 0 siempre,
  # ya que "no significativo" no depende de que fuente se use).
  fake_slope <- -1 * trend$stats$Sm
  direction_slope <- direction_map(trend$stats, fdr_result, slope = fake_slope,
                                   method = "raw", verbose = FALSE)

  vals_stat <- terra::values(direction_stat, mat = FALSE)
  vals_slope <- terra::values(direction_slope, mat = FALSE)
  nonzero <- !is.na(vals_stat) & vals_stat != 0

  expect_true(any(nonzero))
  expect_equal(vals_slope[nonzero], -vals_stat[nonzero])
  # las celdas no significativas (0) deben coincidir exactamente igual,
  # ya que la fuente del signo no afecta a la mascara de significancia
  expect_equal(vals_slope[!nonzero & !is.na(vals_stat)],
               vals_stat[!nonzero & !is.na(vals_stat)])
})

test_that("direction_map() rejects a 'slope' that is not a SpatRaster", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 10, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  fdr_result <- fdr_correction(trend$stats$p, method = "BH", report = FALSE,
                               verbose = FALSE)

  expect_error(
    direction_map(trend$stats, fdr_result, slope = matrix(1:4, 2, 2),
                 method = "BH", verbose = FALSE),
    "'slope' must be a terra SpatRaster"
  )
})
