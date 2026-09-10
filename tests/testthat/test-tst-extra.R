test_that("workflow_tst() errors on non-SpatRaster input", {
  expect_error(workflow_tst(matrix(1:9, 3, 3), verbose = FALSE), "SpatRaster")
})

test_that("workflow_tst() with prewhiten = FALSE skips that step and leaves $prewhiten NULL", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  expect_message(
    result <- workflow_tst(r, prewhiten = FALSE, report = FALSE, verbose = TRUE),
    "prewhitening skipped"
  )
  expect_null(result$prewhiten)
})

test_that("workflow_tst() with export_dw = TRUE includes the prewhitening diagnostics", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, export_dw = TRUE, report = FALSE, verbose = FALSE)
  expect_false(is.null(result$dw_diagnostics))
  expect_true("Modified" %in% names(result$dw_diagnostics))
})

test_that("workflow_tst() with theil_sen = FALSE leaves $theil_sen NULL", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, theil_sen = FALSE, report = FALSE, verbose = FALSE)
  expect_null(result$theil_sen)
})

test_that("workflow_tst() with fdr_method = NULL skips FDR and messages accordingly", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  expect_message(
    result <- workflow_tst(r, fdr_method = NULL, report = FALSE, verbose = TRUE),
    "FDR correction skipped"
  )
  expect_null(result$fdr)
})

test_that("workflow_tst() forwards moran_check through to fdr_correction", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, moran_check = TRUE, report = FALSE, verbose = FALSE)
  expect_false(is.null(result$fdr$moran))
})

test_that("workflow_tst() forwards prewhiten_args, cmk_args, and theil_sen_args correctly", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 20, ar1 = 0.3, seed = 1)$series
  result <- workflow_tst(r,
                 prewhiten_args = list(dw_method = "test"),
                 cmk_args = list(method = "MK"),
                 theil_sen_args = list(max_pairs = 50),
                 report = FALSE, verbose = FALSE)
  # method = "MK" -> classic MK -> trend has "S"/"VarS", not "Sm"/"VarSm"
  expect_true("S" %in% names(result$trend))
  expect_false("Sm" %in% names(result$trend))
})

test_that("workflow_tst() messages every step label when verbose = TRUE", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  msgs <- character(0)
  withCallingHandlers(
    workflow_tst(r, report = FALSE, verbose = TRUE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  full <- paste(msgs, collapse = " ")
  expect_true(grepl("prewhitening", full, fixed = TRUE))
  expect_true(grepl("Contextual Mann-Kendall", full, fixed = TRUE))
  expect_true(grepl("Theil-Sen", full, fixed = TRUE))
  expect_true(grepl("FDR correction", full, fixed = TRUE))
})

test_that("summary.tst runs without error when theil_sen = FALSE", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, theil_sen = FALSE, report = FALSE, verbose = FALSE)
  expect_output(summary(result), "Trend test")
})

test_that("print.tst runs without error when fdr_method = NULL", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, fdr_method = NULL, report = FALSE, verbose = FALSE)
  expect_output(print(result), "True Significant Trends")
})

test_that("plot.tst supports which = 'significance' and which = 'trend'", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, report = FALSE, verbose = FALSE)
  expect_silent(plot(result, which = "significance"))
  expect_silent(plot(result, which = "trend"))
})

test_that("plot.workflow_tst(which = 'significance') also errors without FDR correction", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, fdr_method = NULL, report = FALSE, verbose = FALSE)
  expect_error(plot(result, which = "significance"), "No FDR correction")
})

test_that("print.tst reports the FDR-BH line when workflow_tst() was run with fdr_method = 'BH'", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, fdr_method = "BH", report = FALSE, verbose = FALSE)
  expect_output(print(result), "Significant after FDR-BH")
})

test_that("print.tst reports the Moran's I assessment line when moran_check = TRUE", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, moran_check = TRUE, report = FALSE, verbose = FALSE)
  expect_output(print(result), "Moran's I assessment:")
})

test_that("plot.workflow_tst(method = 'BH') draws the BH direction map when available", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12,
                       trend_shape = "block", signal_size = c(8, 8),
                       trend_fraction = 1, trend_strength = 0.3,
                       noise_sd = 0.3, seed = 1)$series
  result <- workflow_tst(r, fdr_method = "BH", report = FALSE, verbose = FALSE)
  expect_silent(plot(result, method = "BH"))
})

test_that("print.tst reports 'Prewhitening: skipped' when workflow_tst() was run with prewhiten = FALSE", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, prewhiten = FALSE, report = FALSE, verbose = FALSE)
  expect_output(print(result), "Prewhitening: skipped")
})

test_that("plot.workflow_tst(which = 'slope') runs without error, smoothed and unsmoothed", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 12, trend_strength = 0.3,
                        trend_fraction = 1, noise_sd = 0.3, seed = 1)$series
  result <- workflow_tst(r, report = FALSE, verbose = FALSE)
  expect_error(plot(result, which = "slope"), NA)
  expect_error(plot(result, which = "slope", smooth = FALSE), NA)
})

test_that("plot.workflow_tst(which = 'slope') never modifies the underlying tst object", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, report = FALSE, verbose = FALSE)
  theil_sen_before <- terra::values(result$theil_sen)
  plot(result, which = "slope")
  expect_identical(terra::values(result$theil_sen), theil_sen_before)
})

test_that("plot.workflow_tst(which = 'slope') errors without a Theil-Sen slope or FDR correction", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result_no_ts  <- workflow_tst(r, theil_sen = FALSE, report = FALSE, verbose = FALSE)
  result_no_fdr <- workflow_tst(r, fdr_method = NULL, report = FALSE, verbose = FALSE)
  expect_error(plot(result_no_ts, which = "slope"), "No Theil-Sen slope")
  expect_error(plot(result_no_fdr, which = "slope"), "No FDR correction")
})

test_that("plot.workflow_tst(which = 'slope', method = ) errors when that method is unavailable", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, fdr_method = "BH", report = FALSE, verbose = FALSE)
  expect_error(plot(result, which = "slope", method = "BKY"), "not available")
})

test_that("plot.workflow_tst(which = 'direction', method = ) errors when that method is unavailable", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, fdr_method = "BH", report = FALSE, verbose = FALSE)
  expect_error(plot(result, which = "direction", method = "BKY"), "not available")
})

test_that("workflow_tst() does NOT smooth theil_sen by default (aligned with workflow_rta(), neither published method includes this)", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, prewhiten = FALSE, report = FALSE, verbose = FALSE)
  expect_false(result$theil_sen_smoothed)

  slope_raw <- slope_estimator(r, smooth_neighbourhood = FALSE, verbose = FALSE, report = FALSE)$slope
  slope_smoothed <- slope_estimator(r, smooth_neighbourhood = TRUE, verbose = FALSE, report = FALSE)$slope
  expect_equal(terra::values(result$theil_sen), terra::values(slope_raw))
  expect_false(isTRUE(all.equal(terra::values(result$theil_sen), terra::values(slope_smoothed))))
})

test_that("workflow_tst() can be asked for a smoothed theil_sen via theil_sen_args (opt-in, not the default)", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, prewhiten = FALSE,
                 theil_sen_args = list(smooth_neighbourhood = TRUE),
                 report = FALSE, verbose = FALSE)
  expect_true(result$theil_sen_smoothed)
  slope_smoothed <- slope_estimator(r, smooth_neighbourhood = TRUE, verbose = FALSE, report = FALSE)$slope
  expect_equal(terra::values(result$theil_sen), terra::values(slope_smoothed))
})

test_that("plot.workflow_tst(which = 'slope') does not smooth a second time when source smoothing was explicitly requested", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 12, trend_strength = 0.3,
                        trend_fraction = 1, noise_sd = 0.3, seed = 1)$series
  result <- workflow_tst(r, theil_sen_args = list(smooth_neighbourhood = TRUE),
                report = FALSE, verbose = FALSE)
  expect_true(result$theil_sen_smoothed)
  # should run cleanly and not error even though smoothing was already applied
  expect_error(plot(result, which = "slope"), NA)
})

test_that("plot.workflow_tst(which = 'slope', smooth = FALSE) messages when source smoothing was already applied", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, theil_sen_args = list(smooth_neighbourhood = TRUE),
                report = FALSE, verbose = FALSE)
  expect_message(plot(result, which = "slope", smooth = FALSE), "already computed")
})

test_that("plot.workflow_tst(which = 'slope') applies smoothing at plot time by default (workflow_tst() no longer pre-smooths)", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 12, trend_strength = 0.3,
                        trend_fraction = 1, noise_sd = 0.3, seed = 1)$series
  result <- workflow_tst(r, report = FALSE, verbose = FALSE)
  expect_false(result$theil_sen_smoothed)
  expect_error(plot(result, which = "slope"), NA)  # smooth = TRUE by default
})

test_that("plot.workflow_tst(which = 'trend') calls trend_maps() and runs without error", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, report = FALSE, verbose = FALSE)
  expect_error(plot(result, which = "trend"), NA)
})

test_that("plot.workflow_tst() supports the 8 new uncorrected slope_*/pvalue_* views", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 250)$series
  result <- workflow_tst(r, report = FALSE, verbose = FALSE)
  for (w in c("slope_map", "slope_direction", "slope_hist", "slope_bar",
              "pvalue_map", "pvalue_significance", "pvalue_hist", "pvalue_bar")) {
    expect_error(plot(result, which = w), NA, info = w)
  }
})

test_that("plot.workflow_tst(which = 'slope_*') errors clearly when theil_sen = FALSE", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 251)$series
  result <- workflow_tst(r, theil_sen = FALSE, report = FALSE, verbose = FALSE)
  expect_error(plot(result, which = "slope_map"), "theil_sen = TRUE")
  expect_error(plot(result, which = "slope_direction"), "theil_sen = TRUE")
  expect_error(plot(result, which = "slope_hist"), "theil_sen = TRUE")
  expect_error(plot(result, which = "slope_bar"), "theil_sen = TRUE")
})

test_that("slope_direction_map()/slope_direction_barplot() correctly classify positive/negative/zero", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 6, seed = 252)$series[[1]]
  terra::values(r) <- c(1, -1, 0, NA, 2, -2, 0.5, -0.5, 1, 1, -1, -1, 0, 0, NA, 2)
  cls <- sptrends:::slope_direction_map(r)
  vals <- terra::values(cls, mat = FALSE)
  expect_equal(sum(vals == 1, na.rm = TRUE), sum(terra::values(r, mat = FALSE) > 0, na.rm = TRUE))
  expect_equal(sum(vals == -1, na.rm = TRUE), sum(terra::values(r, mat = FALSE) < 0, na.rm = TRUE))

  counts <- sptrends:::slope_direction_barplot(r)
  expect_equal(unname(counts["Positive"]), sum(terra::values(r, mat = FALSE) > 0, na.rm = TRUE))
  expect_equal(unname(counts["Negative"]), sum(terra::values(r, mat = FALSE) < 0, na.rm = TRUE))
})

test_that("slope_direction_barplot()/.plot_pvalue_bar() support probability = TRUE", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 253)$series[[1]]
  expect_error(sptrends:::slope_direction_barplot(r, probability = TRUE), NA)

  trend <- trend_test(
    sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 254)$series,
    report = FALSE, verbose = FALSE
  )
  expect_error(sptrends:::.plot_pvalue_bar(trend$stats, probability = TRUE), NA)
})

test_that(".plot_pvalue_map()/.plot_pvalue_significance()/.plot_pvalue_hist() run without error", {
  trend <- trend_test(
    sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 255)$series,
    report = FALSE, verbose = FALSE
  )
  expect_error(sptrends:::.plot_pvalue_map(trend$stats), NA)
  expect_error(sptrends:::.plot_pvalue_significance(trend$stats), NA)
  expect_error(sptrends:::.plot_pvalue_significance(trend$stats, alpha = 0.1), NA)
  expect_error(sptrends:::.plot_pvalue_hist(trend$stats), NA)
})

test_that("plot.workflow_tst(which = 'direction'/'significance') error clearly when no FDR correction is present", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 40)$series
  result <- workflow_tst(r, prewhiten = FALSE, theil_sen = FALSE,
                          fdr_method = NULL, report = FALSE, verbose = FALSE)
  expect_error(plot(result, which = "direction"), "No FDR correction")
  expect_error(plot(result, which = "significance"), "No FDR correction")
})

test_that("plot.workflow_tst(which = 'slope') errors clearly at each of its 3 missing-piece checks in turn", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 41)$series

  # 1. No Theil-Sen slope at all.
  result_no_slope <- workflow_tst(r, prewhiten = FALSE, theil_sen = FALSE,
                                   report = FALSE, verbose = FALSE)
  expect_error(plot(result_no_slope, which = "slope"), "No Theil-Sen slope")

  # 2. Theil-Sen present, but no FDR correction at all.
  result_no_fdr <- workflow_tst(r, prewhiten = FALSE, theil_sen = TRUE,
                                 fdr_method = NULL, report = FALSE,
                                 verbose = FALSE)
  expect_error(plot(result_no_fdr, which = "slope"), "No FDR correction")

  # 3. FDR present, but only BH was requested -- asking for the
  # (unavailable) BKY-masked slope should name which method is missing.
  result_bh_only <- workflow_tst(r, prewhiten = FALSE, theil_sen = TRUE,
                                  fdr_method = "BH", report = FALSE,
                                  verbose = FALSE)
  expect_error(plot(result_bh_only, which = "slope", method = "BKY"),
               "FDR-BKY is not available")
  # BH itself, on the same object, should work without error.
  expect_error(plot(result_bh_only, which = "slope", method = "BH"), NA)
})

test_that("plot.workflow_tst(which = 'slope_map'/etc.) errors clearly when no Theil-Sen slope is present", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 42)$series
  result <- workflow_tst(r, prewhiten = FALSE, theil_sen = FALSE,
                          report = FALSE, verbose = FALSE)
  expect_error(plot(result, which = "slope_map"), "No Theil-Sen slope")
  expect_error(plot(result, which = "slope_direction"), "No Theil-Sen slope")
  expect_error(plot(result, which = "slope_hist"), "No Theil-Sen slope")
  expect_error(plot(result, which = "slope_bar"), "No Theil-Sen slope")
})

test_that("plot.workflow_tst(which = 'direction') uses the already-smoothed slope as-is when theil_sen_smoothed = TRUE, without re-smoothing", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12,
                       trend_shape = "block", signal_size = c(8, 8),
                       trend_fraction = 1, trend_strength = 0.3,
                       noise_sd = 0.3, seed = 1)$series
  result <- workflow_tst(
    r, fdr_method = "BKY",
    theil_sen_args = list(smooth_neighbourhood = TRUE),
    report = FALSE, verbose = FALSE
  )
  expect_true(isTRUE(result$theil_sen_smoothed))
  expect_error(plot(result, which = "direction"), NA)
})

test_that("direction_map()'s slope input, when smoothed for the direction panel, does not collapse non-significant cells to NA (regression test for a real bug: the grey 'not significant' category silently vanishing)", {
  # Senal solo en la mitad del raster (signal_size cubre menos que el
  # raster completo), garantizando una mezcla real de celdas
  # significativas y no significativas -- necesario para poder
  # detectar si "no significativo" desaparece por error.
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 12,
                       trend_shape = "block", signal_size = c(5, 10),
                       signal_location = "centre",
                       trend_fraction = 1, trend_strength = 0.3,
                       noise_sd = 0.05, seed = 1)$series
  result <- workflow_tst(r, fdr_method = "BKY", report = FALSE,
                         verbose = FALSE)

  reject <- result$fdr$reject_BKY
  n_significativas <- sum(reject, na.rm = TRUE)
  n_no_significativas <- sum(!reject, na.rm = TRUE)

  # Confirmar que el escenario de prueba tiene realmente ambos casos
  # -- si esto fallara, el test en si no probaria nada util
  expect_true(n_significativas > 0)
  expect_true(n_no_significativas > 0)

  slope_suave <- sptrends:::.smooth_slope_for_direction(
    result$theil_sen, TRUE
  )
  direction <- direction_map(result$trend, result$fdr,
                             slope = slope_suave, method = "BKY",
                             verbose = FALSE)
  vals <- terra::values(direction, mat = FALSE)

  # El bug real hacia que TODAS las celdas no significativas
  # colapsaran a NA (desaparecian del mapa) en lugar de mostrarse en 0
  # (gris). No exigimos una coincidencia EXACTA con el recuento crudo
  # de 'reject' -- el suavizado (focal 3x3) puede legitimamente
  # desplazar unas pocas celdas fronterizas de categoria, eso no es
  # el bug que este test debe capturar. Lo que si debe cumplirse
  # siempre: (a) existen celdas en 0 (el gris no desaparece del todo),
  # y (b) el recuento de "no significativo visible" esta en el mismo
  # orden de magnitud que el recuento crudo, no reducido a practicamente
  # cero como pasaria si el bug estuviera presente.
  n_cero <- sum(vals == 0, na.rm = TRUE)
  expect_true(n_cero > 0)
  expect_true(abs(n_cero - n_no_significativas) <= 0.1 * n_no_significativas)
})
