test_that("fdr_correction on a raster returns q-value and significance rasters for both methods", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  res <- fdr_correction(trend$stats$p, method = c("BH", "BKY"), report = FALSE, verbose = FALSE)

  expect_type(res$rasters, "list")
  expect_true(all(c("p_value", "q_BH", "sig_BH", "q_BKY", "sig_BKY", "sig_raw") %in% names(res$rasters)))
  expect_s4_class(res$rasters$sig_BH, "SpatRaster")
  # the significance raster must match the vector version exactly
  expect_equal(
    sum(terra::values(res$rasters$sig_BH, mat = FALSE), na.rm = TRUE),
    sum(res$reject_BH, na.rm = TRUE)
  )
})

test_that("fdr_correction with a single method only returns that method's columns", {
  p <- c(0.001, 0.01, 0.02, 0.5, 0.8)
  res <- fdr_correction(p, method = "BH", report = FALSE, verbose = FALSE)
  expect_true("q_BH" %in% names(res))
  expect_false("q_BKY" %in% names(res))
})

test_that("fdr_correction's moran_check = TRUE runs the diagnostic on a raster p-value input", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  res <- fdr_correction(trend$stats$p, moran_check = TRUE,
                         moran_args = list(nperm = 19, seed = 1),
                         report = FALSE, verbose = FALSE)
  expect_false(is.null(res$moran))
  expect_true(res$moran_assessment %in% c(
    "positive spatial association compatible with BH",
    "inconclusive for FDR-procedure selection"
  ))
  expect_true(is.na(res$moran_recommendation) ||
                identical(res$moran_recommendation, "BH"))
})

test_that("fdr_correction's moran_check = TRUE does not trigger spatial_autocorrelation()'s own full report (avoids double-reporting)", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 12, seed = 23)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  # spatial_autocorrelation_summary()'s own message text ("Null distribution (...permutations)")
  # must NOT appear -- fdr_correction() has its own, separate moran-related
  # messaging, and spatial_autocorrelation()'s new report = TRUE default
  # must be overridden internally to avoid firing both.
  expect_false(
    any(grepl(
      "Null distribution",
      testthat::capture_messages(
        fdr_correction(trend$stats$p, moran_check = TRUE,
                        moran_args = list(nperm = 19, seed = 1),
                        report = FALSE, verbose = TRUE)
      ),
      fixed = TRUE
    ))
  )
})

test_that("fdr_correction's moran_check = TRUE warns (not errors) on a plain numeric vector", {
  p <- c(0.001, 0.01, 0.02, 0.5, 0.8)
  expect_warning(
    fdr_correction(p, moran_check = TRUE, report = FALSE, verbose = FALSE),
    "moran_check=TRUE ignored"
  )
})

test_that("fdr_summary reports a plausible table for raw/BH/BKY", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  res <- fdr_correction(trend$stats$p, method = c("BH", "BKY"), report = FALSE, verbose = FALSE)
  tab <- suppressMessages(fdr_summary(res))

  expect_s3_class(tab, "data.frame")
  expect_setequal(tab$method, c("raw (uncorrected)", "FDR-BH", "FDR-BKY"))
  # BH/BKY should never call MORE cells significant than the uncorrected raw count
  raw_n <- tab$n_significant[tab$method == "raw (uncorrected)"]
  expect_lte(tab$n_significant[tab$method == "FDR-BH"], raw_n)
  expect_lte(tab$n_significant[tab$method == "FDR-BKY"], raw_n)
})

test_that("fdr_pvalue_histogram runs without error", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  expect_silent(fdr_pvalue_histogram(trend$stats$p))
})

test_that("fdr_significance_maps runs without error", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  res <- fdr_correction(trend$stats$p, method = c("BH", "BKY"), report = FALSE, verbose = FALSE)
  expect_silent(fdr_significance_maps(res))
})

test_that("fdr_comparison_barplot runs without error", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  res <- fdr_correction(trend$stats$p, method = c("BH", "BKY"), report = FALSE, verbose = FALSE)
  expect_silent(fdr_comparison_barplot(res))
})

test_that("fdr_threshold_plot runs without error", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  res <- fdr_correction(trend$stats$p, method = c("BH", "BKY"), report = FALSE, verbose = FALSE)
  expect_silent(fdr_threshold_plot(res))
})

test_that("fdr_direction_plot runs without error", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  res <- fdr_correction(trend$stats$p, method = "BH", report = FALSE, verbose = FALSE)
  direction <- direction_map(trend$stats, res, method = "BH", verbose = FALSE)
  expect_silent(fdr_direction_plot(direction))
})

test_that("fdr_bky's pi0_hat is exactly 1 when nothing looks non-null (all large p-values)", {
  p <- c(0.6, 0.7, 0.8, 0.9, 0.95)
  res <- fdr_bky(p, q = 0.05)
  expect_identical(res$pi0_hat, 1)
  expect_identical(sum(res$reject), 0L)
})

test_that("fdr_correction errors on raster values outside [0,1]", {
  r <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3, ymin = 0, ymax = 3)
  r <- terra::setValues(r, c(0.1, 0.5, 1.5, 0.2, 0.3, -0.1, 0.4, 0.6, 0.9))
  expect_error(
    fdr_correction(r, report = FALSE, verbose = FALSE),
    "Values outside \\[0,1\\]"
  )
})

test_that("fdr_correction uses only the first layer when p has more than one, with a message", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 6, seed = 1)$series
  rng <- terra::global(r, "range", na.rm = TRUE)
  global_min <- min(rng$min)
  global_max <- max(rng$max)
  p_multi <- (r - global_min) / (global_max - global_min)  # rescaled to [0, 1], still 6 layers

  expect_message(
    res <- fdr_correction(p_multi, report = FALSE, verbose = TRUE),
    "more than one layer"
  )
  expect_equal(terra::nlyr(res$rasters$p_value), 1)
})

test_that("fdr_correction's report = TRUE runs the full reporting branch without error", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  expect_error(
    fdr_correction(trend$stats$p, method = c("BH", "BKY"), report = TRUE, verbose = FALSE),
    NA
  )
})

test_that("fdr_correction's moran_check = TRUE message reflects the qualified assessment", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 12, ar1 = 0.6, smooth_radius = 2, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  expect_message(
    res <- fdr_correction(trend$stats$p, moran_check = TRUE, moran_args = list(nperm = 19, seed = 1),
                           report = FALSE, verbose = TRUE),
    "Moran's I ="
  )
  expect_true(res$moran_assessment %in% c(
    "positive spatial association compatible with BH",
    "inconclusive for FDR-procedure selection"
  ))
})

test_that("fdr_correction's moran_check message reports the recommend_bh = TRUE branch", {
  # A raster with strong, smooth spatial structure and a clear trend
  # should give a positive, significant Moran's I on the p-values.
  r <- sim_trend_stack(nrow = 12, ncol = 12, n_time = 15, trend_strength = 0.3,
                        trend_fraction = 1, smooth_radius = 3, noise_sd = 0.5, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  res <- fdr_correction(trend$stats$p, moran_check = TRUE, moran_args = list(nperm = 49, seed = 1),
                         report = FALSE, verbose = FALSE)
  skip_if(res$moran$sign != "positive" || res$moran$p >= 0.05,
          "this seed did not give a significant positive Moran's I -- try a different one")
  expect_message(
    fdr_correction(trend$stats$p, moran_check = TRUE, moran_args = list(nperm = 49, seed = 1),
                   report = FALSE, verbose = TRUE),
    "compatible with the dependence structures"
  )
})

test_that("fdr_correction runs report = TRUE and moran_check = TRUE together", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  expect_error(
    fdr_correction(trend$stats$p, method = c("BH", "BKY"), moran_check = TRUE,
                    moran_args = list(nperm = 19, seed = 1), report = TRUE, verbose = FALSE),
    NA
  )
})

test_that("fdr_summary, fdr_pvalue_histogram, fdr_significance_maps, fdr_comparison_barplot, and fdr_threshold_plot can all write PNG/CSV via path", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 15, trend_strength = 0.3,
                        trend_fraction = 1, noise_sd = 0.5, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  res <- fdr_correction(trend$stats$p, method = c("BH", "BKY"), report = FALSE, verbose = FALSE)

  path_csv <- tempfile(fileext = ".csv")
  suppressMessages(fdr_summary(res, path = path_csv))
  expect_true(file.exists(path_csv))
  unlink(path_csv)

  path_hist <- tempfile(fileext = ".png")
  fdr_pvalue_histogram(trend$stats$p, path = path_hist)
  expect_true(file.exists(path_hist))
  unlink(path_hist)

  path_maps <- tempfile(fileext = ".png")
  fdr_significance_maps(res, path = path_maps)
  expect_true(file.exists(path_maps))
  unlink(path_maps)

  path_bar <- tempfile(fileext = ".png")
  fdr_comparison_barplot(res, path = path_bar)
  expect_true(file.exists(path_bar))
  unlink(path_bar)

  path_thresh <- tempfile(fileext = ".png")
  fdr_threshold_plot(res, path = path_thresh)
  expect_true(file.exists(path_thresh))
  unlink(path_thresh)
})

test_that("fdr_significance_maps errors when result has no rasters (non-raster fdr_correction input)", {
  p <- c(0.01, 0.02, 0.5, 0.8, 0.9)
  res <- fdr_correction(p, report = FALSE, verbose = FALSE)
  expect_error(fdr_significance_maps(res), "must be run on a SpatRaster")
})

test_that("fdr_threshold_plot draws the cutoff line when at least one cell is rejected", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 15, trend_strength = 0.3,
                        trend_fraction = 1, noise_sd = 0.3, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  res <- fdr_correction(trend$stats$p, method = c("BH", "BKY"), report = FALSE, verbose = FALSE)
  skip_if(sum(res$reject_BKY, na.rm = TRUE) == 0, "no significant cells with this seed -- try another")
  expect_silent(fdr_threshold_plot(res))
})

test_that("fdr_threshold_plot handles a family with no rejected tests", {
  res <- fdr_correction(
    rep(0.9, 20),
    method = c("BH", "BKY"),
    report = FALSE,
    verbose = FALSE
  )
  expect_false(any(res$reject_BH))
  expect_false(any(res$reject_BKY))
  expect_silent(fdr_threshold_plot(res))
})

test_that("direction_map reports binarised trend counts when verbose = TRUE", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  res <- fdr_correction(trend$stats$p, method = "BH", report = FALSE, verbose = FALSE)
  expect_message(
    direction_map(trend$stats, res, method = "BH", verbose = TRUE),
    "Binarised trend map \\(BH\\)"
  )
})

test_that("fdr_direction_plot can write a PNG via path", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  res <- fdr_correction(trend$stats$p, method = "BH", report = FALSE, verbose = FALSE)
  direction <- direction_map(trend$stats, res, method = "BH", verbose = FALSE)
  path <- tempfile(fileext = ".png")
  fdr_direction_plot(direction, path = path)
  expect_true(file.exists(path))
  unlink(path)
})

test_that("fdr_direction_summary messages which methods are skipped when not all are available", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  res <- fdr_correction(trend$stats$p, method = "BH", report = FALSE, verbose = FALSE)  # no BKY
  expect_message(fdr_direction_summary(trend$stats, res), "Skipping.*BKY")
})

test_that("fdr_summary works when only BH was run (no summary_bky)", {
  p <- c(0.001, 0.01, 0.02, 0.03, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
  res <- fdr_correction(p, method = "BH", report = FALSE, verbose = FALSE)
  tab <- suppressMessages(fdr_summary(res))
  expect_equal(nrow(tab), 2)  # raw + BH rows, no BKY row
})

test_that("fdr_comparison_barplot works when only BH was run (no summary_bky)", {
  p <- c(0.001, 0.01, 0.02, 0.03, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
  res <- fdr_correction(p, method = "BH", report = FALSE, verbose = FALSE)
  expect_error(fdr_comparison_barplot(res), NA)
})

test_that("fdr_threshold_plot errors when BKY was not run (no threshold_data)", {
  p <- c(0.001, 0.01, 0.02, 0.03, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
  res <- fdr_correction(p, method = "BH", report = FALSE, verbose = FALSE)
  expect_error(fdr_threshold_plot(res), "threshold_data")
})

test_that("fdr_by() matches stats::p.adjust(method = 'BY') exactly", {
  p <- c(0.001, 0.01, 0.02, 0.03, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
  res <- sptrends:::fdr_by(p)
  ref <- stats::p.adjust(p, method = "BY")
  expect_equal(res$q_value, ref)
  expect_equal(res$reject, ref <= 0.05)
})

test_that("fdr_correction()'s own default (no method specified) does NOT compute BY -- it remains deliberately opt-in", {
  p <- c(0.001, 0.01, 0.02, 0.03, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
  res <- fdr_correction(p, report = FALSE, verbose = FALSE)
  expect_true(is.null(res$q_BY))
  expect_true(is.null(res$reject_BY))
  expect_false(is.null(res$q_BH))
  expect_false(is.null(res$q_BKY))
})

test_that("fdr_correction(method = 'BY') computes BY alone, not BH/BKY", {
  p <- c(0.001, 0.01, 0.02, 0.03, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
  res <- fdr_correction(p, method = "BY", report = FALSE, verbose = FALSE)
  expect_false(is.null(res$q_BY))
  expect_true(is.null(res$q_BH))
  expect_true(is.null(res$q_BKY))

  ref <- stats::p.adjust(p, method = "BY")
  expect_equal(res$q_BY, ref)
})

test_that("fdr_correction(method = c('BH', 'BY')) computes both, and BY is at least as conservative as BH (fewer or equal significant cells)", {
  set.seed(1101)
  p <- c(runif(80, 0, 1), runif(20, 0, 0.01))
  res <- fdr_correction(p, method = c("BH", "BY"), report = FALSE,
                         verbose = FALSE)
  expect_false(is.null(res$q_BH))
  expect_false(is.null(res$q_BY))
  expect_true(sum(res$reject_BY, na.rm = TRUE) <=
                sum(res$reject_BH, na.rm = TRUE))
})

test_that("fdr_correction() errors clearly on an invalid method", {
  p <- c(0.01, 0.02, 0.5)
  expect_error(
    fdr_correction(p, method = "invalid", report = FALSE, verbose = FALSE),
    "must be one or more of"
  )
  expect_error(
    fdr_correction(p, method = character(0),
                   report = FALSE, verbose = FALSE),
    "must contain at least one"
  )
})

test_that("fdr_correction validates q and p-value inputs", {
  for (bad_q in list(0, 1, -0.1, 1.1, NA_real_, Inf, c(0.05, 0.1))) {
    expect_error(
      fdr_correction(c(0.01, 0.2), q = bad_q,
                     report = FALSE, verbose = FALSE),
      "'q'"
    )
  }
  expect_error(
    fdr_correction(character(), report = FALSE, verbose = FALSE),
    "numeric vector"
  )
  expect_error(
    fdr_correction(numeric(), report = FALSE, verbose = FALSE),
    "at least one"
  )
  expect_error(
    fdr_correction(c(NA_real_, NA_real_), report = FALSE, verbose = FALSE),
    "no valid p-values"
  )
  expect_error(
    fdr_correction(c(0.1, Inf), report = FALSE, verbose = FALSE),
    "finite"
  )
  expect_error(
    fdr_correction(c(0.1, 1.2), report = FALSE, verbose = FALSE),
    "outside \\[0,1\\]"
  )
})

test_that("fdr_correction(verbose = TRUE) prints the BY-specific message when requested", {
  p <- c(0.001, 0.01, 0.02, 0.03, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
  expect_message(
    fdr_correction(p, method = "BY", report = FALSE, verbose = TRUE),
    "FDR-BY:"
  )
})

test_that("fdr_summary() includes a FDR-BY row when BY was run", {
  p <- c(0.001, 0.01, 0.02, 0.03, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
  res <- fdr_correction(p, method = c("BH", "BY"), report = FALSE,
                         verbose = FALSE)
  tab <- suppressMessages(fdr_summary(res))
  expect_true("FDR-BY" %in% tab$method)
})

test_that("fdr_comparison_barplot() runs without error with BY included, using three panels", {
  p <- c(0.001, 0.01, 0.02, 0.03, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
  res <- fdr_correction(p, method = c("BH", "BKY", "BY"), report = FALSE,
                         verbose = FALSE)
  expect_error(fdr_comparison_barplot(res), NA)
})

test_that("fdr_significance_maps() runs without error with BY included, on raster input, using a dynamic panel count", {
  r_sim <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 10, seed = 1102)$series
  trend <- trend_test(r_sim, report = FALSE, verbose = FALSE)
  res <- fdr_correction(trend$stats$p, method = c("BH", "BY"),
                         report = FALSE, verbose = FALSE)
  expect_error(fdr_significance_maps(res), NA)
})

test_that("fdr_correction() on raster input includes BY's own rasters (q_BY, sig_BY) when requested", {
  r_sim <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 10, seed = 1103)$series
  trend <- trend_test(r_sim, report = FALSE, verbose = FALSE)
  res <- fdr_correction(trend$stats$p, method = "BY", report = FALSE,
                         verbose = FALSE)
  expect_true("q_BY" %in% names(res$rasters))
  expect_true("sig_BY" %in% names(res$rasters))
})

test_that("print() shows FDR-BY's own count -- a real bug found by the user: .print_fdr() only ever checked reject_BH/reject_BKY, so a BY-only object's header correctly said 'BY' but the summary itself showed nothing for it", {
  p <- c(0.001, 0.01, 0.02, 0.03, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95)
  res <- fdr_correction(p, method = "BY", report = FALSE, verbose = FALSE)
  expect_output(print(res), "FDR-BY: \\d+ significant")
})
