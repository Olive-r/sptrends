test_that("trend_test returns the expected layers", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 1)$series

  trend_cmk <- trend_test(r, method = "CMK", report = FALSE, verbose = FALSE)
  expect_true(all(c("Sm", "VarSm", "p") %in% names(trend_cmk$stats)))

  trend_mk <- trend_test(r, method = "MK", report = FALSE, verbose = FALSE)
  expect_true(all(c("S", "VarS", "p") %in% names(trend_mk$stats)))

  expect_true(all(terra::values(trend_cmk$stats$p, mat = FALSE) >= 0, na.rm = TRUE))
  expect_true(all(terra::values(trend_cmk$stats$p, mat = FALSE) <= 1, na.rm = TRUE))
})

test_that("trend_test returns a classed object recording which variant was run", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 1)$series

  trend_cmk <- trend_test(r, method = "CMK", report = FALSE, verbose = FALSE)
  expect_identical(class(trend_cmk), c("trend_test", "sptrends"))
  expect_true(trend_cmk$neighbourhood)

  trend_mk <- trend_test(r, method = "MK", report = FALSE, verbose = FALSE)
  expect_identical(class(trend_mk), c("trend_test", "sptrends"))
  expect_false(trend_mk$neighbourhood)
})

test_that("trend_test gives identical S with n_cores = 1 and n_cores = 2", {
  .skip_unless_parallel_tests()
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 2)$series

  seq_result <- trend_test(r, method = "MK", n_cores = 1, report = FALSE, verbose = FALSE)
  par_result <- trend_test(r, method = "MK", n_cores = 2, report = FALSE, verbose = FALSE)

  expect_equal(terra::values(seq_result$stats$S, mat = FALSE), terra::values(par_result$stats$S, mat = FALSE))
})

test_that("prepare_cmk_neighbourhood can be reused across calls", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 6, seed = 3)$series
  X <- terra::values(r, mat = TRUE)
  ok <- stats::complete.cases(X)
  nb <- prepare_cmk_neighbourhood(r, ok)

  trend_a <- trend_test(r, precomputed_neighbourhood = nb, report = FALSE, verbose = FALSE)
  trend_b <- trend_test(r, report = FALSE, verbose = FALSE)
  expect_equal(terra::values(trend_a$stats$Sm, mat = FALSE), terra::values(trend_b$stats$Sm, mat = FALSE))
})

test_that("trend_summary's reference alpha is 0.05, not the strictest value, with the default vector", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)

  # default alpha = c(0.1, 0.05, 0.01) -- the printed reference message
  # must use 0.05, not min(alpha) = 0.01 (regression test for a real bug).
  expect_message(trend_summary(trend$stats), "At alpha=0\\.05")
})

test_that(".reference_alpha falls back to the strictest value when 0.05 is absent", {
  expect_identical(sptrends:::.reference_alpha(c(0.1, 0.2)), 0.1)
  expect_identical(sptrends:::.reference_alpha(c(0.1, 0.05, 0.01)), 0.05)
})

test_that("print.cmk runs without error, distinguishes CMK from classic MK, and returns invisibly", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 10)$series
  trend_cmk <- trend_test(r, method = "CMK", report = FALSE, verbose = FALSE)
  out_cmk <- capture.output(ret <- print(trend_cmk))
  expect_true(any(grepl("Contextual Mann-Kendall", out_cmk, fixed = TRUE)))
  expect_identical(ret, trend_cmk)

  trend_mk <- trend_test(r, method = "MK", report = FALSE, verbose = FALSE)
  out_mk <- capture.output(print(trend_mk))
  expect_true(any(grepl("classic Mann-Kendall", out_mk, fixed = TRUE)))

  trend_5 <- trend_test(r, method = "CMK", window_size = 5L,
                        report = FALSE, verbose = FALSE)
  expect_output(print(trend_5), "Contextual Mann-Kendall \\(5x5\\)")

  legacy <- trend_cmk
  legacy$window_size <- NULL
  expect_output(print(legacy), "Contextual Mann-Kendall \\(3x3\\)")
})

test_that("summary.cmk calls trend_summary() and returns its table", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 11)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  tab <- suppressMessages(summary(trend))
  expect_true("alpha" %in% names(tab))
})

test_that("plot.cmk draws the four trend maps via trend_maps()", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 12)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  expect_error(plot(trend), NA)
  expect_error(plot(trend, alpha = 0.1), NA)
})

test_that("trend_test(method = 'ols') returns beta/se_beta/p layers, class 'trend_test', and matches lm() exactly on a known cell", {
  # A raster built from explicit, hand-picked values, deliberately with
  # realistic (non-negligible) noise around the trend at every cell --
  # sim_trend_stack()'s own random draw can occasionally produce a cell
  # whose fit is essentially perfect (near-zero residuals), a genuinely
  # numerically fragile regime where a closed-form matrix computation
  # and R's own QR-decomposition-based lm() can diverge even though
  # neither is wrong (R's own summary.lm() warns "essentially perfect
  # fit: summary may be unreliable" in exactly that case) -- avoided
  # here entirely by controlling the data directly, rather than hoping
  # a given seed never hits it.
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 10, seed = 20)$series
  terra::values(r)[1, ] <- c(1.2, 2.6, 1.9, 3.4, 2.8, 4.5, 3.6, 5.2, 4.4, 6.1)
  trend_ols <- trend_test(r, method = "OLS", report = FALSE, verbose = FALSE)
  expect_identical(class(trend_ols), c("trend_test", "sptrends"))
  expect_true(all(c("beta", "se_beta", "p") %in% names(trend_ols$stats)))
  expect_false(trend_ols$neighbourhood)

  # Cross-check cell 1's own beta/se_beta/p against a direct lm() fit --
  # the closed-form vectorised computation must agree with the textbook
  # formula, not just "look plausible".
  y <- terra::values(r, mat = TRUE)[1, ]
  fit <- stats::lm(y ~ seq_along(y))
  beta_lm <- unname(coef(fit)[2])
  se_lm <- unname(sqrt(diag(vcov(fit)))[2])
  p_lm <- unname(summary(fit)$coefficients[2, 4])

  beta_ols <- terra::values(trend_ols$stats$beta, mat = FALSE)[1]
  se_ols <- terra::values(trend_ols$stats$se_beta, mat = FALSE)[1]
  p_ols <- terra::values(trend_ols$stats$p, mat = FALSE)[1]

  expect_equal(beta_ols, beta_lm, tolerance = 1e-8)
  expect_equal(se_ols, se_lm, tolerance = 1e-8)
  expect_equal(p_ols, p_lm, tolerance = 1e-8)
})

test_that("trend_test(method = 'ols') respects a custom, non-default t", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 6, seed = 21)$series
  t_custom <- c(2000, 2001, 2003, 2004, 2008, 2009)  # irregular spacing
  trend_ols <- trend_test(r, method = "OLS", t = t_custom, report = FALSE,
                           verbose = FALSE)

  y <- terra::values(r, mat = TRUE)[1, ]
  fit <- stats::lm(y ~ t_custom)
  beta_lm <- unname(coef(fit)[2])
  beta_ols <- terra::values(trend_ols$stats$beta, mat = FALSE)[1]
  expect_equal(beta_ols, beta_lm, tolerance = 1e-8)
})

test_that("trend_test(method = 'ols') validates its own inputs and edge cases", {
  r2 <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 2, seed = 22)$series
  expect_error(trend_test(r2, method = "OLS", report = FALSE, verbose = FALSE),
               "at least 3 layers")

  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 6, seed = 23)$series
  expect_error(trend_test(r, method = "OLS", t = 1:5, report = FALSE,
                           verbose = FALSE),
               "must have length")
  expect_error(trend_test(r, method = "OLS", t = rep(1, 6), report = FALSE,
                           verbose = FALSE),
               "must not contain duplicate")
  expect_error(trend_test(r, method = "OLS", t = c(1, 2, 4, 3, 5, 6),
                           report = FALSE, verbose = FALSE),
               "strictly increasing")
  expect_error(trend_test(r, method = "OLS", t = c(1, 2, 3, 4, 5, Inf),
                           report = FALSE, verbose = FALSE),
               "finite")

  # A perfectly constant cell (zero temporal variance): beta = 0 exactly,
  # se_beta = 0 (zero residual variance too), p should be 1 (definitely
  # not significant), not NaN from a literal 0/0.
  r_const <- r
  terra::values(r_const)[1, ] <- 5
  trend_const <- trend_test(r_const, method = "OLS", report = FALSE,
                             verbose = FALSE)
  beta1 <- terra::values(trend_const$stats$beta, mat = FALSE)[1]
  p1 <- terra::values(trend_const$stats$p, mat = FALSE)[1]
  expect_equal(beta1, 0)
  expect_equal(p1, 1)

  # The other se_beta == 0 edge case: a perfectly linear cell (no noise
  # at all) with a genuinely non-zero slope -- t_stat should come out
  # as +/-Inf, giving p = 0 (definitively significant), not p = 1 or NA
  # from the same 0/0-adjacent branch as the constant-cell case above.
  r_perfect <- r
  terra::values(r_perfect)[2, ] <- seq(1, by = 0.5,
                                        length.out = terra::nlyr(r))
  trend_perfect <- trend_test(r_perfect, method = "OLS", report = FALSE,
                               verbose = FALSE)
  beta2 <- terra::values(trend_perfect$stats$beta, mat = FALSE)[2]
  p2 <- terra::values(trend_perfect$stats$p, mat = FALSE)[2]
  expect_equal(beta2, 0.5, tolerance = 1e-8)
  expect_equal(p2, 0)
})

test_that("print.trend_test correctly labels an OLS result, not 'classic Mann-Kendall'", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 8, seed = 24)$series
  trend_ols <- trend_test(r, method = "OLS", report = FALSE, verbose = FALSE)
  out <- capture.output(print(trend_ols))
  expect_true(any(grepl("OLS", out, fixed = TRUE)))
  expect_false(any(grepl("classic Mann-Kendall", out, fixed = TRUE)))
})

test_that("trend_maps()/trend_summary()/trend_histograms() and direction_map() all recognise the 'beta' layer for an OLS result", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 25)$series
  trend_ols <- trend_test(r, method = "OLS", report = FALSE, verbose = FALSE)
  expect_error(plot(trend_ols), NA)
  expect_error(suppressMessages(summary(trend_ols)), NA)

  fdr_ols <- fdr_correction(trend_ols$stats$p, report = FALSE, verbose = FALSE)
  expect_error(
    sptrends:::direction_map(trend_ols$stats, fdr_ols, verbose = FALSE),
    NA
  )
})

test_that("trend_test(method = 'ols', report = TRUE) exercises the report branch's own trend_summary()/trend_histograms()/trend_maps() calls without error", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 26)$series
  expect_error(
    suppressMessages(trend_test(r, method = "OLS", report = TRUE,
                                 verbose = FALSE)),
    NA
  )
})

test_that("workflow_tst()/workflow_rta() correctly propagate method = 'ols' via cmk_args, exercising rta.R's/tst-methods.R's own 'beta' detection", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 27)$series

  # workflow_tst()/workflow_rta() store $trend as the bare stats
  # SpatRaster directly (out$trend <- cmk_result$stats internally),
  # unlike trend_test()'s own standalone return (a list with its own
  # $stats field) -- verified against workflow_tst.R/rta.R's own source
  # before writing this, not assumed from trend_test()'s own shape.
  result_tst <- workflow_tst(r, prewhiten = FALSE,
                              cmk_args = list(method = "OLS"),
                              report = FALSE, verbose = FALSE)
  expect_true("beta" %in% names(result_tst$trend))
  out_tst <- capture.output(print(result_tst))
  expect_true(any(grepl("beta statistic", out_tst, fixed = TRUE)))
  expect_error(suppressMessages(summary(result_tst)), NA)

  result_rta <- workflow_rta(r, cmk_args = list(method = "OLS"),
                              report = FALSE, verbose = FALSE)
  expect_true("beta" %in% names(result_rta$trend))
  out_rta <- capture.output(print(result_rta))
  expect_true(any(grepl("beta statistic", out_rta, fixed = TRUE)))
  expect_error(suppressMessages(summary(result_rta)), NA)
})

test_that("direction_map() correctly recognises the 'S' layer for a classic MK result (not just 'Sm' or 'beta')", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 28)$series
  trend_mk <- trend_test(r, method = "MK", report = FALSE, verbose = FALSE)
  fdr_mk <- fdr_correction(trend_mk$stats$p, report = FALSE, verbose = FALSE)
  expect_error(
    sptrends:::direction_map(trend_mk$stats, fdr_mk, verbose = FALSE),
    NA
  )
})

test_that("trend_test() errors when x is not a terra SpatRaster", {
  expect_error(trend_test(matrix(1:4, 2, 2), verbose = FALSE),
               "must be a terra SpatRaster")
  expect_error(trend_test(data.frame(a = 1:3), verbose = FALSE),
               "must be a terra SpatRaster")
})

test_that("trend_test() errors when no cell has a complete time series", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 5, seed = 30)$series
  # Force at least one NA into every single cell's own time series --
  # complete.cases() should then be FALSE for all of them.
  X <- terra::values(r, mat = TRUE)
  X[, 1] <- NA
  terra::values(r) <- X
  expect_error(trend_test(r, verbose = FALSE),
               "No cell has a complete time series")
})

test_that("trend_test(method = 'ols', verbose = TRUE) prints its own progress messages", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 50)$series
  expect_message(
    trend_test(r, method = "OLS", report = FALSE, verbose = TRUE),
    "Fitting OLS trend"
  )
})

test_that("print.workflow_tst reports the 'S' statistic name for a classic MK result propagated via cmk_args", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 8, seed = 51)$series
  result <- workflow_tst(r, prewhiten = FALSE,
                          cmk_args = list(method = "MK"),
                          report = FALSE, verbose = FALSE)
  out <- capture.output(print(result))
  expect_true(any(grepl("S statistic", out, fixed = TRUE)))
})
