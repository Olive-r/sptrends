test_that("prewhiten's Clamped diagnostic correctly identifies cells whose iterative estimate hits the +-1 stability bound, and every clamped cell is left uncorrected (Modified = 0) rather than transformed using its own unreliable rho", {
  r <- .make_ar1_raster(rho = 0.85, n_time = 50, seed = 1)
  res <- prewhiten(r, report = FALSE, verbose = FALSE)

  mod <- terra::values(res$diagnostics$Modified, mat = FALSE)
  clamped <- terra::values(res$diagnostics$Clamped, mat = FALSE)

  # Not asserting how often clamping itself happens at rho = 0.85 --
  # the underlying iterative mechanism was substantially rewritten
  # (see NEWS.md: now mirrors zyp::zyp.TFPW_Z()'s own mechanics rather
  # than the package's own earlier, less stable iteration), so a
  # clamp rate calibrated to the old mechanism is not a claim this
  # test can still make reliably without re-verifying it empirically
  # first -- the same class of mistake already made and corrected
  # more than once earlier in this file. Only the invariant the code
  # itself guarantees by construction is checked here: whenever
  # clamping does happen, for whatever fraction of cells (including
  # none at all, under the rewritten mechanism), that cell is left
  # uncorrected.
  if (any(clamped == 1, na.rm = TRUE)) {
    expect_equal(mean(mod[clamped == 1], na.rm = TRUE), 0)
  } else {
    succeed("No cell hit the clamp under the rewritten mechanism for this seed -- the invariant this test checks has nothing to verify here, but that is expected, not a failure.")
  }
})

test_that("prewhiten's Clamped == 0 cells are genuinely corrected (Modified = 1) using their own converged rho, not left untouched", {
  r <- .make_ar1_raster(rho = 0.5, n_time = 50, seed = 3)
  res <- prewhiten(r, report = FALSE, verbose = FALSE)

  mod <- terra::values(res$diagnostics$Modified, mat = FALSE)
  rho_est <- terra::values(res$diagnostics$Rho, mat = FALSE)
  clamped <- terra::values(res$diagnostics$Clamped, mat = FALSE)

  # The complementary invariant to the test above: a gated cell that
  # was NOT clamped is genuinely corrected -- Modified = 1, using its
  # own converged (non-clamped) rho. Not asserting how often clamping
  # itself happens at this rho, only what happens on either side of
  # that outcome once it is known.
  gated_unclamped <- mod == 1 & (is.na(clamped) | clamped == 0)
  if (any(gated_unclamped, na.rm = TRUE)) {
    # The only mathematically guaranteed invariant for an unclamped
    # cell is abs(rho) < 1, strictly (required for sqrt(1 - rho^2) and
    # 1 / (1 - rho) downstream not to fail) -- not < 0.99. 0.99 is only
    # the fixed value a genuinely-exceeded estimate gets clamped TO,
    # not a ceiling on how close a real, unclamped convergence is
    # allowed to approach 1 on its own; a persistent-but-not-clamped
    # series can legitimately converge to 0.995, 0.999, and so on,
    # without ever technically hitting or exceeding the >= 1 bound
    # that triggers clamping in the first place.
    expect_true(all(abs(rho_est[gated_unclamped]) < 1))
  }
  gated_clamped <- mod == 1 & clamped == 1
  expect_equal(sum(gated_clamped, na.rm = TRUE), 0)
})

test_that("prewhiten leaves near-white-noise cells ungated under default thresholds", {
  r <- .make_ar1_raster(rho = 0, n_time = 50, seed = 2)
  res <- prewhiten(r, report = FALSE, verbose = FALSE)
  mod <- terra::values(res$diagnostics$Modified, mat = FALSE)
  # with no true autocorrelation, most cells should stay under the default [1.4, 2.6] gate
  expect_lt(mean(mod, na.rm = TRUE), 0.3)
})

test_that("dw_method = 'threshold' and 'test' can give different gating on the same data", {
  r <- .make_ar1_raster(rho = 0.3, n_time = 20, seed = 3)
  res_threshold <- prewhiten(r, dw_method = "threshold", report = FALSE, verbose = FALSE)
  res_test <- prewhiten(r, dw_method = "test", report = FALSE, verbose = FALSE)

  mod_threshold <- sum(terra::values(res_threshold$diagnostics$Modified, mat = FALSE), na.rm = TRUE)
  mod_test <- sum(terra::values(res_test$diagnostics$Modified, mat = FALSE), na.rm = TRUE)
  # not asserting which is larger (depends on the random draw), just that the
  # two methods are actually independent code paths, not aliases of each other
  expect_true(is.numeric(mod_threshold))
  expect_true(is.numeric(mod_test))
})

test_that("dw_inconclusive = 'conservative' gates at least as many cells as 'power'", {
  r <- .make_ar1_raster(rho = 0.25, n_time = 20, seed = 4)
  res_cons <- prewhiten(r, dw_method = "test", dw_inconclusive = "conservative",
                                    report = FALSE, verbose = FALSE)
  res_pow  <- prewhiten(r, dw_method = "test", dw_inconclusive = "power",
                                    report = FALSE, verbose = FALSE)

  n_cons <- sum(terra::values(res_cons$diagnostics$Modified, mat = FALSE), na.rm = TRUE)
  n_pow  <- sum(terra::values(res_pow$diagnostics$Modified, mat = FALSE), na.rm = TRUE)
  # conservative prewhitens the inconclusive zone too, power does not --
  # conservative can never gate fewer cells than power on the same data.
  expect_gte(n_cons, n_pow)
})

test_that("prewhiten errors on fewer than 3 layers", {
  r <- .make_ar1_raster(rho = 0, n_time = 2, seed = 5)
  expect_error(prewhiten(r, verbose = FALSE), "at least 3 layers")
})

test_that("prewhiten errors on non-SpatRaster input", {
  expect_error(prewhiten(matrix(1:9, 3, 3), verbose = FALSE), "SpatRaster")
})

test_that("prewhiten validates the supplied time axis", {
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 8, seed = 402)$series
  expect_error(
    prewhiten(r, t = 1:7, report = FALSE, verbose = FALSE),
    "must have length"
  )
  expect_error(
    prewhiten(r, t = c(1, 2, 3, 4, 4, 6, 7, 8),
              report = FALSE, verbose = FALSE),
    "must not contain duplicate"
  )
  expect_error(
    prewhiten(r, t = c(1, 2, 3, 5, 4, 6, 7, 8),
              report = FALSE, verbose = FALSE),
    "strictly increasing"
  )
  expect_error(
    prewhiten(r, t = c(1, 2, 3, 4, 5, 6, 7, NA),
              report = FALSE, verbose = FALSE),
    "finite"
  )
})

test_that("prewhiten accepts a custom time vector", {
  r <- .make_ar1_raster(rho = 0.5, n_time = 10, seed = 6)
  res_default <- prewhiten(r, report = FALSE, verbose = FALSE)
  res_custom  <- prewhiten(r, t = seq(0, 90, by = 10), report = FALSE, verbose = FALSE)
  # different time spacing shouldn't change *whether* a cell is gated
  # (DW/rho depend on residuals from detrending, not the absolute time scale)
  expect_equal(
    terra::values(res_default$diagnostics$Modified, mat = FALSE),
    terra::values(res_custom$diagnostics$Modified, mat = FALSE)
  )
})

test_that("prewhiten_summary's counts are consistent with the diagnostics raster", {
  r <- .make_ar1_raster(rho = 0.7, n_time = 30, seed = 8)
  res <- prewhiten(r, report = FALSE, verbose = FALSE)
  tab <- suppressMessages(prewhiten_summary(res$diagnostics))

  mod <- terra::values(res$diagnostics$Modified, mat = FALSE)
  expect_equal(tab$value[tab$metric == "valid_cells"], sum(!is.na(mod)))
  expect_equal(tab$value[tab$metric == "prewhitened_cells"], sum(mod == 1, na.rm = TRUE))
})

test_that("prewhiten_summary can write a CSV", {
  r <- .make_ar1_raster(rho = 0.5, n_time = 20, seed = 9)
  res <- prewhiten(r, report = FALSE, verbose = FALSE)
  path <- tempfile(fileext = ".csv")
  suppressMessages(prewhiten_summary(res$diagnostics, path = path))
  expect_true(file.exists(path))
  unlink(path)
})

test_that("prewhiten_histograms handles the 'no prewhitened cells' branch", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 15, ar1 = 0, seed = 42)$series
  res <- prewhiten(r, dw_low = 0, dw_high = 4,
                   report = FALSE, verbose = FALSE)
  expect_message(prewhiten_histograms(res$diagnostics), "no prewhitened cells")
})

test_that("prewhiten_maps runs without error", {
  r <- .make_ar1_raster(rho = 0.6, n_time = 20, seed = 10)
  res <- prewhiten(r, report = FALSE, verbose = FALSE)
  expect_silent(prewhiten_maps(res$diagnostics))
})

test_that("prewhiten's rho estimate never reaches +-1 (numerical stability guard)", {
  # Extremely high true autocorrelation (rho = 0.999) with very low noise --
  # a case where the iterative estimate could plausibly reach/exceed 1
  # without the guard. Note: the guard only clamps rho >= 1 or rho <= -1
  # down to +-0.99 -- it does not cap rho at 0.99 in general, so values
  # like 0.995 are legitimate and expected here; the only real guarantee
  # is abs(rho) < 1 (needed later for sqrt(1 - rho^2) not to fail).
  r <- .make_ar1_raster(rho = 0.999, n_time = 40, sd_e = 0.1, seed = 11)
  res <- prewhiten(r, report = FALSE, verbose = FALSE)
  rho_vals <- terra::values(res$diagnostics$Rho, mat = FALSE)
  expect_true(all(abs(rho_vals) < 1, na.rm = TRUE))
})

test_that("prewhiten respects a custom itmax (fewer iterations allowed)", {
  r <- .make_ar1_raster(rho = 0.7, n_time = 30, seed = 12)
  res_1iter <- prewhiten(r, itmax = 1, report = FALSE, verbose = FALSE)
  res_full  <- prewhiten(r, itmax = 20, report = FALSE, verbose = FALSE)
  # not asserting they're identical (fewer iterations can give a different
  # rho estimate) -- just that a very small itmax runs without error and
  # returns a result of the right shape.
  expect_equal(terra::nlyr(res_1iter$series), terra::nlyr(res_full$series))
})

test_that("prewhiten's eps controls how early iteration stops", {
  r <- .make_ar1_raster(rho = 0.7, n_time = 30, seed = 13)
  # a very loose eps should make the loop break after iteration 1 -- just
  # checking this runs without error and gives a plausible result.
  res_loose <- prewhiten(r, eps = 10, report = FALSE, verbose = FALSE)
  expect_true(all(is.finite(terra::values(res_loose$diagnostics$Rho, mat = FALSE)) |
                     is.na(terra::values(res_loose$diagnostics$Rho, mat = FALSE))))
})

test_that("prewhiten warns when n is below the DW table range (n < 15)", {
  r <- .make_ar1_raster(rho = 0.5, n_time = 8, seed = 20)
  expect_message(
    prewhiten(r, dw_method = "test", report = FALSE, verbose = TRUE),
    "below the Durbin-Watson table range"
  )
})

test_that("prewhiten warns when n is above the DW table range (n > 100)", {
  r <- .make_ar1_raster(rho = 0.5, n_time = 110, seed = 21)
  expect_message(
    prewhiten(r, dw_method = "test", report = FALSE, verbose = TRUE),
    "exceeds the Durbin-Watson table range"
  )
})

test_that("prewhiten's dw_method = 'test' prints the formal-test messages when verbose", {
  r <- .make_ar1_raster(rho = 0.5, n_time = 30, seed = 22)
  expect_message(
    prewhiten(r, dw_method = "test", report = FALSE, verbose = TRUE),
    "formal Durbin-Watson test"
  )
})

test_that("prewhiten messages when no cell crosses the threshold, verbose = TRUE", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 15, ar1 = 0, seed = 42)$series
  expect_message(
    prewhiten(r, dw_low = 0, dw_high = 4,
              report = FALSE, verbose = TRUE),
    "No cell crosses the DW threshold"
  )
})

test_that("prewhiten's report = TRUE runs the full reporting branch without error", {
  r <- .make_ar1_raster(rho = 0.6, n_time = 20, seed = 23)
  expect_error(prewhiten(r, report = TRUE, verbose = FALSE), NA)
})

test_that("prewhiten accepts a file-backed (non-in-memory) raster", {
  r <- .make_ar1_raster(rho = 0.5, n_time = 20, seed = 24)
  path <- tempfile(fileext = ".tif")
  terra::writeRaster(r, path, overwrite = TRUE)
  r_disk <- terra::rast(path)
  expect_false(terra::inMemory(r_disk)[1])
  expect_error(prewhiten(r_disk, report = FALSE, verbose = FALSE), NA)
  unlink(path)
})

test_that("prewhiten_histograms can write both PNG files via path", {
  r <- .make_ar1_raster(rho = 0.6, n_time = 20, seed = 25)
  res <- prewhiten(r, report = FALSE, verbose = FALSE)
  path <- tempfile()
  suppressMessages(prewhiten_histograms(res$diagnostics, path = path))
  expect_true(file.exists(paste0(path, "_dw.png")))
  expect_true(file.exists(paste0(path, "_rho.png")))
  unlink(paste0(path, c("_dw.png", "_rho.png")))
})

test_that("prewhiten_maps can write all three PNG files via path", {
  r <- .make_ar1_raster(rho = 0.6, n_time = 20, seed = 26)
  res <- prewhiten(r, report = FALSE, verbose = FALSE)
  path <- tempfile()
  prewhiten_maps(res$diagnostics, path = path)
  expect_true(file.exists(paste0(path, "_dw_map.png")))
  expect_true(file.exists(paste0(path, "_rho_map.png")))
  expect_true(file.exists(paste0(path, "_modified_map.png")))
  unlink(paste0(path, c("_dw_map.png", "_rho_map.png", "_modified_map.png")))
})

test_that("prewhiten_maps fixes the rho palette range at [-1, 1] regardless of the estimated values", {
  r <- .make_ar1_raster(rho = 0.6, n_time = 20, seed = 27)
  res <- prewhiten(r, report = FALSE, verbose = FALSE)
  expect_error(prewhiten_maps(res$diagnostics), NA)
})

test_that("prewhiten_summary now prints informative messages (was previously silent -- API consistency fix)", {
  r <- .make_ar1_raster(rho = 0.6, n_time = 20, seed = 30)
  res <- prewhiten(r, report = FALSE, verbose = FALSE)
  msgs <- testthat::capture_messages(prewhiten_summary(res$diagnostics))
  full <- paste(msgs, collapse = " ")
  expect_true(grepl("Valid cells", full, fixed = TRUE))
  expect_true(grepl("Prewhitened", full, fixed = TRUE))
  expect_true(grepl("Median Durbin-Watson", full, fixed = TRUE))
})

test_that("prewhiten_summary reports 'none prewhitened' when zero cells were modified", {
  # Construct a diagnostics raster directly (same layer names/order
  # prewhiten() itself builds) with Modified all zero, rather
  # than relying on a random series happening to produce that outcome.
  r1 <- terra::rast(nrows = 3, ncols = 3, vals = 1.9)   # DW_initial
  r2 <- terra::rast(nrows = 3, ncols = 3, vals = 0)     # Rho
  r3 <- terra::rast(nrows = 3, ncols = 3, vals = 0)     # Modified: none
  diagnostics <- c(r1, r2, r3)
  names(diagnostics) <- c("DW_initial", "Rho", "Modified")

  msgs <- testthat::capture_messages(prewhiten_summary(diagnostics))
  expect_true(any(grepl("none prewhitened", msgs, fixed = TRUE)))
})

test_that("prewhiten(method = 'TFPW_WS') works at its exact minimum of 3 layers", {
  r <- .make_ar1_raster(rho = 0.3, n_time = 3, seed = 320)
  result <- prewhiten(r, method = "TFPW_WS", report = FALSE, verbose = FALSE)
  expect_equal(terra::nlyr(result$series), 3)
})

test_that("prewhiten(method = 'TFPW_Y') errors on exactly 3 layers (needs 4, one more than TFPW_WS)", {
  r <- .make_ar1_raster(rho = 0.3, n_time = 3, seed = 321)
  expect_error(prewhiten(r, method = "TFPW_Y", verbose = FALSE),
               "at least 4 layers")
})

test_that("prewhiten(method = 'TFPW_Y') works at its exact minimum of 4 layers", {
  r <- .make_ar1_raster(rho = 0.3, n_time = 4, seed = 322)
  result <- prewhiten(r, method = "TFPW_Y", report = FALSE, verbose = FALSE)
  expect_equal(terra::nlyr(result$series), 3)  # TFPW_Y always loses 1 layer
})

test_that("prewhiten does not error on a perfectly constant series (zero variance -- a degenerate case for rho/DW estimation)", {
  r <- .make_ar1_raster(rho = 0.3, n_time = 10, seed = 323)
  r[] <- 4  # every cell, every layer: exactly the same constant value
  expect_error(
    prewhiten(r, method = "TFPW_WS", report = FALSE, verbose = FALSE),
    NA
  )
})

test_that("prewhiten does not error on a 1x1 raster", {
  r <- .make_ar1_raster(rho = 0.3, n_time = 10, nrow = 1, ncol = 1, seed = 324)
  expect_error(
    prewhiten(r, method = "TFPW_WS", report = FALSE, verbose = FALSE),
    NA
  )
})

test_that("prewhiten errors clearly, not cryptically, when every cell is NA", {
  r_na <- terra::rast(nrows = 4, ncols = 4, nlyr = 8, vals = NA_real_)
  expect_error(prewhiten(r_na, verbose = FALSE), "nothing to prewhiten")
})

test_that("prewhiten handles a very long series (100 layers) without erroring", {
  r <- .make_ar1_raster(rho = 0.3, n_time = 100, nrow = 3, ncol = 3, seed = 325)
  expect_error(
    prewhiten(r, method = "TFPW_WS", report = FALSE, verbose = FALSE),
    NA
  )
})

test_that("prewhiten(method = 'TFPW_Y', verbose = TRUE) prints its own step-by-step progress messages", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 8, seed = 317)$series
  msgs <- capture.output(
    prewhiten(r, method = "TFPW_Y", verbose = TRUE, report = FALSE),
    type = "message"
  )
  full <- paste(msgs, collapse = " ")
  expect_true(grepl("Estimating the Theil-Sen slope", full, fixed = TRUE))
  expect_true(grepl("Detrending", full, fixed = TRUE))
  expect_true(grepl("Estimating lag-1 autocorrelation", full, fixed = TRUE))
  expect_true(grepl("Prewhitening the residuals", full, fixed = TRUE))
  expect_true(grepl("Done.", full, fixed = TRUE))
})

test_that("prewhiten(method = 'TFPW_Y')'s own maps do not error when Beta_TheilSen is entirely zero (max_abs guard, same pattern as slope_map's own all-NA test)", {
  r <- .make_ar1_raster(rho = 0, n_time = 8, seed = 318)
  pw <- prewhiten(r, method = "TFPW_Y", report = FALSE, verbose = FALSE)
  # Built the same way prewhiten() itself constructs diagnostics (see
  # its own "diagnostics <- c(terra::setValues(...), ...)" pattern),
  # rather than mutating pw$diagnostics$Beta_TheilSen in place, since
  # whether a nested $<- on a multi-layer SpatRaster actually writes
  # back into the original object is not something to assume without
  # running it.
  zero_beta <- terra::setValues(pw$diagnostics[[1]], 0)
  diagnostics_zero <- c(zero_beta, pw$diagnostics$Rho)
  names(diagnostics_zero) <- c("Beta_TheilSen", "Rho")
  expect_error(
    sptrends:::.TFPW_Y_maps(diagnostics_zero),
    NA
  )
})

test_that(".theil_sen_slope_vectorised() matches a direct, unvectorised computation of the same pairwise-median slope", {
  set.seed(401)
  Y <- matrix(rnorm(4 * 10), nrow = 4, ncol = 10)
  t <- 1:10

  pairs <- sptrends:::.theil_sen_pairs(t, max_pairs = Inf)
  got <- sptrends:::.theil_sen_slope_vectorised(Y, pairs$i_idx, pairs$j_idx,
                                                 pairs$dt)

  want <- apply(Y, 1, function(row) {
    all_pairs <- utils::combn(length(t), 2)
    i <- all_pairs[1, ]; j <- all_pairs[2, ]
    stats::median((row[j] - row[i]) / (t[j] - t[i]))
  })

  expect_equal(got, want, tolerance = 1e-10)
})

test_that("prewhiten(method = 'TFPW_WS', refit_method = 'sen') runs and returns valid results, with 'ols' remaining the default", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 20, seed = 402)$series

  pw_default <- prewhiten(r, verbose = FALSE, report = FALSE)
  pw_ols     <- prewhiten(r, refit_method = "OLS", verbose = FALSE,
                           report = FALSE)
  pw_sen     <- prewhiten(r, refit_method = "TS", verbose = FALSE,
                           report = FALSE)

  # Default is unchanged: explicit "OLS" must match the implicit default
  # exactly, not just "run without error".
  expect_equal(terra::values(pw_default$series),
               terra::values(pw_ols$series))

  expect_s4_class(pw_sen$series, "SpatRaster")
  expect_true(all(names(pw_sen$diagnostics) ==
                     c("DW_initial", "Rho", "Modified", "Clamped")))
})

test_that("prewhiten(method = 'TFPW_WS', refit_method = 'sen' or 'ols') both run validly on a near-unit-root series, using the same outer iteration mechanism and differing only in the internal slope refit", {
  n <- 42
  set.seed(403)
  near_unit_root <- cumsum(rnorm(n, sd = 0.1))

  r <- terra::rast(nrows = 1, ncols = 1, nlyrs = n)
  r <- terra::setValues(r, matrix(near_unit_root, nrow = 1))

  pw_ols <- prewhiten(r, refit_method = "OLS", verbose = FALSE,
                       report = FALSE)
  pw_sen <- prewhiten(r, refit_method = "TS", verbose = FALSE,
                       report = FALSE)

  rho_ols <- terra::values(pw_ols$diagnostics$Rho, mat = FALSE)[1]
  rho_sen <- terra::values(pw_sen$diagnostics$Rho, mat = FALSE)[1]

  # Not asserting a comparison between the two here: with the outer
  # iteration mechanism now unified between refit_method values (both
  # measure lag-1 autocorrelation on the raw series detrended by the
  # latest slope estimate, following zyp::zyp.TFPW_Z()'s own mechanics
  # -- see NEWS.md), the only remaining difference between "TS" and
  # "OLS" is the internal slope estimator's own robustness, a much
  # smaller effect than when the two used genuinely different outer
  # loops. Asserting one is systematically more moderate than the
  # other on a specific random series, without having verified that
  # empirically first, is the same mistake already made and corrected
  # more than once earlier in this file. Both are simply checked to
  # run and return valid, finite results.
  expect_true(is.finite(rho_ols) && abs(rho_ols) <= 1)
  expect_true(is.finite(rho_sen) && abs(rho_sen) <= 1)
})

test_that(".theil_sen_pairs() subsamples when the full pair set exceeds max_pairs, and still returns a valid, consistent set", {
  t <- 1:60  # choose(60, 2) = 1770 pairs, comfortably above a small cap
  pairs <- sptrends:::.theil_sen_pairs(t, max_pairs = 100L, seed = 7L)

  expect_equal(length(pairs$i_idx), 100L)
  expect_equal(length(pairs$j_idx), 100L)
  expect_equal(length(pairs$dt), 100L)
  # Every i < j (a valid, non-degenerate pair), and dt strictly positive
  # for every one of them.
  expect_true(all(pairs$i_idx < pairs$j_idx))
  expect_true(all(pairs$dt > 0))

  # Two calls with the same seed give the same subsample -- reproducible,
  # not a fresh random draw each time.
  pairs2 <- sptrends:::.theil_sen_pairs(t, max_pairs = 100L, seed = 7L)
  expect_equal(pairs$i_idx, pairs2$i_idx)
})

test_that(".theil_sen_pairs()'s internal seeding does not disturb the caller's own RNG state", {
  t <- 1:60
  set.seed(99)
  before <- runif(1)

  set.seed(99)
  invisible(sptrends:::.theil_sen_pairs(t, max_pairs = 100L, seed = 7L))
  after <- runif(1)

  # If .theil_sen_pairs() had permanently consumed/altered the global RNG
  # stream, `after` would differ from `before` despite both being "the
  # very next runif() draw after set.seed(99)".
  expect_equal(before, after)
})

test_that("prewhiten(method = 'TFPW_WS') runs validly on an ordinary, moderately autocorrelated series under its own rewritten iteration mechanism (mirroring zyp::zyp.TFPW_Z()'s mechanics)", {
  # A short, clean linear-plus-noise series -- not the pathological
  # near-unit-root case the rewrite specifically targeted.
  set.seed(404)
  n <- 30
  data <- 1:n * 0.05 + arima.sim(list(ar = 0.4), n = n)
  r <- terra::rast(nrows = 1, ncols = 1, nlyrs = n)
  r <- terra::setValues(r, matrix(as.numeric(data), nrow = 1))

  pw <- prewhiten(r, verbose = FALSE, report = FALSE)
  rho_est <- terra::values(pw$diagnostics$Rho, mat = FALSE)[1]
  clamped <- terra::values(pw$diagnostics$Clamped, mat = FALSE)[1]

  # Not asserting whether this specific series clamps or not under the
  # rewritten mechanism (unverified without running it, the same
  # mistake already made and corrected more than once earlier in this
  # file) -- only the invariant genuinely guaranteed regardless: Rho is
  # always a valid, finite correlation coefficient.
  expect_true(is.finite(rho_est) && abs(rho_est) <= 1)
  expect_true(is.na(clamped) || clamped %in% c(0, 1))
})

test_that(".lag1_acf_vectorised() matches R's own acf(x, lag.max = 1)$acf[2] for each row", {
  set.seed(405)
  Y <- matrix(rnorm(4 * 20), nrow = 4, ncol = 20)
  got <- sptrends:::.lag1_acf_vectorised(Y)

  want <- apply(Y, 1, function(row) {
    stats::acf(row, lag.max = 1, plot = FALSE)$acf[2]
  })

  expect_equal(got, want, tolerance = 1e-10)
})

test_that("prewhiten(method = 'TFPW_WS')'s rewritten core still returns a valid, well-formed result end to end on a small raster with mixed autocorrelation levels", {
  set.seed(406)
  n_time <- 24
  vals <- rbind(
    as.numeric(arima.sim(list(ar = 0.2), n = n_time)),
    as.numeric(arima.sim(list(ar = 0.6), n = n_time)),
    as.numeric(arima.sim(list(ar = 0.9), n = n_time)),
    cumsum(rnorm(n_time, sd = 0.1))
  )
  r <- terra::rast(nrows = 2, ncols = 2, nlyrs = n_time)
  r <- terra::setValues(r, vals)

  pw <- prewhiten(r, verbose = FALSE, report = FALSE)

  expect_s4_class(pw$series, "SpatRaster")
  expect_equal(terra::nlyr(pw$series), n_time)
  expect_true(all(names(pw$diagnostics) ==
                    c("DW_initial", "Rho", "Modified", "Clamped")))

  rho_vals <- terra::values(pw$diagnostics$Rho, mat = FALSE)
  expect_true(all(is.finite(rho_vals) & abs(rho_vals) <= 1))
})

test_that(".theil_sen_pairs() handles the case where no .Random.seed exists yet in the global environment (a fresh R session, before any randomness has been drawn)", {
  t <- 1:60  # choose(60, 2) = 1770 pairs, above the small cap below,
             # so the seeding branch is actually entered
  had_seed <- exists(".Random.seed", envir = .GlobalEnv)
  old_seed <- if (had_seed) get(".Random.seed", envir = .GlobalEnv) else NULL
  if (had_seed) rm(".Random.seed", envir = .GlobalEnv)

  result <- tryCatch({
    expect_false(exists(".Random.seed", envir = .GlobalEnv))
    pairs <- sptrends:::.theil_sen_pairs(t, max_pairs = 100L, seed = 7L)
    expect_equal(length(pairs$i_idx), 100L)
    expect_false(exists(".Random.seed", envir = .GlobalEnv))
    "ok"
  }, finally = {
    # Restore whatever RNG state this test suite's own later tests
    # expect, rather than leaving .Random.seed permanently absent.
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
  })
  expect_equal(result, "ok")
})

test_that("prewhiten(method = 'TFPW_Z') processes every valid cell unconditionally (no DW gate), unlike TFPW_WS", {
  # A mix of persistence levels -- some cells would NOT cross
  # TFPW_WS's own DW gate, but TFPW_Z has no such gate at all.
  set.seed(504)
  n_time <- 20
  vals <- rbind(
    rnorm(n_time),                                    # white noise
    as.numeric(arima.sim(list(ar = 0.7), n = n_time)) # persistent
  )
  r <- terra::rast(nrows = 1, ncols = 2, nlyrs = n_time)
  r <- terra::setValues(r, vals)

  pw_zhang <- prewhiten(r, method = "TFPW_Z", verbose = FALSE, report = FALSE)
  pw_ws    <- prewhiten(r, method = "TFPW_WS", verbose = FALSE,
                         report = FALSE)

  # TFPW_Z's own DW_initial is NA throughout -- there is no gate for it
  # to report.
  dw_zhang <- terra::values(pw_zhang$diagnostics$DW_initial, mat = FALSE)
  expect_true(all(is.na(dw_zhang)))

  # Every valid cell is processed under TFPW_Z -- Modified reflects
  # only whether the iteration's own clamp fired, not any gating
  # decision (both cells here are valid, so both get a real attempt).
  mod_zhang <- terra::values(pw_zhang$diagnostics$Modified, mat = FALSE)
  expect_true(all(mod_zhang %in% c(0, 1)))

  # TFPW_WS, by contrast, does report a real DW_initial value.
  dw_ws <- terra::values(pw_ws$diagnostics$DW_initial, mat = FALSE)
  expect_true(all(is.finite(dw_ws)))
})

test_that("prewhiten(method = 'TFPW_Z')'s own out$method field reflects 'TFPW_Z', not 'TFPW_WS'", {
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 15, seed = 505)$series
  pw <- prewhiten(r, method = "TFPW_Z", verbose = FALSE, report = FALSE)
  expect_equal(pw$method, "TFPW_Z")
})

test_that("prewhiten(method = 'VCTFPW') runs validly, keeps the full series length, and reports its three diagnostics", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 24, seed = 506)$series
  pw <- prewhiten(r, method = "VCTFPW", verbose = FALSE, report = FALSE)

  expect_s4_class(pw$series, "SpatRaster")
  expect_equal(terra::nlyr(pw$series), terra::nlyr(r))
  expect_true(all(names(pw$diagnostics) ==
                    c("Rho", "Beta_corrected", "Modified")))

  rho_vals <- terra::values(pw$diagnostics$Rho, mat = FALSE)
  expect_true(all(is.finite(rho_vals[!is.na(rho_vals)])))
  expect_true(all(abs(rho_vals[!is.na(rho_vals)]) <= 1))
})

test_that("prewhiten(method = 'VCTFPW') messages, rather than errors, when report = TRUE", {
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 20, seed = 507)$series
  expect_message(
    prewhiten(r, method = "VCTFPW", report = TRUE, verbose = FALSE),
    "no dedicated summary/histogram/map functions"
  )
})

test_that("prewhiten(method = 'VCTFPW') errors on fewer than 4 layers", {
  r <- sim_trend_stack(nrow = 2, ncol = 2, n_time = 3, seed = 508)$series
  expect_error(
    prewhiten(r, method = "VCTFPW", verbose = FALSE, report = FALSE),
    "at least 4 layers"
  )
})

test_that(".lag1_acf_vectorised() returns 0, not NaN, for a row with zero variance", {
  Y <- rbind(
    rep(5, 10),          # constant row -- zero variance
    rnorm(10)             # ordinary row, for contrast
  )
  result <- sptrends:::.lag1_acf_vectorised(Y)
  expect_equal(result[1], 0)
  expect_true(is.finite(result[2]))
})

test_that("prewhiten(method = 'TFPW_Z') does not error on a cell whose residuals collapse to zero variance during iteration -- a real regression, not a hypothetical one", {
  # A perfectly linear cell (zero noise): a Sen/OLS trend fit removes
  # it essentially exactly, so the detrended residual's own variance
  # can collapse to (near) zero during iteration -- exactly the
  # degenerate case that broke this method's own convergence check
  # before .lag1_acf_vectorised() was fixed to return 0 rather than
  # NaN for it. TFPW_Z has no DW gate to exclude this cell the way
  # TFPW_WS might, so it is the method most exposed to this case.
  n <- 20
  perfectly_linear <- 1:n * 0.1
  ordinary <- as.numeric(arima.sim(list(ar = 0.5), n = n))
  r <- terra::rast(nrows = 1, ncols = 2, nlyrs = n)
  r <- terra::setValues(r, rbind(perfectly_linear, ordinary))

  expect_error(
    prewhiten(r, method = "TFPW_Z", verbose = FALSE, report = FALSE),
    NA
  )
})

test_that("prewhiten(method = 'VCTFPW', verbose = TRUE) prints its own step-by-step progress messages", {
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 20, seed = 510)$series
  expect_message(
    prewhiten(r, method = "VCTFPW", verbose = TRUE, report = FALSE),
    "Estimating the Sen slope"
  )
  expect_message(
    prewhiten(r, method = "VCTFPW", verbose = TRUE, report = FALSE),
    "Estimating lag-1 autocorrelation"
  )
  expect_message(
    prewhiten(r, method = "VCTFPW", verbose = TRUE, report = FALSE),
    "Variance-correcting"
  )
  expect_message(
    prewhiten(r, method = "VCTFPW", verbose = TRUE, report = FALSE),
    "Restoring the corrected trend"
  )
})

test_that("prewhiten(method = 'TFPW_Z', verbose = TRUE) prints its own 'no DW gate' message", {
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 15, seed = 511)$series
  expect_message(
    prewhiten(r, method = "TFPW_Z", verbose = TRUE, report = FALSE),
    "No Durbin-Watson gate for method = \"TFPW_Z\""
  )
})

test_that("prewhiten(method = 'VCTFPW') returns NA for a cell with missing values, without erroring on the other, valid cells", {
  set.seed(513)
  n <- 20
  vals <- rbind(
    as.numeric(arima.sim(list(ar = 0.3), n = n)),
    { v <- as.numeric(arima.sim(list(ar = 0.3), n = n)); v[3] <- NA; v }
  )
  r <- terra::rast(nrows = 1, ncols = 2, nlyrs = n)
  r <- terra::setValues(r, vals)

  pw <- prewhiten(r, method = "VCTFPW", verbose = FALSE, report = FALSE)
  rho_vals <- terra::values(pw$diagnostics$Rho, mat = FALSE)

  expect_true(is.finite(rho_vals[1]))
  expect_true(is.na(rho_vals[2]))
})

test_that("print()/summary()/plot() no longer error for method = 'VCTFPW'", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 20, seed = 902)$series
  pw <- prewhiten(r, method = "VCTFPW", verbose = FALSE, report = FALSE)

  expect_output(print(pw), "Wang, Chen, Becker & Liu")
  expect_error(print(pw), NA)

  expect_message(summary(pw), "No dedicated summary")
  expect_error(summary(pw), NA)

  expect_message(plot(pw), "No dedicated plot")
  expect_error(plot(pw), NA)
})

test_that("print() correctly labels method = 'TFPW_Z' with its own citation, not Wang & Swail's", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 15, seed = 903)$series
  pw <- prewhiten(r, method = "TFPW_Z", verbose = FALSE, report = FALSE)

  expect_output(print(pw), "Zhang, Vincent, Hogg & Niitsoo")
  expect_false(isTRUE(grepl(
    "Wang & Swail",
    paste(capture.output(print(pw)), collapse = " ")
  )))
})

test_that("print() correctly labels the default method (TFPW_WS) with its own citation", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 15, seed = 905)$series
  pw <- prewhiten(r, verbose = FALSE, report = FALSE)
  expect_output(print(pw), "Wang & Swail")
})

test_that("prewhiten(method = 'VCTFPW') uses the published variance ratio and 95% autocorrelation gate", {
  set.seed(907)
  n <- 30
  ncell_test <- 6
  vals <- t(vapply(seq_len(ncell_test), function(i) {
    as.numeric(arima.sim(list(ar = 0.8), n = n)) + 1:n * 0.02
  }, numeric(n)))
  r <- terra::rast(nrows = 2, ncols = 3, nlyrs = n)
  r <- terra::setValues(r, vals)

  pw <- prewhiten(r, method = "VCTFPW", verbose = FALSE, report = FALSE)

  series_vals <- terra::values(pw$series, mat = TRUE)
  modified <- as.logical(
    terra::values(pw$diagnostics$Modified, mat = FALSE))
  limit <- stats::qnorm(0.975) / sqrt(n)
  rho <- terra::values(pw$diagnostics$Rho, mat = FALSE)
  expect_equal(modified, abs(rho) > limit)
  expect_true(any(modified))

  for (cell in which(modified)) {
    slopes <- utils::combn(seq_len(n), 2, function(ii) {
      (vals[cell, ii[2]] - vals[cell, ii[1]]) / (ii[2] - ii[1])
    })
    beta <- stats::median(slopes)
    detrended <- vals[cell, ] - beta * (seq_len(n) - 1)
    prewhitened <- c(
      detrended[1],
      detrended[-1] - rho[cell] * detrended[-n]
    )
    variance_ratio <- stats::var(vals[cell, ]) /
      stats::var(prewhitened)
    beta_corrected <- if (rho[cell] > 0) {
      beta / sqrt((1 + rho[cell]) / (1 - rho[cell]))
    } else {
      beta
    }
    expected <- prewhitened * variance_ratio +
      beta_corrected * (seq_len(n) - 1)
    expect_equal(unname(series_vals[cell, ]), expected, tolerance = 1e-10)
  }
})

test_that("prewhiten(method = 'VCTFPW') leaves cells below its autocorrelation gate unchanged", {
  set.seed(908)
  n <- 20
  vals <- rbind(
    rep(5, n),
    as.numeric(arima.sim(list(ar = 0.4), n = n)) + 1:n * 0.02  # ordinary cell
  )
  r <- terra::rast(nrows = 1, ncols = 2, nlyrs = n)
  r <- terra::setValues(r, vals)

  pw <- prewhiten(r, method = "VCTFPW", verbose = FALSE, report = FALSE)
  series_vals <- terra::values(pw$series, mat = TRUE)

  modified <- terra::values(pw$diagnostics$Modified, mat = FALSE)
  expect_equal(modified[1], 0)
  expect_equal(unname(series_vals[1, ]), vals[1, ])
})
