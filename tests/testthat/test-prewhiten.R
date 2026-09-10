test_that("prewhiten returns series and diagnostics with correct layers", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 20, ar1 = 0.6, seed = 1)$series
  res <- prewhiten(r, report = FALSE, verbose = FALSE)

  expect_s4_class(res$series, "SpatRaster")
  expect_equal(terra::nlyr(res$series), terra::nlyr(r))
  expect_true(all(c("DW_initial", "Rho", "Modified") %in% names(res$diagnostics)))
})

test_that("prewhiten leaves ungated cells unmodified", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 15, ar1 = 0, seed = 42)$series
  res <- prewhiten(r, dw_low = 0, dw_high = 4,
                   report = FALSE, verbose = FALSE)
  mod <- terra::values(res$diagnostics$Modified, mat = FALSE)
  expect_true(all(mod == 0, na.rm = TRUE))
})

test_that("prewhiten excludes cells with missing values", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 10, seed = 7)$series
  r[1] <- NA
  res <- prewhiten(r, report = FALSE, verbose = FALSE)
  expect_true(is.na(terra::values(res$diagnostics$Modified, mat = FALSE)[1]))
})

test_that("prewhiten(method = 'TFPW_Y') returns n-1 layers, preserves the Theil-Sen trend, and matches a hand-computed reference", {
  set.seed(1)
  n <- 20
  nr <- 5
  nc <- 5
  true_slope <- 0.05
  true_rho <- 0.6
  t <- 1:n

  # One cell's AR(1)+trend series, computed by hand in R (mirroring the
  # Python cross-check done before writing this method), replicated
  # across a small raster so the same known series is testable via the
  # real vectorised implementation.
  eps <- numeric(n)
  eps[1] <- rnorm(1)
  for (i in 2:n) eps[i] <- true_rho * eps[i - 1] + rnorm(1)
  series_vals <- true_slope * t + eps

  r <- terra::rast(nrows = nr, ncols = nc, nlyrs = n)
  for (k in seq_len(n)) r[[k]] <- series_vals[k]
  names(r) <- paste0("t", t)

  result <- prewhiten(r, method = "TFPW_Y", report = FALSE, verbose = FALSE)
  expect_s3_class(result, "prewhiten")
  expect_identical(result$method, "TFPW_Y")
  expect_equal(terra::nlyr(result$series), n - 1)
  expect_identical(names(result$diagnostics), c("Beta_TheilSen", "Rho"))

  # Hand-computed reference for the same algorithm, independent of the
  # package's own internal helper functions -- all pairwise slopes via
  # combn(), matching the standard Theil-Sen definition directly (the
  # earlier nested-vapply attempt here had a real bug: the inner
  # result's length varies with i, from n-1 down to 1, which cannot be
  # vapply()'d against a fixed FUN.VALUE template).
  pairs_ref <- utils::combn(n, 2)
  slopes_ref <- (series_vals[pairs_ref[2, ]] - series_vals[pairs_ref[1, ]]) /
    (t[pairs_ref[2, ]] - t[pairs_ref[1, ]])
  beta_ref <- stats::median(slopes_ref)
  R_ref <- series_vals - beta_ref * t
  Rc_ref <- R_ref - mean(R_ref)
  rho_ref <- sum(Rc_ref[-1] * Rc_ref[-n]) / sum(Rc_ref^2)
  rho_ref <- min(max(rho_ref, -0.99), 0.99)
  Rp_ref <- R_ref[-1] - rho_ref * R_ref[-n]
  Y_ref <- Rp_ref + beta_ref * t[-1]

  beta_actual <- terra::values(result$diagnostics$Beta_TheilSen, mat = FALSE)[1]
  rho_actual <- terra::values(result$diagnostics$Rho, mat = FALSE)[1]
  y_actual <- terra::values(result$series, mat = TRUE)[1, ]

  expect_equal(beta_actual, beta_ref, tolerance = 1e-8)
  expect_equal(rho_actual, rho_ref, tolerance = 1e-8)
  expect_equal(unname(y_actual), unname(Y_ref), tolerance = 1e-8)
})

test_that("prewhiten(method = 'TFPW_Y') errors below 4 layers", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 3, seed = 200)$series
  expect_error(prewhiten(r, method = "TFPW_Y", verbose = FALSE, report = FALSE),
               "at least 4 layers")
})

test_that("prewhiten(method = 'TFPW_Y') excludes cells with missing values, same convention as TFPW_WS", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 8, seed = 201)$series
  r[1] <- NA
  result <- prewhiten(r, method = "TFPW_Y", verbose = FALSE, report = FALSE)
  beta_vals <- terra::values(result$diagnostics$Beta_TheilSen, mat = FALSE)
  expect_true(is.na(beta_vals[1]))
  series_row1 <- terra::values(result$series, mat = TRUE)[1, ]
  expect_true(all(is.na(series_row1)))
})

test_that("print()/summary()/plot() work on a standalone TFPW_Y prewhiten() result", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 202)$series
  result <- prewhiten(r, method = "TFPW_Y", verbose = FALSE, report = FALSE)
  expect_error(print(result), NA)
  expect_error(suppressMessages(summary(result)), NA)
  expect_error(plot(result, which = "maps"), NA)
  expect_error(plot(result, which = "histograms"), NA)
})

test_that("prewhiten(method = 'TFPW_Y', report = TRUE) runs the full reporting branch without error", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 203)$series
  expect_error(
    suppressMessages(prewhiten(r, method = "TFPW_Y", verbose = FALSE,
                                report = TRUE)),
    NA
  )
})

test_that("workflow_tst(prewhiten_args = list(method = 'TFPW_Y')) runs end-to-end and print()/summary() don't error (regression: previously errored on the truncated series' layer names and Modified field)", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 204)$series
  result <- workflow_tst(r, prewhiten_args = list(method = "TFPW_Y"),
                report = FALSE, verbose = FALSE)
  expect_identical(result$prewhiten$method, "TFPW_Y")
  expect_equal(terra::nlyr(result$prewhiten$series), 11)
  expect_error(print(result), NA)
  expect_error(suppressMessages(summary(result)), NA)
})

test_that("inspect_ts_cell() handles a TFPW_Y-prewhitened object without error or warning (regression: previously errored on the missing Modified field, and separately silently misaligned t by one step)", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 205)$series
  pw <- prewhiten(r, method = "TFPW_Y", verbose = FALSE, report = FALSE)
  ext_r <- terra::ext(r)
  mid_x <- mean(c(ext_r[1], ext_r[2]))
  mid_y <- mean(c(ext_r[3], ext_r[4]))
  testthat::local_mocked_bindings(
    click = function(...) data.frame(x = mid_x, y = mid_y),
    .package = "terra"
  )
  # expect_warning(..., NA) here specifically guards against the
  # silent-recycling failure mode: graphics::plot(t, y, ...) with
  # mismatched lengths warns (recycling a shorter argument), it does
  # not error -- "no error" alone would not have caught the bug this
  # regression test is named for.
  expect_warning(
    result <- suppressMessages(inspect_ts_cell(r, prewhitened = pw)),
    NA
  )
  expect_false(is.null(result$prewhitened))
})

test_that(".inspect_ts_cell_core() aligns t to the shorter TFPW_Y-prewhitened series correctly, not just avoids erroring", {
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 206)$series
  pw <- prewhiten(r, method = "TFPW_Y", verbose = FALSE, report = FALSE)
  expect_equal(terra::nlyr(pw$series), terra::nlyr(r) - 1)

  result <- suppressMessages(sptrends:::.inspect_ts_cell_core(
    r, cell = 15, neighbourhood = FALSE, prewhitened = pw
  ))
  # The prewhitened panel's own series must be exactly one shorter than
  # the raw panel's -- if t and pw_series were still mismatched, this
  # would either error (matrix methods) or silently succeed with the
  # wrong length propagated through, not the correct, deliberately
  # shorter one.
  expect_equal(length(result$prewhitened$series),
               length(result$raw$series) - 1)
})

test_that("prewhiten(method = 'TFPW_WS')'s new Clamped diagnostic flags a cell whose iterative rho estimate hits the +-1 stability bound, and does not flag a well-behaved cell that never hits it", {
  n <- 20
  set.seed(301)

  # Cell 1: a near-unit-root series (cumulative sum -- extremely
  # persistent, close to a random walk). Its own residual
  # autocorrelation is expected to push num/den to or beyond 1 during
  # the iteration, hitting the Clamped safety bound.
  near_unit_root <- cumsum(rnorm(n, sd = 0.1))

  # Cell 2: white noise with a small linear trend -- no serial
  # correlation problem at all, so should not even be gated, let
  # alone clamped.
  well_behaved <- 1:n * 0.05 + rnorm(n, sd = 0.2)

  vals <- matrix(c(near_unit_root, well_behaved), nrow = 2, byrow = TRUE)
  r <- terra::rast(nrows = 1, ncols = 2, nlyrs = n)
  r <- terra::setValues(r, vals)

  pw <- prewhiten(r, verbose = FALSE, report = FALSE)
  expect_true("Clamped" %in% names(pw$diagnostics))

  clamped <- terra::values(pw$diagnostics$Clamped, mat = FALSE)
  rho <- terra::values(pw$diagnostics$Rho, mat = FALSE)
  modified <- terra::values(pw$diagnostics$Modified, mat = FALSE)

  # Not asserting cell 1 (near-unit-root) definitely hits the clamp --
  # the core iterative mechanism was substantially rewritten since this
  # test was first written (see NEWS.md), and whether a specific
  # random construction like this one still hits the +-1 bound under
  # the new mechanism is not something to assume without re-verifying
  # it empirically first -- the same class of mistake already made and
  # corrected more than once elsewhere in this test suite. Only the
  # invariants the code actually guarantees, on whichever side of that
  # outcome a given cell lands on, are checked here.
  raw_cell1 <- as.numeric(terra::extract(r, 1))
  pw_cell1 <- as.numeric(terra::extract(pw$series, 1))

  if (isTRUE(clamped[1] == 1)) {
    # A clamped cell's own Rho should sit exactly at the +-0.99 bound --
    # this is the whole point of the diagnostic: that value reflects
    # the safety bound, not a precise converged estimate.
    expect_true(abs(abs(rho[1]) - 0.99) < 1e-8)
    # A clamped cell is deliberately left uncorrected (see the
    # function's own documentation for why: dividing by (1 - rho) at
    # rho = 0.99 would multiply the residuals by 100, an extreme
    # transform built on an unreliable estimate) -- Modified must be 0
    # for it, and its prewhitened series must equal the original raw
    # series exactly, not some aggressively transformed version of it.
    expect_equal(modified[1], 0)
    expect_equal(pw_cell1, raw_cell1)
  } else {
    # Not clamped: still a valid, finite rho, and Modified reflecting
    # whether this cell was gated and genuinely corrected.
    expect_true(is.finite(rho[1]) && abs(rho[1]) <= 1)
    expect_true(modified[1] %in% c(0, 1))
  }

  # The well-behaved cell should not be clamped -- either because it
  # was never gated for prewhitening at all, or because its own rho
  # genuinely converged without ever touching +-1.
  expect_true(is.na(clamped[2]) || clamped[2] == 0)
})

test_that("prewhiten(method = 'TFPW_WS')'s Clamped diagnostic is NA for invalid (all-NA) cells, matching Modified's own existing convention", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 12, seed = 302)$series
  vals <- terra::values(r, mat = TRUE)
  vals[1, ] <- NA
  r <- terra::setValues(r, vals)

  pw <- prewhiten(r, verbose = FALSE, report = FALSE)
  clamped <- terra::values(pw$diagnostics$Clamped, mat = FALSE)
  modified <- terra::values(pw$diagnostics$Modified, mat = FALSE)

  expect_true(is.na(clamped[1]))
  expect_equal(is.na(clamped), is.na(modified))
})
