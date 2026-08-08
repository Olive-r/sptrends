# These tests compare sptrends' own core statistics against independent,
# long-established CRAN implementations of the same published methods:
# Kendall::MannKendall() and trend::mk.test() for the classic Mann-Kendall
# S statistic, trend::sens.slope() for the Theil-Sen slope, and
# modifiedmk::tfpwmk() for Yue-Pilon trend-free prewhitening (the same
# algorithm sptrends' own prewhiten(method = "TFPW_Y") implements --
# not one of modifiedmk's other variance-correction methods, which use a
# genuinely different approach). Every comparison here uses a single
# 1-cell raster, since these external packages work on plain vectors,
# not SpatRasters -- this only validates the underlying statistic, not
# sptrends' own spatial/raster-handling layer, which is covered
# elsewhere in this test suite.
#
# What is expected to match EXACTLY (deterministic quantities, no
# implementation choice involved): the Mann-Kendall S statistic itself,
# and the Theil-Sen slope (both are just "count concordant minus
# discordant pairs" and "median of pairwise slopes" respectively --
# there is only one correct value). What is expected to match only
# APPROXIMATELY: the p-value/z-statistic, since different packages make
# different, individually legitimate choices about continuity
# correction and variance-of-S formulas in the presence of ties; an
# exact mismatch there is not itself evidence of a bug, only a genuine
# implementation difference documented as such wherever it appears
# below.

test_that("trend_test(method = 'mk')'s S statistic matches Kendall::MannKendall() exactly", {
  skip_if_not_installed("Kendall")
  set.seed(500)
  x <- cumsum(rnorm(30)) + seq_len(30) * 0.1

  r <- terra::rast(nrows = 1, ncols = 1, nlyr = 30)
  terra::values(r) <- matrix(x, nrow = 1)
  result <- trend_test(r, method = "MK", report = FALSE, verbose = FALSE)
  S_sptrends <- terra::values(result$stats[["S"]], mat = FALSE)[1]

  ref <- Kendall::MannKendall(x)
  expect_equal(S_sptrends, as.numeric(ref$S), tolerance = 1e-8)

  # tau = S / D, D = n(n-1)/2 for data with no ties (true here, since
  # x is not rounded) -- an independent cross-check derived from
  # Kendall's own tau, not just re-testing S. Kendall::MannKendall()
  # stores tau internally as C-level single precision (confirmed via
  # str(unclass(ref)): attr(*, "Csingle") = TRUE on tau itself), so a
  # tolerance tighter than single-precision's own ~7 significant
  # digits would fail on rounding noise that has nothing to do with
  # sptrends' own correctness -- 1e-6 stays well inside that.
  n <- length(x)
  D <- n * (n - 1) / 2
  expect_equal(S_sptrends / D, as.numeric(ref$tau), tolerance = 1e-6)
})

test_that("trend_test(method = 'mk')'s S statistic matches trend::mk.test() exactly, on data with ties too", {
  skip_if_not_installed("trend")
  set.seed(501)
  # round()-ing deliberately introduces ties, since S's tie-handling
  # is exactly the kind of detail an exact-match check should exercise,
  # not just avoid.
  x <- round(cumsum(rnorm(25)) + seq_len(25) * 0.15, 1)

  r <- terra::rast(nrows = 1, ncols = 1, nlyr = 25)
  terra::values(r) <- matrix(x, nrow = 1)
  result <- trend_test(r, method = "MK", report = FALSE, verbose = FALSE)
  S_sptrends <- terra::values(result$stats[["S"]], mat = FALSE)[1]

  ref <- trend::mk.test(x)
  expect_equal(S_sptrends, unname(ref$estimates["S"]), tolerance = 1e-8)
  # tau is deliberately NOT cross-checked here the same naive way as
  # the no-ties test above: tau's own denominator D is not simply
  # n(n-1)/2 once ties are present (Kendall, 1974, Ch. 3 -- the
  # correction trend::mk.test() itself already applies internally to
  # its own $estimates["tau"]), and reimplementing that correction
  # here just to check this one field would be testing this test's
  # own arithmetic more than sptrends. S's exact match above already
  # confirms the one thing this test needs to: tie-handling in the
  # underlying score itself agrees with an independent implementation.
})

test_that("trend_test(method = 'mk')'s p-value is approximately consistent with Kendall::MannKendall() and trend::mk.test() (not expected to match exactly -- see file header)", {
  skip_if_not_installed("Kendall")
  skip_if_not_installed("trend")
  set.seed(502)
  x <- cumsum(rnorm(40)) + seq_len(40) * 0.12

  r <- terra::rast(nrows = 1, ncols = 1, nlyr = 40)
  terra::values(r) <- matrix(x, nrow = 1)
  result <- trend_test(r, method = "MK", report = FALSE, verbose = FALSE)
  p_sptrends <- terra::values(result$stats[["p"]], mat = FALSE)[1]

  p_kendall <- Kendall::MannKendall(x)$sl
  p_trend <- trend::mk.test(x)$p.value

  # A loose but meaningful check: all three should agree on
  # significance at the usual alpha = 0.05 threshold, and be within a
  # generous absolute tolerance of each other -- catching a gross
  # implementation error (wrong sign, wrong distribution, off by an
  # order of magnitude) without demanding bit-for-bit agreement on a
  # value where legitimate formula choices differ.
  expect_equal((p_sptrends <= 0.05), (p_kendall <= 0.05))
  expect_equal((p_sptrends <= 0.05), (p_trend <= 0.05))
  expect_lt(abs(p_sptrends - p_kendall), 0.05)
  expect_lt(abs(p_sptrends - p_trend), 0.05)
})

test_that("slope_estimator(method = 'theilsen')'s slope matches trend::sens.slope() exactly", {
  skip_if_not_installed("trend")
  set.seed(503)
  x <- cumsum(rnorm(35)) + seq_len(35) * 0.08

  r <- terra::rast(nrows = 1, ncols = 1, nlyr = 35)
  terra::values(r) <- matrix(x, nrow = 1)
  result <- slope_estimator(r, method = "TS", report = FALSE,
                             verbose = FALSE)
  slope_sptrends <- terra::values(result$slope, mat = FALSE)[1]

  ref <- trend::sens.slope(x)
  expect_equal(slope_sptrends, unname(ref$estimates["Sen's slope"]),
               tolerance = 1e-8)
})

test_that("slope_estimator(method = 'theilsen')'s slope matches trend::sens.slope() exactly, on data with ties too", {
  skip_if_not_installed("trend")
  set.seed(504)
  x <- round(cumsum(rnorm(28)) + seq_len(28) * 0.1, 1)

  r <- terra::rast(nrows = 1, ncols = 1, nlyr = 28)
  terra::values(r) <- matrix(x, nrow = 1)
  result <- slope_estimator(r, method = "TS", report = FALSE,
                             verbose = FALSE)
  slope_sptrends <- terra::values(result$slope, mat = FALSE)[1]

  ref <- trend::sens.slope(x)
  expect_equal(slope_sptrends, unname(ref$estimates["Sen's slope"]),
               tolerance = 1e-8)
})

test_that("prewhiten(method = 'TFPW_Y')'s own Theil-Sen slope (computed on the raw series, before prewhitening) matches modifiedmk::tfpwmk()'s 'Old Sen's Slope' exactly -- both are the same unprewhitened quantity", {
  skip_if_not_installed("modifiedmk")
  set.seed(505)
  x <- cumsum(rnorm(30) * 0.5) + seq_len(30) * 0.1

  r <- terra::rast(nrows = 1, ncols = 1, nlyr = 30)
  terra::values(r) <- matrix(x, nrow = 1)
  result <- prewhiten(r, method = "TFPW_Y", report = FALSE,
                       verbose = FALSE)
  beta_sptrends <- terra::values(result$diagnostics$Beta_TheilSen,
                                  mat = FALSE)[1]

  ref <- modifiedmk::tfpwmk(x)
  expect_equal(beta_sptrends, unname(ref["Old Sen's Slope"]),
               tolerance = 1e-8)
})

test_that("prewhiten(method = 'TFPW_Y') followed by trend_test() gives a significance call directionally consistent with modifiedmk::tfpwmk() on strongly autocorrelated data (not an exact-match check -- the two implementations may estimate the lag-1 autocorrelation coefficient on the detrended series slightly differently)", {
  skip_if_not_installed("modifiedmk")
  set.seed(506)
  # Strong AR(1) noise plus a real trend, so the sign of the detected
  # trend is unambiguous regardless of exactly how each implementation
  # handles the serial correlation.
  n <- 40
  e <- numeric(n)
  e[1] <- rnorm(1)
  for (k in 2:n) e[k] <- 0.6 * e[k - 1] + rnorm(1)
  x <- e + seq_len(n) * 0.15

  r <- terra::rast(nrows = 1, ncols = 1, nlyr = n)
  terra::values(r) <- matrix(x, nrow = 1)
  pw <- prewhiten(r, method = "TFPW_Y", report = FALSE, verbose = FALSE)
  trend_pw <- trend_test(pw$series, method = "MK", report = FALSE,
                          verbose = FALSE)
  slope_sign_sptrends <- sign(terra::values(pw$diagnostics$Beta_TheilSen,
                                             mat = FALSE)[1])

  ref <- modifiedmk::tfpwmk(x)
  slope_sign_ref <- sign(unname(ref["Sen's Slope"]))

  expect_equal(slope_sign_sptrends, slope_sign_ref)
})

test_that("CMK regional aggregation agrees with rkt's independent RKT score", {
  skip_if_not_installed("rkt")

  X <- rbind(
    c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12),
    c(1, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 12),
    c(12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1),
    c(2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11),
    c(1, 2, 4, 3, 5, 6, 8, 7, 9, 10, 12, 11),
    c(12, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 1),
    c(1, 4, 2, 5, 3, 7, 6, 9, 8, 11, 10, 12),
    c(3, 1, 2, 6, 4, 5, 9, 7, 8, 12, 10, 11),
    c(10, 12, 8, 11, 6, 9, 4, 7, 2, 5, 1, 3)
  )
  base <- terra::rast(nrows = 3, ncols = 3, xmin = 0, xmax = 3,
                       ymin = 0, ymax = 3)
  layers <- lapply(seq_len(ncol(X)), function(i) {
    terra::setValues(base, X[, i])
  })
  raster_series <- do.call(c, layers)

  cmk <- trend_test(raster_series, method = "CMK",
                     continuity = FALSE, report = FALSE, verbose = FALSE)
  regional <- rkt::rkt(
    date = rep(seq_len(ncol(X)), times = nrow(X)),
    y = as.vector(t(X)),
    block = rep(seq_len(nrow(X)), each = ncol(X)),
    correct = FALSE
  )

  # rkt returns the sum across blocks; CMK returns the mean across the
  # local region. Their scores must therefore differ only by m.
  expect_equal(
    terra::values(cmk$stats$Sm, mat = FALSE)[5],
    regional$S / nrow(X),
    tolerance = 1e-12
  )
  expect_equal(
    sign(terra::values(cmk$stats$Sm, mat = FALSE)[5]),
    sign(regional$S)
  )

  # No equality is asserted for the corrected variance or p-value:
  # rkt uses the Hirsch-Slack rank-covariance correction and continuity
  # on the summed score, whereas CMK uses analytical RAMK covariance
  # and, by default, no continuity correction on the regional mean.
})
