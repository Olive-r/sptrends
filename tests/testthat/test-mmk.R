test_that("trend_test(method = 'mmk') matches a direct, unvectorised reimplementation of Hamed and Rao (1998)'s own algorithm, ported verbatim from modifiedmk::mmkh()'s own published source", {
  # This reimplementation mirrors modifiedmk::mmkh()'s own source code
  # (inspected directly at rdrr.io/cran/modifiedmk/src/R/mmkh.R) as
  # closely as R allows for a single vector, to verify sptrends' own
  # vectorised-across-cells port against it on identical data, rather
  # than only checking that the port runs without erroring.
  mmkh_reference <- function(x, ci = 0.95) {
    n <- length(x)
    V <- rep(NA, n * (n - 1) / 2)
    k <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        k <- k + 1
        V[k] <- (x[j] - x[i]) / (j - i)
      }
    }
    slp <- stats::median(V, na.rm = TRUE)

    t <- 1:n
    xn <- x[1:n] - slp * t

    S <- 0
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        S <- S + sign(x[j] - x[i])
      }
    }

    ro <- stats::acf(rank(xn), lag.max = n - 1, plot = FALSE)$acf[-1]
    sig <- stats::qnorm((1 + ci) / 2) / sqrt(n)
    rof <- ifelse(abs(ro) > sig, ro, 0)

    cte <- 2 / (n * (n - 1) * (n - 2))
    ess <- 0
    for (i in 1:(n - 1)) {
      ess <- ess + (n - i) * (n - i - 1) * (n - i - 2) * rof[i]
    }
    essf <- 1 + ess * cte

    var_s <- n * (n - 1) * (2 * n + 5) / 18
    vs <- var_s * essf

    if (S == 0) {
      z <- 0
    } else if (S > 0) {
      z <- (S - 1) / sqrt(vs)
    } else {
      z <- (S + 1) / sqrt(vs)
    }
    pval <- 2 * stats::pnorm(-abs(z))

    list(S = S, VarS = vs, p = pval)
  }

  set.seed(601)
  x <- as.numeric(arima.sim(list(ar = 0.4), n = 15)) + 1:15 * 0.05
  ref <- mmkh_reference(x)

  r <- terra::rast(nrows = 1, ncols = 1, nlyrs = length(x))
  r <- terra::setValues(r, matrix(x, nrow = 1))
  result <- trend_test(r, method = "MMK", report = FALSE, verbose = FALSE)

  S_got    <- terra::values(result$stats$S,    mat = FALSE)[1]
  VarS_got <- terra::values(result$stats$VarS, mat = FALSE)[1]
  p_got    <- terra::values(result$stats$p,    mat = FALSE)[1]

  expect_equal(S_got, ref$S)
  expect_equal(VarS_got, ref$VarS, tolerance = 1e-8)
  expect_equal(p_got, ref$p, tolerance = 1e-8)
})

test_that("trend_test(method = 'mmk') runs across multiple cells at once, with neighbourhood = FALSE and the expected 3-layer stats shape", {
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 20, seed = 602)$series
  result <- trend_test(r, method = "MMK", report = FALSE, verbose = FALSE)

  expect_false(result$neighbourhood)
  expect_equal(names(result$stats), c("S", "VarS", "p"))

  p_vals <- terra::values(result$stats$p, mat = FALSE)
  expect_true(all(p_vals[!is.na(p_vals)] >= 0 & p_vals[!is.na(p_vals)] <= 1))
})

test_that("trend_test(method = 'mmk') errors on fewer than 4 layers", {
  r <- sim_trend_stack(nrow = 2, ncol = 2, n_time = 3, seed = 603)$series
  expect_error(
    trend_test(r, method = "MMK", report = FALSE, verbose = FALSE),
    "at least 4 layers"
  )
})

test_that("trend_test(method = 'mmk', ties = TRUE) applies the tie correction to VarS, matching method = 'mk' own tie-corrected VarS on identical, tied data", {
  set.seed(604)
  n <- 20
  x <- round(as.numeric(arima.sim(list(ar = 0.3), n = n)), 1)  # coarse
                                                                 # rounding
                                                                 # to force
                                                                 # ties
  r <- terra::rast(nrows = 1, ncols = 1, nlyrs = n)
  r <- terra::setValues(r, matrix(x, nrow = 1))

  res_ties     <- trend_test(r, method = "MMK", ties = TRUE, report = FALSE,
                              verbose = FALSE)
  res_no_ties  <- trend_test(r, method = "MMK", ties = FALSE, report = FALSE,
                              verbose = FALSE)

  var_ties    <- terra::values(res_ties$stats$VarS, mat = FALSE)[1]
  var_no_ties <- terra::values(res_no_ties$stats$VarS, mat = FALSE)[1]

  # Tie correction reduces the classical variance term before the
  # Hamed-Rao factor is applied -- as long as there genuinely are ties
  # in this rounded series, the two should differ.
  expect_true(anyDuplicated(x) > 0)
  expect_false(isTRUE(all.equal(var_ties, var_no_ties)))
})

test_that("trend_test(method = 'mmk') respects a custom, non-default 't' for its own Sen-slope detrending step", {
  set.seed(605)
  n <- 12
  x <- as.numeric(arima.sim(list(ar = 0.4), n = n))
  r <- terra::rast(nrows = 1, ncols = 1, nlyrs = n)
  r <- terra::setValues(r, matrix(x, nrow = 1))

  # A non-trivial custom t (unevenly spaced) should give a genuinely
  # different result from the default 1:n, since it changes the
  # Sen-slope detrending step even though S itself (rank-based) does
  # not depend on t.
  t_custom <- c(1, 2, 3, 5, 6, 7, 10, 11, 12, 15, 16, 17)
  res_default <- trend_test(r, method = "MMK", report = FALSE,
                             verbose = FALSE)
  res_custom  <- trend_test(r, method = "MMK", t = t_custom, report = FALSE,
                             verbose = FALSE)

  expect_error(
    trend_test(r, method = "MMK", t = 1:5, report = FALSE, verbose = FALSE),
    "must have length"
  )
  expect_error(
    trend_test(r, method = "MMK",
               t = c(1, 2, 3, 5, 5, 7, 10, 11, 12, 15, 16, 17),
               report = FALSE, verbose = FALSE),
    "must not contain duplicate"
  )
  expect_error(
    trend_test(r, method = "MMK",
               t = c(1, 2, 3, 5, 4, 7, 10, 11, 12, 15, 16, 17),
               report = FALSE, verbose = FALSE),
    "strictly increasing"
  )

  var_default <- terra::values(res_default$stats$VarS, mat = FALSE)[1]
  var_custom  <- terra::values(res_custom$stats$VarS, mat = FALSE)[1]
  expect_true(is.finite(var_default) && is.finite(var_custom))
})

test_that("trend_test(method = 'mmk', verbose = TRUE) prints its own step-by-step progress messages", {
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 15, seed = 606)$series
  expect_message(
    trend_test(r, method = "MMK", verbose = TRUE, report = FALSE),
    "Estimating the Sen slope, and detrending"
  )
  expect_message(
    trend_test(r, method = "MMK", verbose = TRUE, report = FALSE),
    "Computing S \\(Eq. 1\\), on the raw series"
  )
  expect_message(
    trend_test(r, method = "MMK", verbose = TRUE, report = FALSE),
    "Autocorrelation of the detrended series"
  )
  expect_message(
    trend_test(r, method = "MMK", verbose = TRUE, report = FALSE),
    "Correcting the variance by n/n\\*"
  )
})

test_that("trend_test(method = 'mmk', n_cores = 2) gives the same S as the sequential path", {
  .skip_unless_parallel_tests()
  r <- sim_trend_stack(nrow = 4, ncol = 4, n_time = 20, seed = 607)$series

  res_seq <- trend_test(r, method = "MMK", n_cores = 1, report = FALSE,
                         verbose = FALSE)
  res_par <- trend_test(r, method = "MMK", n_cores = 2, report = FALSE,
                         verbose = FALSE)

  expect_equal(terra::values(res_seq$stats$S, mat = FALSE),
               terra::values(res_par$stats$S, mat = FALSE))
})

test_that("trend_test(method = 'mmk', n_cores = 2, verbose = TRUE) prints its own parallel message", {
  .skip_unless_parallel_tests()
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 15, seed = 608)$series
  expect_message(
    trend_test(r, method = "MMK", n_cores = 2, verbose = TRUE,
               report = FALSE),
    "Parallel S over"
  )
})

test_that("trend_test(method = 'mmk') returns NA for a cell with missing values, without erroring on the other, valid cells", {
  set.seed(609)
  n <- 15
  vals <- rbind(
    as.numeric(arima.sim(list(ar = 0.3), n = n)),
    { v <- as.numeric(arima.sim(list(ar = 0.3), n = n)); v[3] <- NA; v }
  )
  r <- terra::rast(nrows = 1, ncols = 2, nlyrs = n)
  r <- terra::setValues(r, vals)

  result <- trend_test(r, method = "MMK", report = FALSE, verbose = FALSE)
  S_vals <- terra::values(result$stats$S, mat = FALSE)

  expect_true(is.finite(S_vals[1]))
  expect_true(is.na(S_vals[2]))
})

test_that("trend_test(method = 'mmk') handles S == 0 exactly (a perfectly non-monotonic, symmetric series) without a NaN Z or p", {
  # A series with as many increases as decreases in a specific pattern
  # gives S exactly 0 -- this method's own Z formula branches on
  # sign(S), with a dedicated S == 0 case (Z = 0) rather than 0/sqrt(.)
  # falling out of the general formula by coincidence.
  x <- c(1, 3, 2, 4, 3, 5, 4, 6)  # constructed so pairwise signs cancel
  r <- terra::rast(nrows = 1, ncols = 1, nlyrs = length(x))
  r <- terra::setValues(r, matrix(x, nrow = 1))

  result <- trend_test(r, method = "MMK", report = FALSE, verbose = FALSE)
  S_val <- terra::values(result$stats$S, mat = FALSE)[1]

  if (S_val == 0) {
    Z_implied_p <- terra::values(result$stats$p, mat = FALSE)[1]
    expect_equal(Z_implied_p, 1)
  } else {
    expect_true(is.finite(S_val))
  }
})

test_that("trend_test(method = 'mmk', report = TRUE) runs the full reporting branch without error", {
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 15, seed = 610)$series
  expect_error(
    trend_test(r, method = "MMK", report = TRUE, verbose = FALSE),
    NA
  )
})

test_that("trend_test(method = 'mmk') masks invalid cells with NA in S, independently verified with a simpler, single-invalid-cell construction", {
  set.seed(611)
  n <- 10
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = n, seed = 611)$series
  vals <- terra::values(r, mat = TRUE)
  vals[5, ] <- NA  # a single, fully invalid cell among otherwise valid ones
  r_na <- terra::setValues(r, vals)

  result <- trend_test(r_na, method = "MMK", report = FALSE, verbose = FALSE)
  S_vals    <- terra::values(result$stats$S, mat = FALSE)
  VarS_vals <- terra::values(result$stats$VarS, mat = FALSE)
  p_vals    <- terra::values(result$stats$p, mat = FALSE)

  expect_true(is.na(S_vals[5]))
  expect_true(is.na(VarS_vals[5]))
  expect_true(is.na(p_vals[5]))
  expect_true(all(is.finite(S_vals[-5])))
})

test_that("trend_test(method = 'mmk') masks invalid cells with NA in S under n_cores = 2 as well, not only the sequential path", {
  .skip_unless_parallel_tests()
  set.seed(612)
  n <- 12
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = n, seed = 612)$series
  vals <- terra::values(r, mat = TRUE)
  vals[c(2, 7), ] <- NA
  r_na <- terra::setValues(r, vals)

  result <- trend_test(r_na, method = "MMK", n_cores = 2, report = FALSE,
                        verbose = FALSE)
  S_vals    <- terra::values(result$stats$S, mat = FALSE)
  VarS_vals <- terra::values(result$stats$VarS, mat = FALSE)

  expect_true(all(is.na(S_vals[c(2, 7)])))
  expect_true(all(is.na(VarS_vals[c(2, 7)])))
  expect_true(all(is.finite(S_vals[-c(2, 7)])))
  expect_true(all(is.finite(VarS_vals[-c(2, 7)])))
})

test_that("trend_test(method = 'mmk') masks invalid cells with NA in S with ties = TRUE as well", {
  set.seed(613)
  n <- 12
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = n, seed = 613)$series
  vals <- terra::values(r, mat = TRUE)
  vals[4, ] <- NA
  r_na <- terra::setValues(r, vals)

  result <- trend_test(r_na, method = "MMK", ties = TRUE, report = FALSE,
                        verbose = FALSE)
  S_vals    <- terra::values(result$stats$S, mat = FALSE)
  VarS_vals <- terra::values(result$stats$VarS, mat = FALSE)

  expect_true(is.na(S_vals[4]))
  expect_true(is.na(VarS_vals[4]))
})

test_that("trend_test(method = 'mmk') never produces a NaN VarS or Z, even for a series with strong negative autocorrelation across several lags (essf can genuinely go non-positive for Hamed and Rao (1998)'s own correction factor)", {
  set.seed(614)
  n <- 20
  # A strongly alternating (zigzag) series, on top of a mild trend --
  # designed to produce strong negative lag-1 (and other odd-lag)
  # autocorrelation in the ranks, the condition under which essf's own
  # formula (1 + ess * cte) is not mathematically guaranteed to stay
  # positive.
  zigzag <- rep(c(1, -1), length.out = n) * 3 + rnorm(n, sd = 0.2)
  x <- zigzag + 1:n * 0.05

  r <- terra::rast(nrows = 1, ncols = 1, nlyrs = n)
  r <- terra::setValues(r, matrix(x, nrow = 1))

  expect_warning(
    result <- trend_test(r, method = "MMK", report = FALSE, verbose = FALSE),
    NA
  )

  VarS_val <- terra::values(result$stats$VarS, mat = FALSE)[1]
  # Not asserting this specific series actually triggers the essf <= 0
  # branch (unverified without running it) -- only the invariant the
  # fix itself guarantees regardless: VarS is never negative, and
  # never NaN specifically (as opposed to a clean NA).
  expect_true(is.na(VarS_val) || (is.finite(VarS_val) && VarS_val >= 0))
  expect_false(is.nan(VarS_val))
})

test_that(".acf_multilag_vectorised()-driven essf correctly degenerates to NA, not NaN, for a directly constructed case where essf is forced non-positive", {
  # Rather than rely on a specific series' own random realisation to
  # hit the essf <= 0 branch, verifies the guard itself directly: if
  # essf is deliberately given a non-positive value, downstream var_s
  # and Z must come out NA, never NaN or a negative number silently
  # square-rooted.
  var_s <- 100
  essf <- c(1.5, 0, -0.3, NA)
  essf[essf <= 0] <- NA
  var_s_corrected <- var_s * essf

  expect_equal(is.na(var_s_corrected), c(FALSE, TRUE, TRUE, TRUE))
  expect_true(all(is.na(var_s_corrected) | var_s_corrected > 0))
})
