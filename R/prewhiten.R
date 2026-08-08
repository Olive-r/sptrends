# Selective AR(1) prewhitening (Wang & Swail 2001 / Prais-Winsten), gated by
# a Durbin-Watson diagnostic. See ?prewhiten for the full method
# description and references.

#' @noRd
.ols_slope_vectorised <- function(Y, t) {
  dt <- t - mean(t)
  sxx <- sum(dt^2)
  ybar <- rowMeans(Y)
  as.numeric((Y - ybar) %*% dt) / sxx
}

# Vectorised lag-1 autocorrelation across every row of Y at once, matching
# R's own acf(x, lag.max = 1)$acf[2] definition (mean-centred, denominator
# the full sum of squares, not the lag-1 cross-product's own count) --
# this is the "raw autocorrelation of the untouched series" measurement
# zyp::zyp.TFPW_Z()'s own initial estimate uses, before any detrending, and
# again at each subsequent iteration on the raw series detrended by the
# latest slope estimate (never on a running transformed series).
#' @noRd
.lag1_acf_vectorised <- function(Y) {
  n <- ncol(Y)
  Yc <- Y - rowMeans(Y)
  num <- rowSums(Yc[, 2:n, drop = FALSE] * Yc[, 1:(n - 1), drop = FALSE])
  den <- rowSums(Yc^2)
  # A row with zero variance (a genuinely constant series, or -- more
  # commonly here -- one whose residuals collapse to (near) zero at some
  # iteration because a fitted trend happens to match it almost exactly)
  # gives den == 0, and num / den would be NaN rather than a number the
  # calling code's own +-1 stability clamp can act on (NaN >= 1 is NA in
  # R, not TRUE, so that clamp silently fails to catch it, and the NaN
  # then poisons every later computation for that cell, including the
  # iteration's own convergence check). A constant series has no serial
  # correlation to speak of, so 0 is the sensible value here, not an
  # error condition to propagate.
  rho <- num / den
  rho[den == 0] <- 0
  rho
}

# Vectorised Theil-Sen slope across every row of Y at once -- the median of
# all pairwise slopes per row, computed as a single matrix operation for the
# pairwise differences (one call, not a per-cell loop), with only the final
# per-row median needing a row-wise apply. Pair indices are the caller's own
# responsibility to build once (see .theil_sen_pairs() below) and reuse
# across repeated calls in an iterative context such as TFPW_WS's own
# refit loop, rather than rebuilding them on every iteration.
#' @noRd
.theil_sen_slope_vectorised <- function(Y, i_idx, j_idx, dt) {
  diffs <- (Y[, j_idx, drop = FALSE] - Y[, i_idx, drop = FALSE]) /
    rep(dt, each = nrow(Y))
  apply(diffs, 1, stats::median)
}

# Builds the pair indices (i_idx, j_idx) and their time differences (dt)
# once, for reuse across repeated .theil_sen_slope_vectorised() calls in an
# iterative loop. Mirrors slope_estimator()'s own max_pairs safeguard
# (random subsample of pairs for long series) but with a fixed, internal
# default rather than a user-facing argument, since this is an
# implementation detail of a single refit step inside prewhiten(), not a
# reported result in its own right.
#' @noRd
.theil_sen_pairs <- function(t, max_pairs = 2000L, seed = 1L) {
  n <- length(t)
  all_pairs <- utils::combn(n, 2)
  n_pairs_total <- ncol(all_pairs)
  if (n_pairs_total > max_pairs) {
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    }
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv,
                        inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
    keep <- sample.int(n_pairs_total, max_pairs)
    all_pairs <- all_pairs[, keep, drop = FALSE]
  }
  list(i_idx = all_pairs[1, ], j_idx = all_pairs[2, ],
       dt = t[all_pairs[2, ]] - t[all_pairs[1, ]])
}

# Durbin-Watson critical value bounds (dL, dU) for k = 1 regressor, alpha =
# 0.05, reproduced (transcription only, not methodology) from Savin & White
# (1977), via the University of York historical statistics archive. Covers
# N = 15 to 100; outside that range the nearest tabulated value is used
# with a warning.
.dw_table_N  <- c(15:40, seq(45, 100, by = 5))
.dw_table_dL <- c(1.077, 1.106, 1.133, 1.158, 1.180, 1.201, 1.221, 1.240, 1.257,
                  1.273, 1.288, 1.302, 1.316, 1.328, 1.341, 1.352, 1.363, 1.373,
                  1.383, 1.393, 1.402, 1.411, 1.419, 1.427, 1.435, 1.442,
                  1.475, 1.503, 1.527, 1.549, 1.567, 1.583, 1.598, 1.611, 1.624,
                  1.635, 1.645, 1.654)
.dw_table_dU <- c(1.361, 1.371, 1.381, 1.392, 1.401, 1.411, 1.420, 1.429, 1.437,
                  1.446, 1.454, 1.461, 1.468, 1.476, 1.483, 1.489, 1.496, 1.502,
                  1.508, 1.514, 1.519, 1.524, 1.530, 1.535, 1.540, 1.544,
                  1.566, 1.585, 1.601, 1.616, 1.629, 1.641, 1.652, 1.662, 1.671,
                  1.679, 1.687, 1.694)

#' @noRd
.dw_critical_values <- function(n, verbose = TRUE) {
  if (verbose && n < min(.dw_table_N)) {
    message(sprintf(
      paste0("Note: n=%d is below the Durbin-Watson table range (N>=15) -- ",
             "using N=15 as an approximation; interpret with caution for ",
             "such short series."),
      n))
  }
  if (verbose && n > max(.dw_table_N)) {
    message(sprintf(
      paste0("Note: n=%d exceeds the Durbin-Watson table range (N<=100) -- ",
             "using N=100 as a conservative approximation."),
      n))
  }
  list(
    dL = stats::approx(.dw_table_N, .dw_table_dL, xout = n, rule = 2)$y,
    dU = stats::approx(.dw_table_N, .dw_table_dU, xout = n, rule = 2)$y
  )
}

#' @noRd
.durbin_watson_ols <- function(X, t) {
  n <- ncol(X)
  ok <- stats::complete.cases(X)

  slope <- .ols_slope_vectorised(X, t)
  intercept <- rowMeans(X) - slope * mean(t)
  pred <- outer(slope, t) + intercept
  resid <- X - pred

  sum_diff2 <- rowSums(
    (resid[, 2:n, drop = FALSE] - resid[, 1:(n - 1), drop = FALSE])^2
  )
  sum_res2 <- rowSums(resid^2)
  dw <- sum_diff2 / sum_res2

  rho1 <- rowSums(resid[, 2:n, drop = FALSE] *
                     resid[, 1:(n - 1), drop = FALSE]) /
    rowSums(resid[, 1:(n - 1), drop = FALSE]^2)

  dw[!ok] <- NA
  rho1[!ok] <- NA
  slope[!ok] <- NA

  list(dw = dw, rho1 = rho1, slope = slope)
}

#' AR(1) prewhitening of raster time series
#'
#' Removes temporal (serial) autocorrelation from each cell's time
#' series while preserving its linear trend, using one of four published
#' methods (see the `method` argument below for the choice). The
#' default, `"TFPW_WS"`: for each cell, fit an OLS trend, compute
#' the Durbin-Watson statistic of the residuals, and -- **only** for
#' cells where DW indicates relevant serial autocorrelation --
#' iteratively estimate the AR(1) coefficient `rho` and apply a
#' trend-preserving transformation (with a Prais-Winsten correction
#' for the first observation). Cells that pass the DW check are left
#' untouched.
#'
#' **Function type:** **Preprocessing function** -- prepares the raw
#' raster time series before trend estimation or significance testing
#' (see [compute_anomalies()] for the other preprocessing step this
#' package offers). Not one of the core trend-analysis pillars itself.
#' This function typically precedes trend estimation ([trend_test()])
#' and slope estimation ([slope_estimator()]) within the standard
#' sptrends workflow.
#'
#' @section Typical use:
#' This function is the preprocessing step between reading data and
#' testing it for a trend:
#' [read_ordered_stack()]/[read_netcdf_stack()] -> [compute_anomalies()]
#' (optional) -> `prewhiten()` -> [trend_test()] -> [slope_estimator()]
#' -> [fdr_correction()]. [workflow_tst()] runs this whole chain in one
#' call.
#'
#' ```
#' raster time series
#'     |
#' prewhiten()
#'     |
#' transformed time series (`result$series`) + diagnostics
#'     |
#' trend_test() and slope_estimator()
#' ```
#' If the input has a seasonal cycle, remove it first with
#' [compute_anomalies()]. Pass `result$series` to later stages; the
#' original input object is not modified and remains available under
#' the name supplied by the caller.
#'
#' @section Methodological details:
#' **Methods and method selection**
#'
#' - **Wang & Swail (2001), `method = "TFPW_WS"`**: trend-free
#'   prewhitening (a Prais-Winsten-style correction that preserves the
#'   linear trend, unlike simple differencing).
#' - **Yue, Pilon & Cavadias (2002), `method = "TFPW_Y"`**: also a
#'   trend-preserving prewhitening method, but a different mechanism:
#'   `rho` is estimated on the detrended residuals directly, every
#'   valid cell is processed unconditionally (no DW-based gate), and
#'   the classic transform loses the first time step. See the `method`
#'   argument below for the practical difference this makes to the
#'   output.
#' - **Zhang et al. (2000), `method = "TFPW_Z"`**: uses the same
#'   iterative mechanism as `TFPW_WS`, but applies it to every valid
#'   cell without a Durbin-Watson gate.
#' - **Wang et al. (2015), `method = "VCTFPW"`**: applies the published
#'   variance and slope corrections only where lag-1 autocorrelation
#'   crosses its two-sided 95% gate.
#' - **Main references**: Wang & Swail (2001) for the `TFPW_WS`
#'   transformation itself; Durbin & Watson (1950, 1951) for the
#'   statistic used as its selective gate; Yue & Wang (2002) for why
#'   gating selectively, rather than prewhitening every cell, matters
#'   for `TFPW_WS` specifically; Yue, Pilon & Cavadias (2002) for
#'   `TFPW_Y`; Zhang et al. (2000) for `TFPW_Z`; and Wang et al. (2015)
#'   for `VCTFPW`. Full citations appear under "References" below.
#' - **Typical applications**: hydrology and climatology time series
#'   with suspected serial autocorrelation, applied before a
#'   Mann-Kendall-family trend test to avoid inflating its false
#'   positive rate.
#'
#' Classical prewhitening is deliberately not provided as a fifth
#' method because filtering the observed series directly can remove or
#' attenuate part of the trend and reduce the power of the subsequent
#' test. The four implemented procedures explicitly preserve or restore
#' the trend, or otherwise correct the transformation.
#'
#' **How it works**
#'
#' Prewhitening every cell indiscriminately (including cells that
#' already behave like white noise) reduces the statistical power of
#' any trend test applied afterwards (Yue, Pilon & Cavadias, 2002);
#' hence the selective gate.
#'
#' The method preserves the linear trend, unlike simple differencing,
#' by solving `Y_t = a + b*t + X_t` with `X_t = rho * X_(t-1) + e_t`
#' for `W_t = (Y_t - rho * Y_(t-1)) / (1 - rho)` (Wang & Swail, 2001,
#' "trend-free prewhitening").
#'
#' This selective, DW-gated approach to prewhitening follows the same
#' overall methodology as the Earth Trends Modeler (ETM) module of
#' TerrSet (Clark Labs) -- an independent implementation, not a port of
#' that software's code.
#'
#' **Implementation notes**
#'
#' **On the core iterative mechanism**: this function's own iteration
#' was substantially rewritten in this version to mirror
#' `zyp::zyp.TFPW_Z()`'s own mechanics, after empirically comparing this
#' function's own `rho` estimate against two independent
#' implementations of essentially the same published method (`zyp`'s
#' `"TFPW_Z"`, and `MannKendallTrends`'s `nanprewhite.AR()`/`prewhite()`)
#' on a real near-unit-root cell: this function's own earlier iteration
#' hit the `+-1` stability clamp (`0.99`), while `zyp` converged to
#' `0.776` and `MannKendallTrends` to `0.860` -- both substantially more
#' moderate, and reasonably close to each other despite differing
#' implementations. Tracing both algorithms step by step on the
#' identical series (see `NEWS.md` for the specific numbers at each
#' iteration) found two mechanical differences, both now adopted:
#' first, each iteration measures lag-1 autocorrelation on the *raw*
#' series detrended by the latest slope estimate, never on a running
#' transformed series -- so what changes between iterations is only
#' which slope estimate is used to detrend the same raw data, not the
#' data being detrended itself; second, the very first estimate, before
#' any detrending, is the raw series' own lag-1 autocorrelation
#' directly, matching `zyp`'s own initial step, rather than the
#' residuals of an initial trend fit (which, on a near-unit-root
#' series, can themselves carry *more* apparent autocorrelation than
#' the raw series did, pushing the very first estimate toward the
#' clamp before any subsequent iteration gets a chance to recover from
#' it).
#'
#' **On `refit_method` and the evidence behind its default**: TerrSet's
#' own Earth Trends Modeler documentation for this method states only
#' that its iteration is "determined exactly as described by Wang and
#' Swail (2001)", without specifying which slope estimator that
#' refitting step itself uses -- this function's own primary source does
#' not settle the question directly. Independently, several papers
#' describing the original Wang and Swail (2001) procedure in more
#' methodological detail (including Zhang and Zwiers (2004), a direct
#' comment on this literature, and Collaud Coen et al. (2020)) describe
#' its own iterative refit as using the Sen slope specifically, and
#' `MannKendallTrends`'s own `prewhite()` source, inspected directly,
#' confirms this: every refit inside its own `TFPW.WS` branch calls its
#' package's own `sen.slope()`, never an OLS fit. Zhang and Zwiers
#' (2004) describe substituting an OLS refit as their own added variant
#' for comparison, not as what Wang and Swail (2001) themselves used.
#' This function's own default (`refit_method = "OLS"`) was kept
#' despite this evidence, both because the primary source available to
#' this package does not confirm it directly, and because `"OLS"` was
#' this function's own behaviour in earlier package versions -- with
#' the outer iteration mechanism now unified between the two values
#' (see above), the choice of refit estimator itself is a materially
#' smaller effect than it was before that rewrite; `"TS"` remains
#' available for a user who prefers it or wants to compare both
#' directly, as recommended above for a `Clamped = 1` cell specifically.
#'
#' **Computational considerations**
#'
#' **On TerrSet's own iteration cap**: TerrSet's documentation states
#' its own maximum of 5 iterations "to avoid the rare cases that fail
#' to converge" -- tested directly against a real near-unit-root cell
#' under this function's own earlier iteration, raising `itmax` from 20
#' to 100 left the estimate unchanged (still at the clamp), showing the
#' earlier instability was a genuine fixed point of that iteration, not
#' merely an iteration budget cut short -- so `itmax` was not lowered to
#' match TerrSet's own cap; this function's default (20) is kept.
#'
#' `TFPW_WS` and `TFPW_Z` are iterative. Within either procedure,
#' `refit_method = "OLS"` is computationally lighter than `"TS"`;
#' Theil-Sen refitting gains robustness at the cost of evaluating many
#' pairwise slopes. The Theil-Sen steps in `TFPW_Y` and `VCTFPW` also
#' become more expensive as the number of time points increases.
#'
#' **Statistical assumptions**
#'
#' All four methods assume an AR(1) noise process (no higher orders)
#' and an approximately linear trend.
#'
#' **Limitations**
#'
#' This is a purely temporal (per-cell) preprocessing step with no
#' spatial component. Cells with any missing value in the time series
#' are excluded entirely from the output.
#'
#' **`method = "TFPW_WS"`**: the `"threshold"` method uses a
#' widely-used but non-universal empirical DW cutoff (`[1.4, 2.6]`);
#' `rho` is assumed constant across the whole series.
#'
#' **`method = "TFPW_Y"`**: loses one observation to lag-1
#' differencing; unlike `TFPW_WS`, has no selective gate, so a cell
#' with genuinely no autocorrelation still gets a (small) `rho`
#' estimated and applied from that specific sample -- see
#' `TFPW_WS` above for the alternative that only touches cells that
#' need it.
#'
#' **`method = "TFPW_Z"`**: processes every valid cell without a
#' selective gate. As with `TFPW_WS`, an iteration that reaches the
#' stability clamp is reported and that cell is left uncorrected.
#'
#' **`method = "VCTFPW"`**: uses the method's published two-sided 95%
#' lag-1-autocorrelation gate and variance correction. Its slope
#' correction follows the published positive-autocorrelation rule.
#'
#' **Quality assurance**
#'
#' Automated tests compare the transformed series and diagnostics with
#' hand-calculated references, verify method-specific gates and
#' first-year retention, exercise constant/invalid series and boundary
#' cases, and require identical sequential and parallel results.
#' Yue-Pilon trend-free prewhitening is additionally compared with
#' `modifiedmk::tfpwmk()` for the quantities both implementations define
#' identically. See `?sptrends` for the package-wide release-check
#' protocol; current check results belong only in `cran-comments.md`.
#'
#' @param x A `terra::SpatRaster`; each layer is one time step, in
#'   increasing chronological order.
#' @param method Which prewhitening method to use, sharing one
#'   interface between four related procedures (see "Methodological
#'   details" below for the full comparison and citations).
#'
#'   `"TFPW_WS"` (default; TFPW-WS, also seen as "Wang-Swail
#'   prewhitening" in the literature):
#'   - **Description**: the selective AR(1) prewhitening described
#'     above, gated by a Durbin-Watson diagnostic.
#'   - **Characteristics**: only cells that need it get modified;
#'     `rho` is estimated on the raw trend residuals.
#'   - **Advantage**: avoids the power loss of prewhitening cells that
#'     already behave like white noise (see "Statistical rationale"
#'     above).
#'
#'   `"TFPW_Y"` (TFPW-Y, also seen simply as "TFPW" -- this is the
#'   best-known formulation of trend-free prewhitening in the
#'   literature):
#'   - **Description**: trend-free pre-whitening (Yue, Pilon, Cavadias,
#'     2002).
#'   - **Characteristics**: estimates a Theil-Sen slope per cell first,
#'     removes it, estimates lag-1 autocorrelation on the *detrended*
#'     residuals (not the raw series, which would over-estimate rho
#'     when a real trend is present), prewhitens those residuals, and
#'     adds the slope back; every valid cell is processed
#'     unconditionally, there is no DW gate.
#'   - **Advantage**: a single, simple rule applied uniformly, with no
#'     threshold or test choice to make.
#'
#'   `"TFPW_Z"` (TFPW-Z here -- not a universally established acronym
#'   in the literature the way the others are, so defined explicitly
#'   wherever used):
#'   - **Description**: Zhang, Vincent, Hogg and Niitsoo (2000), as
#'     later refined into `TFPW_WS` above -- this method is the
#'     earlier, unrefined form the same iterative mechanism (see
#'     "Implementation notes" above) is traced against and mirrors.
#'   - **Characteristics**: identical iteration to `TFPW_WS`, but
#'     with no Durbin-Watson gate -- every valid cell is processed
#'     unconditionally, matching `zyp::zyp.TFPW_Z()`'s own published
#'     behaviour.
#'   - **Advantage**: no gating decision to make, at the same cost
#'     `TFPW_Y` has for the same reason.
#'
#'   `"VCTFPW"` (VCTFPW):
#'   - **Description**: variance-corrected trend-free prewhitening
#'     (Wang, Chen, Becker and Liu, 2015).
#'   - **Characteristics**: removes a Sen slope, estimates lag-1
#'     autocorrelation on the detrended residuals and applies the
#'     published transformation only where that autocorrelation is
#'     significant at the two-sided 95% level. It then applies the
#'     published variance-ratio correction and corrects the restored
#'     slope for positive autocorrelation. Cells that do not cross
#'     the gate retain their original series.
#'   - **Advantage**: addresses a specific, published criticism of
#'     simpler trend-free methods -- their own reduced variance biases
#'     the slope and its significance -- at the cost of a
#'     substantially more involved procedure.
#'
#'   **Comparison**: `TFPW_WS` uses a Durbin-Watson gate, `TFPW_Z`
#'   and `TFPW_Y` process every valid cell unconditionally, and
#'   `VCTFPW` uses its published 95% lag-1-autocorrelation gate before
#'   applying its variance and slope corrections.
#'   See "Methodological details" above for the full citations and
#'   practical implications of each.
#' @param t Numeric vector of time points, one per layer. Defaults to
#'   `1:nlyr(x)`. Irregular spacing is supported, but values must be finite,
#'   unique and strictly increasing.
#' @param dw_low,dw_high Only used when `method = "TFPW_WS"`. Fixed
#'   Durbin-Watson gate thresholds, used only when `dw_method =
#'   "threshold"`.
#' @param dw_method Only used when `method = "TFPW_WS"`.
#'   `"threshold"` (default): use the fixed `dw_low`/`dw_high`
#'   cutoffs. `"test"`: use the formal Durbin-Watson test with critical
#'   values (dL, dU) tabulated for the actual sample size, at the 5% level.
#' @param dw_inconclusive Only used when `method = "TFPW_WS"` and
#'   `dw_method = "test"`. The classic
#'   DW test has an inconclusive zone (`dL <= DW <= dU`) where it cannot
#'   decide -- this is *not* the same as "autocorrelation present".
#'   `"conservative"` (default) prewhitens that zone too (favours not
#'   leaving autocorrelation uncorrected); `"power"` prewhitens only when
#'   the test actually rejects H0 (`DW < dL`), favouring statistical power
#'   downstream. Neither is "the correct" choice in absolute terms.
#' @param eps Only used when `method = "TFPW_WS"` or `"TFPW_Z"`.
#'   Convergence threshold for the iterative estimation of `rho`.
#' @param itmax Only used when `method = "TFPW_WS"` or `"TFPW_Z"`.
#'   Maximum number
#'   of iterations. TerrSet's own Earth Trends Modeler, which implements
#'   this same method, caps this at 5 "to avoid the rare cases that fail
#'   to converge" (its own documentation's wording) -- `sptrends` uses a
#'   higher default (20) instead, relying on `Clamped` (see "Value" below)
#'   to flag a cell whose own iteration did not settle within that budget,
#'   rather than capping the budget itself as tightly as TerrSet does.
#' @param refit_method Only used when `method = "TFPW_WS"` or
#'   `"TFPW_Z"`. Which
#'   slope estimator re-fits the trend at each iteration, before removing
#'   it to isolate the residual autocorrelation for the next `rho`
#'   estimate: `"OLS"` (default, matching this function's own behaviour
#'   in earlier package versions) or `"TS"` (Theil-Sen, matching the
#'   original iterative procedure as multiple independent sources
#'   describe Wang and Swail (2001) actually publishing it -- see
#'   "Implementation notes" below for the specific evidence and why the
#'   default was not simply changed to match it outright). `"TS"`'s own
#'   robustness to outliers can matter specifically for a cell whose
#'   iteration is already unstable (see `Clamped` below): an extreme
#'   value in one iteration's own transformed series distorts an
#'   `"OLS"` refit more than a `"TS"` one, which can itself feed into a
#'   worse `rho` estimate the following iteration. Comparing both directly
#'   on a specific `Clamped = 1` cell, the same way comparing
#'   `TFPW_WS` against `TFPW_Y` is already recommended for that
#'   case (see "Value" below), is worth doing before trusting either
#'   estimate in isolation.
#' @param report Logical. If `TRUE` (default), automatically print a
#'   summary and draw diagnostic histograms and maps after computing
#'   (which functions, exactly, depends on `method` -- see "Value"
#'   below). Set to `FALSE` for programmatic use (e.g. inside a
#'   loop, or when called from [workflow_tst()]) where you don't want
#'   console output or plots as a side effect.
#' @param verbose Logical. Print progress messages and elapsed time.
#'
#' @return Returns a list of class `"prewhiten"`, with:
#'   \item{series}{A `SpatRaster`. For `method = "TFPW_WS"`: same
#'     structure as `x`, with prewhitened cells replaced and the rest
#'     unchanged; layer names get a `"_prewhitened"` suffix. For
#'     `method = "TFPW_Y"`: **one fewer layer than `x`** (the
#'     classic Yue-Pilon transform loses the first time step to lag-1
#'     differencing); every valid cell is prewhitened, with its
#'     Theil-Sen trend preserved; layer names come from `x`'s own
#'     layers 2 through n, with a `"_prewhitened"` suffix. For
#'     `method = "TFPW_Z"`: the same number of layers as `x`; every
#'     valid cell enters the same iterative transformation used by
#'     `TFPW_WS`, but without the Durbin-Watson gate. For
#'     `method = "VCTFPW"`: the same number of layers as `x`;
#'     significantly autocorrelated cells contain the published VCTFPW
#'     transformation and other valid cells remain unchanged.}
#'   \item{diagnostics}{A `SpatRaster`. For `method = "TFPW_WS"`:
#'     4 layers, `DW_initial`, `Rho`, `Modified` (0/1), and `Clamped`
#'     (0/1) -- `Clamped` is `1` for a gated cell whose iterative
#'     `rho` estimate hit the +-1 stability bound (see "Implementation
#'     notes" below) at any point during the iteration, even if a
#'     later iteration's own `rho - rho_old` difference happened to
#'     fall under `eps` immediately afterwards, since a `Rho` reaching
#'     that bound reflects the bound itself, not necessarily a
#'     precise, well-converged estimate. A clamped cell is
#'     deliberately left uncorrected (`series` keeps its original
#'     observed values there, `Modified` is `0` for it) rather than
#'     transformed using that unreliable `rho`: the correction divides
#'     by `(1 - rho)`, so at `rho = 0.99` the residuals would be
#'     multiplied by 100 -- an extreme, disproportionate transform for
#'     what was, underneath, an unstable rather than a trustworthy
#'     estimate. `Rho` still records the clamped value itself, so the
#'     reason a cell went uncorrected remains visible in the
#'     diagnostics rather than indistinguishable from a cell that
#'     never needed correcting in the first place. Reachable via
#'     `summary()`/`plot()` as
#'     [prewhiten_summary()]/[prewhiten_histograms()]/
#'     [prewhiten_maps()]. For `method =
#'     "TFPW_Y"`: 2 layers, `Beta_TheilSen` (the slope removed and
#'     restored) and `Rho` (lag-1 autocorrelation of the detrended
#'     residuals); reachable via internal reporting functions specific
#'     to this method, not the `TFPW_WS` ones above (their
#'     structure differs -- there is no DW gate, `Modified`, or
#'     `Clamped` field for `TFPW_Y`, which does not iterate to
#'     estimate `rho` and so cannot hit this same failure mode --
#'     precisely why comparing both methods on the same data is worth
#'     doing for a cell flagged `Clamped = 1`, rather than trusting
#'     either one in isolation: agreement between them is reassuring,
#'     and disagreement is itself informative about how sensitive that
#'     cell's own result is to the prewhitening method chosen.).
#'     For `method = "TFPW_Z"`: the same 4 layers as `TFPW_WS`, but
#'     `DW_initial` is `NA` because this method has no Durbin-Watson
#'     gate; every valid cell enters the iteration, while `Modified`
#'     and `Clamped` retain the meanings described above.
#'     For `method = "VCTFPW"`: 3 layers, `Rho` (lag-1
#'     autocorrelation of the detrended series), `Beta_corrected`, and
#'     `Modified` (0/1, indicating whether the 95% autocorrelation gate
#'     was crossed).}
#'   \item{method}{The selected method, recorded as one of `"TFPW_WS"`,
#'     `"TFPW_Y"`, `"TFPW_Z"`, or `"VCTFPW"`.}
#'
#' @references
#' Primary method reference (`method = "TFPW_WS"`):
#' - Wang, X.L. and Swail, V.R. (2001) Changes of Extreme Wave Heights in
#'   Northern Hemisphere Oceans and Related Atmospheric Circulation
#'   Regimes. Journal of Climate, 14(10), 2204-2221.
#'
#' On the evidence behind `refit_method`, and the core iterative
#' mechanism this version's rewrite mirrors (see "Implementation notes"
#' above): the original iterative method this function's own
#' `method = "TFPW_WS"` implements, as refined by Wang and Swail
#' (2001) from an earlier procedure, and the R package whose own
#' mechanics this version's rewrite was traced against and now mirrors:
#' - Zhang, X., Vincent, L.A., Hogg, W.D. and Niitsoo, A. (2000)
#'   Temperature and precipitation trends in Canada during the 20th
#'   century. Atmosphere-Ocean, 38(3), 395-429.
#'   \doi{10.1080/07055900.2000.9649654}
#' - Bronaugh, D. and Werner, A. (2013) zyp: Zhang + Yue-Pilon Trends
#'   Package. R package.
#'   \url{https://CRAN.R-project.org/package=zyp}
#' - Zhang, X. and Zwiers, F.W. (2004) Comment on "Applicability of
#'   prewhitening to eliminate the influence of serial correlation on
#'   the Mann-Kendall test" by Sheng Yue and Chun Yuan Wang. Water
#'   Resources Research, 40(3), W03805. \doi{10.1029/2003WR002073}
#' - Collaud Coen, M., Andrews, E., Bigi, A., Martucci, G., Romanens,
#'   G., Vogt, F.P.A. and Vuilleumier, L. (2020) Effects of the
#'   prewhitening method, the time granularity, and the time
#'   segmentation on the Mann-Kendall trend detection and the
#'   associated Sen's slope. Atmospheric Measurement Techniques, 13(12),
#'   6945-6964. \doi{10.5194/amt-13-6945-2020}
#'
#' Primary method reference (`method = "TFPW_Y"`):
#' - Yue, S., Pilon, P., Phinney, B. and Cavadias, G. (2002) The
#'   influence of autocorrelation on the ability to detect trend in
#'   hydrological series. Hydrological Processes, 16(9), 1807-1829.
#'   \doi{10.1002/hyp.1095}
#'
#' (Corrected from an earlier, different Yue, Pilon and Cavadias (2002)
#' paper this citation previously pointed to -- Power of the
#' Mann-Kendall and Spearman's rho tests for detecting monotonic
#' trends in hydrological series, Journal of Hydrology 259, a related
#' but distinct paper by essentially the same authors, published the
#' same year, that does not itself specify the TFPW procedure. Found
#' by tracing which paper `zyp`'s own documentation and the
#' `mannkendall` project's own "Spirit of mannkendall" reference for
#' this same method -- both cite the paper now cited here.)
#'
#' Primary method reference (`method = "VCTFPW"`):
#' - Wang, W., Chen, Y., Becker, S. and Liu, B. (2015) Variance
#'   Correction Prewhitening Method for Trend Detection in
#'   Autocorrelated Data. Journal of Hydrologic Engineering, 20(12),
#'   04015033. \doi{10.1061/(ASCE)HE.1943-5584.0001234}
#'
#' `VCTFPW`'s own implementation here was adapted from the logic of
#' (not copied verbatim from -- see "Implementation notes" above for
#' this package's own vectorised helpers used instead) the R package
#' whose scientific article is already cited above for `refit_method`:
#' - Collaud Coen, M., Andrews, E., Bigi, A., Martucci, G., Romanens,
#'   G., Vogt, F.P.A. and Vuilleumier, L. (2020) Effects of the
#'   prewhitening method, the time granularity, and the time
#'   segmentation on the Mann-Kendall trend detection and the
#'   associated Sen's slope. Atmospheric Measurement Techniques, 13(12),
#'   6945-6964. \doi{10.5194/amt-13-6945-2020}
#'
#' Source of the Durbin-Watson statistic used as the selective gate
#' (`method = "TFPW_WS"` only):
#' - Durbin, J. and Watson, G.S. (1950) Testing for Serial Correlation in
#'   Least Squares Regression, I. Biometrika, 37(3-4), 409-428.
#'   \doi{10.1093/biomet/37.3-4.409}
#' - Durbin, J. and Watson, G.S. (1951) Testing for Serial Correlation in
#'   Least Squares Regression, II. Biometrika, 38(1-2), 159-178.
#'   \doi{10.1093/biomet/38.1-2.159}
#'
#' Basis for gating prewhitening selectively rather than applying it to
#' every cell (see "Methodological details" above):
#' - Yue, S. and Wang, C.Y. (2002) Applicability of prewhitening to
#'   eliminate the influence of serial correlation on the Mann-Kendall
#'   test. Water Resources Research, 38(6), 4-1. \doi{10.1029/2001WR000861}
#'
#' Official software implementation (independent re-implementation of the
#' published method, not a port of this module's code):
#' - Eastman, J.R. (2016) TerrSet Geospatial Monitoring and Modeling
#'   System: Earth Trends Modeler. Clark Labs, Clark University,
#'   Worcester, MA.
#'
#' This function is used (not authored) by the following study:
#' - Gutiérrez-Hernández, O. and García, L.V. (2025) Uncovering true
#'   significant trends in global greening. Remote Sensing Applications:
#'   Society and Environment, 37, 101377. \doi{10.1016/j.rsase.2024.101377}
#'
#' @examples
#' \donttest{
#' # Annual mean NDVI from the bundled environmental dataset.
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#'
#' # Remove serial autocorrelation before testing for a trend -- only
#' # cells that actually show relevant autocorrelation are modified;
#' # the rest pass through unchanged.
#' result <- prewhiten(r, report = FALSE, verbose = FALSE)
#'
#' # result$series is the prewhitened stack, same number of layers as the
#' # input -- feed this into trend_test()/slope_estimator()
#' # next, not the raw r.
#' terra::nlyr(result$series)
#' summary(result)
#' plot(result)   # Rho map and the DW-before/after comparison, real data
#'
#' # Trend-free pre-whitening (Yue-Pilon): every valid cell is
#' # processed, its Theil-Sen trend preserved -- result_yp$series has
#' # one fewer layer than r (the classic transform loses the first
#' # time step to lag-1 differencing).
#' result_yp <- prewhiten(r, method = "TFPW_Y", report = FALSE,
#'                         verbose = FALSE)
#' terra::nlyr(result_yp$series)
#' terra::plot(result_yp$diagnostics$Beta_TheilSen,
#'             main = "Preserved Theil-Sen slope (Yue-Pilon)")
#' }
#'
#' @family Prewhitening functions
#' @export
prewhiten <- function(x, method = c("TFPW_WS", "TFPW_Y", "TFPW_Z",
                                     "VCTFPW"), t,
                       dw_low = 1.4,
                       dw_high = 2.6,
                       dw_method = c("threshold", "test"),
                       dw_inconclusive = c("conservative", "power"),
                       eps = 1e-4, itmax = 20,
                       refit_method = c("OLS", "TS"), report = TRUE,
                       verbose = TRUE) {
  finish_timer <- .sptrends_elapsed_timer("prewhiten()", verbose)
  on.exit(finish_timer(), add = TRUE)
  method <- match.arg(method)
  dw_method <- match.arg(dw_method)
  dw_inconclusive <- match.arg(dw_inconclusive)
  refit_method <- match.arg(refit_method)

  if (!inherits(x, "SpatRaster")) stop("'x' must be a terra SpatRaster.")
  if (!is.numeric(dw_low) || length(dw_low) != 1L || is.na(dw_low) ||
      !is.finite(dw_low) || !is.numeric(dw_high) ||
      length(dw_high) != 1L || is.na(dw_high) || !is.finite(dw_high) ||
      dw_low < 0 || dw_high > 4 || dw_low >= dw_high) {
    stop("'dw_low' and 'dw_high' must be finite scalars satisfying ",
         "0 <= dw_low < dw_high <= 4.")
  }
  eps <- .validate_positive_numeric(eps, "eps")
  itmax <- .validate_positive_integer(itmax, "itmax")
  if (!all(terra::inMemory(x))) x <- x + 0

  X <- terra::values(x, mat = TRUE)
  n <- ncol(X)
  ncell_ <- nrow(X)
  if (missing(t)) t <- seq_len(n)
  t <- .validate_time_axis(t, n)
  ok <- stats::complete.cases(X)
  if (!any(ok)) {
    stop("No cell has a complete time series (every cell has at least ",
         "one NA layer) -- there is nothing to prewhiten.")
  }

  if (method == "VCTFPW") {
    if (n < 4) {
      stop("'x' must have at least 4 layers (time steps) for VCTFPW ",
           "prewhitening.")
    }
    if (verbose) {
      message("[1/4] Estimating the Sen slope on the raw series and ",
              "removing it...")
    }
    ts_pairs_full <- .theil_sen_pairs(t)
    b0_or <- .theil_sen_slope_vectorised(
      X, ts_pairs_full$i_idx, ts_pairs_full$j_idx, ts_pairs_full$dt)
    detrend_or <- X - outer(b0_or, t - t[1])

    if (verbose) {
      message("[2/4] Estimating lag-1 autocorrelation on the ",
              "detrended series and applying the published 95% ",
              "significance gate...")
    }
    c_VCTFPW <- .lag1_acf_vectorised(detrend_or)
    c_VCTFPW[!ok] <- NA_real_
    acf_limit <- stats::qnorm(0.975) / sqrt(n)
    modified <- ok & is.finite(c_VCTFPW) &
      abs(c_VCTFPW) > acf_limit

    # Wang et al. (2015), Eqs. 7-8. Retaining the first detrended
    # observation gives every raster cell the same n-layer geometry;
    # cells that do not cross the autocorrelation gate remain equal to
    # the original series, as the published procedure prescribes.
    resid_or <- detrend_or
    if (any(modified)) {
      resid_or[modified, 2:n] <-
        detrend_or[modified, 2:n, drop = FALSE] -
        c_VCTFPW[modified] *
        detrend_or[modified, 1:(n - 1), drop = FALSE]
    }

    if (verbose) {
      message(sprintf(
        "[3/4] Variance-correcting %d of %d valid cells...",
        sum(modified), sum(ok)))
    }
    var_data <- apply(X, 1, stats::var, na.rm = TRUE)
    var_resid <- apply(resid_or, 1, stats::var, na.rm = TRUE)
    # Published VCTFPW factor (Wang et al., 2015, Eq. 9; Collaud Coen
    # et al., 2020, Eq. 9). Although a square-root ratio would force
    # exact equality of the sample variances algebraically, that would
    # define a different transformation from the published method.
    scale_factor <- var_data / var_resid
    invalid_scale <- modified &
      (!is.finite(scale_factor) | var_data <= 0 | var_resid <= 0)
    resid_rescaled <- resid_or
    valid_modified <- modified & !invalid_scale
    resid_rescaled[valid_modified, ] <-
      resid_or[valid_modified, , drop = FALSE] *
      scale_factor[valid_modified]

    b_vc <- b0_or
    positive_modified <- valid_modified & c_VCTFPW > 0
    b_vc[positive_modified] <- b0_or[positive_modified] /
      sqrt((1 + c_VCTFPW[positive_modified]) /
             (1 - c_VCTFPW[positive_modified]))
    b_vc[invalid_scale] <- NA_real_

    if (verbose) {
      message("[4/4] Restoring the corrected trend...")
    }
    Y <- X
    Y[valid_modified, ] <-
      resid_rescaled[valid_modified, , drop = FALSE] +
      outer(b_vc[valid_modified], t - t[1])
    Y[!ok, ] <- NA_real_
    Y[invalid_scale, ] <- NA_real_

    r <- x[[1]]
    Y_layers <- lapply(seq_len(n), function(k) terra::setValues(r, Y[, k]))
    Y_out <- do.call(c, Y_layers)
    names(Y_out) <- paste0(names(x), "_prewhitened")

    modified_diagnostic <- as.integer(modified)
    modified_diagnostic[!ok] <- NA_integer_
    diagnostics <- c(terra::setValues(r, c_VCTFPW),
                      terra::setValues(r, b_vc),
                      terra::setValues(r, modified_diagnostic))
    names(diagnostics) <- c("Rho", "Beta_corrected", "Modified")

    out <- list(series = Y_out, diagnostics = diagnostics,
                method = "VCTFPW")
    class(out) <- c("prewhiten", "sptrends")
    if (isTRUE(report)) {
      message("report = TRUE has no dedicated summary/histogram/map ",
              "functions yet for method = \"VCTFPW\" (its own ",
              "Rho/Beta_corrected/Modified diagnostics do not match the shape ",
              "prewhiten_summary()/prewhiten_maps() expect). Inspect ",
              "`diagnostics` directly, e.g. summary(out$diagnostics) ",
              "or terra::plot(out$diagnostics).")
    }
    if (verbose) message("Done.")
    return(out)
  }

  if (method == "TFPW_Y") {
    if (n < 4) {
      stop("'x' must have at least 4 layers (time steps) for Yue-Pilon ",
           "prewhitening -- it loses one observation to lag-1 ",
           "differencing, and needs enough points left after that to ",
           "test for a trend.")
    }
    if (verbose) {
      message("[1/4] Estimating the Theil-Sen slope per cell ",
              "(the trend Yue-Pilon preserves)...")
    }
    beta_result <- slope_estimator(x, method = "TS", t = t,
                                    report = FALSE, verbose = FALSE)
    beta <- terra::values(beta_result$slope, mat = FALSE)

    if (verbose) message("[2/4] Detrending...")
    R <- X - outer(beta, t)

    if (verbose) {
      message("[3/4] Estimating lag-1 autocorrelation of the ",
              "detrended residuals...")
    }
    R_mean <- rowMeans(R)
    Rc <- R - R_mean
    num <- rowSums(Rc[, 2:n, drop = FALSE] * Rc[, 1:(n - 1), drop = FALSE])
    den <- rowSums(Rc^2)
    rho <- num / den
    rho[!ok] <- NA_real_
    # Same safety cap as the TFPW_WS branch above: a rho estimate at
    # or beyond +-1 would make the prewhitening transform below blow up
    # or divide by zero.
    rho_capped <- pmin(pmax(rho, -0.99), 0.99)

    if (verbose) {
      message("[4/4] Prewhitening the residuals and restoring the ",
              "trend...")
    }
    Rp <- R[, 2:n, drop = FALSE] -
      rho_capped * R[, 1:(n - 1), drop = FALSE]
    Y <- Rp + outer(beta, t[-1])
    Y[!ok, ] <- NA_real_

    r <- x[[1]]
    Y_layers <- lapply(seq_len(n - 1), function(k) {
      terra::setValues(r, Y[, k])
    })
    Y_out <- do.call(c, Y_layers)
    names(Y_out) <- paste0(names(x)[-1], "_prewhitened")

    diagnostics <- c(terra::setValues(r, beta),
                      terra::setValues(r, rho_capped))
    names(diagnostics) <- c("Beta_TheilSen", "Rho")

    out <- list(series = Y_out, diagnostics = diagnostics, method = "TFPW_Y")
    class(out) <- c("prewhiten", "sptrends")

    if (isTRUE(report)) {
      .TFPW_Y_summary(diagnostics)
      .TFPW_Y_histograms(diagnostics)
      .TFPW_Y_maps(diagnostics)
    }
    if (verbose) message("Done.")
    return(out)
  }

  if (n < 3) {
    stop("'x' must have at least 3 layers (time steps) for ",
         "Durbin-Watson/AR(1) prewhitening to be meaningful.")
  }

  if (method == "TFPW_Z") {
    # Zhang et al. (2000), as implemented in zyp::zyp.TFPW_Z(): every
    # valid cell is processed unconditionally, with no Durbin-Watson
    # pre-filter -- unlike TFPW_WS, this method's own iteration
    # (identical below) is not gated by any diagnostic. dw_init is
    # left NA throughout: there is no DW-based decision here for it to
    # report.
    if (verbose) {
      message("[1/3] No Durbin-Watson gate for method = \"TFPW_Z\" -- ",
              "every valid cell is processed.")
    }
    dw_init <- rep(NA_real_, ncell_)
    gated <- ok
  } else if (dw_method == "threshold") {
    if (verbose) {
      message("[1/3] Computing initial Durbin-Watson (closed-form OLS, ",
              "vectorised)...")
      message(sprintf(
        "      Method: fixed threshold (dw_low=%.2f, dw_high=%.2f).",
        dw_low, dw_high))
    }
    diag0 <- .durbin_watson_ols(X, t)
    dw_init <- diag0$dw
    gated <- ok & (dw_init < dw_low | dw_init > dw_high)
  } else {
    if (verbose) {
      message("[1/3] Computing initial Durbin-Watson (closed-form OLS, ",
              "vectorised)...")
    }
    diag0 <- .durbin_watson_ols(X, t)
    dw_init <- diag0$dw
    cv <- .dw_critical_values(n, verbose = verbose)
    limit <- if (dw_inconclusive == "conservative") cv$dU else cv$dL
    if (verbose) {
      message(sprintf(
        paste0("      Method: formal Durbin-Watson test (n=%d, alpha=0.05): ",
               "dL=%.3f, dU=%.3f"),
        n, cv$dL, cv$dU))
      message(sprintf(
        "      Inconclusive zone (%.3f to %.3f) treated, by design, as: %s.",
        cv$dL, cv$dU, dw_inconclusive))
    }
    gated <- ok & (dw_init < limit | dw_init > (4 - limit))
  }

  if (verbose) {
    message(sprintf(
      "      %d of %d valid cells cross the criterion and will be prewhitened.",
      sum(gated, na.rm = TRUE), sum(ok)))
  }

  W <- X
  rho_final <- rep(0, ncell_)
  modified <- rep(0L, ncell_)
  clamped_final <- rep(0L, ncell_)
  rho_final[!ok] <- NA
  modified[!ok] <- NA
  clamped_final[!ok] <- NA

  if (any(gated, na.rm = TRUE)) {
    idx_g <- which(gated)
    Xg <- X[idx_g, , drop = FALSE]
    # Tracks, per cell, whether the iteration's own num/den estimate ever
    # hit the +-1 stability clamp below at any iteration -- not just the
    # last one. A cell whose true rho is near a unit root can have
    # num/den repeatedly land at or beyond +-1, so the clamped value
    # (0.99/-0.99) stops changing between iterations and the convergence
    # check below (max(abs(rho - rho_old)) < eps) is satisfied immediately
    # -- looking exactly like a clean, well-behaved convergence from the
    # output alone, when it is actually the safety bound repeatedly
    # firing, not a precise, trustworthy estimate of rho. Recorded here
    # (once TRUE, always TRUE for that cell, across all iterations) so a
    # user can tell the two cases apart in the returned diagnostics,
    # rather than treating every rho near +-1 with the same confidence.
    clamped <- rep(FALSE, length(idx_g))

    if (verbose) {
      message("[2/3] Iteratively estimating rho over the gated cells...")
    }

    # Initial estimate: lag-1 autocorrelation of the untouched series
    # itself, exactly as zyp::zyp.TFPW_Z()'s own first step (a raw acf()
    # call, before any detrending) -- not the residuals of an initial
    # trend fit, which on a near-unit-root series can themselves carry
    # more apparent persistence than the raw series did (see NEWS.md
    # for the real near-unit-root cell this was traced against).
    c_val <- .lag1_acf_vectorised(Xg)
    yt <- t[1:(n - 1)]
    ts_pairs <- if (refit_method == "TS") .theil_sen_pairs(yt) else NULL
    refit <- function(Y) {
      if (refit_method == "TS") {
        .theil_sen_slope_vectorised(Y, ts_pairs$i_idx, ts_pairs$j_idx,
                                     ts_pairs$dt)
      } else {
        .ols_slope_vectorised(Y, yt)
      }
    }

    y <- (Xg[, 2:n, drop = FALSE] - c_val * Xg[, 1:(n - 1), drop = FALSE]) /
      (1 - c_val)
    trend <- refit(y)

    pb <- .sptrends_progress(itmax, "Estimating rho (Wang & Swail)", verbose)

    for (iter in seq_len(itmax)) {
      # Detrend the RAW series (never the running transformed y) using
      # the latest trend estimate, then re-measure lag-1 autocorrelation
      # on those residuals -- mirroring zyp::zyp.TFPW_Z()'s own
      # `x <- data - trend * t; c <- acf(x, ...)` step exactly (named
      # `resid` here, not `x`, to avoid shadowing this function's own
      # `x` parameter -- the original input raster, still needed after
      # this loop ends to assemble the final result).
      resid <- Xg - outer(trend, t)
      c_new <- .lag1_acf_vectorised(resid)
      c_old <- c_val
      clamped <- clamped | (c_new >= 1) | (c_new <= -1)
      c_new[c_new >= 1] <- 0.99
      c_new[c_new <= -1] <- -0.99
      c_val <- c_new

      y <- (Xg[, 2:n, drop = FALSE] - c_val * Xg[, 1:(n - 1), drop = FALSE]) /
        (1 - c_val)
      trend <- refit(y)

      .sptrends_progress_step(pb, iter)
      # na.rm = TRUE as a second layer of defence: .lag1_acf_vectorised()
      # itself now guards its own den == 0 case (see its own comment),
      # but should any other, currently unforeseen degenerate
      # configuration produce a stray NaN here, max() with na.rm = FALSE
      # would turn the whole vectorised convergence check itself into
      # NaN -- silently swallowing every OTHER cell's genuine
      # convergence and crashing this if() with exactly the "missing
      # value where TRUE/FALSE needed" error this comment is here to
      # prevent recurring.
      if (max(abs(c_val - c_old), na.rm = TRUE) < eps) break
    }
    .sptrends_progress_close(pb)
    rho <- c_val

    Wg <- matrix(NA_real_, nrow(Xg), n)
    Wg[, 1] <- Xg[, 1] * sqrt(1 - rho^2)
    Wg[, 2:n] <- (Xg[, 2:n, drop = FALSE] -
                    rho * Xg[, 1:(n - 1), drop = FALSE]) / (1 - rho)

    # A clamped cell's rho did not genuinely converge (see the comment
    # above the clamped tracker) -- the correction built from it divides
    # by (1 - rho), which at rho = 0.99 means multiplying the residuals
    # by 100, an extreme, disproportionate transformation for what was,
    # underneath, an unstable rather than a trustworthy estimate. Rather
    # than apply that transformation anyway, a clamped cell's own series
    # is left as originally observed (Xg, not Wg) -- Modified = 0 for it,
    # same as a cell that was never gated in the first place, and Rho
    # still records the clamped estimate itself so the reason a cell
    # went uncorrected remains visible in the diagnostics rather than
    # silently indistinguishable from "never needed correction".
    Wg[clamped, ] <- Xg[clamped, ]

    W[idx_g, ] <- Wg
    rho_final[idx_g] <- rho
    modified[idx_g] <- as.integer(!clamped)
    clamped_final[idx_g] <- as.integer(clamped)
  } else if (verbose) {
    message("[2/3] No cell crosses the DW threshold -- nothing is prewhitened.")
  }

  if (verbose) message("[3/3] Assembling results...")
  r <- x[[1]]
  W_layers <- lapply(seq_len(n), function(k) terra::setValues(r, W[, k]))
  W_out <- do.call(c, W_layers)
  names(W_out) <- paste0(names(x), "_prewhitened")

  diagnostics <- c(terra::setValues(r, dw_init),
                    terra::setValues(r, rho_final),
                    terra::setValues(r, modified),
                    terra::setValues(r, clamped_final))
  names(diagnostics) <- c("DW_initial", "Rho", "Modified", "Clamped")

  out <- list(series = W_out, diagnostics = diagnostics, method = method)
  class(out) <- c("prewhiten", "sptrends")

  if (isTRUE(report)) {
    prewhiten_summary(diagnostics)
    prewhiten_histograms(diagnostics)
    prewhiten_maps(diagnostics)
  }

  out
}

#' Summarise a prewhitening diagnostics raster
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- called internally by `report = TRUE`, and reachable
#' from outside the package via `summary()`.
#'
#' @param diagnostics The `diagnostics` element returned by
#'   [prewhiten()].
#' @param path Character or `NULL`. If supplied, write the summary table to
#'   this CSV path.
#'
#' @return Invisibly, a data frame with summary metrics.
#' @family Prewhitening functions
#' @references
#' - Wang, X.L. and Swail, V.R. (2001) Changes of Extreme Wave Heights in
#'   Northern Hemisphere Oceans and Related Atmospheric Circulation
#'   Regimes. Journal of Climate, 14(10), 2204-2221.
#' @examples
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' pw <- prewhiten(r, report = FALSE, verbose = FALSE)
#'
#' # How many cells needed correcting, and by how much (Durbin-Watson
#' # before/after, estimated AR(1) coefficient rho).
#' sptrends:::prewhiten_summary(pw$diagnostics)
#'
#' @keywords internal
prewhiten_summary <- function(diagnostics, path = NULL) {
  dw <- terra::values(diagnostics$DW_initial, mat = FALSE)
  rho <- terra::values(diagnostics$Rho, mat = FALSE)
  mod <- terra::values(diagnostics$Modified, mat = FALSE)

  n_total <- sum(!is.na(mod))
  n_mod <- sum(mod == 1, na.rm = TRUE)
  mean_rho_mod <- if (n_mod > 0) {
    round(mean(rho[mod == 1], na.rm = TRUE), 4)
  } else {
    NA
  }
  median_dw <- round(stats::median(dw, na.rm = TRUE), 4)

  message(sprintf("Valid cells: %d", n_total))
  message(sprintf("Prewhitened: %d (%.1f%%)", n_mod, 100 * n_mod / n_total))
  if (n_mod > 0) {
    message(sprintf("Mean rho among prewhitened cells: %.4f", mean_rho_mod))
  } else {
    message("Mean rho among prewhitened cells: none prewhitened.")
  }
  message(sprintf("Median Durbin-Watson (all valid cells): %.4f", median_dw))

  tab <- data.frame(
    metric = c("valid_cells", "prewhitened_cells", "pct_prewhitened",
               "mean_rho_prewhitened", "median_dw_global"),
    value = c(n_total, n_mod, round(100 * n_mod / n_total, 2),
              mean_rho_mod, median_dw)
  )

  if (!is.null(path)) {
    utils::write.csv(tab, path, row.names = FALSE)
    message(sprintf("Table written to: %s", path))
  }

  invisible(tab)
}

#' Histograms of prewhitening diagnostics
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- called internally by `report = TRUE`, and reachable
#' from outside the package via `plot(x, which = "histograms")`.
#'
#' @inheritParams prewhiten_summary
#' @param path Character or `NULL`. If supplied, PNGs are written using
#'   this as a filename prefix: setting `path` to `out` writes files named
#'   `out_dw.png` and `out_rho.png`.
#' @return `NULL`, invisibly. Called for its plotting side effect.
#' @family Prewhitening functions
#' @references
#' - Wang, X.L. and Swail, V.R. (2001) Changes of Extreme Wave Heights in
#'   Northern Hemisphere Oceans and Related Atmospheric Circulation
#'   Regimes. Journal of Climate, 14(10), 2204-2221.
#' @examples
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' pw <- prewhiten(r, report = FALSE, verbose = FALSE)
#'
#' # Distributions of the Durbin-Watson statistic and the estimated rho
#' # across all cells, before vs. after prewhitening.
#' sptrends:::prewhiten_histograms(pw$diagnostics)
#'
#' @keywords internal
prewhiten_histograms <- function(diagnostics, path = NULL) {
  dw <- terra::values(diagnostics$DW_initial, mat = FALSE)
  rho <- terra::values(diagnostics$Rho, mat = FALSE)
  mod <- terra::values(diagnostics$Modified, mat = FALSE)

  graphics::hist(dw, breaks = 40, col = "steelblue", border = "white",
       main = "Distribution of initial Durbin-Watson", xlab = "DW")
  graphics::abline(v = c(1.4, 2.6), col = "red", lwd = 2, lty = 2)
  if (!is.null(path)) .save_current_plot(paste0(path, "_dw.png"))

  rho_mod <- rho[mod == 1]
  if (length(rho_mod) > 0) {
    graphics::hist(rho_mod, breaks = 40, col = "orange", border = "white",
         main = "Distribution of rho (prewhitened cells only)", xlab = "rho")
    if (!is.null(path)) .save_current_plot(paste0(path, "_rho.png"))
  } else {
    message("(no prewhitened cells -- rho histogram skipped)")
  }
  invisible(NULL)
}

#' Maps of prewhitening diagnostics
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- called internally by `report = TRUE`, and reachable
#' from outside the package via `plot()`.
#'
#' @inheritParams prewhiten_histograms
#' @return `NULL`, invisibly. Called for its plotting side effect.
#' @family Prewhitening functions
#' @references
#' - Wang, X.L. and Swail, V.R. (2001) Changes of Extreme Wave Heights in
#'   Northern Hemisphere Oceans and Related Atmospheric Circulation
#'   Regimes. Journal of Climate, 14(10), 2204-2221.
#' @examples
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' pw <- prewhiten(r, report = FALSE, verbose = FALSE)
#'
#' # Where, spatially, cells were modified by prewhitening.
#' sptrends:::prewhiten_maps(pw$diagnostics)
#'
#' @keywords internal
prewhiten_maps <- function(diagnostics, path = NULL) {
  pal_seq <- grDevices::hcl.colors(50, "Viridis")
  pal_div <- .sptrends_diverging_palette(50)

  terra::plot(diagnostics$DW_initial, col = pal_seq,
              main = "Initial Durbin-Watson")
  if (!is.null(path)) .save_current_plot(paste0(path, "_dw_map.png"))

  # rho is theoretically bounded in [-1, 1] -- fixing the range there
  # (rather than letting terra::plot() auto-range to whatever rho values
  # happen to occur) keeps zero exactly at the palette's white midpoint,
  # so blue always means genuinely negative autocorrelation, not merely
  # "less positive than the largest rho in this dataset".
  terra::plot(diagnostics$Rho, col = pal_div, range = c(-1, 1),
              main = "Estimated rho (AR1)")
  if (!is.null(path)) .save_current_plot(paste0(path, "_rho_map.png"))

  mod <- diagnostics$Modified == 1
  .safe_categorical_plot(mod, values = c(0, 1),
                          colours = c("grey85", "firebrick"),
                          labels = c("Not prewhitened", "Prewhitened"),
                          main = "Prewhitened cells")
  if (!is.null(path)) .save_current_plot(paste0(path, "_modified_map.png"))

  invisible(NULL)
}

#' @noRd
.TFPW_Y_summary <- function(diagnostics) {
  beta <- terra::values(diagnostics$Beta_TheilSen, mat = FALSE)
  rho <- terra::values(diagnostics$Rho, mat = FALSE)
  n_total <- sum(!is.na(rho))

  message(sprintf("Valid cells: %d", n_total))
  message("Prewhitened: all valid cells (Yue-Pilon has no DW gate --",
          " every cell is detrended, decorrelated, and restored).")
  message(sprintf("Mean rho: %.4f", round(mean(rho, na.rm = TRUE), 4)))
  message(sprintf("Median Theil-Sen slope preserved: %.4g",
                   stats::median(beta, na.rm = TRUE)))
  invisible(NULL)
}

#' @noRd
.TFPW_Y_histograms <- function(diagnostics) {
  beta <- terra::values(diagnostics$Beta_TheilSen, mat = FALSE)
  rho <- terra::values(diagnostics$Rho, mat = FALSE)

  graphics::hist(rho, breaks = 40, col = "orange", border = "white",
                  main = "Distribution of rho (lag-1, on detrended residuals)",
                  xlab = "rho")

  graphics::hist(beta, breaks = 40, col = "steelblue", border = "white",
                  main = "Distribution of the preserved Theil-Sen slope",
                  xlab = "Slope")
  invisible(NULL)
}

#' @noRd
.TFPW_Y_maps <- function(diagnostics) {
  pal_div <- .sptrends_diverging_palette(50)

  # Same reasoning as prewhiten_maps()'s own Rho map: fixing the range
  # to rho's theoretical bounds keeps zero exactly at the palette's
  # white midpoint, not merely "less positive than this dataset's max".
  terra::plot(diagnostics$Rho, col = pal_div, range = c(-1, 1),
              main = "Estimated rho (lag-1, on detrended residuals)")

  range_lim <- .robust_diverging_range(diagnostics$Beta_TheilSen)
  # fill_range = TRUE: without it, cells beyond the robust range would
  # render as blank (NA) instead of saturating to the extreme colour --
  # see slope_map()'s own comment on this for the full reasoning.
  terra::plot(diagnostics$Beta_TheilSen, col = pal_div,
              range = range_lim, fill_range = TRUE,
              main = "Preserved Theil-Sen slope")
  invisible(NULL)
}
