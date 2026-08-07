# Pixel-wise false discovery rate correction (BH, BKY and BY). See
# ?fdr_correction for method notes and references.

#' Benjamini-Hochberg (1995) false discovery rate correction
#'
#' Thin wrapper around [stats::p.adjust()] with `method = "BH"` -- base R
#' already implements this correctly, so it is not reimplemented.
#'
#' **Function type:** **Support function** -- computes the BH procedure
#' used internally by [fdr_correction()]. It is not exported; call
#' `fdr_correction(p, method = "BH")` for a BH-only result.
#'
#' @section Typical use:
#' Supply one family of raw p-values to [fdr_correction()] with
#' `method = "BH"`; this internal helper returns the corresponding
#' adjusted values and rejection decisions.
#'
#' @section Methodological details:
#' **Methods and method selection**
#'
#' - **Original publication**: Benjamini & Hochberg (1995), the paper
#'   that introduced false discovery rate control itself.
#' - **Main references**: Benjamini & Hochberg (1995) for the procedure;
#'   Benjamini & Yekutieli (2001) for why it remains valid under the
#'   positive dependence typical of gridded spatial data. Full citations
#'   appear under "References" below.
#' - **Typical applications**: correcting for multiple testing when many
#'   hypotheses are tested at once (e.g. one Mann-Kendall test per pixel
#'   in a raster) and a fixed, non-adaptive guarantee is preferred over
#'   [fdr_bky()]'s adaptive one -- see [workflow_rta()] for a workflow that
#'   defaults to this method specifically for that reason.
#'
#' **Statistical assumptions and limitations**
#'
#' BH controls FDR under independence and recognised positive-dependence
#' conditions such as PRDS. Spatial autocorrelation diagnostics can be
#' compatible with those conditions but do not prove them. Use [fdr_by()]
#' when control under arbitrary dependence is required.
#'
#' **Quality assurance**
#'
#' Adjusted values are generated directly by `stats::p.adjust()` and are
#' also checked through the package's automated FDR tests.
#'
#' @param p Numeric vector of raw p-values in `[0, 1]` (may contain `NA`).
#' @param q Numeric. Target FDR level. It is not a renamed or adjusted
#'   `alpha`: `alpha` is a per-test Type I error threshold, whereas `q`
#'   bounds the expected proportion of false discoveries among all
#'   discoveries -- a property of the complete rejection set, not of any
#'   one test. A `q` of `0.05` means that, over repeated applications under
#'   the procedure's assumptions, the expected false-discovery proportion
#'   is controlled at 5%; it is not a guarantee that every realised set
#'   contains at most 5% false discoveries.
#'
#' @return A list with `q_value` (BH-adjusted p-values) and `reject`
#'   (logical).
#' @references
#' Primary method reference:
#' - Benjamini, Y., & Hochberg, Y. (1995). Controlling the False
#'   Discovery Rate: A Practical and Powerful Approach to Multiple
#'   Testing. Journal of the Royal Statistical Society: Series B, 57,
#'   289-300. \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Theoretical justification for why BH remains valid under the positive
#' spatial dependence typical of gridded data (for which
#' [spatial_autocorrelation()] or `moran_check` can provide a diagnostic):
#' - Benjamini, Y., & Yekutieli, D. (2001). The control of the false
#'   discovery rate in multiple testing under dependency. Annals of
#'   Statistics, 29(4), 1165-1188. \doi{10.1214/aos/1013699998}
#'
#' On why multiple testing must be addressed at all in gridded remote
#' sensing data (general problem statement):
#' - Gutiérrez-Hernández, O. and García, L.V. (2025, September 17)
#'   Multiple Testing in Remote Sensing: Addressing the Elephant in the
#'   Room. Available at SSRN: https://ssrn.com/abstract=4891512.
#'   \doi{10.2139/ssrn.4891512}
#'
#' On FDR estimation and control specifically under the spatial dependence
#' structure of gridded data (directly motivates this function):
#' - Gutiérrez-Hernández, O., & García, L.V. (2025). False discovery rate
#'   estimation and control in remote sensing: reliable statistical
#'   significance in spatially dependent gridded data. Remote Sensing
#'   Letters, 16(5), 537-548. \doi{10.1080/2150704X.2025.2478664}
#'
#' This function is used (not authored) by [workflow_rta()], this package's own
#' non-prewhitened workflow ([workflow_tst()] instead defaults to the adaptive
#' [fdr_bky()]):
#' - Gutiérrez-Hernández, O. and García, L.V. (2024) Robust Trend Analysis
#'   in Environmental Remote Sensing: A Case Study of Cork Oak Forest
#'   Decline. Remote Sensing, 16(20), 3886. \doi{10.3390/rs16203886}
#' @examples
#' # Five p-values, ranked from most to least significant -- fdr_bh()
#' # tells you how many survive correction at q = 0.05 (the default).
#' sptrends:::fdr_bh(c(0.001, 0.01, 0.02, 0.5, 0.8))
#' @family FDR correction functions
#' @keywords internal
fdr_bh <- function(p, q = 0.05) {
  q_value <- stats::p.adjust(p, method = "BH")
  list(q_value = q_value, reject = q_value <= q)
}

#' Benjamini-Yekutieli (2001) false discovery rate correction
#'
#' Thin wrapper around [stats::p.adjust()] with `method = "BY"` -- base R
#' already implements this correctly, so it is not reimplemented.
#'
#' **Function type:** **Support function** -- computes the BY safeguard
#' used internally by [fdr_correction()]. It is not exported; call
#' `fdr_correction(p, method = "BY")` for a BY-only result.
#'
#' @section Typical use:
#' Supply one family of raw p-values to [fdr_correction()] with
#' `method = "BY"` when arbitrary dependence is a material concern, or
#' compare it with BH as a sensitivity analysis.
#'
#' @section Methodological details:
#' **Methods and method selection**
#'
#' [fdr_bh()] controls FDR under independence and specified forms of
#' positive dependence. A positive Moran statistic can be compatible
#' with that setting, but it does not prove the formal PRDS condition.
#' BY remains valid under arbitrary dependence, at the cost of being
#' more conservative (usually fewer rejected hypotheses) than BH for
#' the same data. It is therefore available as an explicit safeguard
#' when arbitrary dependence is a scientifically material concern, not
#' as the package default. Agreement between BH and BY is a useful
#' sensitivity result, but it does not by itself establish BH's
#' dependence assumptions.
#'
#' **Statistical assumptions and limitations**
#'
#' BY controls FDR without requiring the independence or PRDS conditions
#' used by BH. Its harmonic correction can be substantially conservative,
#' so this broader guarantee may considerably reduce power.
#'
#' **Quality assurance**
#'
#' Adjusted values are generated directly by `stats::p.adjust()` and are
#' also checked through the package's automated FDR tests.
#'
#' @param p Numeric vector of raw p-values in `[0, 1]` (may contain `NA`).
#' @param q Numeric. Target FDR level -- see [fdr_bh()]'s own `q`
#'   documentation for the distinction between this and a per-test
#'   `alpha`, which applies identically here.
#'
#' @return A list with `q_value` (BY-adjusted p-values) and `reject`
#'   (logical).
#' @references
#' Primary method reference:
#' - Benjamini, Y., & Yekutieli, D. (2001). The control of the false
#'   discovery rate in multiple testing under dependency. Annals of
#'   Statistics, 29(4), 1165-1188. \doi{10.1214/aos/1013699998}
#' @examples
#' # The same five p-values fdr_bh()'s own example uses -- compare the
#' # two directly on identical data.
#' sptrends:::fdr_by(c(0.001, 0.01, 0.02, 0.5, 0.8))
#' @family FDR correction functions
#' @keywords internal
fdr_by <- function(p, q = 0.05) {
  q_value <- stats::p.adjust(p, method = "BY")
  list(q_value = q_value, reject = q_value <= q)
}

#' Benjamini-Krieger-Yekutieli (2006) adaptive two-stage FDR correction ("TSBH")
#'
#' Independent implementation (no dependency on `multtest`/`cp4p`) that
#' follows the *code* of `multtest::mt.rawp2adjp(proc = "TSBH")` rather than
#' Definition 6 of the 2006 paper literally. Both share Stage 1 (same
#' `m0_hat` estimate, same `q' = q / (1 + q)`); they differ in Stage 2:
#' Definition 6 recomputes a full second step-up pass at level `q* = q' *
#' m / m0_hat`, while the `multtest` code instead rescales the standard
#' BH-adjusted p-value directly by `m0_hat / m`. The `multtest` threshold is
#' always at least as permissive as Definition 6's (never stricter).
#' Validated bit-for-bit against `cp4p::adjust.p(pi0.method = "bky")` (which
#' calls the real `multtest`) on three real p-value rasters of different
#' sizes and `pi0`, with 100% agreement in all three.
#'
#' **Function type:** **Support function** -- computes the adaptive BKY
#' procedure used internally by [fdr_correction()]. It is not exported;
#' call `fdr_correction(p, method = "BKY")` for a BKY-only result.
#'
#' @section Typical use:
#' Supply one family of raw p-values to [fdr_correction()] with
#' `method = "BKY"`. Choose `bky_implementation = "original"` there only
#' when the literal Definition 6 procedure is required.
#'
#' @section Methodological details:
#' **Methods and method selection**
#'
#' - **Original publication**: Benjamini, Krieger & Yekutieli (2006), the
#'   adaptive two-stage extension of the original BH procedure.
#' - **Main references**: Benjamini, Krieger & Yekutieli (2006) for the
#'   procedure itself; Benjamini & Hochberg (1995) for the underlying FDR
#'   framework it adapts. Full citations appear under "References" below.
#' - **Typical applications**: correcting for multiple testing when a lot
#'   of real signal is expected to be present (the common case in
#'   gridded environmental data) and the extra statistical power over
#'   [fdr_bh()]'s fixed threshold is worth the added complexity -- see
#'   [workflow_tst()] for a workflow that defaults to this method for
#'   that reason.
#'
#' **Statistical assumptions**
#'
#' BKY is adaptive, not a safeguard against arbitrary dependence. Its
#' FDR interpretation requires the assumptions of the selected two-stage
#' procedure; estimating `pi0` does not itself remove dependence among
#' tests. Use [fdr_by()] when an arbitrary-dependence guarantee is needed.
#'
#' **Computational considerations**
#'
#' Both implementations are dominated by ordering or BH adjustment of
#' the valid p-values and are lightweight relative to raster trend tests.
#'
#' **Limitations**
#'
#' The `"multtest"` and `"original"` branches are genuine but distinct
#' second-stage conventions. Under `"original"`, the function returns a
#' rejection decision but `q_value` is `NA`, because that branch does not
#' define the monotone adjusted p-values returned by the rescaling branch.
#'
#' **Quality assurance**
#'
#' The default branch is validated against the behaviour reached through
#' `cp4p` and `multtest`; the literal-paper branch is checked against
#' direct step-up calculations. Automated tests cover missing values,
#' degenerate stages, monotonicity, thresholds, and returned diagnostics.
#'
#' @inheritParams fdr_bh
#' @return A list with `q_value` (or `NA` under `implementation =
#'   "original"`), `reject`, `pi0_hat` (estimated proportion
#'   of true null hypotheses), `m0_hat`, `m`, `r1` (Stage 1 figures), and
#'   `p_sorted`, `thresh_bh`, `thresh_bky` (for [fdr_threshold_plot()]).
#' @references
#' Primary method reference:
#' - Benjamini, Y., Krieger, A. M., & Yekutieli, D. (2006). Adaptive
#'   Linear Step-Up Procedures that Control the False Discovery Rate.
#'   Biometrika, 93(3), 491-507. \doi{10.1093/biomet/93.3.491}
#'
#' Source of the `multtest::mt.rawp2adjp(proc = "TSBH")` behaviour this
#' implementation follows (see "Methodological details" above):
#' - Pollard, K. S., Dudoit, S., & van der Laan, M. J. (2005). Multiple
#'   Testing Procedures: the multtest Package and Applications to
#'   Genomics. In R. Gentleman, V. Carey, W. Huber, R. Irizarry, & S.
#'   Dudoit (eds.), Bioinformatics and Computational Biology Solutions
#'   Using R and Bioconductor, Chapter 15, pp. 249-271. Springer, New
#'   York. \doi{10.1007/0-387-29362-0_15}
#'
#' Theoretical justification for FDR control under positive dependence,
#' relevant background for the adaptive two-stage procedure:
#' - Benjamini, Y., & Yekutieli, D. (2001). The control of the false
#'   discovery rate in multiple testing under dependency. Annals of
#'   Statistics, 29(4), 1165-1188. \doi{10.1214/aos/1013699998}
#'
#' On why multiple testing must be addressed at all in gridded remote
#' sensing data (general problem statement):
#' - Gutiérrez-Hernández, O. and García, L.V. (2025, September 17)
#'   Multiple Testing in Remote Sensing: Addressing the Elephant in the
#'   Room. Available at SSRN: https://ssrn.com/abstract=4891512.
#'   \doi{10.2139/ssrn.4891512}
#'
#' On implementing this specific adaptive procedure for spatiotemporal
#' trend testing (directly motivates this function):
#' - Gutiérrez-Hernández, O., & García, L.V. (2025). Implementing the
#'   Linear Adaptive False Discovery Rate Procedure for Spatiotemporal
#'   Trend Testing. Mathematics, 13(22), 3630. \doi{10.3390/math13223630}
#' @param implementation `"multtest"` (default, unchanged from previous
#'   versions) or `"original"`. Both share Stage 1 (the same `m0_hat`
#'   estimate, the same `q' = q / (1 + q)`); they differ in Stage 2.
#'   `"multtest"` reproduces the implementation used by the CRAN
#'   package `multtest` (specifically `multtest::mt.rawp2adjp(proc =
#'   "TSBH")`, and `cp4p::adjust.p(pi0.method = "bky")`, which calls
#'   it): it rescales the standard BH-adjusted p-value directly by
#'   `m0_hat / m`. `"original"` instead follows the procedure
#'   described in the original Benjamini, Krieger & Yekutieli (2006)
#'   paper literally (Definition 6, Step 3): a full second linear
#'   step-up pass, run on the original p-values, at level `q* = q' *
#'   m / m0_hat`. Neither is more "correct" than the other -- they are
#'   two distinct, real implementations of the same published
#'   procedure, not a canonical version and a variant of it.
#'   `"multtest"`'s threshold is always at least as permissive as
#'   `"original"`'s (never stricter) -- see the "References" section
#'   below for how `"multtest"` was validated against `cp4p`.
#' @examples
#' # 15 p-values, most already close to 0 (a case where several cells
#' # likely have a real trend) -- the adaptive BKY threshold can reject
#' # more hypotheses than BH while targeting FDR control under its
#' # assumptions.
#' p <- c(0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298, 0.0344,
#'        0.0459, 0.3240, 0.4262, 0.5719, 0.6528, 0.7590, 1.0000)
#' sptrends:::fdr_bky(p, q = 0.05)
#'
#' # The original-paper implementation can differ (it is never more
#' # permissive than the default "multtest" implementation, only ever
#' # equal or stricter).
#' sptrends:::fdr_bky(p, q = 0.05, implementation = "original")
#' @family FDR correction functions
#' @keywords internal
fdr_bky <- function(p, q = 0.05, implementation = c("multtest", "original")) {
  implementation <- match.arg(implementation)
  ok <- !is.na(p)
  p_ok <- p[ok]
  m <- length(p_ok)

  bh_adjp <- stats::p.adjust(p_ok, method = "BH")
  q_prime <- q / (1 + q)
  r1 <- sum(bh_adjp < q_prime)
  m0_hat <- m - r1
  pi0_hat <- m0_hat / m

  if (implementation == "multtest") {
    q_value_ok <- bh_adjp * pi0_hat
    reject_ok <- q_value_ok <= q
  } else {
    # Definition 6, Step 3 (Benjamini, Krieger & Yekutieli, 2006): a full
    # second linear step-up pass on the *original* p-values, at level
    # q* = q' * m / m0_hat -- not a rescaling of the BH-adjusted p-value.
    if (r1 == 0) {
      reject_ok <- rep(FALSE, m)
    } else if (r1 == m) {
      reject_ok <- rep(TRUE, m)
    } else {
      q_star <- q_prime * m / m0_hat
      ord <- order(p_ok)
      p_sorted_tmp <- p_ok[ord]
      thresh_tmp <- (seq_len(m) / m) * q_star
      below <- which(p_sorted_tmp <= thresh_tmp)
      k <- if (length(below) > 0) max(below) else 0
      reject_ok <- logical(m)
      reject_ok[ord] <- seq_len(m) <= k
    }
    # No per-cell q-value under this implementation's step-up logic (it
    # is a reject/do-not-reject decision, not a monotone-adjusted
    # p-value) -- report NA rather than inventing one.
    q_value_ok <- rep(NA_real_, m)
  }

  q_value <- rep(NA_real_, length(p))
  q_value[ok] <- q_value_ok
  reject <- rep(NA, length(p))
  reject[ok] <- reject_ok

  p_sorted <- sort(p_ok)
  rank_ <- seq_len(m)
  thresh_bh <- (rank_ / m) * q
  thresh_bky <- (rank_ / m) * q / pi0_hat

  list(q_value = q_value, reject = reject, pi0_hat = pi0_hat, m0_hat = m0_hat,
       m = m, r1 = r1, p_sorted = p_sorted, thresh_bh = thresh_bh,
       thresh_bky = thresh_bky, implementation = implementation)
}

#' Apply false discovery rate (FDR) correction to multiple p-values
#'
#' Applies one or more false discovery rate (FDR) correction procedures
#' to a previously computed set of p-values. Input may be either a
#' single-layer `SpatRaster` or a numeric vector. This function does
#' not compute p-values itself -- it is the final inferential step of
#' the standard sptrends workflow, applied once one raw p-value has
#' already been obtained for every spatial cell (e.g. with
#' [trend_test()] or local [spatial_autocorrelation()]), and it controls
#' the expected false discovery rate
#' across the complete set of simultaneous tests.
#'
#' **Function type:** **Core function** -- one of the core building
#' blocks of TST and RTA. Typically follows [trend_test()] or another
#' inferential function, once one raw p-value exists per spatial cell.
#'
#' @section Typical use:
#' ```
#' raster time series or one spatial field
#'     |
#' trend_test() or spatial_autocorrelation(scope = "local")
#'     |
#' raw *p*-values (`trend$stats$p` or `local$p`)
#'     |
#' fdr_correction()
#'     |
#' adjusted *p*-values + rejection maps at target *q*
#' ```
#' This function does not accept the time series itself. It operates on
#' the family of raw *p*-values produced by a preceding inferential method,
#' whether that method comes from sptrends or elsewhere.
#'
#' @section Methodological details:
#' **Why multiple-testing correction matters for raster data**
#'
#' A raster trend analysis runs one significance test per cell -- often
#' thousands at once. Comparing every raw *p*-value with \eqn{\alpha} = 0.05
#' treats each cell as though it were the only test. That unadjusted
#' interpretation is normally unsuitable as a final raster-wide result:
#' chance alone can produce many apparently significant cells. FDR
#' correction instead evaluates the complete family of *p*-values.
#' Its target *q* limits the expected proportion of false discoveries
#' among the cells declared significant; *q* is not another name for the
#' per-test level \eqn{\alpha}.
#'
#' An FDR rejection is therefore a statement about membership in a
#' family of discoveries whose expected false-discovery proportion is
#' controlled. It is not a separate local error guarantee for that cell.
#' Calling rejected cells "significant pixels" is convenient shorthand,
#' but does not mean that each cell has error probability *q* or satisfies
#' an individually adjusted \eqn{\alpha}. The testing domain must also be
#' defined before inspecting the results: changing it afterwards changes
#' the family, ranks, thresholds, and inferential interpretation.
#'
#' **Methods and method selection**
#'
#' **Why BH and BKY are recommended**: Benjamini & Hochberg (1995) is
#' the classic, widely adopted general-purpose procedure. Across the
#' empirical validation and applications considered during sptrends
#' development, BH has behaved stably, without anomalous inflation of
#' discoveries being observed. Benjamini, Krieger & Yekutieli (2006)
#' extends it
#' adaptively -- estimating how much of the map is likely non-null
#' first, then adjusting the threshold accordingly -- gaining power in
#' exactly the case common to gridded environmental data, where a real
#' trend affects a sizeable share of the map. Neither supersedes the
#' other; see the `method` argument below for when each is the more
#' defensible choice. BY is available for justified arbitrary-dependence
#' settings, but is not recommended routinely because its safeguard can
#' impose a substantial loss of power.
#'
#' **Statistical assumptions**
#'
#' **Role and limit of the spatial autocorrelation diagnostic**: BH's
#' own guarantee formally holds under independence or recognised forms
#' of positive dependence, commonly expressed through PRDS (Benjamini &
#' Yekutieli, 2001). Moran's I,
#' computed via `moran_check = TRUE`, is a *diagnostic* for that
#' assumption -- not a mathematical precondition the function checks
#' or enforces. A positive, significant Moran's I is evidence
#' consistent with positive dependence; it is not proof, and its
#' magnitude depends on the neighbourhood definition used. This
#' function never runs that diagnostic automatically (`moran_check`
#' defaults to `FALSE`) -- the analyst decides whether to spend the
#' extra computation checking it.
#'
#' **Computational considerations**
#'
#' FDR correction is extraordinarily efficient relative to the preceding
#' cell-wise trend analysis. BH and BY are dominated by sorting the
#' *p*-values; BKY adds a lightweight adaptive stage. With large raster
#' families, raster input/output may cost more than the correction itself.
#'
#' **BKY is still an FDR procedure**: its adaptive first stage
#' estimates a proportion of true nulls to sharpen the threshold, but
#' the underlying guarantee it targets is the same false discovery
#' rate BH targets -- it is a refinement of the same family, not a
#' conceptually different kind of correction.
#'
#' **Common interface**
#'
#' Rather than exposing BH and BKY as completely separate workflows,
#' this function provides a common interface returning comparable
#' outputs regardless of the chosen correction -- the same
#' one-interface-many-methods philosophy this package uses throughout
#' (see, e.g., [prewhiten()]'s `method` argument). This function wraps
#' three distinct, individually citable methods -- see [fdr_bh()],
#' [fdr_bky()] and [fdr_by()] for their methodological background,
#' original publications and full reference lists. In short:
#' Benjamini & Hochberg (1995) for BH, Benjamini, Krieger & Yekutieli
#' (2006) for the adaptive BKY extension. Typical applications:
#' correcting for multiple testing in gridded environmental data,
#' whichever method's guarantee (fixed vs. adaptive) fits the analysis
#' at hand.
#'
#' **Limitations and interpretation**
#'
#' FDR control concerns the expected proportion of false discoveries
#' among rejected hypotheses; it does not control the probability of
#' making even one false rejection. The guarantee also depends on the
#' assumptions of the selected procedure. `moran_check` provides only a
#' spatial diagnostic and cannot prove PRDS or choose a valid procedure
#' automatically. The unadjusted comparison stored as `reject_raw` uses
#' the numerical value of `q` only as a visual reference threshold; it
#' is not an FDR-controlled result.
#'
#' **Quality assurance**
#'
#' BH and BY adjusted values are required to match
#' `stats::p.adjust()` exactly. BKY stages, thresholds, rejection
#' monotonicity, missing values, empty and degenerate inputs, and the
#' equivalence of numeric-vector and `SpatRaster` interfaces are covered
#' by automated tests. Integration tests verify that workflow printing,
#' maps, and summaries report only methods actually requested. See
#' `?sptrends` for the package-wide release-check protocol.
#'
#' @param p A single-layer `terra::SpatRaster` of *p*-values, or a numeric
#'   vector of raw *p*-values in `[0, 1]`.
#' @param method Character vector: which correction(s) to compute.
#'   The Usage lists all supported values, `c("BH", "BKY", "BY")`;
#'   when `method` is omitted, the established default remains
#'   `c("BH", "BKY")`. BH is the classic,
#'   fixed-threshold procedure -- see [fdr_bh()] for the full detail.
#'   BKY is the adaptive method originally used in the TST methodology
#'   (Gutiérrez-Hernández & García, 2025): it adapts to the estimated
#'   proportion of true nulls, gaining statistical power over BH while
#'   still controlling the false discovery rate -- see [fdr_bky()] for
#'   the full distinction. `"BY"` (Benjamini and Yekutieli, 2001) is a
#'   third, deliberately opt-in option -- it is never computed unless
#'   explicitly requested
#'   (e.g. `method = "BY"` or `method = c("BH", "BY")`) -- offered as a
#'   comparison point and theoretical safeguard against arbitrary
#'   dependence structures, not as a routine recommendation; see
#'   [fdr_by()] for why.
#' @param q Numeric. Target FDR level. This is not a renamed or adjusted
#'   \eqn{\alpha}: \eqn{\alpha} is a per-test Type I error threshold,
#'   whereas *q* bounds the expected proportion of false discoveries
#'   among all discoveries. A cell is declared significant when its
#'   FDR-adjusted *p*-value is at or below *q*. See [fdr_bh()] for the
#'   full distinction.
#' @param moran_check Logical. If `TRUE`, and `p` is a `SpatRaster`,
#'   automatically run [spatial_autocorrelation()] on the p-value raster.
#'
#'   **What it does**: computes Moran's I on the p-value raster as a
#'   diagnostic of spatial association relevant to the positive-dependence
#'   assumption behind FDR-BH (Benjamini & Yekutieli, 2001), and reports
#'   a plain-language assessment.
#'
#'   **Why**: BH's guarantee formally relies on independence or positive
#'   dependence between tests. A significant positive Moran's I is
#'   compatible with that setting, but does not establish the full PRDS
#'   condition. A non-significant or negative result is inconclusive and
#'   does not make adaptive BKY a dependence safeguard; use BY when control
#'   under arbitrary dependence is scientifically required.
#'
#'   **Limitations**: adds computation time (a permutation loop);
#'   ignored with a warning if `p` is a plain numeric vector (no
#'   spatial structure to test); a diagnostic, not a guarantee -- see
#'   "Methodological details" above. Default `FALSE`.
#'
#'   **More control**: for full control over the test (custom
#'   `nperm`, `connectivity`, `n_cores`, `alternative`), call
#'   [spatial_autocorrelation()] directly instead -- this argument only
#'   covers the common case.
#' @param moran_args A named list of extra arguments, passed unchanged
#'   to [spatial_autocorrelation()] when `moran_check = TRUE` (e.g.
#'   `list(nperm = 999, n_cores = 4)`) -- no intermediate processing.
#' @param report Logical. If `TRUE` (default), automatically print the
#'   summary table ([fdr_summary()]) and draw the diagnostic maps/plots
#'   ([fdr_significance_maps()], [fdr_comparison_barplot()], and, if `BKY`
#'   was requested, [fdr_threshold_plot()]) after computing. Set to `FALSE`
#'   for programmatic use (e.g. when called from [workflow_tst()]).
#' @param verbose Logical. Print progress messages and elapsed time.
#'
#' @return Returns an object of class `c("fdr", "sptrends")`: a list
#'   with `q` (the target FDR level), `method` (the requested
#'   procedures), `p` (the raw input p-values, so `plot()`/`summary()` are
#'   self-contained and don't need the original input kept around
#'   separately), and, for each requested method (`BH`, `BKY`, and/or
#'   `BY` -- `q_BH`/`reject_BH`, `q_BKY`/`reject_BKY`,
#'   `q_BY`/`reject_BY` respectively), numeric vectors of FDR-adjusted
#'   p-values and rejection decisions, plus
#'   `summary_bky` (Stage-1 figures, if BKY was requested) and
#'   `threshold_data` (for [fdr_threshold_plot()]), and `reject_raw`
#'   (the unadjusted comparison at the same numerical threshold). When
#'   the input is a raster, `rasters` additionally contains the raw
#'   p-value raster and corresponding adjusted-value and significance
#'   rasters for the requested methods. If
#'   `moran_check = TRUE`, the list also contains `moran`, the diagnostic
#'   result; `moran_assessment`, a deliberately qualified interpretation;
#'   and the compatibility field `moran_recommendation`, which is `"BH"`
#'   only for significant positive association and `NA` otherwise.
#'   Neither field establishes the positive-regression-dependence
#'   condition required by BH. Use
#'   `print()`/`summary()`/`plot()` -- see [print.sptrends()],
#'   [summary.sptrends()], and [plot.sptrends()].
#' @family FDR correction functions
#' @references
#' Combines the three methods below -- see [fdr_bh()], [fdr_bky()] and
#' [fdr_by()] for the full reference list and the reasoning behind each
#' citation, including why `BY` is offered but not run by default.
#' - Benjamini, Y., & Hochberg, Y. (1995) Controlling the False Discovery
#'   Rate: A Practical and Powerful Approach to Multiple Testing. Journal
#'   of the Royal Statistical Society: Series B, 57, 289-300.
#'   \doi{10.1111/j.2517-6161.1995.tb02031.x}
#' - Benjamini, Y., Krieger, A. M., & Yekutieli, D. (2006) Adaptive
#'   Linear Step-Up Procedures that Control the False Discovery Rate.
#'   Biometrika, 93(3), 491-507. \doi{10.1093/biomet/93.3.491}
#' - Benjamini, Y., & Yekutieli, D. (2001) The control of the false
#'   discovery rate in multiple testing under dependency. Annals of
#'   Statistics, 29(4), 1165-1188. \doi{10.1214/aos/1013699998}
#' @examples
#' # Annual mean NDVI from the bundled environmental dataset.
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' trend <- trend_test(r, report = FALSE, verbose = FALSE)
#'
#' # Corrects trend$stats$p for the fact that every cell was tested at once
#' # (see the "Warning" section of ?trend_test).
#' fdr_result <- fdr_correction(
#'   trend$stats$p,
#'   method = c("BH", "BKY", "BY"),
#'   report = FALSE,
#'   verbose = FALSE
#' )
#' fdr_result$q_BH   # BH-adjusted p-values, one per cell
#' fdr_result$q_BKY  # BKY-adjusted p-values, one per cell
#' fdr_result$q_BY   # BY-adjusted p-values, one per cell
#' summary(fdr_result)
#' # Significance maps and a comparison barplot.
#' plot(fdr_result)
#'
#' # Each method can also be requested separately -- useful when only
#' # one is actually needed downstream.
#' fdr_bh_only  <- fdr_correction(trend$stats$p, method = "BH",
#'                                 report = FALSE, verbose = FALSE)
#' fdr_bky_only <- fdr_correction(trend$stats$p, method = "BKY",
#'                                 report = FALSE, verbose = FALSE)
#' fdr_by_only  <- fdr_correction(trend$stats$p, method = "BY",
#'                                 report = FALSE, verbose = FALSE)
#'
#' \donttest{
#' # The same FDR implementation accepts local spatial p-values. Global
#' # spatial_autocorrelation() performs one test, so it needs no FDR.
#' field <- r[[1]]
#' local <- spatial_autocorrelation(
#'   field, scope = "local", nperm = 99, seed = 1,
#'   report = FALSE, verbose = FALSE
#' )
#' local_fdr <- fdr_correction(
#'   local$p, method = c("BH", "BKY", "BY"), q = 0.05,
#'   report = FALSE, verbose = FALSE
#' )
#' plot(local_fdr)
#' }
#'
#' @param bky_implementation `"multtest"` (default, unchanged from previous
#'   versions) or `"original"`, passed straight to [fdr_bky()] --
#'   see its own documentation for the difference. Ignored if `"BKY"`
#'   is not in `method`.
#' @export
fdr_correction <- function(p, method = c("BH", "BKY", "BY"), q = 0.05,
                            bky_implementation = c("multtest", "original"),
                            moran_check = FALSE, moran_args = list(),
                            report = TRUE, verbose = TRUE) {
  finish_timer <- .sptrends_elapsed_timer("fdr_correction()", verbose)
  on.exit(finish_timer(), add = TRUE)
  # Usage lists every valid procedure, while omission deliberately keeps
  # the established BH+BKY default. BY remains opt-in.
  if (missing(method)) method <- c("BH", "BKY")
  valid_methods <- c("BH", "BKY", "BY")
  if (length(method) == 0) {
    stop("'method' must contain at least one FDR procedure.")
  }
  if (!all(method %in% valid_methods)) {
    stop("'method' must be one or more of: ",
         paste(sprintf('"%s"', valid_methods), collapse = ", "))
  }
  if (!is.numeric(q) || length(q) != 1L || is.na(q) || !is.finite(q) ||
      q <= 0 || q >= 1) {
    stop("'q' must be one finite numeric value strictly between 0 and 1.")
  }
  bky_implementation <- match.arg(bky_implementation)

  is_raster <- inherits(p, "SpatRaster")
  if (is_raster) {
    if (terra::nlyr(p) > 1) {
      if (verbose) message("Note: more than one layer -- using only the first.")
      p <- p[[1]]
    }
    p_vals <- terra::values(p, mat = FALSE)
  } else {
    if (!is.numeric(p)) {
      stop("'p' must be a numeric vector or a single-layer terra SpatRaster.")
    }
    p_vals <- p
  }
  if (length(p_vals) == 0L) {
    stop("'p' must contain at least one value.")
  }
  if (all(is.na(p_vals))) {
    stop("'p' contains no valid p-values (every value is NA).")
  }
  if (any(!is.na(p_vals) & !is.finite(p_vals))) {
    stop("'p' must contain only finite values or NA.")
  }
  rng <- range(p_vals, na.rm = TRUE)
  if (rng[1] < 0 || rng[2] > 1) {
    stop(sprintf(
      paste0("Values outside [0,1] (range: [%.6g, %.6g]) -- this does ",
             "not look like p-value input."),
      rng[1], rng[2]))
  }

  out <- list(q = q, method = method, p = p_vals)

  if ("BH" %in% method) {
    res_bh <- .with_timer("FDR-BH", fdr_bh(p_vals, q = q), verbose)
    out$q_BH <- res_bh$q_value
    out$reject_BH <- res_bh$reject
  }
  if ("BKY" %in% method) {
    res_bky <- .with_timer("FDR-BKY",
                            fdr_bky(p_vals, q = q,
                                    implementation = bky_implementation),
                            verbose)
    out$q_BKY <- res_bky$q_value
    out$reject_BKY <- res_bky$reject
    out$summary_bky <- list(pi0_hat = res_bky$pi0_hat, m0_hat = res_bky$m0_hat,
                             m = res_bky$m, r1 = res_bky$r1)
    out$threshold_data <- list(p_sorted = res_bky$p_sorted,
                                thresh_bh = res_bky$thresh_bh,
                                thresh_bky = res_bky$thresh_bky, m = res_bky$m)
  }
  if ("BY" %in% method) {
    res_by <- .with_timer("FDR-BY", fdr_by(p_vals, q = q), verbose)
    out$q_BY <- res_by$q_value
    out$reject_BY <- res_by$reject
  }

  out$reject_raw <- p_vals <= q

  if (is_raster) {
    t1 <- p[[1]]
    rasterise <- function(v) terra::setValues(t1, as.numeric(v))
    out$rasters <- list(p_value = p)
    if ("BH" %in% method) {
      out$rasters$q_BH <- rasterise(out$q_BH)
      out$rasters$sig_BH <- rasterise(out$reject_BH)
    }
    if ("BKY" %in% method) {
      out$rasters$q_BKY <- rasterise(out$q_BKY)
      out$rasters$sig_BKY <- rasterise(out$reject_BKY)
    }
    if ("BY" %in% method) {
      out$rasters$q_BY <- rasterise(out$q_BY)
      out$rasters$sig_BY <- rasterise(out$reject_BY)
    }
    out$rasters$sig_raw <- rasterise(out$reject_raw)
    for (nm in names(out$rasters)) names(out$rasters[[nm]]) <- nm
  }

  n_valid <- sum(!is.na(p_vals))
  if (verbose) {
    message(sprintf("Valid cells: %d | raw: %d significant",
                     n_valid, sum(out$reject_raw, na.rm = TRUE)))
    if ("BH" %in% method) {
      message(sprintf("FDR-BH: %d significant",
                       sum(out$reject_BH, na.rm = TRUE)))
    }
    if ("BKY" %in% method) {
      message(sprintf("FDR-BKY: %d significant",
                       sum(out$reject_BKY, na.rm = TRUE)))
    }
    if ("BY" %in% method) {
      message(sprintf("FDR-BY: %d significant",
                       sum(out$reject_BY, na.rm = TRUE)))
    }
  }

  if (isTRUE(moran_check)) {
    if (is_raster) {
      if (verbose) {
        message("Running Moran's I diagnostic on the p-value raster ",
                "(moran_check = TRUE)...")
      }
      moran_full_args <- utils::modifyList(
        list(x = p, method = "moran", report = FALSE, verbose = verbose),
        moran_args)
      out$moran <- do.call(spatial_autocorrelation, moran_full_args)
      bh_compatible <- out$moran$sign == "positive" && out$moran$p < 0.05
      out$moran_recommendation <- if (bh_compatible) "BH" else NA_character_
      out$moran_assessment <- if (bh_compatible) {
        "positive spatial association compatible with BH"
      } else {
        "inconclusive for FDR-procedure selection"
      }
      if (verbose) {
        if (bh_compatible) {
          message(sprintf(
            paste0("Moran's I = %.4f (positive, p = %.4f): the observed ",
                   "positive spatial dependence is compatible with the ",
                   "dependence structures under which BH is commonly ",
                   "applied (a significant, positive Moran's I does not ",
                   "formally establish the full PRDS condition BH's own ",
                   "validity relies on)."),
            out$moran$statistic, out$moran$p))
        } else {
          message(sprintf(
            paste0("Moran's I = %.4f (%s, p = %.4f): positive spatial ",
                   "association not confirmed. This diagnostic is ",
                   "inconclusive for choosing an FDR procedure: BKY is ",
                   "adaptive, not a safeguard against arbitrary dependence; ",
                   "consider BY only when such a safeguard is scientifically ",
                   "required."),
            out$moran$statistic, out$moran$sign, out$moran$p))
        }
      }
    } else {
      warning("moran_check=TRUE ignored: 'p' is a plain numeric vector ",
              "with no spatial structure to test.")
    }
  }

  class(out) <- c("fdr", "sptrends")

  if (isTRUE(report)) {
    fdr_summary(out)
    if (is_raster) fdr_significance_maps(out)
    fdr_pvalue_histogram(p)
    fdr_comparison_barplot(out)
    if ("BKY" %in% method) fdr_threshold_plot(out)
    if (!is.null(out$moran)) spatial_autocorrelation_null_plot(out$moran)
  }

  out
}

#' Summary table for an fdr_correction() result
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- called internally by `report = TRUE`, and reachable
#' from outside the package via `summary()`.
#'
#' @param result Output of [fdr_correction()].
#' @param path Character or `NULL`. If supplied, write the table to this
#'   CSV path.
#' @return Invisibly, a data frame.
#' @family FDR correction functions
#' @references
#' See [fdr_bh()] and [fdr_bky()] for the full reference list and the
#' reasoning behind each citation.
#' - Benjamini, Y., & Yekutieli, D. (2001) The control of the false
#'   discovery rate in multiple testing under dependency. Annals of
#'   Statistics, 29(4), 1165-1188. \doi{10.1214/aos/1013699998}
#' @examples
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' trend <- trend_test(r, report = FALSE, verbose = FALSE)
#' fdr_result <- fdr_correction(trend$stats$p, report = FALSE, verbose = FALSE)
#'
#' # A table comparing raw and the selected FDR methods: how many
#' # cells each one calls significant, side by side.
#' sptrends:::fdr_summary(fdr_result)
#'
#' @keywords internal
fdr_summary <- function(result, path = NULL) {
  m <- if (!is.null(result$summary_bky)) {
    result$summary_bky$m
  } else {
    sum(!is.na(result$reject_raw))
  }
  n_raw <- sum(result$reject_raw, na.rm = TRUE)
  tab_rows <- list(data.frame(
    method = "raw (uncorrected)", n_significant = n_raw,
    n_not_significant = m - n_raw,
    pct_significant = round(100 * n_raw / m, 3)))
  if (!is.null(result$reject_BH)) {
    n_bh <- sum(result$reject_BH, na.rm = TRUE)
    tab_rows <- c(tab_rows, list(data.frame(
      method = "FDR-BH", n_significant = n_bh,
      n_not_significant = m - n_bh,
      pct_significant = round(100 * n_bh / m, 3))))
  }
  if (!is.null(result$reject_BKY)) {
    n_bky <- sum(result$reject_BKY, na.rm = TRUE)
    tab_rows <- c(tab_rows, list(data.frame(
      method = "FDR-BKY", n_significant = n_bky,
      n_not_significant = m - n_bky,
      pct_significant = round(100 * n_bky / m, 3))))
  }
  if (!is.null(result$reject_BY)) {
    n_by <- sum(result$reject_BY, na.rm = TRUE)
    tab_rows <- c(tab_rows, list(data.frame(
      method = "FDR-BY", n_significant = n_by,
      n_not_significant = m - n_by,
      pct_significant = round(100 * n_by / m, 3))))
  }
  tab <- do.call(rbind, tab_rows)

  message(sprintf("Valid cells (m): %d | target q: %.3g", m, result$q))
  if (!is.null(result$summary_bky)) {
    message(sprintf("BKY -- pi0_hat: %.6f | m0_hat: %.1f | r1 (stage 1): %d",
                     result$summary_bky$pi0_hat, result$summary_bky$m0_hat,
                     result$summary_bky$r1))
  }

  if (!is.null(path)) {
    utils::write.csv(tab, path, row.names = FALSE)
    message(sprintf("Table written to: %s", path))
  }

  invisible(tab)
}

#' Histogram of input p-values
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- called internally by `report = TRUE`, and reachable
#' from outside the package via `plot(x, which = "pvalue_histogram")`.
#'
#' @param p Numeric vector or single-layer `SpatRaster` of p-values.
#' @param path Character or `NULL`. If supplied, a PNG is written there.
#' @return `NULL`, invisibly.
#' @family FDR correction functions
#' @references
#' See [fdr_bh()] and [fdr_bky()] for the full reference list and the
#' reasoning behind each citation.
#' - Benjamini, Y., & Yekutieli, D. (2001) The control of the false
#'   discovery rate in multiple testing under dependency. Annals of
#'   Statistics, 29(4), 1165-1188. \doi{10.1214/aos/1013699998}
#' @examples
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' trend <- trend_test(r, report = FALSE, verbose = FALSE)
#' fdr_result <- fdr_correction(trend$stats$p, report = FALSE, verbose = FALSE)
#'
#' # A spike near 0 suggests real trends are present; a flat histogram
#' # is what you would expect under the null (no trend anywhere).
#' sptrends:::fdr_pvalue_histogram(trend$stats$p)
#'
#' @keywords internal
fdr_pvalue_histogram <- function(p, path = NULL) {
  p_vals <- if (inherits(p, "SpatRaster")) terra::values(p, mat = FALSE) else p
  graphics::hist(p_vals, breaks = seq(0, 1, by = 0.05), col = "lightblue",
       border = "white",
       xlab = "p-values", main = "Histogram of input p-values")
  if (!is.null(path)) .save_current_plot(path)
  invisible(NULL)
}

#' Significance maps for the selected FDR procedures
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- called internally by `report = TRUE`, and reachable
#' from outside the package via `plot()`.
#'
#' @param result Output of [fdr_correction()] (must have been run on a
#'   raster, so `result$rasters` is populated).
#' @param path Character or `NULL`. If supplied, a PNG is written there.
#' @return `NULL`, invisibly.
#' @family FDR correction functions
#' @references
#' See [fdr_bh()] and [fdr_bky()] for the full reference list and the
#' reasoning behind each citation.
#' - Benjamini, Y., & Yekutieli, D. (2001) The control of the false
#'   discovery rate in multiple testing under dependency. Annals of
#'   Statistics, 29(4), 1165-1188. \doi{10.1214/aos/1013699998}
#' @examples
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' trend <- trend_test(r, report = FALSE, verbose = FALSE)
#' fdr_result <- fdr_correction(trend$stats$p, report = FALSE, verbose = FALSE)
#'
#' # Side-by-side maps: which cells are significant under BH vs. BKY.
#' sptrends:::fdr_significance_maps(fdr_result)
#'
#' @keywords internal
fdr_significance_maps <- function(result, path = NULL) {
  if (is.null(result$rasters)) {
    stop("'result' has no rasters -- fdr_correction() must be run on a ",
         "SpatRaster.")
  }
  # This package's own official brand navy (see .sptrends_brand,
  # sourced from the exact palette used across every version of
  # man/figures/logo.png -- not re-derived from the image itself),
  # used here deliberately rather than a generic dark blue: significant
  # cells get this package's own signature colour, non-significant
  # ones a neutral grey -- a visual identity a reader learns to
  # recognise across every "is this significant" map sptrends draws,
  # not a one-off choice for this function alone.
  colours <- c("grey95", .sptrends_brand$navy)
  labels <- c("Not significant", "Significant")

  n_panels <- sum(!is.null(result$rasters$sig_BH),
                   !is.null(result$rasters$sig_BKY),
                   !is.null(result$rasters$sig_BY))
  op <- graphics::par(mfrow = c(1, max(n_panels, 1)), mar = c(3, 3, 3, 5))
  on.exit(graphics::par(op), add = TRUE)
  if (!is.null(result$rasters$sig_BH)) {
    .safe_categorical_plot(result$rasters$sig_BH, values = c(0, 1),
                            colours = colours, labels = labels,
                            main = sprintf("FDR-BH (q = %.2f)", result$q))
  }
  if (!is.null(result$rasters$sig_BKY)) {
    .safe_categorical_plot(result$rasters$sig_BKY, values = c(0, 1),
                            colours = colours, labels = labels,
                            main = sprintf("FDR-BKY (q = %.2f)", result$q))
  }
  if (!is.null(result$rasters$sig_BY)) {
    .safe_categorical_plot(result$rasters$sig_BY, values = c(0, 1),
                            colours = colours, labels = labels,
                            main = sprintf("FDR-BY (q = %.2f)", result$q))
  }
  if (!is.null(path)) {
    .save_current_plot(path, width = 1400, height = 900, res = 130)
  }
  invisible(NULL)
}

#' Comparison bar chart: raw vs. the selected FDR procedures
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- called internally by `report = TRUE`, and reachable from
#' outside the package via `plot(x, which = "comparison")`.
#'
#' @inheritParams fdr_significance_maps
#' @return `NULL`, invisibly.
#' @family FDR correction functions
#' @references
#' See [fdr_bh()] and [fdr_bky()] for the full reference list and the
#' reasoning behind each citation.
#' - Benjamini, Y., & Yekutieli, D. (2001) The control of the false
#'   discovery rate in multiple testing under dependency. Annals of
#'   Statistics, 29(4), 1165-1188. \doi{10.1214/aos/1013699998}
#' @examples
#' # The plotting function needs an FDR result, not a complete raster trend
#' # analysis. Using a short p-value vector keeps the example focused and
#' # fast while exercising the same plotting code.
#' p <- c(0.001, 0.008, 0.02, 0.04, 0.3, 0.8)
#' fdr_result <- fdr_correction(p, report = FALSE, verbose = FALSE)
#'
#' # A bar chart version of fdr_summary()'s table -- raw vs. BH vs. BKY,
#' # visually.
#' sptrends:::fdr_comparison_barplot(fdr_result)
#'
#' @keywords internal
fdr_comparison_barplot <- function(result, path = NULL) {
  m <- if (!is.null(result$summary_bky)) {
    result$summary_bky$m
  } else {
    sum(!is.na(result$reject_raw))
  }

  panel <- function(sig, title) {
    n_sig <- sum(sig, na.rm = TRUE)
    vals <- rev(c(n_sig, m - n_sig))
    labels <- rev(c("Significant", "Not significant"))
    pct <- round(100 * vals / m, 1)
    txt <- paste0(vals, " (", pct, "%)")
    bp <- graphics::barplot(vals, names.arg = labels,
                             col = c("grey95", .sptrends_brand$navy),
                             main = title, xlab = "Number of cells",
                             horiz = TRUE, las = 1,
                             xlim = c(0, m), xaxs = "i")
    graphics::text(x = vals + m * 0.02, y = bp, labels = txt, pos = 4,
                    cex = 1, col = "black")
  }

  n_panels <- 1 + !is.null(result$reject_BH) + !is.null(result$reject_BKY) +
    !is.null(result$reject_BY)
  op <- graphics::par(mfrow = c(1, n_panels), mar = c(4, 7, 3, 3))
  on.exit(graphics::par(op), add = TRUE)
  panel(result$reject_raw, sprintf("Raw, uncorrected (p <= %.2f)", result$q))
  if (!is.null(result$reject_BH)) {
    panel(result$reject_BH, sprintf("FDR-BH (q = %.2f)", result$q))
  }
  if (!is.null(result$reject_BKY)) {
    panel(result$reject_BKY, sprintf("FDR-BKY (q = %.2f)", result$q))
  }
  if (!is.null(result$reject_BY)) {
    panel(result$reject_BY, sprintf("FDR-BY (q = %.2f)", result$q))
  }

  if (!is.null(path)) {
    .save_current_plot(path, width = 1400, height = 900, res = 130)
  }
  invisible(NULL)
}

#' Rejection-threshold plot for FDR-BH and FDR-BKY
#'
#' Two side-by-side panels, one for FDR-BH and one for FDR-BKY, each
#' showing the sorted p-values against their step-up rejection threshold.
#' Following Benjamini & Hochberg (1995), p-values are shown in increasing
#' order as `p_(i)`: the x-axis is the rank `i` (from 1 to the total
#' number of valid cells `m`), and the y-axis is the ordered p-value
#' `p_(i)` itself.
#'
#' In the left panel, the line is the BH linear step-up threshold,
#' `p_(i) = (i / m) * q`; in the right panel, it is the adaptive BKY
#' threshold, `p_(i) = (i / m) * q_star`, which uses `q_star` -- rescaled
#' from `q` using the estimated proportion of true nulls, `pi0_hat =
#' m0_hat / m` -- rather than a fixed `q`. Points below their panel's
#' Every rank at or before the final cutoff is drawn in blue (rejected),
#' whereas ranks after it are drawn in grey (not rejected). The
#' dashed vertical line marks the cutoff rank `k`: the step-up rule
#' rejects every hypothesis from rank 1 up to `k`, not just the
#' individual points that happen to fall under the line -- `k` is the
#' *last* point (in increasing p order) still below the threshold, and
#' everything at or before it is rejected even if a point in between sits
#' slightly above the line by chance. Because BKY's threshold adapts to
#' `pi0_hat`, it typically sits above BH's fixed-slope line whenever
#' `pi0_hat < 1` -- i.e. whenever some cells are estimated to have a real
#' trend -- which is why the right panel usually shows more rejections
#' than the left one for the same nominal `q`.
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- called internally by `report = TRUE`, and reachable
#' from outside the package via `plot(x, which = "threshold")`.
#'
#' @param result Output of [fdr_correction()] (must include BKY, since
#'   `threshold_data` is only populated then).
#' @param path Character or `NULL`. If supplied, a PNG is written there.
#' @return `NULL`, invisibly.
#' @family FDR correction functions
#' @references
#' See [fdr_bh()] and [fdr_bky()] for the full reference list and the
#' reasoning behind each citation.
#' - Benjamini, Y., & Hochberg, Y. (1995) Controlling the False Discovery
#'   Rate: A Practical and Powerful Approach to Multiple Testing. Journal
#'   of the Royal Statistical Society: Series B, 57, 289-300.
#'   \doi{10.1111/j.2517-6161.1995.tb02031.x}
#' - Benjamini, Y., & Yekutieli, D. (2001) The control of the false
#'   discovery rate in multiple testing under dependency. Annals of
#'   Statistics, 29(4), 1165-1188. \doi{10.1214/aos/1013699998}
#'
#' Detailed graphical interpretation of this exact figure, including a
#' simulation study of the stability of `pi0_hat` across resamples:
#' - Gutiérrez-Hernández, O., & García, L.V. (2025) Implementing the
#'   Linear Adaptive False Discovery Rate Procedure for Spatiotemporal
#'   Trend Testing. Mathematics, 13(22), 3630. \doi{10.3390/math13223630}
#' @examples
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' trend <- trend_test(r, report = FALSE, verbose = FALSE)
#' fdr_result <- fdr_correction(trend$stats$p, report = FALSE, verbose = FALSE)
#'
#' # See ?fdr_threshold_plot for what each axis/line/colour means --
#' # the ordered p-values against the BH/BKY rejection thresholds.
#' sptrends:::fdr_threshold_plot(fdr_result)
#'
#' @keywords internal
fdr_threshold_plot <- function(result, path = NULL) {
  if (is.null(result$threshold_data)) {
    stop("'result' has no threshold_data -- run fdr_correction() with ",
         "method including 'BKY'.")
  }
  du <- result$threshold_data
  rng <- seq_len(du$m)
  q <- result$q

  cutoff <- function(thresh) {
    idx <- which(du$p_sorted <= thresh)
    if (length(idx)) max(idx) else NA_integer_
  }
  cut_bh <- cutoff(du$thresh_bh)
  cut_bky <- cutoff(du$thresh_bky)

  panel <- function(thresh, k, title, ylim) {
    # This package's own official brand colours (.sptrends_brand):
    # navy marks rejected (significant) p-values -- the same colour
    # fdr_significance_maps() and its own barplot use for
    # "Significant", so a reader learns to recognise it as this
    # package's one consistent "yes, significant" signal wherever it
    # appears. Cyan marks the threshold line itself, distinct from
    # both navy and "firebrick" (the cutoff), so all three stay
    # visually separable.
    rejected <- if (is.na(k)) rep(FALSE, length(rng)) else rng <= k
    graphics::plot(rng, du$p_sorted, pch = 19,
         col = ifelse(rejected, .sptrends_brand$navy, "grey"),
         xlab = "p-value rank", ylab = "p-value", main = title, ylim = ylim)
    graphics::lines(rng, thresh, col = .sptrends_brand$cyan, lwd = 2)
    if (!is.na(k)) {
      graphics::abline(v = k, col = "firebrick", lty = 2, lwd = 1.5)
      graphics::text(k, 0.95 * ylim[2], sprintf(" cutoff: k=%d", k),
                      col = "firebrick", pos = 4, cex = 0.9)
    }
    graphics::legend("bottomright",
           legend = c("Rejected", "Not rejected", "Threshold", "Cutoff (k)"),
           col = c(.sptrends_brand$navy, "grey", .sptrends_brand$cyan,
                   "firebrick"),
           pch = c(19, 19, NA, NA),
           lty = c(NA, NA, 1, 2), lwd = c(NA, NA, 2, 1.5), bty = "n",
           cex = 0.85)
    graphics::grid()
  }

  op <- graphics::par(mfrow = c(1, 2), mar = c(5, 5, 4, 2) + 0.1,
                       cex.lab = 1.2, cex.axis = 1.1, cex.main = 1.3)
  on.exit(graphics::par(op), add = TRUE)
  panel(du$thresh_bh, cut_bh, sprintf("FDR-BH threshold (q = %.2f)", q),
        c(0, 1))
  panel(du$thresh_bky, cut_bky,
        sprintf("Adaptive FDR-BKY threshold (q = %.2f)", q),
        c(0, max(du$p_sorted, na.rm = TRUE)))

  if (!is.null(path)) {
    .save_current_plot(path, width = 1400, height = 900, res = 130)
  }
  invisible(NULL)
}

#' Plot a binarised trend map
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported. Unlike the other functions internalised alongside it, this
#' one has no direct single-object S3 wrapper -- it takes a `direction`
#' raster (from `direction_map()`, also not exported) rather
#' than a `"trend_test"`/`"fdr"` object, so there is nothing for a method to
#' dispatch on. Reachable from outside the package via
#' `plot(tst_result, which = "direction")`/`plot(rta_result, which =
#' "direction")` for `workflow_tst()`/`workflow_rta()` results
#' specifically; both this
#' function and `direction_map()` are reachable with `:::` for any
#' other case.
#'
#' @param direction Output of [direction_map()].
#' @param path Character or `NULL`. If supplied, a PNG is written there.
#' @return `NULL`, invisibly.
#' @family FDR correction functions
#' @references
#' This combination of FDR-corrected significance with trend direction is
#' this package's own contribution, not from an external method paper;
#' cited here as the source of the overall TST workflow it belongs to:
#' - Gutiérrez-Hernández, O. and García, L.V. (2025) Uncovering true
#'   significant trends in global greening. Remote Sensing Applications:
#'   Society and Environment, 37, 101377. \doi{10.1016/j.rsase.2024.101377}
#'
#' Underlying theoretical justification for the FDR-BH assumption this
#' significance is based on:
#' - Benjamini, Y., & Yekutieli, D. (2001) The control of the false
#'   discovery rate in multiple testing under dependency. Annals of
#'   Statistics, 29(4), 1165-1188. \doi{10.1214/aos/1013699998}
#' @examples
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' trend <- trend_test(r, report = FALSE, verbose = FALSE)
#' fdr_result <- fdr_correction(trend$stats$p, report = FALSE, verbose = FALSE)
#' direction <- sptrends:::direction_map(trend$stats, fdr_result)
#'
#' # Draws the map computed above.
#' sptrends:::fdr_direction_plot(direction)
#'
#' @keywords internal
fdr_direction_plot <- function(direction, path = NULL) {
  .safe_categorical_plot(direction, values = c(-1, 0, 1),
                          colours = c("firebrick", "grey85", "forestgreen"),
                          labels = c("Decrease", "Not significant", "Increase"),
                          main = "Binarised trend map")
  if (!is.null(path)) .save_current_plot(path)
  invisible(NULL)
}

#' Compare direction of change across raw, FDR-BH, and FDR-BKY
#'
#' A single table with one row per correction method, so you can see how
#' the significant-increase/decrease counts shrink (or don't) as the
#' correction gets stricter.
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic.
#'
#' @param trend The `$stats` field of [trend_test()]'s output.
#' @param fdr_result Output of [fdr_correction()], run on `trend$stats$p`.
#' @param slope Optional single-layer `SpatRaster` (e.g.
#'   [slope_estimator()]'s own `$slope`) whose sign determines direction
#'   instead of the trend test's own statistic. `NULL` (default): use
#'   `trend`'s own `Sm`/`S`/`beta`, matching [direction_map()]'s own
#'   default.
#' @param methods Character vector of methods to include, matching
#'   whichever rejection vectors are present in `fdr_result`
#'   (`c("raw", "BH", "BKY")` by default -- methods not present in
#'   `fdr_result` are skipped with a message, not an error).
#' @param path Character or `NULL`. If supplied, write the table to this
#'   CSV path.
#'
#' @return Invisibly, a data frame with one row per method:
#'   `n_increase`, `n_decrease`, `n_not_significant`, `pct_increase`,
#'   `pct_decrease`.
#'
#' @examples
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' trend <- trend_test(r, report = FALSE, verbose = FALSE)
#' fdr_result <- fdr_correction(trend$stats$p, report = FALSE, verbose = FALSE)
#'
#' # A table version of direction_map(): cell counts and percentages
#' # for increase/decrease/not-significant results, raw and the selected FDR
#' # methods side by side (BH and BKY by default; BY when explicitly
#' # requested).
#' sptrends:::fdr_direction_summary(trend$stats, fdr_result)
#'
#' @family FDR correction functions
#' @references
#' This combination of FDR-corrected significance with trend direction is
#' this package's own contribution, not from an external method paper;
#' cited here as the source of the overall TST workflow it belongs to:
#' - Gutiérrez-Hernández, O. and García, L.V. (2025) Uncovering true
#'   significant trends in global greening. Remote Sensing Applications:
#'   Society and Environment, 37, 101377. \doi{10.1016/j.rsase.2024.101377}
#'
#' Underlying theoretical justification for the FDR-BH assumption this
#' significance is based on:
#' - Benjamini, Y., & Yekutieli, D. (2001) The control of the false
#'   discovery rate in multiple testing under dependency. Annals of
#'   Statistics, 29(4), 1165-1188. \doi{10.1214/aos/1013699998}
#' Not exported. Same reasoning as [fdr_direction_plot()]: it takes
#' `trend`/`fdr_result` as two separate objects rather than one classed
#' object, so there is nothing for a single S3 method to dispatch on.
#' Unlike every other function internalised alongside it, no S3 method
#' calls this one either -- it has no reporting/derived function
#' equivalent left reachable at all from outside the package (other than
#' `:::`). `direction_map()` (also not exported) plus your own
#' tabulation is the closest standalone alternative.
#'
#' @keywords internal
fdr_direction_summary <- function(trend, fdr_result, slope = NULL,
                                   methods = c("raw", "BH", "BKY"),
                                   path = NULL) {
  available <- c(
    raw = !is.null(fdr_result$reject_raw),
    BH  = !is.null(fdr_result$reject_BH),
    BKY = !is.null(fdr_result$reject_BKY)
  )
  methods <- methods[methods %in% names(available)[available]]
  skipped <- setdiff(c("raw", "BH", "BKY"), methods)
  if (length(skipped) > 0) {
    message(sprintf("Skipping (not present in fdr_result): %s",
                     paste(skipped, collapse = ", ")))
  }

  rows <- lapply(methods, function(m) {
    direction <- direction_map(trend, fdr_result, slope = slope,
                                    method = m,
                                    verbose = FALSE)
    vals <- terra::values(direction, mat = FALSE)
    n_tot <- sum(!is.na(vals))
    n_inc <- sum(vals == 1, na.rm = TRUE)
    n_dec <- sum(vals == -1, na.rm = TRUE)
    n_ns <- n_tot - n_inc - n_dec
    data.frame(method = m, n_increase = n_inc, n_decrease = n_dec,
               n_not_significant = n_ns,
               pct_increase = round(100 * n_inc / n_tot, 2),
               pct_decrease = round(100 * n_dec / n_tot, 2))
  })
  tab <- do.call(rbind, rows)

  message(paste(utils::capture.output(print(tab, row.names = FALSE)),
                collapse = "\n"))

  if (!is.null(path)) {
    utils::write.csv(tab, path, row.names = FALSE)
    message(sprintf("Table written to: %s", path))
  }

  invisible(tab)
}
