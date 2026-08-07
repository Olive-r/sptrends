# Benjamini-Krieger-Yekutieli (2006) adaptive two-stage FDR correction ("TSBH")

Independent implementation (no dependency on `multtest`/`cp4p`) that
follows the *code* of `multtest::mt.rawp2adjp(proc = "TSBH")` rather
than Definition 6 of the 2006 paper literally. Both share Stage 1 (same
`m0_hat` estimate, same `q' = q / (1 + q)`); they differ in Stage 2:
Definition 6 recomputes a full second step-up pass at level
`q* = q' * m / m0_hat`, while the `multtest` code instead rescales the
standard BH-adjusted p-value directly by `m0_hat / m`. The `multtest`
threshold is always at least as permissive as Definition 6's (never
stricter). Validated bit-for-bit against
`cp4p::adjust.p(pi0.method = "bky")` (which calls the real `multtest`)
on three real p-value rasters of different sizes and `pi0`, with 100%
agreement in all three.

## Usage

``` r
fdr_bky(p, q = 0.05, implementation = c("multtest", "original"))
```

## Arguments

- p:

  Numeric vector of raw p-values in `[0, 1]` (may contain `NA`).

- q:

  Numeric. Target FDR level. It is not a renamed or adjusted `alpha`:
  `alpha` is a per-test Type I error threshold, whereas `q` bounds the
  expected proportion of false discoveries among all discoveries – a
  property of the complete rejection set, not of any one test. A `q` of
  `0.05` means that, over repeated applications under the procedure's
  assumptions, the expected false-discovery proportion is controlled at
  5%; it is not a guarantee that every realised set contains at most 5%
  false discoveries.

- implementation:

  `"multtest"` (default, unchanged from previous versions) or
  `"original"`. Both share Stage 1 (the same `m0_hat` estimate, the same
  `q' = q / (1 + q)`); they differ in Stage 2. `"multtest"` reproduces
  the implementation used by the CRAN package `multtest` (specifically
  `multtest::mt.rawp2adjp(proc = "TSBH")`, and
  `cp4p::adjust.p(pi0.method = "bky")`, which calls it): it rescales the
  standard BH-adjusted p-value directly by `m0_hat / m`. `"original"`
  instead follows the procedure described in the original Benjamini,
  Krieger & Yekutieli (2006) paper literally (Definition 6, Step 3): a
  full second linear step-up pass, run on the original p-values, at
  level `q* = q' * m / m0_hat`. Neither is more "correct" than the other
  – they are two distinct, real implementations of the same published
  procedure, not a canonical version and a variant of it. `"multtest"`'s
  threshold is always at least as permissive as `"original"`'s (never
  stricter) – see the "References" section below for how `"multtest"`
  was validated against `cp4p`.

## Value

A list with `q_value` (or `NA` under `implementation = "original"`),
`reject`, `pi0_hat` (estimated proportion of true null hypotheses),
`m0_hat`, `m`, `r1` (Stage 1 figures), and `p_sorted`, `thresh_bh`,
`thresh_bky` (for
[`fdr_threshold_plot()`](https://olivergh.github.io/sptrends/reference/fdr_threshold_plot.md)).

## Details

**Function type:** **Support function** – computes the adaptive BKY
procedure used internally by
[`fdr_correction()`](https://olivergh.github.io/sptrends/reference/fdr_correction.md).
It is not exported; call `fdr_correction(p, method = "BKY")` for a
BKY-only result.

## Typical use

Supply one family of raw p-values to
[`fdr_correction()`](https://olivergh.github.io/sptrends/reference/fdr_correction.md)
with `method = "BKY"`. Choose `bky_implementation = "original"` there
only when the literal Definition 6 procedure is required.

## Methodological details

**Methods and method selection**

- **Original publication**: Benjamini, Krieger & Yekutieli (2006), the
  adaptive two-stage extension of the original BH procedure.

- **Main references**: Benjamini, Krieger & Yekutieli (2006) for the
  procedure itself; Benjamini & Hochberg (1995) for the underlying FDR
  framework it adapts. Full citations appear under "References" below.

- **Typical applications**: correcting for multiple testing when a lot
  of real signal is expected to be present (the common case in gridded
  environmental data) and the extra statistical power over
  [`fdr_bh()`](https://olivergh.github.io/sptrends/reference/fdr_bh.md)'s
  fixed threshold is worth the added complexity – see
  [`workflow_tst()`](https://olivergh.github.io/sptrends/reference/workflow_tst.md)
  for a workflow that defaults to this method for that reason.

**Statistical assumptions**

BKY is adaptive, not a safeguard against arbitrary dependence. Its FDR
interpretation requires the assumptions of the selected two-stage
procedure; estimating `pi0` does not itself remove dependence among
tests. Use
[`fdr_by()`](https://olivergh.github.io/sptrends/reference/fdr_by.md)
when an arbitrary-dependence guarantee is needed.

**Computational considerations**

Both implementations are dominated by ordering or BH adjustment of the
valid p-values and are lightweight relative to raster trend tests.

**Limitations**

The `"multtest"` and `"original"` branches are genuine but distinct
second-stage conventions. Under `"original"`, the function returns a
rejection decision but `q_value` is `NA`, because that branch does not
define the monotone adjusted p-values returned by the rescaling branch.

**Quality assurance**

The default branch is validated against the behaviour reached through
`cp4p` and `multtest`; the literal-paper branch is checked against
direct step-up calculations. Automated tests cover missing values,
degenerate stages, monotonicity, thresholds, and returned diagnostics.

## References

Primary method reference:

- Benjamini, Y., Krieger, A. M., & Yekutieli, D. (2006). Adaptive Linear
  Step-Up Procedures that Control the False Discovery Rate. Biometrika,
  93(3), 491-507.
  [doi:10.1093/biomet/93.3.491](https://doi.org/10.1093/biomet/93.3.491)

Source of the `multtest::mt.rawp2adjp(proc = "TSBH")` behaviour this
implementation follows (see "Methodological details" above):

- Pollard, K. S., Dudoit, S., & van der Laan, M. J. (2005). Multiple
  Testing Procedures: the multtest Package and Applications to Genomics.
  In R. Gentleman, V. Carey, W. Huber, R. Irizarry, & S. Dudoit (eds.),
  Bioinformatics and Computational Biology Solutions Using R and
  Bioconductor, Chapter 15, pp. 249-271. Springer, New York.
  [doi:10.1007/0-387-29362-0_15](https://doi.org/10.1007/0-387-29362-0_15)

Theoretical justification for FDR control under positive dependence,
relevant background for the adaptive two-stage procedure:

- Benjamini, Y., & Yekutieli, D. (2001). The control of the false
  discovery rate in multiple testing under dependency. Annals of
  Statistics, 29(4), 1165-1188.
  [doi:10.1214/aos/1013699998](https://doi.org/10.1214/aos/1013699998)

On why multiple testing must be addressed at all in gridded remote
sensing data (general problem statement):

- Gutiérrez-Hernández, O. and García, L.V. (2025, September 17) Multiple
  Testing in Remote Sensing: Addressing the Elephant in the Room.
  Available at SSRN: https://ssrn.com/abstract=4891512.
  [doi:10.2139/ssrn.4891512](https://doi.org/10.2139/ssrn.4891512)

On implementing this specific adaptive procedure for spatiotemporal
trend testing (directly motivates this function):

- Gutiérrez-Hernández, O., & García, L.V. (2025). Implementing the
  Linear Adaptive False Discovery Rate Procedure for Spatiotemporal
  Trend Testing. Mathematics, 13(22), 3630.
  [doi:10.3390/math13223630](https://doi.org/10.3390/math13223630)

## See also

Other FDR correction functions:
[`fdr_bh()`](https://olivergh.github.io/sptrends/reference/fdr_bh.md),
[`fdr_by()`](https://olivergh.github.io/sptrends/reference/fdr_by.md),
[`fdr_comparison_barplot()`](https://olivergh.github.io/sptrends/reference/fdr_comparison_barplot.md),
[`fdr_correction()`](https://olivergh.github.io/sptrends/reference/fdr_correction.md),
[`fdr_direction_plot()`](https://olivergh.github.io/sptrends/reference/fdr_direction_plot.md),
[`fdr_direction_summary()`](https://olivergh.github.io/sptrends/reference/fdr_direction_summary.md),
[`fdr_pvalue_histogram()`](https://olivergh.github.io/sptrends/reference/fdr_pvalue_histogram.md),
[`fdr_significance_maps()`](https://olivergh.github.io/sptrends/reference/fdr_significance_maps.md),
[`fdr_summary()`](https://olivergh.github.io/sptrends/reference/fdr_summary.md),
[`fdr_threshold_plot()`](https://olivergh.github.io/sptrends/reference/fdr_threshold_plot.md)

## Examples

``` r
# 15 p-values, most already close to 0 (a case where several cells
# likely have a real trend) -- the adaptive BKY threshold can reject
# more hypotheses than BH while targeting FDR control under its
# assumptions.
p <- c(0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298, 0.0344,
       0.0459, 0.3240, 0.4262, 0.5719, 0.6528, 0.7590, 1.0000)
sptrends:::fdr_bky(p, q = 0.05)
#> $q_value
#>  [1] 0.001100000 0.002200000 0.006966667 0.026125000 0.044220000 0.046828571
#>  [7] 0.046828571 0.047300000 0.056100000 0.356400000 0.426200000 0.524241667
#> [13] 0.552369231 0.596357143 0.733333333
#> 
#> $reject
#>  [1]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE
#> [13] FALSE FALSE FALSE
#> 
#> $pi0_hat
#> [1] 0.7333333
#> 
#> $m0_hat
#> [1] 11
#> 
#> $m
#> [1] 15
#> 
#> $r1
#> [1] 4
#> 
#> $p_sorted
#>  [1] 0.0001 0.0004 0.0019 0.0095 0.0201 0.0278 0.0298 0.0344 0.0459 0.3240
#> [11] 0.4262 0.5719 0.6528 0.7590 1.0000
#> 
#> $thresh_bh
#>  [1] 0.003333333 0.006666667 0.010000000 0.013333333 0.016666667 0.020000000
#>  [7] 0.023333333 0.026666667 0.030000000 0.033333333 0.036666667 0.040000000
#> [13] 0.043333333 0.046666667 0.050000000
#> 
#> $thresh_bky
#>  [1] 0.004545455 0.009090909 0.013636364 0.018181818 0.022727273 0.027272727
#>  [7] 0.031818182 0.036363636 0.040909091 0.045454545 0.050000000 0.054545455
#> [13] 0.059090909 0.063636364 0.068181818
#> 
#> $implementation
#> [1] "multtest"
#> 

# The original-paper implementation can differ (it is never more
# permissive than the default "multtest" implementation, only ever
# equal or stricter).
sptrends:::fdr_bky(p, q = 0.05, implementation = "original")
#> $q_value
#>  [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 
#> $reject
#>  [1]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE
#> [13] FALSE FALSE FALSE
#> 
#> $pi0_hat
#> [1] 0.7333333
#> 
#> $m0_hat
#> [1] 11
#> 
#> $m
#> [1] 15
#> 
#> $r1
#> [1] 4
#> 
#> $p_sorted
#>  [1] 0.0001 0.0004 0.0019 0.0095 0.0201 0.0278 0.0298 0.0344 0.0459 0.3240
#> [11] 0.4262 0.5719 0.6528 0.7590 1.0000
#> 
#> $thresh_bh
#>  [1] 0.003333333 0.006666667 0.010000000 0.013333333 0.016666667 0.020000000
#>  [7] 0.023333333 0.026666667 0.030000000 0.033333333 0.036666667 0.040000000
#> [13] 0.043333333 0.046666667 0.050000000
#> 
#> $thresh_bky
#>  [1] 0.004545455 0.009090909 0.013636364 0.018181818 0.022727273 0.027272727
#>  [7] 0.031818182 0.036363636 0.040909091 0.045454545 0.050000000 0.054545455
#> [13] 0.059090909 0.063636364 0.068181818
#> 
#> $implementation
#> [1] "original"
#> 
```
