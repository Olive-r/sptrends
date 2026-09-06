# Benjamini-Hochberg (1995) false discovery rate correction

Thin wrapper around
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) with
`method = "BH"` – base R already implements this correctly, so it is not
reimplemented.

## Usage

``` r
fdr_bh(p, q = 0.05)
```

## Arguments

- p:

  Numeric vector of raw p-values in `[0, 1]` (may contain `NA`).

- q:

  Numeric. Target FDR level.

## Value

A list with `q_value` (BH-adjusted p-values) and `reject` (logical).

## Details

**Function type:** **Support function** – computes the BH procedure used
internally by
[`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md).
Exported as a standalone convenience for callers who only need a plain
vector of p-values corrected, without a raster or the rest of
[`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md)'s
own output structure; call `fdr_correction(p, method = "BH")` instead
for the full result (rasters, when `p` is one, plus the other methods).

## Typical use

Supply one family of raw p-values to
[`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md)
with `method = "BH"`; this internal helper returns the corresponding
adjusted values and rejection decisions.

## Methodological details

**Methods and method selection**

- **Original publication**: Benjamini & Hochberg (1995), the paper that
  introduced false discovery rate control itself.

- **Main references**: Benjamini & Hochberg (1995) for the procedure;
  Benjamini & Yekutieli (2001) for why it remains valid under the
  positive dependence typical of gridded spatial data. Full citations
  appear under "References" below.

- **Typical applications**: correcting for multiple testing when many
  hypotheses are tested at once (e.g. one Mann-Kendall test per pixel in
  a raster) and a fixed, non-adaptive guarantee is preferred over
  [`fdr_bky()`](https://olive-r.github.io/sptrends/reference/fdr_bky.md)'s
  adaptive one – see
  [`workflow_rta()`](https://olive-r.github.io/sptrends/reference/workflow_rta.md)
  for a workflow that defaults to this method specifically for that
  reason.

**Statistical assumptions and limitations**

BH controls FDR under independence and recognised positive-dependence
conditions such as PRDS. Spatial autocorrelation diagnostics can be
compatible with those conditions but do not prove them. Use
[`fdr_by()`](https://olive-r.github.io/sptrends/reference/fdr_by.md)
when control under arbitrary dependence is required.

**Quality assurance**

Adjusted values are generated directly by
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) and are
also checked through the package's automated FDR tests.

## References

Primary method reference:

- Benjamini, Y., & Hochberg, Y. (1995). Controlling the False Discovery
  Rate: A Practical and Powerful Approach to Multiple Testing. Journal
  of the Royal Statistical Society: Series B, 57, 289-300.
  [doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

Theoretical justification for why BH remains valid under the positive
spatial dependence typical of gridded data (for which
[`spatial_autocorrelation()`](https://olive-r.github.io/sptrends/reference/spatial_autocorrelation.md)
or `moran_check` can provide a diagnostic):

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

On FDR estimation and control specifically under the spatial dependence
structure of gridded data (directly motivates this function):

- Gutiérrez-Hernández, O., & García, L.V. (2025). False discovery rate
  estimation and control in remote sensing: reliable statistical
  significance in spatially dependent gridded data. Remote Sensing
  Letters, 16(5), 537-548.
  [doi:10.1080/2150704X.2025.2478664](https://doi.org/10.1080/2150704X.2025.2478664)

This function is used (not authored) by
[`workflow_rta()`](https://olive-r.github.io/sptrends/reference/workflow_rta.md),
this package's own non-prewhitened workflow
([`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)
instead defaults to the adaptive
[`fdr_bky()`](https://olive-r.github.io/sptrends/reference/fdr_bky.md)):

- Gutiérrez-Hernández, O. and García, L.V. (2024) Robust Trend Analysis
  in Environmental Remote Sensing: A Case Study of Cork Oak Forest
  Decline. Remote Sensing, 16(20), 3886.
  [doi:10.3390/rs16203886](https://doi.org/10.3390/rs16203886)

## See also

Other FDR correction functions:
[`fdr_bky()`](https://olive-r.github.io/sptrends/reference/fdr_bky.md),
[`fdr_by()`](https://olive-r.github.io/sptrends/reference/fdr_by.md),
[`fdr_comparison_barplot()`](https://olive-r.github.io/sptrends/reference/fdr_comparison_barplot.md),
[`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md),
[`fdr_direction_plot()`](https://olive-r.github.io/sptrends/reference/fdr_direction_plot.md),
[`fdr_direction_summary()`](https://olive-r.github.io/sptrends/reference/fdr_direction_summary.md),
[`fdr_pvalue_histogram()`](https://olive-r.github.io/sptrends/reference/fdr_pvalue_histogram.md),
[`fdr_significance_maps()`](https://olive-r.github.io/sptrends/reference/fdr_significance_maps.md),
[`fdr_summary()`](https://olive-r.github.io/sptrends/reference/fdr_summary.md),
[`fdr_threshold_plot()`](https://olive-r.github.io/sptrends/reference/fdr_threshold_plot.md)

## Examples

``` r
# Five p-values, ranked from most to least significant -- fdr_bh()
# tells you how many survive correction at q = 0.05 (the default).
fdr_bh(c(0.001, 0.01, 0.02, 0.5, 0.8))
#> $q_value
#> [1] 0.00500000 0.02500000 0.03333333 0.62500000 0.80000000
#> 
#> $reject
#> [1]  TRUE  TRUE  TRUE FALSE FALSE
#> 
```
