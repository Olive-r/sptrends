# Benjamini-Yekutieli (2001) false discovery rate correction

Thin wrapper around
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) with
`method = "BY"` – base R already implements this correctly, so it is not
reimplemented.

## Usage

``` r
fdr_by(p, q = 0.05)
```

## Arguments

- p:

  Numeric vector of raw p-values in `[0, 1]` (may contain `NA`).

- q:

  Numeric. Target FDR level – see
  [`fdr_bh()`](https://olivergh.github.io/sptrends/reference/fdr_bh.md)'s
  own `q` documentation for the distinction between this and a per-test
  `alpha`, which applies identically here.

## Value

A list with `q_value` (BY-adjusted p-values) and `reject` (logical).

## Details

**Function type:** **Support function** – computes the BY safeguard used
internally by
[`fdr_correction()`](https://olivergh.github.io/sptrends/reference/fdr_correction.md).
It is not exported; call `fdr_correction(p, method = "BY")` for a
BY-only result.

## Typical use

Supply one family of raw p-values to
[`fdr_correction()`](https://olivergh.github.io/sptrends/reference/fdr_correction.md)
with `method = "BY"` when arbitrary dependence is a material concern, or
compare it with BH as a sensitivity analysis.

## Methodological details

**Methods and method selection**

[`fdr_bh()`](https://olivergh.github.io/sptrends/reference/fdr_bh.md)
controls FDR under independence and specified forms of positive
dependence. A positive Moran statistic can be compatible with that
setting, but it does not prove the formal PRDS condition. BY remains
valid under arbitrary dependence, at the cost of being more conservative
(usually fewer rejected hypotheses) than BH for the same data. It is
therefore available as an explicit safeguard when arbitrary dependence
is a scientifically material concern, not as the package default.
Agreement between BH and BY is a useful sensitivity result, but it does
not by itself establish BH's dependence assumptions.

**Statistical assumptions and limitations**

BY controls FDR without requiring the independence or PRDS conditions
used by BH. Its harmonic correction can be substantially conservative,
so this broader guarantee may considerably reduce power.

**Quality assurance**

Adjusted values are generated directly by
[`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html) and are
also checked through the package's automated FDR tests.

## References

Primary method reference:

- Benjamini, Y., & Yekutieli, D. (2001). The control of the false
  discovery rate in multiple testing under dependency. Annals of
  Statistics, 29(4), 1165-1188.
  [doi:10.1214/aos/1013699998](https://doi.org/10.1214/aos/1013699998)

## See also

Other FDR correction functions:
[`fdr_bh()`](https://olivergh.github.io/sptrends/reference/fdr_bh.md),
[`fdr_bky()`](https://olivergh.github.io/sptrends/reference/fdr_bky.md),
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
# The same five p-values fdr_bh()'s own example uses -- compare the
# two directly on identical data.
sptrends:::fdr_by(c(0.001, 0.01, 0.02, 0.5, 0.8))
#> $q_value
#> [1] 0.01141667 0.05708333 0.07611111 1.00000000 1.00000000
#> 
#> $reject
#> [1]  TRUE FALSE FALSE FALSE FALSE
#> 
```
