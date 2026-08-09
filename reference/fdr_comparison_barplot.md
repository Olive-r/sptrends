# Comparison bar chart: raw vs. the selected FDR procedures

**Function type:** **Reporting/derived function** – summarises or plots
the output of another function; it does not compute any new statistic.
Not exported – called internally by `report = TRUE`, and reachable from
outside the package via `plot(x, which = "comparison")`.

## Usage

``` r
fdr_comparison_barplot(result, path = NULL)
```

## Arguments

- result:

  Output of
  [`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md)
  (must have been run on a raster, so `result$rasters` is populated).

- path:

  Character or `NULL`. If supplied, a PNG is written there.

## Value

`NULL`, invisibly.

## References

See [`fdr_bh()`](https://olive-r.github.io/sptrends/reference/fdr_bh.md)
and
[`fdr_bky()`](https://olive-r.github.io/sptrends/reference/fdr_bky.md)
for the full reference list and the reasoning behind each citation.

- Benjamini, Y., & Yekutieli, D. (2001) The control of the false
  discovery rate in multiple testing under dependency. Annals of
  Statistics, 29(4), 1165-1188.
  [doi:10.1214/aos/1013699998](https://doi.org/10.1214/aos/1013699998)

## See also

Other FDR correction functions:
[`fdr_bh()`](https://olive-r.github.io/sptrends/reference/fdr_bh.md),
[`fdr_bky()`](https://olive-r.github.io/sptrends/reference/fdr_bky.md),
[`fdr_by()`](https://olive-r.github.io/sptrends/reference/fdr_by.md),
[`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md),
[`fdr_direction_plot()`](https://olive-r.github.io/sptrends/reference/fdr_direction_plot.md),
[`fdr_direction_summary()`](https://olive-r.github.io/sptrends/reference/fdr_direction_summary.md),
[`fdr_pvalue_histogram()`](https://olive-r.github.io/sptrends/reference/fdr_pvalue_histogram.md),
[`fdr_significance_maps()`](https://olive-r.github.io/sptrends/reference/fdr_significance_maps.md),
[`fdr_summary()`](https://olive-r.github.io/sptrends/reference/fdr_summary.md),
[`fdr_threshold_plot()`](https://olive-r.github.io/sptrends/reference/fdr_threshold_plot.md)

## Examples

``` r
# The plotting function needs an FDR result, not a complete raster trend
# analysis. Using a short p-value vector keeps the example focused and
# fast while exercising the same plotting code.
p <- c(0.001, 0.008, 0.02, 0.04, 0.3, 0.8)
fdr_result <- fdr_correction(p, report = FALSE, verbose = FALSE)

# A bar chart version of fdr_summary()'s table -- raw vs. BH vs. BKY,
# visually.
sptrends:::fdr_comparison_barplot(fdr_result)



```
