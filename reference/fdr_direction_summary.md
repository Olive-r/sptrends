# Compare direction of change across raw, FDR-BH, and FDR-BKY

A single table with one row per correction method, so you can see how
the significant-increase/decrease counts shrink (or don't) as the
correction gets stricter.

## Usage

``` r
fdr_direction_summary(
  trend,
  fdr_result,
  slope = NULL,
  methods = c("raw", "BH", "BKY", "BY"),
  path = NULL
)
```

## Arguments

- trend:

  The `$stats` field of
  [`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md)'s
  output.

- fdr_result:

  Output of
  [`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md),
  run on `trend$stats$p`.

- slope:

  Optional single-layer `SpatRaster` (e.g.
  [`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md)'s
  own `$slope`) whose sign determines direction instead of the trend
  test's own statistic. `NULL` (default): use `trend`'s own
  `Sm`/`S`/`beta`, matching
  [`direction_map()`](https://olive-r.github.io/sptrends/reference/direction_map.md)'s
  own default.

- methods:

  Character vector of methods to include, matching whichever rejection
  vectors are present in `fdr_result` (`c("raw", "BH", "BKY", "BY")` by
  default – methods not present in `fdr_result` are skipped with a
  message, not an error).

- path:

  Character or `NULL`. If supplied, write the table to this CSV path.

## Value

Invisibly, a data frame with one row per method: `n_increase`,
`n_decrease`, `n_not_significant`, `pct_increase`, `pct_decrease`.

## Details

**Function type:** **Reporting/derived function** – summarises or plots
the output of another function; it does not compute any new statistic.

## References

This combination of FDR-corrected significance with trend direction is
this package's own contribution, not from an external method paper;
cited here as the source of the overall TST workflow it belongs to:

- Gutiérrez-Hernández, O. and García, L.V. (2025) Uncovering true
  significant trends in global greening. Remote Sensing Applications:
  Society and Environment, 37, 101377.
  [doi:10.1016/j.rsase.2024.101377](https://doi.org/10.1016/j.rsase.2024.101377)

Underlying theoretical justification for the FDR-BH assumption this
significance is based on:

- Benjamini, Y., & Yekutieli, D. (2001) The control of the false
  discovery rate in multiple testing under dependency. Annals of
  Statistics, 29(4), 1165-1188.
  [doi:10.1214/aos/1013699998](https://doi.org/10.1214/aos/1013699998)
  Not exported. Same reasoning as
  [`fdr_direction_plot()`](https://olive-r.github.io/sptrends/reference/fdr_direction_plot.md):
  it takes `trend`/`fdr_result` as two separate objects rather than one
  classed object, so there is nothing for a single S3 method to dispatch
  on. Unlike every other function internalised alongside it, no S3
  method calls this one either – it has no reporting/derived function
  equivalent left reachable at all from outside the package (other than
  `:::`).
  [`direction_map()`](https://olive-r.github.io/sptrends/reference/direction_map.md)
  (also not exported) plus your own tabulation is the closest standalone
  alternative.

## See also

Other FDR correction functions:
[`fdr_bh()`](https://olive-r.github.io/sptrends/reference/fdr_bh.md),
[`fdr_bky()`](https://olive-r.github.io/sptrends/reference/fdr_bky.md),
[`fdr_by()`](https://olive-r.github.io/sptrends/reference/fdr_by.md),
[`fdr_comparison_barplot()`](https://olive-r.github.io/sptrends/reference/fdr_comparison_barplot.md),
[`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md),
[`fdr_direction_plot()`](https://olive-r.github.io/sptrends/reference/fdr_direction_plot.md),
[`fdr_pvalue_histogram()`](https://olive-r.github.io/sptrends/reference/fdr_pvalue_histogram.md),
[`fdr_significance_maps()`](https://olive-r.github.io/sptrends/reference/fdr_significance_maps.md),
[`fdr_summary()`](https://olive-r.github.io/sptrends/reference/fdr_summary.md),
[`fdr_threshold_plot()`](https://olive-r.github.io/sptrends/reference/fdr_threshold_plot.md)

## Examples

``` r
r <- read_ordered_stack(example_data("vhp_ndvi"))
#> Temporal order auto-detected with pattern '(19[0-9]{2}|20[0-9]{2})'.
#> Automatic mode: order detected from file names. For higher reliability -- especially if the series is not annual -- supplying 'files' explicitly (with 'time' or 'cycle_type') is recommended. See ?read_ordered_stack.
#> Temporal order verification (mandatory, cannot be skipped):
#>  stack_position detected_number                         file
#>               1            1982 VHP_SMN_annual_ndvi_1982.tif
#>               2            1983 VHP_SMN_annual_ndvi_1983.tif
#>               3            1984 VHP_SMN_annual_ndvi_1984.tif
#>               4            1985 VHP_SMN_annual_ndvi_1985.tif
#>               5            1986 VHP_SMN_annual_ndvi_1986.tif
#>               6            1987 VHP_SMN_annual_ndvi_1987.tif
#>               7            1988 VHP_SMN_annual_ndvi_1988.tif
#>               8            1989 VHP_SMN_annual_ndvi_1989.tif
#>               9            1990 VHP_SMN_annual_ndvi_1990.tif
#>              10            1991 VHP_SMN_annual_ndvi_1991.tif
#>              11            1992 VHP_SMN_annual_ndvi_1992.tif
#>              12            1993 VHP_SMN_annual_ndvi_1993.tif
#>              13            1994 VHP_SMN_annual_ndvi_1994.tif
#>              14            1995 VHP_SMN_annual_ndvi_1995.tif
#>              15            1996 VHP_SMN_annual_ndvi_1996.tif
#>              16            1997 VHP_SMN_annual_ndvi_1997.tif
#>              17            1998 VHP_SMN_annual_ndvi_1998.tif
#>              18            1999 VHP_SMN_annual_ndvi_1999.tif
#>              19            2000 VHP_SMN_annual_ndvi_2000.tif
#>              20            2001 VHP_SMN_annual_ndvi_2001.tif
#>              21            2002 VHP_SMN_annual_ndvi_2002.tif
#>              22            2003 VHP_SMN_annual_ndvi_2003.tif
#>              23            2004 VHP_SMN_annual_ndvi_2004.tif
#>              24            2005 VHP_SMN_annual_ndvi_2005.tif
#>              25            2006 VHP_SMN_annual_ndvi_2006.tif
#>              26            2007 VHP_SMN_annual_ndvi_2007.tif
#>              27            2008 VHP_SMN_annual_ndvi_2008.tif
#>              28            2009 VHP_SMN_annual_ndvi_2009.tif
#>              29            2010 VHP_SMN_annual_ndvi_2010.tif
#>              30            2011 VHP_SMN_annual_ndvi_2011.tif
#>              31            2012 VHP_SMN_annual_ndvi_2012.tif
#>              32            2013 VHP_SMN_annual_ndvi_2013.tif
#>              33            2014 VHP_SMN_annual_ndvi_2014.tif
#>              34            2015 VHP_SMN_annual_ndvi_2015.tif
#>              35            2016 VHP_SMN_annual_ndvi_2016.tif
#>              36            2017 VHP_SMN_annual_ndvi_2017.tif
#>              37            2018 VHP_SMN_annual_ndvi_2018.tif
#>              38            2019 VHP_SMN_annual_ndvi_2019.tif
#>              39            2020 VHP_SMN_annual_ndvi_2020.tif
#>              40            2021 VHP_SMN_annual_ndvi_2021.tif
#>              41            2022 VHP_SMN_annual_ndvi_2022.tif
#>              42            2023 VHP_SMN_annual_ndvi_2023.tif

#> Stack built: 42 layers, 146 x 338 cells.
#> >> [read_ordered_stack()] elapsed: 0.12 s
trend <- trend_test(r, report = FALSE, verbose = FALSE)
fdr_result <- fdr_correction(trend$stats$p, report = FALSE, verbose = FALSE)

# A table version of direction_map(): cell counts and percentages
# for increase/decrease/not-significant results, raw and the selected FDR
# methods side by side (BH and BKY by default; BY when explicitly
# requested).
sptrends:::fdr_direction_summary(trend$stats, fdr_result)
#> Skipping (not present in fdr_result): BY
#>  method n_increase n_decrease n_not_significant pct_increase pct_decrease
#>     raw       7167       1837              6671        45.72        11.72
#>      BH       6428       1537              7710        41.01         9.81
#>     BKY       7368       1929              6378        47.00        12.31
```
