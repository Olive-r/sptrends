# Summarise CMK trend significance

**Function type:** **Reporting/derived function** – summarises or plots
the output of another function; it does not compute any new statistic.
Not exported – called internally by `report = TRUE`, and reachable from
outside the package via
[`summary()`](https://rdrr.io/r/base/summary.html).

## Usage

``` r
trend_summary(trend, alpha = c(0.1, 0.05, 0.01), path = NULL, verbose = TRUE)
```

## Arguments

- trend:

  The `$stats` field of
  [`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md)'s
  output (a 3-layer `SpatRaster` with `p` and `Sm`/`S`).

- alpha:

  Numeric vector of significance thresholds to report. **Uncorrected**
  for multiple testing across cells – run
  [`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md)
  on `trend$p` before treating any of this as a final significance
  result (see
  [`?trend_test`](https://olive-r.github.io/sptrends/reference/trend_test.md),
  section "Warning"). The returned table has one row per value in
  `alpha`; the printed increase/decrease/no-change message uses a single
  reference threshold from among them – `0.05` if present (the default
  vector includes it), otherwise the strictest (smallest) value
  supplied. These three default values are not interchangeable: `0.05`
  is the conventional standard; `0.1` is a more liberal threshold not
  unusual in exploratory trend studies; `0.01` is markedly more
  conservative.

- path:

  Character or `NULL`. If supplied, write the summary table to this CSV
  path.

- verbose:

  Logical. Print the narrative messages (cell count,
  increase/decrease/no-change breakdown). Default `TRUE`. Set to `FALSE`
  to get the returned table silently – used internally by
  [`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)/[`workflow_rta()`](https://olive-r.github.io/sptrends/reference/workflow_rta.md)
  when they call this function only to populate their own
  `trend_summary_table` field, independent of their own `report`
  argument (which already governs whether
  [`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md)
  printed this same summary once, earlier in the same call).

## Value

Invisibly, a data frame with one row per `alpha`.

## References

- Neeti, N. and Eastman, J.R. (2011) A Contextual Mann-Kendall Approach
  for the Assessment of Trend Significance in Image Time Series.
  Transactions in GIS, 15(5), 599-611.
  [doi:10.1111/j.1467-9671.2011.01280.x](https://doi.org/10.1111/j.1467-9671.2011.01280.x)

## See also

Other Contextual Mann-Kendall functions:
[`prepare_cmk_neighbourhood()`](https://olive-r.github.io/sptrends/reference/prepare_cmk_neighbourhood.md),
[`trend_histograms()`](https://olive-r.github.io/sptrends/reference/trend_histograms.md),
[`trend_maps()`](https://olive-r.github.io/sptrends/reference/trend_maps.md)

## Examples

``` r
r <- read_ordered_stack(example_data("vhp_ndvi"))
#> Temporal order auto-detected with pattern '(19[0-9]{2}|20[0-9]{2})'.
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
#> >> [read_ordered_stack()] elapsed: 0.06 s
trend <- trend_test(r, report = FALSE, verbose = FALSE)

# A table with one row per alpha threshold, plus a printed one-line
# summary at the "reference" threshold (0.05 by default -- see the
# alpha argument above).
sptrends:::trend_summary(trend$stats)
#> Cells with complete time series: 15675
#> At alpha=0.05 -- increase: 7167 (45.7%) | decrease: 1837 (11.7%) | no change: 6671 (42.6%)
```
