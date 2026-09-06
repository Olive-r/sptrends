# Informal descriptive label for a Moran's I value

**Explicit warning**: this is a descriptive convention specific to this
package, not a recognised disciplinary standard – there is no universal,
citable threshold for Moran's I magnitude (it depends on the weight
matrix used). Use only as an informal aid, never as citable
justification in a formal write-up. `method = "moran"` results only –
Getis-Ord General G's natural range and expected value under H0 differ
enough from Moran's I that these same thresholds would not transfer
meaningfully (see
[`?spatial_autocorrelation`](https://olive-r.github.io/sptrends/reference/spatial_autocorrelation.md)'s
"Limitations" section).

## Usage

``` r
classify_moran(I)
```

## Arguments

- I:

  Numeric. A Moran's I value (e.g. `result$statistic` from a
  `method = "moran"` result).

## Value

Invisibly, a character label (`"low"`, `"moderate"`, or `"strong"`).

## Details

**Function type:** **Reporting/derived function** – summarises or plots
the output of another function; it does not compute any new statistic.
Not exported – folded into
[`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)
of a `"spatial_autocorrelation"` object (see
[`?print.sptrends`](https://olive-r.github.io/sptrends/reference/print.sptrends.md)),
the category now shown there for every `method = "moran"`
[`spatial_autocorrelation()`](https://olive-r.github.io/sptrends/reference/spatial_autocorrelation.md)
result automatically, without a separate call.

## References

- Moran, P.A.P. (1950) Notes on Continuous Stochastic Phenomena.
  Biometrika, 37(1-2), 17-23.
  [doi:10.1093/biomet/37.1-2.17](https://doi.org/10.1093/biomet/37.1-2.17)

## See also

Other Spatial autocorrelation diagnostic functions:
[`spatial_autocorrelation()`](https://olive-r.github.io/sptrends/reference/spatial_autocorrelation.md),
[`spatial_autocorrelation_null_plot()`](https://olive-r.github.io/sptrends/reference/spatial_autocorrelation_null_plot.md),
[`spatial_autocorrelation_summary()`](https://olive-r.github.io/sptrends/reference/spatial_autocorrelation_summary.md)

## Examples

``` r
# A moderately positive Moran's I would be labelled by this
# function's own convention as follows (see the "Warning" message it
# prints about that label not being a disciplinary standard) --
# called internally wherever a categorical Moran's I label is
# shown, e.g. spatial_autocorrelation()'s own summary output:
# classify_moran(0.32)
```
