# Inspect a single cell's (or area's) raw time series interactively

Click a cell (or draw a polygon) on **whatever map is currently
displayed** – a trend map, a significance map, a raw data map, it does
not matter what it shows, only that it shares the same spatial extent as
`x`.

## Usage

``` r
inspect_ts_cell(
  x,
  prewhitened = NULL,
  selection_type = c("point", "polygon"),
  neighbourhood = TRUE,
  connectivity = c("queen", "rook"),
  conf_level = 0.95,
  t = NULL,
  show_neighbours = FALSE,
  slope_method = c("TS", "OLS", "RM"),
  compare_slopes = FALSE,
  verbose = TRUE,
  ...
)
```

## Arguments

- x:

  The full time series stack (a
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  one layer per time step) – this is always where the plotted raw series
  comes from, regardless of what is currently displayed on screen.

- prewhitened:

  Optional. The full list returned by
  [`prewhiten()`](https://olive-r.github.io/sptrends/reference/prewhiten.md)
  (not just its `$series`) – if supplied, a second panel shows the same
  location's prewhitened series, its own Theil-Sen fit and confidence
  interval, and whether that location was actually modified by
  prewhitening (many cells are not, if their own Durbin-Watson statistic
  never crossed the gating threshold; see
  [`?prewhiten`](https://olive-r.github.io/sptrends/reference/prewhiten.md)).
  Default `NULL`: only the raw-data panel is drawn.

- selection_type:

  `"point"` (default): click a single cell. `"polygon"`: draw a polygon
  (left-click to add vertices, press `Esc` or right-click to finish).

- neighbourhood:

  Logical, only used when `selection_type = "point"`. If `TRUE`
  (default), the clicked cell and its queen (or rook) neighbours are
  combined – see the "How each mode aggregates its series" section below
  – before estimating anything, borrowing spatial context the same way
  [`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md)
  does for significance. If `FALSE`, only the clicked cell's own series
  is used. Ignored when `selection_type = "polygon"` (a polygon is
  already an explicit choice of area). Intended for exploratory
  visualisation, at this one location, rather than a substitute for
  formal inference – the slope shown here does not carry the same
  significance guarantee CMK's own variance-adjusted statistic does.

- connectivity:

  `"queen"` (default, 8 neighbours) or `"rook"` (4 neighbours), passed
  to
  [`terra::adjacent()`](https://rspatial.github.io/terra/reference/adjacent.html).
  Only relevant when `neighbourhood = TRUE`.

- conf_level:

  Confidence level for the fitted slope's confidence interval, reported
  in the legend as text (not drawn as a shaded band – see "Confidence
  interval" below). Default `0.95`. Ignored when `slope_method = "RM"`
  or `compare_slopes = TRUE` (see both below).

- t:

  Finite numeric vector of unique, strictly increasing time points, one
  per layer. Defaults to `1:nlyr(x)`.

- show_neighbours:

  Logical, only used when `neighbourhood = TRUE` and
  `selection_type = "point"`. Answers, at a glance, a question any user
  looking at an unusual pixel asks: **is this cell's trend
  representative of its neighbourhood, or an outlier the aggregation is
  smoothing over?** If `TRUE`, draws a second figure after the main
  panel(s): a small-multiples grid, one mini-panel per cell actually
  aggregated (the clicked cell, highlighted, plus each of its individual
  queen/rook neighbours), each with its own raw series and fitted line
  (same `slope_method` as the main panel) – not the median-aggregated
  series the main panel shows, but each contributing cell on its own.
  Ignored (with a message) if the clicked cell has no neighbours with
  complete data (e.g. a corner cell with all-NA neighbours), since there
  would be nothing to compare against. Default `FALSE`.

- slope_method:

  Which estimator
  [`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md)
  uses for the fitted line(s) shown here. `"TS"` (default): robust to
  outliers, matches this package's own default trend workflow. `"OLS"`:
  ordinary least squares – faster, but sensitive to outliers the way a
  single anomalous time step in the plot above can be. `"RM"`: Siegel's
  repeated median – more robust than `"TS"`, but no confidence interval
  is shown (none implemented for it in this package). Ignored when
  `compare_slopes = TRUE`. This only affects this function's own quick
  single-cell fit; it has no bearing on which method a
  [`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)/[`workflow_rta()`](https://olive-r.github.io/sptrends/reference/workflow_rta.md)
  run itself used.

- compare_slopes:

  Logical. If `TRUE`, ignores `slope_method` and draws all three
  estimators (`"TS"`, `"OLS"`, `"RM"`) as separate lines on the same
  panel(s), point estimates only – no confidence intervals in this mode,
  since the three use genuinely different inferential frameworks (or
  none, for `"RM"`) and no single interval formula would honestly apply
  to all three. Default `FALSE`.

- verbose:

  Logical. If `TRUE` (default), print the click/draw instructions and
  the `show_neighbours` no-effect note; matches the `verbose` convention
  used throughout this package, though the interactive click/draw step
  itself still happens either way – this only silences the accompanying
  messages and elapsed time, not the interaction.

- ...:

  Ignored.

## Value

Returns invisibly, a list with `cell` (the clicked/representative cell
number) and `raw` (a list with `series`, `slope`, `ci_lower`,
`ci_upper`, `conf_level` for the raw-data panel). If `prewhitened` was
supplied, also `prewhitened` (the same fields for that panel) and
`n_modified`/`n_total` (how many of the aggregated cells were actually
modified by prewhitening). If `show_neighbours = TRUE` and at least one
neighbour had a complete series, also `neighbours`: a list, one element
per plotted cell (including the clicked one), each with `cell`, `slope`,
and `is_centre`.

## Details

**Why inspect a single cell?**: raster trend maps summarise thousands of
time series into one image. Inspecting an individual location helps
determine whether an apparently unusual pixel reflects a genuine
temporal pattern, an isolated outlier, or a behaviour representative of
its surrounding neighbourhood – a question no summary map, on its own,
can answer.

**What it shows**: the raw time series behind the clicked location, with
a single fitted line overlaid (Theil-Sen or OLS – see `slope_method`),
together with its confidence interval.

**Prewhitening comparison**: optionally, a second panel shows the same
location's *prewhitened* series side by side (see `prewhitened` below),
so you can see the effect of that step at exactly the location you are
looking at.

**When to use it**: to get oriented early on (click around right after
loading data) and to investigate a specific cell that looks anomalous on
a map you have already produced.

**Built from the package's own pieces**: this function is not a separate
implementation of its own – the fitted line calls
[`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md)'s
own estimator, the neighbourhood aggregation follows
[`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md)'s
own logic, and the prewhitened panel reads
[`prewhiten()`](https://olive-r.github.io/sptrends/reference/prewhiten.md)'s
own output directly. It is, in that sense, less a standalone utility
than a visual demonstration of how the rest of sptrends' pieces fit
together, applied to one location at a time.

**Function type:** **Interactive exploration function** – combines
existing sptrends estimators and neighbourhood logic at one selected
location; it introduces no new inferential method.

## Typical use

    raster time series + a displayed map with matching geometry
        |
    inspect_ts_cell()
        |
    selected cell or area -> temporal plot + fitted slope

Optionally supply the complete
[`prewhiten()`](https://olive-r.github.io/sptrends/reference/prewhiten.md)
result to compare the raw and transformed series at the same location.

## Methodological details

**Methods and method selection**

`slope_method` selects Theil-Sen, OLS, or repeated median for the
exploratory fit. `compare_slopes = TRUE` displays all three point
estimates together; confidence intervals are then omitted because no
single inferential framework applies to all three estimators.

**How each mode aggregates its series**

All three modes – point-only, point-with-neighbourhood, and polygon –
follow the *same* two-step procedure, in the *same* order, so a single
confidence interval formula applies rigorously to all of them: first
take the per-time-step **median of raw values** across whichever cells
are included (just the one clicked cell, the clicked cell plus its
neighbours, or every cell inside the drawn polygon), producing one
aggregated series; *then* fit a single Theil-Sen slope to that series.
The same aggregation (same cells) is applied to the prewhitened series
too, when `prewhitened` is supplied, so the two panels are directly
comparable. Aggregating values first and estimating second is what makes
the Sen/Gilbert confidence interval below valid – it is defined for the
pairwise slopes of a single series, not for a median of several
already-computed slope estimates (the reverse order), which has no
standard interval formula. This also means `neighbourhood = TRUE`'s
result is **not** the same number as
`slope_estimator(x, smooth_neighbourhood = TRUE)` at that same cell in a
full-map run – that other mechanism deliberately aggregates in the
opposite order (median of independently-computed per-cell slopes); see
its own documentation for why. The two exist for different purposes and
are not meant to match numerically.

**Statistical assumptions and confidence intervals**

**What it represents**: uncertainty in the *rate of change* itself, not
a prediction interval around the raw data points. Reported as text in
the legend (`ci_lower`/`ci_upper` in the returned value too), not drawn
as a shaded band on the plot – only the single fitted line and its
slope/intercept are drawn. Not shown when `slope_method = "RM"` or
`compare_slopes = TRUE`.

**How it is computed** (the standard rank-based method for the Theil-Sen
slope; Sen, 1968; Gilbert, 1987, pp. 217-219): given the \\N\\ pairwise
slopes of the (already-aggregated) series, sorted, the interval is
\\\[\hat\beta\_{(M1)}, \hat\beta\_{(M2+1)}\]\\, where \\M1 = (N -
C\_\alpha)/2\\, \\M2 = (N + C\_\alpha)/2\\, and \\C\_\alpha =
z\_{1-\alpha/2}\sqrt{Var(S)}\\ (the same Mann-Kendall variance of `S`
used throughout this package, tie-corrected). `M1`/ `M2` are rounded to
the nearest integer rather than interpolated between adjacent ranked
slopes, unlike Gilbert's own recommendation – a deliberate
simplification, documented here rather than left silent. This formula
applies when `slope_method = "TS"` (the default); for
`slope_method = "OLS"`, the interval instead uses the standard
parametric simple-linear-regression interval (which assumes
normally-distributed, homoscedastic residuals) – a genuinely different
formula, not the same computation reused, since it rests on different
assumptions matching its own point estimate.

**Limitations**

This is an exploratory, interactive view of selected locations, not a
raster-wide significance procedure. Neighbourhood or polygon medians
change the series being fitted, and repeated-median or multi-estimator
displays do not include confidence intervals.

**Quality assurance**

Tests cover point and polygon selection, neighbourhood aggregation,
raw/prewhitened comparisons, all slope choices, confidence intervals,
time metadata, incomplete series, interactive selection failures and
silent operation. Core numerical results are compared with direct
calculations and
[`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md).
See
[`?sptrends`](https://olive-r.github.io/sptrends/reference/sptrends-package.md)
for the common release-check protocol.

## References

Confidence interval method:

- Sen, P.K. (1968) Estimates of the regression coefficient based on
  Kendall's tau. Journal of the American Statistical Association, 63,
  1379-1389.
  [doi:10.1080/01621459.1968.10480934](https://doi.org/10.1080/01621459.1968.10480934)

- Gilbert, R.O. (1987) Statistical Methods for Environmental Pollution
  Monitoring. Van Nostrand Reinhold, New York, pp. 217-219.

Origin of the Var(S) formula reused for the confidence interval (the
same one
[`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md)
builds on for CMK's adjusted variance; see its own references for that
extension):

- Mann, H.B. (1945) Nonparametric tests against trend. Econometrica,
  13(3), 245-259. [doi:10.2307/1907187](https://doi.org/10.2307/1907187)

- Kendall, M.G. (1975) Rank Correlation Methods (4th edn). Charles
  Griffin, London. No DOI available (pre-DOI-era publication).

## See also

[`prewhiten()`](https://olive-r.github.io/sptrends/reference/prewhiten.md)
for transformed time series,
[`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md)
for trend significance, and
[`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md)
for raster-wide magnitude estimation.

## Examples

``` r
if (FALSE) { # \dontrun{
# Interactive -- run this yourself, requires clicking on a plot.

# Annual mean NDVI from the bundled environmental dataset.
r <- read_ordered_stack(example_data("vhp_ndvi"))
terra::plot(r[[1]])         # any single-layer map with the same extent as r

# That cell's own raw series, nothing borrowed from its surroundings.
inspect_ts_cell(r, neighbourhood = FALSE)

# Compare against the function's own default -- the same cell combined
# with its queen neighbourhood -- to see how much borrowing spatial
# context changes the fit.
inspect_ts_cell(r)

inspect_ts_cell(r, selection_type = "polygon")       # draw an area instead

# Raw vs. prewhitened, side by side, at the clicked location
pw <- prewhiten(r, report = FALSE, verbose = FALSE)
inspect_ts_cell(r, prewhitened = pw)

# Is the clicked cell's trend representative of its neighbourhood,
# or an outlier the median aggregation is smoothing over? Draws a
# second figure: one mini-panel per neighbour, plus the clicked cell
# itself (highlighted), each with its own Theil-Sen fit.
inspect_ts_cell(r, show_neighbours = TRUE)

# Faster, but not robust to the effect a single anomalous time step
# can have on the fitted line -- see the "slope_method" argument.
inspect_ts_cell(r, slope_method = "OLS")

# Siegel's repeated median -- more robust than Theil-Sen, but no
# confidence interval is shown for it (none is implemented here).
inspect_ts_cell(r, slope_method = "RM")

# All three estimators at once, as three lines with no intervals.
inspect_ts_cell(r, compare_slopes = TRUE)
} # }
```
