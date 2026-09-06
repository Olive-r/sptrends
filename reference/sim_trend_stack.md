# Generate a synthetic gridded time series with known true trends

Generates benchmark datasets with a known ground truth, for evaluating
spatiotemporal trend-detection methods. Concretely: a synthetic gridded
time series, with the true, known slope at every cell kept alongside the
data itself, for runnable `@examples`, vignettes, and unit tests, rather
than as a realistic environmental simulation. Because the ground truth
is known exactly (something no real dataset offers), this is the tool
behind sptrends' own worked examples of what a trend-detection method
actually gets right or wrong: compare `true_slope` against a fitted
method's estimated slope, and against which cells it calls significant,
to see estimation error and Type I/Type II error directly instead of
having to take a method's output on faith – the mechanism that lets
power, false discovery rate, sensitivity, specificity, and robustness be
demonstrated *objectively* for any method, not just described in prose.

## Usage

``` r
sim_trend_stack(
  nrow = 15,
  ncol = 15,
  n_time = 10,
  trend_strength = 0.15,
  trend_shape = c("radial", "gradient", "block", "square", "rectangle", "ellipse",
    "custom"),
  trend_fraction = 0.9,
  ar1 = 0.3,
  noise_sd = 1,
  noise_dist = c("gaussian", "t"),
  t_df = 4,
  smooth_radius = 1L,
  spatial_rho = 1,
  spatial_model = c("legacy", "independent", "gaussian", "exponential", "matern"),
  spatial_range = 1,
  spatial_smoothness = 0.5,
  signal_size = NULL,
  signal_location = c("centre", "random"),
  signal_angle = 0,
  signal_axis_ratio = 1,
  custom_mask = NULL,
  constant_block = TRUE,
  break_type = c("none", "mean", "slope"),
  break_time = NULL,
  break_fraction = 0.3,
  break_magnitude = 2,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- nrow, ncol:

  Integer. Raster dimensions in cells.

- n_time:

  Integer. Number of time steps (layers).

- trend_strength:

  Numeric. Maximum magnitude of the true slope (see "How the trend is
  generated – trend field" above for how it varies spatially according
  to `trend_shape`).

- trend_shape:

  Spatial pattern of the true slope field: `"radial"`, `"gradient"`,
  `"block"`, `"square"`, `"rectangle"`, `"ellipse"`, or `"custom"`; see
  "How the trend is generated – trend field" above.

- trend_fraction:

  Numeric in `[0, 1]`. Proportion of cells that keep a non-zero true
  slope; the rest are forced to exactly `0`. Default `0.9`. Set to `0`
  for a complete null field (no trend anywhere), or `1` for a trend at
  every cell `trend_shape` assigns one to.

- ar1:

  Numeric in `(-1, 1)`. AR(1) autocorrelation coefficient of the noise
  term across time, applied independently within each cell.

- noise_sd:

  Numeric `> 0`. Standard deviation of the noise term – the main control
  over how easy or hard the true trend is to detect against the
  background noise (signal-to-noise ratio). Default `1`.

- noise_dist:

  `"gaussian"` (default) or `"t"` – the distribution of the noise term,
  always rescaled to the standard deviation `noise_sd` regardless of
  shape. This matters for comparing methods: Mann-Kendall-based tests
  (and Theil-Sen) are rank-based and robust to heavy-tailed,
  outlier-prone noise, while a parametric method like OLS is the more
  *efficient* choice specifically when noise is Gaussian – the classic
  case for preferring a non-parametric method only shows up with
  non-Gaussian noise. `"t"` draws from a Student's t distribution with
  `t_df` degrees of freedom instead, giving heavier tails (more extreme
  outliers) at the same standard deviation.

- t_df:

  Numeric `> 2`. Degrees of freedom for the noise distribution when
  `noise_dist = "t"`; ignored otherwise. Lower values give heavier tails
  (more frequent/extreme outliers); values just above `2` are extreme,
  values above `30` are practically indistinguishable from Gaussian.
  Must exceed `2` for the distribution to have finite variance (required
  to rescale it to `noise_sd`).

- smooth_radius:

  Integer `>= 0`. Radius, in cells, of the focal mean window used to
  introduce spatial autocorrelation into the noise (see "How spatial
  autocorrelation is generated" above). `0` disables this: cells get
  independent noise with no spatial structure beyond the trend field
  itself. `1` (default) smooths over a 3x3 window; larger values give
  smoother, more strongly autocorrelated fields, at the cost of more
  computation.

- spatial_rho:

  Numeric in `[0, 1]`. For the legacy model, the weight assigned to
  focal-smoothed noise; this is not a target Moran's I. For Gaussian and
  exponential covariance models, the adjacent-cell correlation of the
  latent Gaussian field. It must be positive for these formal models. It
  is ignored by the independent and Matérn models.

- spatial_model:

  Spatial noise model. `"legacy"` preserves the focal smoother used by
  earlier versions. `"independent"`, `"gaussian"`, `"exponential"`, and
  `"matern"` use a formal stationary covariance model simulated by
  circulant embedding and FFT.

- spatial_range:

  Positive covariance scale for the Matérn model.

- spatial_smoothness:

  Positive Matérn smoothness parameter.

- signal_size:

  One number or a two-number vector giving the height and width, in
  cells, of block, square, rectangular, or elliptical signal regions.
  Defaults to half the raster in each dimension for `"block"` (its
  historical default), or a third of the raster for `"square"`,
  `"rectangle"` and `"ellipse"`.

- signal_location:

  `"centre"` or `"random"`; location of an exact geometric signal region
  (`"block"`, `"square"`, `"rectangle"`, `"ellipse"` or `"custom"`). Has
  no effect on `"radial"` or `"gradient"`, which are not exact geometric
  regions.

- signal_angle:

  Rotation of a geometric signal region in degrees, or `"random"` to
  draw a new orientation in each realisation.

- signal_axis_ratio:

  Positive width multiplier for an ellipse.

- custom_mask:

  Logical/numeric matrix, vector, or single-layer `SpatRaster` defining
  the signal when `trend_shape = "custom"`.

- constant_block:

  Logical. If `TRUE`, sets a small block of cells to a constant value
  (zero temporal variance) to exercise degenerate-case handling.
  Independent of `trend_shape = "block"` above, which is about the
  *trend* pattern, not a zero-variance edge case.

- break_type:

  `"none"` (default): no structural break, this function's original
  behaviour, unchanged. `"mean"`: a step change in level at `break_time`
  – ground truth for a change-point method like Pettitt's test.
  `"slope"`: a change in slope at `break_time` (continuous in level, a
  kink not a jump) – ground truth for distinguishing "did the trend
  change partway through" from "is there a trend at all". A break and a
  monotonic trend (`trend_strength`/`trend_fraction` above) compose
  independently on any given cell – a cell can have a trend only, a
  break only, both, or neither; the two are not mutually exclusive, nor
  does one replace the other.

- break_time:

  Integer or `NULL`. Time step *after* which the break takes effect (so
  `break_time = 5` means steps 1-5 are unaffected, the break shows in
  step 6 onward). `NULL` (default): the middle of the series,
  `round(n_time / 2)`. Must be in `[1, n_time - 1]` so there is at least
  one time step on each side. Ignored when `break_type = "none"`.

- break_fraction:

  Fraction of the same spatially coherent blocks used for
  `trend_fraction` above that get a break, independently of which blocks
  (if any) got a trend. Default `0.3`. Ignored when
  `break_type = "none"`.

- break_magnitude:

  Numeric. For `break_type = "mean"`: the size of the step added to the
  level from `break_time` onward. For `break_type = "slope"`: the size
  of the *additional* slope (per time step) that kicks in from
  `break_time` onward, on top of whatever slope (possibly zero) that
  cell already had. Default `2`. Ignored when `break_type = "none"`.

- seed:

  Integer or `NULL`. Random seed for reproducibility.

- verbose:

  Logical. If `TRUE`, report simulation progress, elapsed duration and
  estimated time remaining. Set to `FALSE` inside large Monte Carlo
  experiments when
  [`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md)
  provides the outer progress display.

## Value

An object of class `"sptrends_simulation"` and `"sptrends"`, supporting
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html), with:

- series:

  A
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  with `n_time` layers named `"t1"`, `"t2"`, ...; each layer represents
  one time step. Pass this to
  [`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)
  or the individual trend functions.

- true_slope:

  A single-layer
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html):
  the exact true slope used to generate each cell's series (`0` where
  `trend_fraction` assigned no trend). Compare directly against a fitted
  [`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md)
  to see estimation error, or against
  [`direction_map()`](https://olive-r.github.io/sptrends/reference/direction_map.md)
  to see which cells a method calls significant vs. which cells truly
  have a trend.

- true_signal:

  A binary `SpatRaster` equal to `1` where the exact true slope is
  non-zero.

- true_direction:

  A `SpatRaster` containing `-1`, `0`, or `1` for the known trend
  direction.

- true_break:

  A single-layer
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html),
  `1` where a structural break was actually applied (per
  `break_fraction`), `0` elsewhere – ground truth for evaluating a
  change-point detection method the same way `true_slope` already lets
  you evaluate a trend test, via
  [`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md).
  Always present, even when `break_type = "none"` (all `0` in that
  case), so downstream code doesn't need to branch on whether the field
  exists.

- break_time:

  The actual time step used as the break point (see `break_time` above),
  or `NULL` if `break_type = "none"` or `break_fraction = 0` (no cell
  actually got a break).

- parameters:

  The data-generating parameters retained with the realisation for
  reproducible benchmarking.

- diagnostics:

  Simulation metadata, including the covariance model and target
  unit-lag correlation where applicable.

## Details

**Function type:** **Benchmarking function** – generates known-truth
data for
[`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md)
and external methods. It is not a component of the inferential workflows
themselves.

## Typical use

Single-run benchmarking:

    sim_trend_stack()
        |
    run one or more methods on sim$series
    (trend_test(), workflow_tst(), workflow_rta(), or your own)
        |
    compare_detections()

For replicated benchmarking, repeat simulation and detection across
seeds, collect the results, and aggregate them:

    list of detections + list of ground truths
        |
    compare_detections(replicates = TRUE)

See
[`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md)'s
examples for both routes worked through in full. One call creates one
replicate; Monte Carlo repetition remains outside this function so the
same simulator can benchmark sptrends methods, future procedures, or
external software. Use `ar1` for temporal dependence and a formal
`spatial_model` with its covariance parameters for spatial dependence.

## Methodological details

**Scope and limitations**

This function is intended for methodological benchmarking rather than
realistic environmental simulation. Its objective is to generate
controlled datasets with a known ground truth, not to reproduce the
statistical properties of a particular environmental variable such as
NDVI, temperature or precipitation. Formal Gaussian, exponential and
Matérn covariance models and exact geometric signal regions provide
controlled temporal and spatial structures for evaluating detection
methods; they are not a calibrated environmental process model.

**How the trend is generated – trend field**

`trend_shape` sets the *pattern* of the true slope field:

- `"radial"` (default): slope is `trend_strength` at the centre cell and
  decays smoothly (exponentially) with distance from it, reaching close
  to zero towards the edges.

- `"gradient"`: slope varies linearly from `-trend_strength` at the left
  edge to `+trend_strength` at the right edge – decreasing on one side,
  increasing on the other, with no radial symmetry.

- `"block"`: `trend_strength` everywhere inside a centred square block
  covering about half the raster's area, `0` everywhere outside it – a
  sharp-edged region of trend against a flat background, rather than a
  smooth spatial gradient.

**How the trend is generated – trend masking**

`trend_fraction` then decides *how much of the raster actually has a
trend at all*: the raster is partitioned into coarse, contiguous spatial
blocks (not individual cells), a random `trend_fraction` proportion of
*blocks* keep the slope values `trend_shape` assigned them, and the rest
are forced to an exact true slope of `0`, regardless of shape. Acting on
whole spatial patches rather than scattering individual cells at random
matters here:
[`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md)'s
whole rationale is borrowing statistical strength from a cell's
neighbours, which only helps when neighbouring cells are plausibly
trending together – a cell-by-cell (salt-and-pepper) random mask would
quietly defeat that assumption and penalise CMK for no good reason.
`trend_fraction = 0` gives a complete null field (every cell's true
slope is exactly `0` – there is no trend anywhere to find, the sharpest
way to check a method's false-positive rate). `trend_fraction = 1` gives
a trend everywhere `trend_shape` says there should be one. Among the
blocks that do keep a trend, a random 10% have their sign flipped as a
whole block, so the field is not purely "everything increases" even
within the trending region, while still keeping each flipped patch
spatially coherent.

**How spatial autocorrelation is generated – the idea**

For formal benchmarking, choose `spatial_model = "gaussian"`,
`"exponential"`, or `"matern"`. Independent stationary fields are drawn
at each time step by circulant embedding and FFT. If an admissible
embedding is unavailable, modest grids use an exact covariance
eigendecomposition. For Gaussian and exponential models, `spatial_rho`
is the theoretical correlation between horizontally or vertically
adjacent cells. For Matérn fields, `spatial_range` and
`spatial_smoothness` define the covariance. The `"independent"` model
provides the exact zero-dependence baseline.

`spatial_model = "legacy"` retains the earlier didactic smoother.
Independent noise, with standard deviation `noise_sd`, is drawn per cell
and per time step (with AR(1) correlation across time, controlled by
`ar1`, but no spatial structure yet). If `smooth_radius > 0` and
`spatial_rho > 0`, each time step's noise layer is smoothed with a focal
mean filter. `smooth_radius` controls spatial scale, whereas
`spatial_rho` controls how strongly the original and smoothed fields are
blended. This is a pragmatic moving-average smoother, not a fitted
CAR/SAR model or a Gaussian random field with a specified covariance
function – it is meant to give
[`spatial_autocorrelation()`](https://olive-r.github.io/sptrends/reference/spatial_autocorrelation.md)
something genuine to detect in examples, not to match any particular
real-world spatial process (see "What this simulator is not" above).

**Computational considerations**

**How spatial autocorrelation is generated – implementation**

The smoothing filter is a square window of side `2 * smooth_radius + 1`
cells (via
[`terra::focal()`](https://rspatial.github.io/terra/reference/focal.html)).
The smoothed and blended noise are rescaled so spatial controls do not
intentionally change `noise_sd`. `spatial_rho = 0` retains the
independent field; `spatial_rho = 1` retains the fully smoothed field
and reproduces the behaviour used before this parameter was added.
Intermediate values are blend weights, not target Moran's I values. All
`n_time` layers are smoothed in a single batched
[`terra::focal()`](https://rspatial.github.io/terra/reference/focal.html)
call on the whole multi-layer noise stack, rather than one call per
layer –
[`terra::focal()`](https://rspatial.github.io/terra/reference/focal.html)
already smooths each layer of a multi-layer input independently (there
is no cross-layer mixing), so this produces identical output while
avoiding `n_time - 1` redundant round trips through `terra`'s internals.
Matters mainly for a large `n_time`; for the modest-sized rasters this
function is typically used to build (examples, tests), both versions are
effectively instant either way.
[`terra::focal()`](https://rspatial.github.io/terra/reference/focal.html)
itself additionally requires the focal window (side
`2 * smooth_radius + 1`) to be no more than twice the raster's own size
in each direction – a small raster with a proportionally large
`smooth_radius` (most sharply, any raster with `nrow` or `ncol` of 1,
since even the smallest window, `smooth_radius = 1`, already has side 3)
will not fit this. Rather than letting
[`terra::focal()`](https://rspatial.github.io/terra/reference/focal.html)'s
own internal error propagate up, this is checked beforehand and, if
violated, a warning is issued and unsmoothed noise is used for that call
instead.

**Statistical assumptions and interpretation**

**Noise distribution and comparing methods**

`noise_dist` controls the *shape* of the noise, independently of its
standard deviation (`noise_sd`) or its temporal (`ar1`) and spatial
(`smooth_radius`) correlation. This matters specifically for comparing a
rank-based method (Mann-Kendall, Contextual Mann-Kendall, Theil-Sen)
against a parametric one (e.g. ordinary least squares): under Gaussian
noise (`noise_dist = "gaussian"`, the default), OLS is the more
*efficient* estimator, so a comparison run only under Gaussian noise
will not show the robustness advantage rank-based methods are typically
chosen for. Set `noise_dist = "t"` with a low `t_df` (e.g. `3` or `4`)
to generate heavy-tailed, outlier-prone noise instead, at the same
`noise_sd` – this is the condition under which rank-based methods are
expected to hold up better than OLS.

The simulation parameters control a synthetic data-generating process;
they are not fitted environmental parameters. Under
`spatial_model = "legacy"`, `spatial_rho` is a blend weight rather than
a target correlation. Under Gaussian or exponential covariance with
Gaussian noise, it is the target correlation between horizontally or
vertically adjacent cells. With Student-t noise, it controls the latent
Gaussian dependence used by the copula, so realised Pearson correlation
need not equal it exactly. One call produces one realisation; use
[`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md)
for a Monte Carlo experiment.

**Quality assurance**

Tests verify known slope and break fields, reproducibility, null and
complete-signal cases, noise scale and distribution, degenerate grid
sizes, focal-window safeguards, temporal/spatial controls and direct
compatibility with
[`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md).
The `spatial_rho` tests specifically protect both endpoint semantics and
the pre-0.89 default. Independent validation additionally compared the
analytical Matérn correlation with
[`fields::Matern()`](https://rdrr.io/pkg/fields/man/Exponential.html),
evaluated 1,000 Gaussian fields per spatial model, and checked temporal
AR(1), marginal distributions, exact truth fields, detection metrics and
reproducibility. All 33 external validation controls passed in the
recorded 0.96.3 full run. The retained script and numerical summaries
are described under `inst/validation`. See
[`?sptrends`](https://olive-r.github.io/sptrends/reference/sptrends-package.md)
for the package-wide release-check protocol.

## References

Dietrich, C.R. and Newsam, G.N. (1997). Fast and exact simulation of
stationary Gaussian processes through circulant embedding of the
covariance matrix. *SIAM Journal on Scientific Computing*, 18(4),
1088-1107.
[doi:10.1137/S1064827592240555](https://doi.org/10.1137/S1064827592240555)

## See also

Other example data functions:
[`example_data()`](https://olive-r.github.io/sptrends/reference/example_data.md)

## Examples

``` r
# A small synthetic dataset with a known trend -- 8 time steps on a
# 12x12 grid. terra::global(..., "range") shows the true slope's
# minimum and maximum across all cells.
sim <- sim_trend_stack(nrow = 12, ncol = 12, n_time = 8, seed = 42)
#> >> [sim_trend_stack()] elapsed: 0.03 s
terra::nlyr(sim$series)
#> [1] 8
terra::global(sim$true_slope, "range", na.rm = TRUE)
#>                    min       max
#> true_slope -0.06043355 0.1250629

# A complete null field: every true slope is exactly zero.
sim_null <- sim_trend_stack(nrow = 12, ncol = 12, n_time = 8,
                             trend_fraction = 0, seed = 1)
#> >> [sim_trend_stack()] elapsed: 0.03 s
terra::global(sim_null$true_slope, "range", na.rm = TRUE)
#>            min max
#> true_slope   0   0

# Ground truth for a change-point method (e.g. Pettitt's test): a
# mean-shift break in 30% of the map, no monotonic trend at all.
sim_break <- sim_trend_stack(nrow = 12, ncol = 12, n_time = 10,
                              trend_fraction = 0, break_type = "mean",
                              break_fraction = 0.3, seed = 2)
#> >> [sim_trend_stack()] elapsed: 0.04 s
sim_break$break_time
#> [1] 5
terra::global(sim_break$true_break, "sum", na.rm = TRUE)
#>            sum
#> true_break  44

# Spatial scale and intensity are separate: smooth_radius defines the
# focal scale, while spatial_rho blends independent and smoothed noise.
# Compare Moran's I for one time step at the two intensity extremes.
# \donttest{
r0 <- sim_trend_stack(nrow = 20, ncol = 20, n_time = 1,
                       smooth_radius = 3, spatial_rho = 0,
                       seed = 1)$series[[1]]
#> >> [sim_trend_stack()] elapsed: 0.01 s
r1 <- sim_trend_stack(nrow = 20, ncol = 20, n_time = 1,
                       smooth_radius = 3, spatial_rho = 1,
                       seed = 1)$series[[1]]
#> >> [sim_trend_stack()] elapsed: 0.02 s
spatial_autocorrelation(r0, nperm = 99, seed = 1, verbose = FALSE,
                         report = FALSE)$statistic
#> [1] 0.07156812
spatial_autocorrelation(r1, nperm = 99, seed = 1, verbose = FALSE,
                         report = FALSE)$statistic
#> [1] 0.6934035
# }

# Heavy-tailed (outlier-prone) noise instead of Gaussian, at the same
# standard deviation -- for comparing a rank-based method against OLS.
sim_heavy <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 10,
                              noise_dist = "t", t_df = 3, seed = 1)
#> >> [sim_trend_stack()] elapsed: 0.04 s
terra::nlyr(sim_heavy$series)
#> [1] 10
```
