#' Generate a synthetic gridded time series with known true trends
#'
#' Generates benchmark datasets with a known ground truth, for
#' evaluating spatiotemporal trend-detection methods. Concretely: a
#' synthetic gridded time series, with the true, known slope at every
#' cell kept alongside the data itself, for runnable `@examples`,
#' vignettes, and unit tests, rather than as a realistic environmental
#' simulation. Because the ground truth is known exactly (something no
#' real dataset offers), this is the tool behind sptrends' own worked
#' examples of what a trend-detection method actually gets right or
#' wrong: compare `true_slope` against a fitted method's estimated
#' slope, and against which cells it calls significant, to see
#' estimation error and Type I/Type II error directly instead of having
#' to take a method's output on faith -- the mechanism that lets power,
#' false discovery rate, sensitivity, specificity, and robustness be
#' demonstrated *objectively* for any method, not just described in
#' prose.
#'
#' **Function type:** **Benchmarking function** -- generates known-truth
#' data for [compare_detections()] and external methods. It is not a
#' component of the inferential workflows themselves.
#'
#' @section Typical use:
#' Single-run benchmarking:
#' ```
#' sim_trend_stack()
#'     |
#' run one or more methods on sim$series
#' (trend_test(), workflow_tst(), workflow_rta(), or your own)
#'     |
#' compare_detections()
#' ```
#' For replicated benchmarking, repeat simulation and detection across
#' seeds, collect the results, and aggregate them:
#' ```
#' list of detections + list of ground truths
#'     |
#' compare_detections(replicates = TRUE)
#' ```
#' See [compare_detections()]'s examples for both routes worked through
#' in full. One call creates one replicate; Monte Carlo repetition
#' remains outside this function so the same simulator can benchmark
#' sptrends methods, future procedures, or external software. Use `ar1` for
#' temporal dependence and a formal `spatial_model` with its covariance
#' parameters for spatial dependence.
#'
#' @section Methodological details:
#' **Scope and limitations**
#'
#' This function is intended for methodological benchmarking rather than
#' realistic environmental simulation. Its objective is to generate
#' controlled datasets with a known ground truth, not to reproduce the
#' statistical properties of a particular environmental variable such
#' as NDVI, temperature or precipitation. Formal Gaussian, exponential and
#' Matérn covariance models and exact geometric signal regions provide
#' controlled temporal and spatial structures for evaluating detection
#' methods; they are not a calibrated environmental process model.
#'
#' **How the trend is generated -- trend field**
#'
#' `trend_shape` sets the *pattern* of the true slope field:
#' - `"radial"` (default): slope is `trend_strength` at the centre cell
#'   and decays smoothly (exponentially) with distance from it, reaching
#'   close to zero towards the edges.
#' - `"gradient"`: slope varies linearly from `-trend_strength` at the
#'   left edge to `+trend_strength` at the right edge -- decreasing on
#'   one side, increasing on the other, with no radial symmetry.
#' - `"block"`: `trend_strength` everywhere inside a centred square block
#'   covering about half the raster's area, `0` everywhere outside it --
#'   a sharp-edged region of trend against a flat background, rather than
#'   a smooth spatial gradient.
#'
#' **How the trend is generated -- trend masking**
#'
#' `trend_fraction` then decides *how much of the raster actually has a
#' trend at all*: the raster is partitioned into coarse, contiguous
#' spatial blocks (not individual cells), a random `trend_fraction`
#' proportion of *blocks* keep the slope values `trend_shape` assigned
#' them, and the rest are forced to an exact true slope of `0`, regardless
#' of shape. Acting on whole spatial patches rather than scattering
#' individual cells at random matters here: [trend_test()]'s
#' whole rationale is borrowing statistical strength from a cell's
#' neighbours, which only helps when neighbouring cells are plausibly
#' trending together -- a cell-by-cell (salt-and-pepper) random mask would
#' quietly defeat that assumption and penalise CMK for no good reason.
#' `trend_fraction = 0` gives a complete null field (every cell's true
#' slope is exactly `0` -- there is no trend anywhere to find, the
#' sharpest way to check a method's false-positive rate). `trend_fraction
#' = 1` gives a trend everywhere `trend_shape` says there should be one.
#' Among the blocks that do keep a trend, a random 10% have their sign
#' flipped as a whole block, so the field is not purely "everything
#' increases" even within the trending region, while still keeping each
#' flipped patch spatially coherent.
#'
#' **How spatial autocorrelation is generated -- the idea**
#'
#' For formal benchmarking, choose `spatial_model = "gaussian"`,
#' `"exponential"`, or `"matern"`. Independent stationary fields are drawn
#' at each time step by circulant embedding and FFT. If an admissible
#' embedding is unavailable, modest grids use an exact covariance
#' eigendecomposition. For Gaussian and
#' exponential models, `spatial_rho` is the theoretical correlation between
#' horizontally or vertically adjacent cells. For Matérn fields,
#' `spatial_range` and `spatial_smoothness` define the covariance. The
#' `"independent"` model provides the exact zero-dependence baseline.
#'
#' `spatial_model = "legacy"` retains the earlier didactic smoother.
#' Independent noise, with standard deviation `noise_sd`, is drawn per
#' cell and per time step (with AR(1) correlation across time, controlled
#' by `ar1`, but no spatial structure yet). If `smooth_radius > 0` and
#' `spatial_rho > 0`, each time step's noise layer is smoothed with a
#' focal mean filter. `smooth_radius` controls spatial scale, whereas
#' `spatial_rho` controls how strongly the original and smoothed fields
#' are blended. This is a
#' pragmatic moving-average smoother, not a fitted CAR/SAR model or a
#' Gaussian random field with a specified covariance function -- it is
#' meant to give [spatial_autocorrelation()] something genuine to
#' detect in examples, not to match any particular real-world spatial
#' process (see "What this simulator is not" above).
#'
#' **Computational considerations**
#'
#' **How spatial autocorrelation is generated -- implementation**
#'
#' The smoothing filter is a square window of side `2 * smooth_radius +
#' 1` cells (via `terra::focal()`). The smoothed and blended noise are
#' rescaled so spatial controls do not intentionally change `noise_sd`.
#' `spatial_rho = 0` retains the independent field; `spatial_rho = 1`
#' retains the fully smoothed field and reproduces the behaviour used
#' before this parameter was added. Intermediate values are blend
#' weights, not target Moran's I values. All `n_time` layers are
#' smoothed in a single batched
#' `terra::focal()` call on the whole multi-layer noise stack, rather
#' than one call per layer -- `terra::focal()` already smooths each
#' layer of a multi-layer input independently (there is no cross-layer
#' mixing), so this produces identical output while avoiding
#' `n_time - 1` redundant round trips through `terra`'s internals.
#' Matters mainly for a large `n_time`; for the modest-sized rasters
#' this function is typically used to build (examples, tests), both
#' versions are effectively instant either way. `terra::focal()`
#' itself additionally requires the focal window (side
#' `2 * smooth_radius + 1`) to be no more than twice the raster's own
#' size in each direction -- a small raster with a proportionally
#' large `smooth_radius` (most sharply, any raster with `nrow` or
#' `ncol` of 1, since even the smallest window, `smooth_radius = 1`,
#' already has side 3) will not fit this. Rather than letting
#' `terra::focal()`'s own internal error propagate up, this is checked
#' beforehand and, if violated, a warning is issued and unsmoothed
#' noise is used for that call instead.
#'
#' **Statistical assumptions and interpretation**
#'
#' **Noise distribution and comparing methods**
#'
#' `noise_dist` controls the *shape* of the noise, independently of its
#' standard deviation (`noise_sd`) or its temporal (`ar1`) and spatial
#' (`smooth_radius`) correlation. This matters specifically for comparing
#' a rank-based method (Mann-Kendall, Contextual Mann-Kendall, Theil-Sen)
#' against a parametric one (e.g. ordinary least squares): under Gaussian
#' noise (`noise_dist = "gaussian"`, the default), OLS is the more
#' *efficient* estimator, so a comparison run only under Gaussian noise
#' will not show the robustness advantage rank-based methods are
#' typically chosen for. Set `noise_dist = "t"` with a low `t_df` (e.g.
#' `3` or `4`) to generate heavy-tailed, outlier-prone noise instead, at
#' the same `noise_sd` -- this is the condition under which rank-based
#' methods are expected to hold up better than OLS.
#'
#' The simulation parameters control a synthetic data-generating process;
#' they are not fitted environmental parameters. Under `spatial_model =
#' "legacy"`, `spatial_rho` is a blend weight rather than a target
#' correlation. Under Gaussian or exponential covariance with Gaussian
#' noise, it is the target correlation between horizontally or vertically
#' adjacent cells. With Student-t noise, it controls the latent Gaussian
#' dependence used by the copula, so realised Pearson correlation need not
#' equal it exactly. One call produces one realisation; use
#' [benchmark_methods()] for a Monte Carlo experiment.
#'
#' **Quality assurance**
#'
#' Tests verify known slope and break fields, reproducibility, null and
#' complete-signal cases, noise scale and distribution, degenerate grid
#' sizes, focal-window safeguards, temporal/spatial controls and direct
#' compatibility with [compare_detections()]. The `spatial_rho` tests
#' specifically protect both endpoint semantics and the pre-0.89 default.
#' Independent validation additionally compared the analytical Matérn
#' correlation with `fields::Matern()`, evaluated 1,000 Gaussian fields per
#' spatial model, and checked temporal AR(1), marginal distributions, exact
#' truth fields, detection metrics and reproducibility. All 33 external
#' validation controls passed in the recorded 0.96.3 full run. The retained
#' script and numerical summaries are described under `inst/validation`.
#' See `?sptrends` for the package-wide release-check protocol.
#'
#' @param nrow,ncol Integer. Raster dimensions in cells.
#' @param n_time Integer. Number of time steps (layers).
#' @param trend_strength Numeric. Maximum magnitude of the true slope (see
#'   "How the trend is generated -- trend field" above for how it
#'   varies spatially according to `trend_shape`).
#' @param trend_shape Spatial pattern of the true slope field: `"radial"`,
#'   `"gradient"`, `"block"`, `"square"`, `"rectangle"`, `"ellipse"`,
#'   or `"custom"`; see "How the trend is generated -- trend field" above.
#' @param trend_fraction Numeric in `[0, 1]`. Proportion of cells that
#'   keep a non-zero true slope; the rest are forced to exactly `0`.
#'   Default `0.9`. Set to `0` for a complete null field (no trend
#'   anywhere), or `1` for a trend at every cell `trend_shape` assigns one
#'   to.
#' @param ar1 Numeric in `(-1, 1)`. AR(1) autocorrelation coefficient of the
#'   noise term across time, applied independently within each cell.
#' @param noise_sd Numeric `> 0`. Standard deviation of the noise term --
#'   the main control over how easy or hard the true trend is to detect
#'   against the background noise (signal-to-noise ratio). Default `1`.
#' @param noise_dist `"gaussian"` (default) or `"t"` -- the distribution
#'   of the noise term, always rescaled to the standard deviation
#'   `noise_sd` regardless of shape. This matters for comparing methods:
#'   Mann-Kendall-based tests (and Theil-Sen) are rank-based and robust to
#'   heavy-tailed, outlier-prone noise, while a parametric method like OLS
#'   is the more *efficient* choice specifically when noise is Gaussian --
#'   the classic case for preferring a non-parametric method only shows up
#'   with non-Gaussian noise. `"t"` draws from a Student's t distribution
#'   with `t_df` degrees of freedom instead, giving heavier tails (more
#'   extreme outliers) at the same standard deviation.
#' @param t_df Numeric `> 2`. Degrees of freedom for the noise
#'   distribution when `noise_dist = "t"`; ignored otherwise. Lower values
#'   give heavier tails (more frequent/extreme outliers); values just
#'   above `2` are extreme, values above `30` are practically
#'   indistinguishable from Gaussian. Must exceed `2` for the
#'   distribution to have finite variance (required to rescale it to
#'   `noise_sd`).
#' @param smooth_radius Integer `>= 0`. Radius, in cells, of the focal mean
#'   window used to introduce spatial autocorrelation into the noise (see
#'   "How spatial autocorrelation is generated" above). `0` disables
#'   this: cells get independent noise with no spatial structure beyond
#'   the trend field itself. `1` (default) smooths over a 3x3 window;
#'   larger values give smoother, more strongly autocorrelated fields, at
#'   the cost of more computation.
#' @param spatial_rho Numeric in `[0, 1]`. For the legacy model, the weight
#'   assigned to focal-smoothed noise; this is not a target Moran's I. For
#'   Gaussian and exponential covariance models, the adjacent-cell
#'   correlation of the latent Gaussian field. It must be positive for these
#'   formal models. It is ignored by the independent and Matérn models.
#' @param spatial_model Spatial noise model. `"legacy"` preserves the focal
#'   smoother used by earlier versions. `"independent"`, `"gaussian"`,
#'   `"exponential"`, and `"matern"` use a formal stationary covariance
#'   model simulated by circulant embedding and FFT.
#' @param spatial_range Positive covariance scale for the Matérn model.
#' @param spatial_smoothness Positive Matérn smoothness parameter.
#' @param signal_size One number or a two-number vector giving the height and
#'   width, in cells, of block, square, rectangular, or elliptical signal
#'   regions. Defaults to half the raster in each dimension for `"block"`
#'   (its historical default), or a third of the raster for `"square"`,
#'   `"rectangle"` and `"ellipse"`.
#' @param signal_location `"centre"` or `"random"`; location of an exact
#'   geometric signal region (`"block"`, `"square"`, `"rectangle"`,
#'   `"ellipse"` or `"custom"`). Has no effect on `"radial"` or `"gradient"`,
#'   which are not exact geometric regions.
#' @param signal_angle Rotation of a geometric signal region in degrees, or
#'   `"random"` to draw a new orientation in each realisation.
#' @param signal_axis_ratio Positive width multiplier for an ellipse.
#' @param custom_mask Logical/numeric matrix, vector, or single-layer
#'   `SpatRaster` defining the signal when `trend_shape = "custom"`.
#' @param constant_block Logical. If `TRUE`, sets a small block of cells to
#'   a constant value (zero temporal variance) to exercise degenerate-case
#'   handling. Independent of `trend_shape = "block"` above, which is
#'   about the *trend* pattern, not a zero-variance edge case.
#' @param break_type `"none"` (default): no structural break, this
#'   function's original behaviour, unchanged. `"mean"`: a step change
#'   in level at `break_time` -- ground truth for a change-point method
#'   like Pettitt's test. `"slope"`: a change in slope at `break_time`
#'   (continuous in level, a kink not a jump) -- ground truth for
#'   distinguishing "did the trend change partway through" from "is
#'   there a trend at all". A break and a monotonic trend
#'   (`trend_strength`/`trend_fraction` above) compose independently on
#'   any given cell -- a cell can have a trend only, a break only, both,
#'   or neither; the two are not mutually exclusive, nor does one
#'   replace the other.
#' @param break_time Integer or `NULL`. Time step *after* which the
#'   break takes effect (so `break_time = 5` means steps 1-5 are
#'   unaffected, the break shows in step 6 onward). `NULL` (default):
#'   the middle of the series, `round(n_time / 2)`. Must be in
#'   `[1, n_time - 1]` so there is at least one time step on each side.
#'   Ignored when `break_type = "none"`.
#' @param break_fraction Fraction of the same spatially coherent blocks
#'   used for `trend_fraction` above that get a break, independently of
#'   which blocks (if any) got a trend. Default `0.3`. Ignored when
#'   `break_type = "none"`.
#' @param break_magnitude Numeric. For `break_type = "mean"`: the size
#'   of the step added to the level from `break_time` onward. For
#'   `break_type = "slope"`: the size of the *additional* slope (per
#'   time step) that kicks in from `break_time` onward, on top of
#'   whatever slope (possibly zero) that cell already had. Default `2`.
#'   Ignored when `break_type = "none"`.
#' @param seed Integer or `NULL`. Random seed for reproducibility.
#' @param verbose Logical. If `TRUE`, report simulation progress, elapsed
#'   duration and estimated time remaining. Set to `FALSE` inside large
#'   Monte Carlo experiments when [benchmark_methods()] provides the outer
#'   progress display.
#'
#' @return An object of class `"sptrends_simulation"` and `"sptrends"`,
#'   supporting `print()`, `summary()`, and `plot()`, with:
#'   \item{series}{A `terra::SpatRaster` with `n_time` layers named
#'     `"t1"`, `"t2"`, ...; each layer represents one time step. Pass this
#'     to [workflow_tst()] or the individual trend functions.}
#'   \item{true_slope}{A single-layer `terra::SpatRaster`: the exact true
#'     slope used to generate each cell's series (`0` where
#'     `trend_fraction` assigned no trend). Compare directly against a
#'     fitted [slope_estimator()] to see estimation error, or against
#'     [direction_map()] to see which cells a method calls
#'     significant vs. which cells truly have a trend.}
#'   \item{true_signal}{A binary `SpatRaster` equal to `1` where the exact
#'     true slope is non-zero.}
#'   \item{true_direction}{A `SpatRaster` containing `-1`, `0`, or `1` for
#'     the known trend direction.}
#'   \item{true_break}{A single-layer `terra::SpatRaster`, `1` where a
#'     structural break was actually applied (per `break_fraction`),
#'     `0` elsewhere -- ground truth for evaluating a change-point
#'     detection method the same way `true_slope` already lets you
#'     evaluate a trend test, via [compare_detections()]. Always
#'     present, even when `break_type = "none"` (all `0` in that
#'     case), so downstream code doesn't need to branch on whether the
#'     field exists.}
#'   \item{break_time}{The actual time step used as the break point (see
#'     `break_time` above), or `NULL` if `break_type = "none"` or
#'     `break_fraction = 0` (no cell actually got a break).}
#'   \item{parameters}{The data-generating parameters retained with the
#'     realisation for reproducible benchmarking.}
#'   \item{diagnostics}{Simulation metadata, including the covariance model
#'     and target unit-lag correlation where applicable.}
#'
#' @examples
#' # A small synthetic dataset with a known trend -- 8 time steps on a
#' # 12x12 grid. terra::global(..., "range") shows the true slope's
#' # minimum and maximum across all cells.
#' sim <- sim_trend_stack(nrow = 12, ncol = 12, n_time = 8, seed = 42)
#' terra::nlyr(sim$series)
#' terra::global(sim$true_slope, "range", na.rm = TRUE)
#'
#' # A complete null field: every true slope is exactly zero.
#' sim_null <- sim_trend_stack(nrow = 12, ncol = 12, n_time = 8,
#'                              trend_fraction = 0, seed = 1)
#' terra::global(sim_null$true_slope, "range", na.rm = TRUE)
#'
#' # Ground truth for a change-point method (e.g. Pettitt's test): a
#' # mean-shift break in 30% of the map, no monotonic trend at all.
#' sim_break <- sim_trend_stack(nrow = 12, ncol = 12, n_time = 10,
#'                               trend_fraction = 0, break_type = "mean",
#'                               break_fraction = 0.3, seed = 2)
#' sim_break$break_time
#' terra::global(sim_break$true_break, "sum", na.rm = TRUE)
#'
#' # Spatial scale and intensity are separate: smooth_radius defines the
#' # focal scale, while spatial_rho blends independent and smoothed noise.
#' # Compare Moran's I for one time step at the two intensity extremes.
#' \donttest{
#' r0 <- sim_trend_stack(nrow = 20, ncol = 20, n_time = 1,
#'                        smooth_radius = 3, spatial_rho = 0,
#'                        seed = 1)$series[[1]]
#' r1 <- sim_trend_stack(nrow = 20, ncol = 20, n_time = 1,
#'                        smooth_radius = 3, spatial_rho = 1,
#'                        seed = 1)$series[[1]]
#' spatial_autocorrelation(r0, nperm = 99, seed = 1, verbose = FALSE,
#'                          report = FALSE)$statistic
#' spatial_autocorrelation(r1, nperm = 99, seed = 1, verbose = FALSE,
#'                          report = FALSE)$statistic
#' }
#'
#' # Heavy-tailed (outlier-prone) noise instead of Gaussian, at the same
#' # standard deviation -- for comparing a rank-based method against OLS.
#' sim_heavy <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 10,
#'                               noise_dist = "t", t_df = 3, seed = 1)
#' terra::nlyr(sim_heavy$series)
#'
#' @references
#' Dietrich, C.R. and Newsam, G.N. (1997). Fast and exact simulation of
#' stationary Gaussian processes through circulant embedding of the
#' covariance matrix. *SIAM Journal on Scientific Computing*, 18(4),
#' 1088-1107. \doi{10.1137/S1064827592240555}
#'
#' @family example data functions
#' @export
sim_trend_stack <- function(nrow = 15, ncol = 15, n_time = 10,
                             trend_strength = 0.15,
                             trend_shape = c("radial", "gradient", "block",
                                             "square", "rectangle", "ellipse",
                                             "custom"),
                             trend_fraction = 0.9, ar1 = 0.3, noise_sd = 1,
                             noise_dist = c("gaussian", "t"), t_df = 4,
                             smooth_radius = 1L, spatial_rho = 1,
                             spatial_model = c("legacy", "independent",
                                               "gaussian", "exponential",
                                               "matern"),
                             spatial_range = 1,
                             spatial_smoothness = 0.5,
                             signal_size = NULL,
                             signal_location = c("centre", "random"),
                             signal_angle = 0, signal_axis_ratio = 1,
                             custom_mask = NULL,
                             constant_block = TRUE,
                             break_type = c("none", "mean", "slope"),
                             break_time = NULL, break_fraction = 0.3,
                             break_magnitude = 2,
                             seed = NULL, verbose = TRUE) {
  finish_timer <- .sptrends_elapsed_timer("sim_trend_stack()", verbose)
  on.exit(finish_timer(), add = TRUE)
  trend_shape <- match.arg(trend_shape)
  noise_dist <- match.arg(noise_dist)
  break_type <- match.arg(break_type)
  spatial_model <- match.arg(spatial_model)
  signal_location <- match.arg(signal_location)
  if (!is.numeric(smooth_radius) || length(smooth_radius) != 1L ||
      is.na(smooth_radius) || !is.finite(smooth_radius) ||
      smooth_radius < 0 || smooth_radius != floor(smooth_radius)) {
    stop("'smooth_radius' must be one non-negative integer.")
  }
  smooth_radius <- as.integer(smooth_radius)
  if (trend_fraction < 0 || trend_fraction > 1) {
    stop("'trend_fraction' must be in [0, 1].")
  }
  if (break_fraction < 0 || break_fraction > 1) {
    stop("'break_fraction' must be in [0, 1].")
  }
  if (!is.numeric(spatial_rho) || length(spatial_rho) != 1L ||
      is.na(spatial_rho) || !is.finite(spatial_rho) ||
      spatial_rho < 0 || spatial_rho > 1) {
    stop("'spatial_rho' must be one finite numeric value in [0, 1].")
  }
  if (spatial_model %in% c("gaussian", "exponential") &&
      spatial_rho <= 0) {
    stop("For a formal spatial model, 'spatial_rho' must be in (0, 1].")
  }
  if (!is.numeric(spatial_range) || length(spatial_range) != 1L ||
      !is.finite(spatial_range) || spatial_range <= 0) {
    stop("'spatial_range' must be one positive finite number.")
  }
  if (!is.numeric(spatial_smoothness) ||
      length(spatial_smoothness) != 1L ||
      !is.finite(spatial_smoothness) || spatial_smoothness <= 0) {
    stop("'spatial_smoothness' must be one positive finite number.")
  }
  if (!is.null(seed)) set.seed(seed)
  if (is.character(signal_angle)) {
    if (length(signal_angle) != 1L || signal_angle != "random") {
      stop("Character 'signal_angle' must be 'random'.")
    }
    signal_angle <- stats::runif(1, 0, 180)
  }
  if (!is.numeric(signal_angle) || length(signal_angle) != 1L ||
      !is.finite(signal_angle)) {
    stop("'signal_angle' must be one finite number or 'random'.")
  }
  if (!is.numeric(signal_axis_ratio) || length(signal_axis_ratio) != 1L ||
      !is.finite(signal_axis_ratio) || signal_axis_ratio <= 0) {
    stop("'signal_axis_ratio' must be one positive finite number.")
  }
  if (!is.null(signal_size) &&
      (!is.numeric(signal_size) || !length(signal_size) %in% c(1L, 2L) ||
       any(!is.finite(signal_size)) || any(signal_size <= 0))) {
    stop("'signal_size' must contain one or two positive finite numbers.")
  }

  ncell_ <- nrow * ncol
  rows <- rep(seq_len(nrow), each = ncol)
  cols <- rep(seq_len(ncol), times = nrow)

  slope <- switch(trend_shape,
    radial = {
      centre_r <- (nrow + 1) / 2
      centre_c <- (ncol + 1) / 2
      dist_centre <- sqrt((rows - centre_r)^2 + (cols - centre_c)^2)
      max_dist <- max(dist_centre)
      if (max_dist == 0) {
        # A 1x1 grid: every cell (there is exactly one) IS the centre,
        # so there is no "distance from the centre" to decay over --
        # the un-guarded formula below would divide by zero here.
        rep(trend_strength, ncell_)
      } else {
        trend_strength * exp(-dist_centre / (max_dist / 2))
      }
    },
    gradient = {
      centre_c <- (ncol + 1) / 2
      half_width <- (ncol - 1) / 2
      if (half_width == 0) {
        # ncol = 1: a single column has no left-to-right span to form a
        # gradient across -- same degenerate case as "radial" above.
        rep(trend_strength, ncell_)
      } else {
        trend_strength * (cols - centre_c) / half_width
      }
    },
    block = rep(trend_strength, ncell_),
    square = rep(trend_strength, ncell_),
    rectangle = rep(trend_strength, ncell_),
    ellipse = rep(trend_strength, ncell_),
    custom = rep(trend_strength, ncell_)
  )

  if (trend_shape %in% c("block", "square", "rectangle", "ellipse", "custom")) {
    if (is.null(signal_size)) {
      # "block" keeps its own historical default (half the raster in each
      # dimension) so existing calls that never touched signal_size are
      # unaffected by this fix -- only signal_location/signal_angle/
      # signal_size, previously silently ignored for "block", now work.
      signal_size <- if (trend_shape == "block") {
        c(max(1, round(nrow * 0.5)), max(1, round(ncol * 0.5)))
      } else {
        c(max(1, round(nrow / 3)), max(1, round(ncol / 3)))
      }
    }
    exact_signal <- .simulation_signal_mask(
      nrow, ncol, trend_shape, signal_size, signal_location,
      signal_angle, signal_axis_ratio, custom_mask)
    slope[!exact_signal] <- 0
  }

  # trend_fraction and the sign flip both act on spatially coherent PATCHES
  # (coarse contiguous blocks of cells), not on individual cells picked at
  # random -- CMK's whole rationale is borrowing strength from a cell's
  # neighbours, which only works if a "trending" cell's neighbours are
  # plausibly trending too. A cell-by-cell (salt-and-pepper) random mask
  # would silently defeat that assumption and unfairly penalise CMK.
  block_size <- max(2L, round(min(nrow, ncol) / 5))
  block_row <- ((rows - 1) %/% block_size)
  block_col <- ((cols - 1) %/% block_size)
  block_id <- block_row * (((ncol - 1) %/% block_size) + 1) + block_col
  unique_blocks <- unique(block_id)

  if (trend_fraction < 1) {
    n_active <- round(trend_fraction * length(unique_blocks))
    active_blocks <- sample(unique_blocks, n_active)
    slope[!(block_id %in% active_blocks)] <- 0
  }
  exact_shape <- trend_shape %in%
    c("square", "rectangle", "ellipse", "custom")
  if (any(slope != 0) && !exact_shape) {
    nonzero_blocks <- unique(block_id[slope != 0])
    n_flip <- floor(0.1 * length(nonzero_blocks))
    if (n_flip > 0) {
      flip_blocks <- sample(nonzero_blocks, n_flip)
      slope[block_id %in% flip_blocks] <- -slope[block_id %in% flip_blocks]
    }
  }

  # Break-affected cells are a SEPARATE, independent block draw from the
  # trend-affected ones above -- a cell can have a trend only, a break
  # only, both, or neither (confirmed as the intended design: these are
  # two independent phenomena, not variants of the same one). Same
  # spatially-coherent block reasoning as trend_fraction applies here:
  # a genuine structural break should affect a contiguous patch of
  # cells, not a cell-by-cell salt-and-pepper pattern, or a spatial
  # change-detection method's own neighbourhood-borrowing assumption
  # (mirroring CMK's) would be unfairly penalised the same way.
  has_break <- rep(FALSE, ncell_)
  break_time_used <- NULL
  if (break_type != "none" && break_fraction > 0) {
    break_time_used <- if (is.null(break_time)) {
      round(n_time / 2)
    } else {
      break_time
    }
    if (break_time_used < 1 || break_time_used >= n_time) {
      stop("'break_time' must be in [1, n_time - 1] so there is at ",
           "least one time step on each side of the break.")
    }
    n_break_active <- round(break_fraction * length(unique_blocks))
    if (n_break_active > 0) {
      break_blocks <- sample(unique_blocks, n_break_active)
      has_break <- block_id %in% break_blocks
    }
  }

  noise <- matrix(0, ncell_, n_time)
  progress <- .sptrends_progress(
    n_time, "Simulating temporal layers", verbose
  )
  on.exit(.sptrends_progress_close(progress), add = TRUE)
  if (spatial_model == "legacy") {
    noise[, 1] <- .sim_noise(ncell_, noise_sd, noise_dist, t_df)
    .sptrends_progress_step(progress, 1L)
    if (n_time > 1) {
      for (k in 2:n_time) {
        noise[, k] <- ar1 * noise[, k - 1] +
          .sim_noise(ncell_, noise_sd * sqrt(1 - ar1^2), noise_dist, t_df)
        .sptrends_progress_step(progress, k)
      }
    }
  } else {
    innovation_sd <- noise_sd * sqrt(1 - ar1^2)
    for (k in seq_len(n_time)) {
      current_sd <- if (k == 1L) noise_sd else innovation_sd
      innovation <- .simulate_spatial_innovation(
        nrow, ncol, spatial_model, spatial_rho, spatial_range,
        spatial_smoothness, current_sd, noise_dist, t_df)
      innovation <- as.vector(t(innovation))
      noise[, k] <- if (k == 1L) {
        innovation
      } else {
        ar1 * noise[, k - 1L] + innovation
      }
      .sptrends_progress_step(progress, k)
    }
  }
  .sptrends_progress_close(progress)

  r <- terra::rast(nrows = nrow, ncols = ncol, xmin = 0, xmax = ncol,
                    ymin = 0, ymax = nrow)

  if (spatial_model == "legacy" && isTRUE(smooth_radius > 0) &&
      spatial_rho > 0) {
    window_size <- 2L * as.integer(smooth_radius) + 1L
    # terra::focal() itself requires the focal window to be no more
    # than 2x the raster's own dimension in each direction -- a real
    # constraint of terra's own C++ implementation, not something this
    # package can smooth over (so to speak). A raster too small for
    # the requested smooth_radius (most sharply, a 1x1 grid, where
    # *any* smooth_radius >= 1 already violates this) has no genuine
    # neighbourhood to average over anyway, so smoothing is skipped
    # for that raster with a warning, rather than letting terra's own
    # much less informative internal error propagate up.
    if (window_size > 2L * nrow || window_size > 2L * ncol) {
      warning(sprintf(
        paste0("smooth_radius = %d skipped: the resulting %dx%d focal ",
               "window does not fit within this %dx%d raster (terra::",
               "focal() requires a window no more than twice the ",
               "raster's own size in each direction). Returning ",
               "unsmoothed noise instead."),
        smooth_radius, window_size, window_size, nrow, ncol))
    } else {
    # Batched across all n_time layers in a single terra::focal() call,
    # rather than looping setValues()/focal()/values() once per layer:
    # terra::focal() already accepts (and is designed for) multi-layer
    # input, smoothing each layer independently either way, so this
    # produces identical results while avoiding n_time - 1 redundant
    # R-to-C round trips. Matters for large n_time -- for the small
    # rasters this function is typically used to build (demos, tests),
    # both versions are effectively instant.
    orig_sd <- apply(noise, 2, stats::sd)
    layers_r <- terra::rast(r, nlyr = n_time)
    terra::values(layers_r) <- noise
    smoothed <- terra::focal(layers_r, w = window_size, fun = "mean",
                              na.rm = TRUE)
    smoothed_vals <- terra::values(smoothed, mat = TRUE)
    new_sd <- apply(smoothed_vals, 2, stats::sd)
    rescale <- ifelse(new_sd > 0, orig_sd / new_sd, 1)
    smoothed_scaled <- sweep(smoothed_vals, 2, rescale, `*`)
    if (spatial_rho == 1) {
      # Preserve the pre-0.89 default exactly, without an unnecessary
      # second floating-point rescaling of the fully smoothed field.
      noise <- smoothed_scaled
    } else {
      blended <- (1 - spatial_rho) * noise +
        spatial_rho * smoothed_scaled
      blended_sd <- apply(blended, 2, stats::sd)
      blend_rescale <- ifelse(blended_sd > 0, orig_sd / blended_sd, 1)
      noise <- sweep(blended, 2, blend_rescale, `*`)
    }
    }
  }

  t <- seq_len(n_time)
  X <- outer(slope, t) + noise + 10

  if (any(has_break)) {
    after_break <- t > break_time_used
    if (break_type == "mean") {
      # A step up in level after break_time_used -- ground truth for a
      # Pettitt-style mean-shift changepoint.
      X[has_break, after_break] <- X[has_break, after_break] +
        break_magnitude
    } else if (break_type == "slope") {
      # An additional slope kicks in from break_time_used onward, on
      # top of (and independent from) any existing trend slope for
      # that cell -- continuous in level at the break itself (no jump),
      # a kink in slope, not a step in level. offset*(t - break_time)
      # for t > break_time is the standard piecewise-linear structural
      # break parametrisation.
      offset <- outer(rep(break_magnitude, sum(has_break)),
                       t[after_break] - break_time_used)
      X[has_break, after_break] <- X[has_break, after_break] + offset
    }
  }

  if (isTRUE(constant_block)) {
    block <- which(rows <= 2 & cols <= 2)
    X[block, ] <- 5
    slope[block] <- 0
    has_break[block] <- FALSE
  }
  # NULL specifically when no cell ends up with a break -- either
  # break_fraction rounded down to 0 active blocks, or constant_block
  # happened to remove every break-affected cell above.
  break_time_returned <- if (any(has_break)) break_time_used else NULL

  layers <- lapply(seq_len(n_time), function(k) terra::setValues(r, X[, k]))
  series <- do.call(c, layers)
  names(series) <- paste0("t", seq_len(n_time))

  true_slope <- terra::setValues(r, slope)
  names(true_slope) <- "true_slope"

  true_break <- terra::setValues(r, as.numeric(has_break))
  names(true_break) <- "true_break"

  true_signal <- terra::setValues(r, as.numeric(slope != 0))
  names(true_signal) <- "true_signal"
  true_direction <- terra::setValues(r, sign(slope))
  names(true_direction) <- "true_direction"

  parameters <- list(
    nrow = nrow, ncol = ncol, n_time = n_time,
    trend_strength = trend_strength, trend_shape = trend_shape,
    trend_fraction = trend_fraction, ar1 = ar1, noise_sd = noise_sd,
    noise_dist = noise_dist, t_df = t_df,
    smooth_radius = smooth_radius, spatial_model = spatial_model,
    spatial_rho = spatial_rho, spatial_range = spatial_range,
    spatial_smoothness = spatial_smoothness,
    signal_size = signal_size, signal_location = signal_location,
    signal_angle = signal_angle, signal_axis_ratio = signal_axis_ratio,
    custom_mask = custom_mask, constant_block = constant_block,
    break_type = break_type, break_time = break_time_returned,
    break_fraction = break_fraction, break_magnitude = break_magnitude,
    seed = seed)
  diagnostics <- list(
    target_lag1_correlation = if (spatial_model %in%
      c("gaussian", "exponential")) spatial_rho else NA_real_,
    covariance_model = spatial_model,
    one_call_is_one_replicate = TRUE)

  result <- list(
    series = series, true_slope = true_slope,
    true_signal = true_signal, true_direction = true_direction,
    true_break = true_break, break_time = break_time_returned,
    parameters = parameters, diagnostics = diagnostics)
  class(result) <- c("sptrends_simulation", "sptrends", "list")
  result
}

#' @noRd
.sim_noise <- function(n, sd, dist, df) {
  if (dist == "gaussian") {
    stats::rnorm(n, sd = sd)
  } else {
    if (df <= 2) {
      stop("'t_df' must be > 2 for the t-distribution to have finite variance.")
    }
    scale_factor <- sqrt((df - 2) / df)
    sd * scale_factor * stats::rt(n, df = df)
  }
}
