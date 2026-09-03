#' Slope estimators for raster time series
#'
#' Estimates the magnitude of a monotonic trend, independently of its
#' statistical significance. [trend_test()] primarily tells you
#' *whether* there is evidence of a trend and its direction; although
#' its `"OLS"` branch necessarily returns a fitted coefficient, this
#' function provides the dedicated and method-independent estimation of
#' the *rate* of change, per cell, using one of three published
#' estimators (see the `method` argument below).
#'
#' **Function type:** **Companion function to `trend_test()`** -- one
#' of the three quantities a standard trend analysis normally reports
#' (alongside significance and multiple-testing correction), not a mere
#' supporting utility, though not one of `trend_test()` or
#' `fdr_correction()`'s own inferential building blocks either.
#'
#' @section Typical use:
#' ```
#' raster time series
#'     |
#' slope_estimator()
#'     |
#' change per time unit (`result$slope`)
#' ```
#' Run this beside [trend_test()] on the same analytical series: the
#' slope measures magnitude, whereas the trend test assesses evidence.
#' If the test uses a prewhitened series, choose deliberately whether the
#' scientific estimand is the slope of that transformed series or of the
#' original series.
#'
#' @section Methodological details:
#' **Why estimate the slope separately from the trend test?**
#'
#' Statistical significance and effect size answer different
#' questions. A statistically significant trend may be negligible in
#' magnitude, whereas a large estimated slope may fail to reach
#' significance when uncertainty is high. For this reason, sptrends
#' separates trend testing ([trend_test()]) from slope estimation
#' (this function), allowing each quantity to be interpreted
#' independently -- neither substitutes for the other.
#'
#' **Methods and method selection**
#'
#' - **Original publication**: Theil (1950), later generalised by Sen
#'   (1968) into the median-of-pairwise-slopes form used here (linking it
#'   to Kendall's tau).
#' - **Main references**: Theil (1950) and Sen (1968), the two papers
#'   this estimator is directly named after. Full citations in
#'   "References" below.
#' - **Typical applications**: estimating the magnitude of a monotonic
#'   trend when the data may contain outliers or heavy-tailed noise --
#'   robust up to a ~29% breakdown point, unlike ordinary least squares.
#' - **Theil-Sen vs. OLS**: both quantify a linear rate of change, but
#'   they define and estimate that rate differently, under different
#'   statistical assumptions -- Theil-Sen makes no distributional
#'   assumption about the noise and tolerates a substantial fraction
#'   of outliers before breaking down; OLS is highly sensitive to
#'   influential observations and heavy-tailed errors. Normality is
#'   not required to calculate the OLS slope itself, although
#'   conventional small-sample inference for it (its own standard
#'   errors and significance test) relies on additional distributional
#'   assumptions the point estimate does not need. See the `method`
#'   argument below for when each is the more defensible choice.
#'
#' **Implementation notes (`method = "RM"` specifically)**
#'
#' `"RM"` implements Siegel's (1982) repeated median estimator using
#' the upper-median convention `robslopes::RepeatedMedian()` itself
#' uses: for a vector of length `k`, the upper median is
#' `sort(v)[floor((k + 2) / 2)]`; this differs from `stats::median()`,
#' which averages the two central observations when `k` is even.
#'
#' This implementation computes the estimator directly in `O(n^2)`
#' time. It therefore reproduces the estimator itself, but not the
#' faster, quasilinear-time algorithm (Matousek, Mount and Netanyahu,
#' 1998) `robslopes` uses internally in C++.
#'
#' The intercept follows the "direct" (also called "separate")
#' repeated-median convention, applying the same nested upper-median
#' operation to the pairwise intercepts, rather than the simpler
#' "hierarchical" convention (`median(y - slope * t)`, using the
#' already-computed overall slope) -- both are legitimate conventions
#' Siegel (1982) himself describes, but only the former matches
#' `robslopes::RepeatedMedian()`'s own intercept output.
#'
#' **Statistical assumptions and interpretation**
#'
#' All methods require finite, unique and strictly increasing time
#' values. Duplicates would produce undefined pairwise slopes for the
#' robust estimators, while unordered values would no longer describe
#' the chronological layer order; both are rejected before calculation.
#' Valid cells must have a complete time series. Each estimator
#' summarises change as one overall linear rate: `"OLS"` targets the
#' least-squares slope, `"TS"` the median pairwise slope, and `"RM"` the
#' repeated-median slope. The robust estimators do not require normally
#' distributed errors; OLS can still be calculated without normality,
#' but is more sensitive to outliers and influential observations.
#'
#' **Linear long-term change and seasonality**
#'
#' These estimators quantify an overall linear rate of change; they do
#' not model seasonal cycles or nonlinear temporal structure. For
#' strongly seasonal series, consider estimating the slope from
#' anomalies (see [compute_anomalies()]) or using an appropriate
#' seasonal model instead.
#'
#' **Computational considerations**
#'
#' **Performance -- read this before using long time series**
#'
#' The estimators have substantially different computational profiles.
#' OLS is the fastest because it uses one vectorised closed-form
#' calculation. Exact Theil-Sen is slower because it evaluates pairwise
#' slopes, but it offers a strong practical balance between robustness,
#' interpretation and computational cost and is therefore the
#' recommended general-purpose estimator. The directly implemented
#' repeated median is usually much slower than both because it calculates
#' nested medians for every observation in every cell; reserve RM for
#' cases in which its higher breakdown point is scientifically needed
#' and its additional cost is acceptable.
#'
#' The exact estimator needs **every** pairwise slope, `n*(n-1)/2` of them
#' per cell -- unlike the Mann-Kendall S statistic, this cannot be
#' accumulated as a running sum, so it does not benefit from the same
#' O(n^2)-time-but-O(1)-memory trick. For a modest `n` (a few dozen time
#' steps) this is fast. For long series (hundreds to thousands of steps --
#' e.g. multi-decade monthly data) the number of pairs per cell can reach
#' the hundreds of thousands, and computing an exact median that many times
#' per cell, over every cell, gets slow. Two ways to cope:
#' \itemize{
#'   \item `n_cores > 1`: splits cells across a `parallel::makeCluster()`
#'     PSOCK cluster (each cell's pairs still computed exactly).
#'   \item `max_pairs`: if `n*(n-1)/2` exceeds this, a random sample of
#'     `max_pairs` pairs is used per cell instead of the full set (a
#'     standard approximation for Theil-Sen on long series). Random
#'     pair sampling approximates the exact Theil-Sen slope;
#'     increasing `max_pairs` generally reduces Monte Carlo
#'     variability between runs and improves agreement with the exact
#'     estimate, though it is not guaranteed to be exactly unbiased in
#'     every finite sample. `seed` below controls the reproducibility
#'     of that approximation across runs. Set to `Inf` to force the
#'     exact computation regardless of `n`.
#' }
#' If speed genuinely matters more than robustness to outliers -- very
#' large rasters, very long series, and residuals not expected to be
#' heavy-tailed or outlier-prone -- `method = "OLS"` sidesteps this
#' entire performance question: a single closed-form matrix operation,
#' with no pairwise sampling and no per-cell iteration at all.
#'
#' **Optional queen-neighbourhood smoothing -- read the caveats first**
#'
#' This is an optional visual/post-processing step, not part of the
#' published Theil-Sen estimator itself.
#'
#' **What it does**: `smooth_neighbourhood = TRUE` replaces each cell's
#' slope with the median of its own slope and its 8 queen neighbours (a
#' 3x3 focal median), computed *after* the per-cell Theil-Sen estimation
#' above -- it does not change how any individual cell's own slope is
#' estimated.
#'
#' **What it does not inherit**: this is **not the same thing as**
#' [trend_test()]'s `neighbourhood` argument, and does not carry the
#' same justification. CMK's neighbourhood-adjusted statistic follows a
#' published, validated method (Neeti & Eastman, 2011) for the specific
#' question "is there a trend", where borrowing spatial evidence for a
#' yes/no decision is on solid ground. Smoothing the *magnitude* this
#' way has no equivalent literature backing here -- it assumes
#' neighbouring cells share similar true slopes, which is not always
#' true (e.g. a valley cell next to a sunlit slope can have a genuinely
#' different real trend from its neighbours), and in that case this
#' would blend two different real signals into one, less accurate,
#' number for both.
#'
#' **When it might make sense**: purely as a display/visual-smoothing
#' aid for a map that will be read at a glance, where averaging out
#' cell-to-cell estimation noise is more valuable than preserving every
#' individual cell's own independent estimate -- not as a way to
#' "improve" the underlying statistic.
#'
#' **Warnings**: off by default for the reasons above. If you use it,
#' validate it first for your own data and resolution using
#' [sim_trend_stack()] and [compare_detections()] against known ground
#' truth, rather than judging it only by whether the smoothed map looks
#' visually more coherent -- a smoother-looking map is not the same as
#' a more accurate one.
#'
#' One thing this smoothing does **not** do: extend the footprint of
#' "has data" beyond where data actually exists. A cell with no complete
#' time series of its own (its own slope is `NA`) stays `NA` after
#' smoothing, even if all 8 of its neighbours have valid slopes --
#' `terra::focal()`'s own default behaviour would otherwise fill such a
#' cell in from its neighbours, which would be presenting an estimate
#' for a location that was never actually observed.
#'
#' **Limitations**
#'
#' These estimators describe one overall linear rate and do not identify
#' breakpoints, nonlinear trajectories, or seasonal components. Cells
#' with incomplete time series are returned as `NA`. Optional spatial
#' smoothing changes the reported local magnitude and must not be
#' interpreted as part of any of the three published estimators.
#'
#' **Quality assurance**
#'
#' Theil-Sen slopes are compared exactly with
#' `trend::sens.slope()`, including tied data; OLS results are checked
#' against `stats::lm()`; and repeated-median slopes and intercepts are
#' checked against direct hand implementations. Automated tests also
#' cover irregular time coordinates, missing and constant series,
#' optional neighbourhood smoothing, raster return types, and
#' sequential/parallel equivalence. See `?sptrends` for the common
#' release-check protocol.
#'
#' @param x A `terra::SpatRaster`; each layer is one time step, in
#'   increasing chronological order.
#' @param method Which slope estimator to use.
#'
#'   `"TS"` (default): the classic Theil-Sen estimator -- the
#'   median of all pairwise slopes per cell, in the same units as `x`
#'   per unit of `t`. Robust (outlier-resistant) up to a ~29% breakdown
#'   point, at the cost of the performance profile described in
#'   "Computational considerations" below. This package's own default --
#'   see "Methodological details" below for why.
#'
#'   `"OLS"`: ordinary least squares -- the closed-form linear
#'   regression slope, `cov(t, y) / var(t)` per cell, computed directly
#'   with a single vectorised matrix operation (no pairwise sampling, no
#'   iteration -- `max_pairs`, `seed`, and `n_cores` below are all
#'   ignored). Faster, and the more familiar choice to most users and to
#'   readers used to classical linear regression. Worth
#'   choosing over the default specifically when: speed matters more
#'   than robustness (very large rasters or very long series -- see
#'   "Computational considerations" below); the residuals are genuinely
#'   expected to be
#'   close to normally distributed with no substantial outlier risk (in
#'   which case OLS is also the more statistically efficient of the
#'   two, with lower variance for the same data); or the result needs
#'   to be directly comparable to other analyses reported as classical
#'   linear-regression slopes, rather than a median-of-pairwise-slopes
#'   estimate.
#'
#'   `"RM"`: Siegel's (1982) repeated median estimator -- for each
#'   cell, the median slope from every observation to every other
#'   observation is computed first (one median per observation), then
#'   the median of those per-observation medians is taken as the final
#'   slope. Like `"TS"`, median-based and highly robust -- but with a
#'   50% breakdown point rather than Theil-Sen's own ~29%, at a real
#'   computational cost: this is a naive O(n^2) port (see
#'   "Implementation notes" below), not the quasilinear-time algorithm
#'   the reference implementation this was ported from and verified
#'   against, `robslopes::RepeatedMedian()`, itself uses. Choose this
#'   over `"TS"` specifically when a dataset is suspected to have
#'   enough outliers or leverage points that Theil-Sen's own lower
#'   breakdown point might not fully resist them; `max_pairs`, `seed`
#'   and `n_cores` below are not used by this method (unlike `"TS"`,
#'   it does not sample pairs, and is not currently parallelised).
#' @param t Numeric vector of time points, one per layer. Defaults to
#'   `1:nlyr(x)`. Irregular spacing is supported, but values must be finite,
#'   unique and strictly increasing.
#' @param max_pairs Integer or `Inf`. Only used when `method =
#'   "TS"`. If the number of possible pairs
#'   (`n*(n-1)/2`) exceeds this, a random sample of `max_pairs` pairs is
#'   used per cell (see "Computational considerations" below). Default
#'   `100000`; set to
#'   `Inf` for the exact computation regardless of series length.
#' @param seed Integer or `NULL`. Only used when `method = "TS"`
#'   and subsampling actually happens (random seed for the pair
#'   sampling).
#' @param n_cores Integer. Only used when `method = "TS"`. Number
#'   of cores for the cell loop. `1`
#'   (default): sequential. `> 1`: uses a `parallel::makeCluster()` PSOCK
#'   cluster.
#' @param smooth_neighbourhood Logical. Default `FALSE` -- the
#'   estimated slope (whichever `method` was used) is left as computed,
#'   per cell independently. If `TRUE`, each cell's
#'   slope is replaced by the **median** of the already-estimated slopes
#'   in its queen 3x3 neighbourhood (itself and its 8 surrounding cells)
#'   -- a queen-neighbourhood median filter applied to the per-cell
#'   slopes already computed, not a different way of estimating
#'   any individual cell's own slope. See the dedicated section below;
#'   read it before setting this to `TRUE`.
#' @param report Logical. If `TRUE` (default), automatically print the
#'   summary ([slope_summary()]) and draw the map ([slope_map()])
#'   once the slope finishes computing.
#' @param verbose Logical. Print progress messages and elapsed time.
#' @param shared_cluster Advanced; most users never need this directly.
#'   An already-running `parallel::makeCluster()` PSOCK cluster (e.g.
#'   from [workflow_tst()]'s own `n_cores`) to reuse instead of
#'   building a new one from `n_cores` above. When supplied, `n_cores`
#'   is ignored -- the shared cluster's own size was already decided
#'   by whoever built it. `NULL` (default): builds and tears down its
#'   own cluster from `n_cores`, exactly as before this argument
#'   existed.
#'
#' @return Returns an object of class `c("slope", "sptrends")`: a list with
#'   \item{slope}{A single-layer `terra::SpatRaster`, named
#'     `"theilsen_slope"`, `"ols_slope"` or `"rm_slope"` depending on
#'     `method` (units: `x` per unit of `t`, e.g. per year if `t` is in
#'     years).}
#'   \item{intercept}{Only for `method = "RM"`: a single-layer
#'     `terra::SpatRaster` (named `"rm_intercept"`, `NA` for invalid
#'     cells), the repeated median estimator's own intercept -- see
#'     "Implementation notes" below for how it is computed. Smoothed
#'     the same way `slope` itself is when `smooth_neighbourhood =
#'     TRUE`, for the same reason: a `SpatRaster`, matching `slope`'s
#'     own representation, rather than a plain numeric vector.}
#'   \item{method}{Character: which `method` produced this result
#'     (`"TS"`, `"OLS"` or `"RM"`).}
#'   \item{smoothed}{Logical: whether post-estimation neighbourhood
#'     smoothing (`smooth_neighbourhood = TRUE`) was applied to
#'     `slope` after estimation, not whether smoothing was itself part
#'     of how `slope` was computed.}
#'
#'   Use `print()`/`summary()`/`plot()` -- see [print.sptrends()],
#'   [summary.sptrends()], and [plot.sptrends()].
#'
#' @references
#' `method = "TS"`, original estimator:
#' - Theil, H. (1950) A rank-invariant method of linear and polynomial
#'   regression analysis. Indagationes Mathematicae, 12, 85-91 (Part I;
#'   published in three parts). No DOI available (pre-DOI-era publication).
#'
#' Generalisation into the median-of-pairwise-slopes estimator used by
#' this function, linking it to Kendall's tau:
#' - Sen, P.K. (1968) Estimates of the regression coefficient based on
#'   Kendall's tau. Journal of the American Statistical Association, 63,
#'   1379-1389. \doi{10.1080/01621459.1968.10480934}
#'
#' `method = "OLS"`: the classical least-squares regression slope,
#' independently attributed to both of the following (no single
#' original paper; both pre-DOI-era):
#' - Legendre, A.M. (1805) Nouvelles méthodes pour la détermination des
#'   orbites des comètes. Firmin Didot, Paris.
#' - Gauss, C.F. (1809) Theoria motus corporum coelestium in sectionibus
#'   conicis solem ambientium. Perthes und Besser, Hamburg.
#'
#' `method = "RM"`, original estimator:
#' - Siegel, A.F. (1982) Robust regression using repeated medians.
#'   Biometrika, 69(1), 242-244. \doi{10.1093/biomet/69.1.242}
#'
#' The quasilinear-time algorithm the reference implementation this
#' was ported from and verified against, `robslopes::RepeatedMedian()`,
#' itself uses (not this package's own, deliberately simpler `O(n^2)`
#' port -- see "Implementation notes" above):
#' - Matoušek, J., Mount, D.M. and Netanyahu, N.S. (1998) Efficient
#'   Randomized Algorithms for the Repeated Median Line Estimator.
#'   Algorithmica, 20(2), 136-150. \doi{10.1007/PL00009190}
#'
#' Documents `robslopes` itself, including the exact formula and upper-
#' median convention this port was verified against:
#' - Raymaekers, J. (2023) robslopes: Efficient Computation of the
#'   (Repeated) Median Slope. The R Journal, 15(1), 249-260.
#'   \doi{10.32614/RJ-2023-012}
#'
#' This function is used (not authored) by both of this package's own
#' integrated workflows, [workflow_tst()] and [workflow_rta()], and by
#' the studies behind
#' them:
#' - Gutiérrez-Hernández, O. and García, L.V. (2025) Uncovering true
#'   significant trends in global greening. Remote Sensing Applications:
#'   Society and Environment, 37, 101377. \doi{10.1016/j.rsase.2024.101377}
#' - Gutiérrez-Hernández, O. and García, L.V. (2024) Robust Trend Analysis
#'   in Environmental Remote Sensing: A Case Study of Cork Oak Forest
#'   Decline. Remote Sensing, 16(20), 3886. \doi{10.3390/rs16203886}
#'
#' @seealso [trend_test()] for statistical evidence of temporal change;
#'   [workflow_trends()], [workflow_tst()], and [workflow_rta()] for
#'   workflows combining significance and magnitude.
#'
#' @examples
#' \donttest{
#' # Annual mean NDVI from the bundled environmental dataset.
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#'
#' # The rate of change per cell, in NDVI units per year (the raster's
#' # own time unit) -- this is "how fast", not "is it significant" (that
#' # is what trend_test()/workflow_tst() answer instead).
#' result <- slope_estimator(r, verbose = FALSE, report = FALSE)
#'
#' # A "slope" object: result$slope is the SpatRaster itself; plot() and
#' # summary() give the diverging map and the descriptive statistics
#' # without reconstructing either by hand.
#' plot(result)
#' summary(result)
#'
#' # Ordinary least squares: much faster (a single closed-form matrix
#' # operation, no pairwise sampling) -- see "Methodological details"
#' # for when this efficiency is, and is not, worth
#' # its own trade-offs relative to the two rank-based methods.
#' result_ols <- slope_estimator(r, method = "OLS", verbose = FALSE,
#'                                report = FALSE)
#' }
#'
#' @family Slope estimation functions
#' @export
slope_estimator <- function(x, method = c("TS", "OLS", "RM"), t,
                             max_pairs = 100000,
                             seed = NULL, n_cores = 1,
                             smooth_neighbourhood = FALSE, report = TRUE,
                             verbose = TRUE,
                             shared_cluster = NULL) {
  finish_timer <- .sptrends_elapsed_timer("slope_estimator()", verbose)
  on.exit(finish_timer(), add = TRUE)
  method <- match.arg(method)
  if (!inherits(x, "SpatRaster")) stop("'x' must be a terra SpatRaster.")
  max_pairs <- .validate_max_pairs(max_pairs)
  seed <- .validate_seed(seed)
  n_cores <- .validate_positive_integer(n_cores, "n_cores")

  X <- terra::values(x, mat = TRUE)
  n <- ncol(X)
  ncell_ <- nrow(X)
  if (n < 2) {
    stop("'x' must have at least 2 layers (time steps) to estimate a slope.")
  }
  if (missing(t)) t <- seq_len(n)
  t <- .validate_time_axis(t, n)
  ok <- stats::complete.cases(X)
  if (!any(ok)) {
    stop("No cell has a complete time series (every cell has at least ",
         "one NA layer) -- there is nothing to estimate a slope on.")
  }

  if (method == "OLS") {
    if (verbose) {
      message("Ordinary least squares slope (closed-form, vectorised)...")
    }
    t_centred <- t - mean(t)
    denom <- sum(t_centred^2)
    X_centred <- X - rowMeans(X)
    slope <- as.numeric(X_centred %*% t_centred) / denom
    slope[!ok] <- NA_real_

    r <- x[[1]]
    out <- terra::setValues(r, slope)
    names(out) <- "ols_slope"

    if (isTRUE(smooth_neighbourhood)) {
      if (verbose) {
        message("Applying queen-3x3 median smoothing ",
                "(smooth_neighbourhood = TRUE) -- ",
                "see ?slope_estimator for why this is off by default.")
      }
      w <- matrix(1, nrow = 3, ncol = 3)
      out <- terra::focal(out, w = w, fun = "median", na.rm = TRUE,
                           na.policy = "omit")
      names(out) <- "ols_slope"
    }

    if (verbose) message("Done.")
    if (isTRUE(report)) {
      slope_summary(out)
      slope_map(out)
      slope_direction_map(out)
      slope_histogram(out)
      slope_direction_barplot(out)
    }
    result <- list(slope = out, smoothed = isTRUE(smooth_neighbourhood),
                   method = "OLS")
    class(result) <- c("slope", "sptrends")
    return(result)
  }

  if (method == "RM") {
    if (n < 3) {
      stop("'x' must have at least 3 layers (time steps) for the ",
           "repeated median estimator to be meaningful.")
    }
    if (verbose) {
      message("Siegel (1982) repeated median slope (naive O(n^2) port, ",
              "vectorised across cells; see \"Implementation notes\" in ",
              "?slope_estimator for why this is not the quasilinear-time ",
              "algorithm robslopes::RepeatedMedian() itself uses)...")
    }
    # The specific order statistic Siegel's own estimator uses for both
    # the inner (per-observation) and outer medians is not R's own
    # median() (which averages the two middle values for an even-length
    # input) -- robslopes' own R Journal paper and source define both as
    # the *upper* median specifically, floor((k + 2) / 2) for a
    # sorted vector of length k, confirmed directly against robslopes'
    # own R-level wrapper source (`medind0 <- floor((n + 2) / 2)`, the
    # same formula) at the user's request.
    .upper_median_rows <- function(M) {
      k <- ncol(M)
      idx <- floor((k + 2) / 2)
      apply(M, 1, function(row) sort(row)[idx])
    }

    inner_medians <- matrix(NA_real_, ncell_, n)
    inner_medians_intercept <- matrix(NA_real_, ncell_, n)
    pb <- .sptrends_progress(n, "Repeated median slope (inner)", verbose)
    for (i in seq_len(n)) {
      j_idx <- setdiff(seq_len(n), i)
      dt_i <- t[j_idx] - t[i]
      slopes_i <- (X[, j_idx, drop = FALSE] - X[, i]) /
        rep(dt_i, each = ncell_)
      # The intercept of the line through observation i and each other
      # observation j is exactly X[, i] - slopes_i * t[i] (that line
      # passes through (t[i], X[, i]) by construction) -- computed
      # alongside the slope's own inner median with the same nested
      # median-of-medians structure and the same upper-median order
      # statistic, matching the "direct"/"separate" intercept
      # convention (Wikipedia's own terminology for Siegel's (1982)
      # own alternative to the simpler "hierarchical" one, median(y -
      # slope * t) using the already-computed overall slope) --
      # confirmed empirically against robslopes::RepeatedMedian()'s
      # own output at the user's request, after the "hierarchical"
      # convention this section used in an earlier version of this
      # port did not match it (a real, meaningful discrepancy, not
      # floating-point noise: 0.119 on a real comparison).
      intercepts_i <- X[, i] - slopes_i * t[i]
      inner_medians[, i] <- .upper_median_rows(slopes_i)
      # nocov start
      # covr reports these three lines (this assignment, the
      # progress-step call right after it, and the loop's own closing
      # brace) as uncovered, despite sitting inside the exact same
      # loop iteration as the line immediately above (which covr does
      # report as covered) -- every test in test-rm-slope.R that
      # exercises method = "RM" genuinely executes this line, verified
      # by the intercept's own correctness tests (which could not
      # otherwise pass; inner_medians_intercept unpopulated would make
      # every intercept NA, not merely uncovered by coverage
      # tooling). Attempted, and abandoned as unsafe, a restructuring
      # to combine this call with the one above via cbind() -- that
      # would have silently computed one median mixing slope and
      # intercept values together, a real correctness bug, not a
      # coverage fix -- confirmed and reverted before ever being
      # packaged. Left as a documented covr/byte-compilation reporting
      # artifact (R's own default byte-compilation on install is the
      # most likely cause of two structurally near-identical
      # consecutive single-line statements losing independent
      # per-line tracking) rather than contorting correct,
      # already-tested code to chase a coverage percentage.
      inner_medians_intercept[, i] <- .upper_median_rows(intercepts_i)
      .sptrends_progress_step(pb, i)
    }
    # nocov end
    .sptrends_progress_close(pb)

    slope <- .upper_median_rows(inner_medians)
    slope[!ok] <- NA_real_

    intercept <- .upper_median_rows(inner_medians_intercept)
    intercept[!ok] <- NA_real_

    r <- x[[1]]
    out <- terra::setValues(r, slope)
    names(out) <- "rm_slope"
    out_intercept <- terra::setValues(r, intercept)
    names(out_intercept) <- "rm_intercept"

    if (isTRUE(smooth_neighbourhood)) {
      if (verbose) {
        message("Applying queen-3x3 median smoothing ",
                "(smooth_neighbourhood = TRUE) -- ",
                "see ?slope_estimator for why this is off by default.")
      }
      w <- matrix(1, nrow = 3, ncol = 3)
      out <- terra::focal(out, w = w, fun = "median", na.rm = TRUE,
                           na.policy = "omit")
      names(out) <- "rm_slope"
      # The intercept is smoothed the same way, for the same reason
      # smoothing slope itself is already a display/interpretation
      # choice rather than a literal recomputation from the
      # underlying pairwise data (unavailable at this point): a
      # smoothed slope and an independently-smoothed intercept do not
      # jointly define a single fitted line reconstructed exactly
      # from either alone, the same caveat that already applies to
      # the smoothed slope by itself -- see "Optional queen-
      # neighbourhood smoothing" above.
      out_intercept <- terra::focal(out_intercept, w = w, fun = "median",
                                     na.rm = TRUE, na.policy = "omit")
      names(out_intercept) <- "rm_intercept"
    }

    if (verbose) message("Done.")
    if (isTRUE(report)) {
      slope_summary(out)
      slope_map(out)
    }
    result <- list(slope = out, intercept = out_intercept,
                   smoothed = isTRUE(smooth_neighbourhood), method = "RM")
    class(result) <- c("slope", "sptrends")
    return(result)
  }

  all_pairs <- utils::combn(n, 2)
  n_pairs_total <- ncol(all_pairs)

  if (n_pairs_total > max_pairs) {
    if (verbose) {
      message(sprintf(
        paste0("n=%d gives %d possible pairs, above max_pairs=%d -- using a ",
               "random sample of %d pairs per cell (approximate Theil-Sen). ",
               "Set max_pairs=Inf for the exact estimator."),
        n, n_pairs_total, max_pairs, max_pairs))
    }
    if (!is.null(seed)) set.seed(seed)
    keep <- sample.int(n_pairs_total, max_pairs)
    all_pairs <- all_pairs[, keep, drop = FALSE]
  } else if (verbose) {
    message(sprintf("Exact Theil-Sen: %d pairs per cell.", n_pairs_total))
  }

  i_idx <- all_pairs[1, ]
  j_idx <- all_pairs[2, ]
  dt <- t[j_idx] - t[i_idx]

  compute_cell <- function(row) stats::median((row[j_idx] - row[i_idx]) / dt)

  if (n_cores <= 1) {
    slope <- rep(NA_real_, ncell_)
    idx_valid <- which(ok)
    pb <- .sptrends_progress(length(idx_valid), "Theil-Sen slope", verbose)
    for (k in seq_along(idx_valid)) {
      cell <- idx_valid[k]
      slope[cell] <- compute_cell(X[cell, ])
      .sptrends_progress_step(pb, k)
    }
    .sptrends_progress_close(pb)
  } else {
    if (verbose) {
      message(sprintf("      Parallel over %d cores (split by cell)...",
                       n_cores))
    }
    idx_valid <- which(ok)
    n_chunks <- min(n_cores, length(idx_valid))
    chunks <- lapply(parallel::splitIndices(length(idx_valid), n_chunks),
                      function(ii) idx_valid[ii])
    process_chunk <- function(idx_chunk) {
      vapply(idx_chunk, function(k) compute_cell(X[k, ]), numeric(1))
    }
    results <- .sptrends_parallel_lapply(
      chunks, process_chunk, n_cores = n_cores,
      export_vars = c("X", "i_idx", "j_idx", "dt", "compute_cell"),
      export_env = environment(), shared_cluster = shared_cluster)
    slope <- rep(NA_real_, ncell_)
    for (m in seq_along(chunks)) slope[chunks[[m]]] <- results[[m]]
  }

  r <- x[[1]]
  out <- terra::setValues(r, slope)
  names(out) <- "theilsen_slope"

  if (isTRUE(smooth_neighbourhood)) {
    if (verbose) {
      message("Applying queen-3x3 median smoothing ",
              "(smooth_neighbourhood = TRUE) -- ",
              "see ?slope_estimator for why this is off by default.")
    }
    w <- matrix(1, nrow = 3, ncol = 3)
    # na.policy = "omit": cells that were NA (no data at that location)
    # stay NA -- terra::focal()'s own default ("all") would otherwise
    # assign a smoothed value to a no-data cell just because it has a
    # valid neighbour, silently extending the footprint of "has data"
    # beyond where any data actually exists.
    out <- terra::focal(out, w = w, fun = "median", na.rm = TRUE,
                         na.policy = "omit")
    names(out) <- "theilsen_slope"
  }

  if (verbose) message("Done.")
  if (isTRUE(report)) {
    slope_summary(out)
    slope_map(out)
    slope_direction_map(out)
    slope_histogram(out)
    slope_direction_barplot(out)
  }
  result <- list(slope = out, smoothed = isTRUE(smooth_neighbourhood),
                 method = "TS")
  class(result) <- c("slope", "sptrends")
  result
}

#' Summarise a slope_estimator() result
#'
#' Descriptive statistics of the slope raster: how many cells have a
#' valid estimate, the range and central tendency of the slope values,
#' and what fraction are increasing, decreasing, or exactly flat.
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- called internally by `report = TRUE`, and reachable from
#' outside the package via `summary()`.
#'
#' @param slope Output of [slope_estimator()] (a single-layer
#'   `terra::SpatRaster`).
#' @param path Character or `NULL`. If supplied, the table is written to
#'   this CSV path.
#' @return Invisibly, a data frame with `metric`/`value` columns.
#' @family Theil-Sen functions
#' @references
#' - Theil, H. (1950) A rank-invariant method of linear and polynomial
#'   regression analysis. Indagationes Mathematicae, 12, 85-91 (Part I;
#'   published in three parts). No DOI available (pre-DOI-era publication).
#' - Sen, P.K. (1968) Estimates of the regression coefficient based on
#'   Kendall's tau. Journal of the American Statistical Association, 63,
#'   1379-1389. \doi{10.1080/01621459.1968.10480934}
#' @examples
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' slope_result <- slope_estimator(r, report = FALSE, verbose = FALSE)
#' # Called internally by summary() on a slope_estimator() result --
#' # the public entry point is:
#' summary(slope_result)
#'
#' @keywords internal
slope_summary <- function(slope, path = NULL) {
  vals <- terra::values(slope, mat = FALSE)
  ok <- !is.na(vals)
  n_valid <- sum(ok)
  v <- vals[ok]

  n_pos <- sum(v > 0)
  n_neg <- sum(v < 0)
  n_zero <- sum(v == 0)

  message(sprintf("Valid cells: %d", n_valid))
  message(sprintf("Slope range: [%.4g, %.4g]", min(v), max(v)))
  message(sprintf("Median slope: %.4g | Mean slope: %.4g",
                   stats::median(v), mean(v)))
  message(sprintf(
    "Increasing: %d (%.1f%%) | Decreasing: %d (%.1f%%) | Flat: %d (%.1f%%)",
    n_pos, 100 * n_pos / n_valid, n_neg, 100 * n_neg / n_valid,
    n_zero, 100 * n_zero / n_valid))

  tab <- data.frame(
    metric = c("valid_cells", "min_slope", "median_slope", "mean_slope",
               "max_slope", "n_increasing", "pct_increasing",
               "n_decreasing", "pct_decreasing", "n_flat", "pct_flat"),
    value = c(n_valid, round(min(v), 6), round(stats::median(v), 6),
              round(mean(v), 6), round(max(v), 6),
              n_pos, round(100 * n_pos / n_valid, 2),
              n_neg, round(100 * n_neg / n_valid, 2),
              n_zero, round(100 * n_zero / n_valid, 2))
  )

  if (!is.null(path)) {
    utils::write.csv(tab, path, row.names = FALSE)
    message(sprintf("Table written to: %s", path))
  }
  invisible(tab)
}

#' Map of a slope_estimator() result
#'
#' A single diverging map of the slope raster, zero-centred so the
#' palette's midpoint always sits at "no change" -- see `trend_maps()`'s
#' own internal comment for why this cannot be left to
#' `terra::plot()`'s automatic range.
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- called internally by `report = TRUE`, and reachable from
#' outside the package via `plot()`.
#'
#' @param slope Output of [slope_estimator()] (a single-layer
#'   `terra::SpatRaster`).
#' @param path Character or `NULL`. If supplied, the map is saved as a
#'   PNG at this path.
#' @return `NULL`, invisibly. Called for its plotting side effect.
#' @family Theil-Sen functions
#' @references
#' - Theil, H. (1950) A rank-invariant method of linear and polynomial
#'   regression analysis. Indagationes Mathematicae, 12, 85-91 (Part I;
#'   published in three parts). No DOI available (pre-DOI-era publication).
#' - Sen, P.K. (1968) Estimates of the regression coefficient based on
#'   Kendall's tau. Journal of the American Statistical Association, 63,
#'   1379-1389. \doi{10.1080/01621459.1968.10480934}
#' @examples
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' slope_result <- slope_estimator(r, report = FALSE, verbose = FALSE)
#' # Called internally by plot() on a slope_estimator() result -- the
#' # public entry point is:
#' plot(slope_result, which = "map")
#'
#' @keywords internal
slope_map <- function(slope, path = NULL) {
  range_lim <- .robust_diverging_range(slope)

  # fill_range = TRUE is essential here, not cosmetic: terra::plot()'s
  # own default for values outside 'range' is to colour them NA (blank),
  # not saturate them to the palette's extreme colour -- without this,
  # any cell whose slope exceeds the robust range (the exact case that
  # range exists to handle) would silently vanish from the map instead
  # of showing up as strongly red/blue.
  plot_map <- function() {
    terra::plot(slope, col = .sptrends_diverging_palette(50),
                range = range_lim, fill_range = TRUE,
                main = "Theil-Sen slope")
  }
  plot_map()
  if (!is.null(path)) {
    .save_current_plot(paste0(path, "_theilsen_map.png"), plot_map)
  }
  invisible(NULL)
}

#' Direction-of-change map for a Theil-Sen slope raster
#'
#' A categorical map classifying every valid cell as `"Positive"`,
#' `"Negative"`, or exactly `"Zero"` -- the binary positive/negative
#' question a reader usually wants, with the (rare, for continuous
#' Theil-Sen estimates) exact-zero case shown separately rather than
#' folded into either side.
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- called internally by `report = TRUE`, and reachable
#' from outside the package via `plot()`.
#'
#' @param slope A `SpatRaster` of Theil-Sen slopes, e.g.
#'   `slope_estimator()`'s own `$slope`.
#' @param ... Passed on to [terra::plot()].
#'
#' @return Invisibly, a `SpatRaster` with the classification (`-1`
#'   negative, `0` zero, `1` positive).
#' @family Theil-Sen slope functions
#' @examples
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' slope <- slope_estimator(r, report = FALSE, verbose = FALSE)
#' plot(slope, which = "direction")
#'
#' @keywords internal
slope_direction_map <- function(slope, ...) {
  cls <- terra::ifel(slope > 0, 1, terra::ifel(slope < 0, -1, 0))
  .safe_categorical_plot(cls, values = c(-1, 0, 1),
                          colours = c("#1f77b4", "#D9D9D9", "#d62728"),
                          labels = c("Negative", "Zero", "Positive"),
                          main = "Direction of the Theil-Sen slope",
                          reverse = TRUE, ...)
  invisible(cls)
}

#' Histogram of a Theil-Sen slope raster's own values
#'
#' A histogram (with a kernel density overlay) of the raw, continuous
#' slope values across every valid cell -- the distribution the
#' diverging map ([slope_map()]) shows spatially, seen instead as a
#' single one-dimensional summary.
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- called internally by `report = TRUE`, and reachable
#' from outside the package via `plot()`.
#'
#' @param slope A `SpatRaster` of Theil-Sen slopes, e.g.
#'   `slope_estimator()`'s own `$slope`.
#' @param breaks Passed to [graphics::hist()]'s own `breaks` argument.
#' @param ... Passed on to [graphics::hist()].
#'
#' @return Invisibly, the `histogram` object [graphics::hist()] itself
#'   returns.
#' @family Theil-Sen slope functions
#' @examples
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' slope <- slope_estimator(r, report = FALSE, verbose = FALSE)
#' plot(slope, which = "histogram")
#'
#' @keywords internal
slope_histogram <- function(slope, breaks = 40, ...) {
  vals <- terra::values(slope, mat = FALSE)
  vals <- vals[!is.na(vals)]
  h <- graphics::hist(vals, breaks = breaks, col = "steelblue",
                       border = "white", freq = FALSE,
                       main = "Distribution of the Theil-Sen slope",
                       xlab = "Slope", ...)
  d <- stats::density(vals)
  graphics::lines(d, col = "firebrick", lwd = 2)
  graphics::abline(v = 0, col = "grey30", lwd = 1, lty = 2)
  invisible(h)
}

#' Bar chart of positive/negative/zero counts for a Theil-Sen slope raster
#'
#' A bar chart of how many valid cells have a positive, negative, or
#' exactly zero Theil-Sen slope -- the same classification
#' [slope_direction_map()] draws spatially, summarised here as raw
#' counts (or percentages).
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- called internally by `report = TRUE`, and reachable
#' from outside the package via `plot()`.
#'
#' @param slope A `SpatRaster` of Theil-Sen slopes, e.g.
#'   `slope_estimator()`'s own `$slope`.
#' @param probability Logical. If `TRUE`, bars show percentage of valid
#'   cells rather than raw counts. Default `FALSE`.
#' @param ... Passed on to [graphics::barplot()].
#'
#' @return Invisibly, a named numeric vector with the raw counts
#'   (`Positive`, `Negative`, `Zero`), regardless of `probability`.
#' @family Theil-Sen slope functions
#' @examples
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' slope <- slope_estimator(r, report = FALSE, verbose = FALSE)
#' plot(slope, which = "bar")
#'
#' @keywords internal
slope_direction_barplot <- function(slope, probability = FALSE, ...) {
  vals <- terra::values(slope, mat = FALSE)
  vals <- vals[!is.na(vals)]
  n_valid <- length(vals)
  n_pos <- sum(vals > 0)
  n_neg <- sum(vals < 0)
  n_zero <- n_valid - n_pos - n_neg
  counts <- c("Positive" = n_pos, "Negative" = n_neg, "Zero" = n_zero)
  vals_bar <- if (isTRUE(probability)) 100 * counts / n_valid else counts
  ylab <- if (isTRUE(probability)) "% of valid cells" else "Number of cells"
  graphics::barplot(vals_bar, col = c("#d62728", "#1f77b4", "#D9D9D9"),
                     main = "Direction of the Theil-Sen slope",
                     ylab = ylab, ...)
  invisible(counts)
}
