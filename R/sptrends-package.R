#' sptrends: Statistical Inference for Spatiotemporal Trends in Gridded
#' Data
#'
#' sptrends is a framework for statistical inference of spatiotemporal
#' trends in gridded environmental data. Trend analysis over
#' thousands of spatially structured time series raises three
#' fundamental statistical challenges: serial correlation, spatial
#' dependence, and multiple testing. sptrends addresses these through a
#' coherent set of methods and integrated workflows for preprocessing,
#' trend estimation, statistical inference, effect-size estimation, and
#' multiple-testing correction. Originally developed as the reference
#' implementation of the published True Significant Trends (TST)
#' framework, it has evolved into a broader platform for developing,
#' comparing and applying methods for spatiotemporal trend analysis.
#'
#' Testing thousands of grid cells at once, each with its own
#' autocorrelated time series and its own spatial neighbours, breaks the
#' independence assumptions behind a plain, cell-by-cell trend test in
#' three separate ways -- serial correlation, spatial autocorrelation, and
#' multiple testing. This package includes two published, integrated
#' workflows that address these explicitly rather than quietly assuming
#' them away, differing in exactly how:
#'
#' \strong{True Significant Trends (TST)}, via [workflow_tst()], of
#' Gutiérrez-Hernández & García (2025): selective AR(1) prewhitening,
#' Contextual Mann-Kendall testing, Theil-Sen slope estimation, and
#' adaptive (BKY) false discovery rate control.
#'
#' \strong{Robust Trend Analysis (RTA)}, via [workflow_rta()], of
#' Gutiérrez-Hernández & García (2024): Theil-Sen slope estimation,
#' Contextual Mann-Kendall testing, and standard (BH) false discovery
#' rate control, without prewhitening.
#'
#' See `?workflow_tst` and `?workflow_rta` -- especially the latter's
#' "Comparison with TST" subsection under "Methodological details" --
#' for how the two differ and when you might
#' prefer one over the other. Both are offered as genuinely different,
#' published, citable methods; this package does not treat one as
#' superseding the other.
#'
#' A third workflow function, [workflow_trends()], is not a published
#' method in its own right -- it lets you assemble your own choice of
#' prewhitening, trend testing, slope estimation and multiple-testing
#' correction method, for the case where neither TST nor RTA matches
#' what a given analysis needs, or for comparing how sensitive a
#' result is to that choice. See `vignette("g-workflow-trends")`.
#'
#' @section TST's four steps:
#' \enumerate{
#'   \item [prewhiten()] -- selective AR(1) prewhitening,
#'     removing serial autocorrelation that would otherwise inflate false
#'     positives in the trend test.
#'   \item [trend_test()] -- the Contextual Mann-Kendall trend
#'     test (TST's own published choice; [trend_test()] itself also
#'     offers classic, non-contextual Mann-Kendall via `method = "MK"`,
#'     and a classical OLS-based trend test via `method = "OLS"`, both
#'     for use outside this specific workflow), which borrows
#'     statistical strength from each cell's spatial
#'     neighbourhood. Its default 3 by 3 region follows Neeti and Eastman
#'     (2011), as implemented in TerrSet's Kendall module; larger odd
#'     `window_size` values change
#'     the contextual scale using the same equations. [slope_estimator()]
#'     gives the accompanying rate of change.
#'   \item [fdr_correction()] -- false discovery rate correction
#'     (Benjamini-Hochberg, adaptive two-stage
#'     Benjamini-Krieger-Yekutieli, or opt-in Benjamini-Yekutieli) of
#'     the resulting p-values,
#'     controlling for the fact that many cells were tested at once.
#'     [spatial_autocorrelation()] (directly, or via `moran_check = TRUE`)
#'     provides a qualified diagnostic of spatial dependence relevant to
#'     interpreting this step; it does not prove its complete assumptions.
#' }
#'
#' Run them one at a time, or all together with [workflow_tst()]:
#' `result <- workflow_tst(x); plot(result)`. See
#' `vignette("a-getting-started")` for a
#' quick start, or `vignette("b-prewhitening")`,
#' `vignette("c-trend-test")`,
#' `vignette("d-slope-estimation")`,
#' `vignette("e-fdr-correction")`, and
#' `vignette("g-workflow-trends")` for concise introductory guides.
#' `workflow_rta()`
#' runs its own,
#' shorter three-step workflow the same way:
#' `result <- workflow_rta(x); print(result)` -- see `?workflow_rta`.
#'
#' @section Core, preprocessing, support, and reporting functions:
#' Every exported function's help page states its type under "Function
#' type", but as an index:
#'
#' **Core** -- the building blocks of TST and RTA:
#' [read_ordered_stack()]/[read_netcdf_stack()],
#' [trend_test()] (with its own
#' `print()`/`summary()`/`plot()` -- see [print.sptrends()]),
#' [fdr_correction()] (with its own
#' `print()`/`summary()`/`plot()` -- see [print.sptrends()]),
#' [workflow_tst()] (which
#' chains preprocessing, the previous two, and [slope_estimator()]) and
#' [workflow_rta()] (which chains [trend_test()], [slope_estimator()],
#' and [fdr_correction()], without preprocessing).
#'
#' **Preprocessing** -- prepare the raw raster time series before trend
#' estimation or significance testing, not one of the core pillars
#' themselves: [prewhiten()] (with its own
#' `print()`/`summary()`/`plot()` -- see [print.sptrends()]) and
#' [compute_anomalies()].
#'
#' **Support** -- compute something real, but are not one of the core
#' building blocks above (used internally by, or as a standalone
#' diagnostic alongside, a core function): [slope_estimator()] and
#' [spatial_autocorrelation()] (both with their own
#' `print()`/`summary()`/`plot()` -- see [print.sptrends()]),
#' whose local p-value raster can be passed to [fdr_correction()] for
#' BH, BKY or BY multiple-testing control,
#' [sim_trend_stack()], [simulation_design()], [compare_detections()],
#' [benchmark_methods()], and [benchmark_summary()]
#' (simulation, design, comparison and benchmark results also have unified
#' `print()`/`summary()`/`plot()` methods; tabular results remain ordinary
#' data frames underneath -- including `compare_detections()`'s
#' `replicates = TRUE` mode for
#' aggregating scores across several simulated runs, folded in rather
#' than kept as a separate function). `fdr_bh()`,
#' `fdr_bky()`, and `prepare_cmk_neighbourhood()` also belong here
#' conceptually, but are not exported -- see [fdr_correction()]'s
#' `method` argument (computes both by default) and
#' [trend_test()]'s `precomputed_neighbourhood` argument
#' respectively. [example_data()] gives the path to the package's
#' bundled real-world example dataset.
#'
#' **Reporting/derived** -- what remains standalone after the redesign
#' above gave `workflow_tst()`, `workflow_rta()`, [trend_test()],
#' [slope_estimator()], [prewhiten()], [fdr_correction()],
#' [spatial_autocorrelation()], and [compare_detections()] their own
#' `print()`/`summary()`/`plot()` methods:
#' [inspect_ts_cell()] (an interactive, click-to-inspect single-cell/
#' polygon time series viewer). `.moran_category()` (an internal
#' helper, distinct from the separate, still-`:::`-reachable
#' `classify_moran()` diagnostic) and
#' `direction_map()` also belong here conceptually, but are not
#' exported -- the former is folded into `print()`/`summary()` of a
#' `"spatial_autocorrelation"` object automatically; the latter
#' computes a genuine new
#' raster (unlike the report/plot wrappers around it, which were
#' folded into [plot.sptrends()]) but is reachable via `plot(x, which
#' = "direction")` for `workflow_tst()`/`workflow_rta()` results, since
#' the corrected and
#' uncorrected direction are the same underlying computation with an
#' optional significance filter, not two genuinely different results.
#'
#' @section Software quality assurance:
#' Statistical validation and software checks are treated as separate,
#' complementary layers. Numerical agreement with another package does
#' not replace tests of raster handling, and high code coverage does not
#' establish statistical correctness.
#'
#' **Internal automated tests** cover published equations, hand-worked
#' examples, edge cases, invalid inputs, missing and constant series,
#' serial and parallel equivalence, raster geometry, stable return
#' structures, and user-facing reporting. Regression tests protect
#' previously corrected behaviour, including exact messages or object
#' shapes where these form part of the public interface.
#'
#' **Release-check protocol**:
#'
#' - regenerate `NAMESPACE` and all `.Rd` files with
#'   `devtools::document()`;
#' - run `devtools::check(args = "--as-cran")`, with compact vignettes;
#' - inspect global and per-file coverage with
#'   `covr::package_coverage()` and use `covr::zero_coverage()` to locate
#'   individual unexercised expressions; each gap is reviewed, because
#'   coverage is a diagnostic rather than a substitute for numerical
#'   validation;
#' - check documentation spelling, URLs, code style, and broader package
#'   practice with `spelling`, `urlchecker`, `lintr`, and `goodpractice`;
#' - run win-builder for R-devel, R-release, and R-oldrelease, and use
#'   r-hub for additional operating systems and toolchains;
#' - build the pkgdown site and manual for visual review, and inspect
#'   packaged data with `tools::checkRdaFiles()`;
#' - record only real, recent results from those runs in
#'   `cran-comments.md`; no success result is inferred from the existence
#'   of this protocol.
#'
#' **External numerical controls** compare independent implementations
#' only where they calculate the same estimand:
#'
#' - classic Mann-Kendall against `Kendall` and `trend`;
#' - CMK against its published RAMK equations and the open-source
#'   `ConMK` implementation;
#' - regional score aggregation against `rkt`, without requiring its
#'   Hirsch-Slack variance to equal CMK's analytical RAMK variance;
#' - Theil-Sen against `trend`, repeated median against the convention
#'   implemented by `robslopes`, and OLS against `stats`;
#' - prewhitening against `modifiedmk`, with additional source- and
#'   output-level checks involving `zyp` and `MannKendallTrends`;
#' - BH/BY against `stats::p.adjust()` and the two-stage BKY behaviour
#'   against the implementation documented by `multtest`;
#' - raster values, geometry, time metadata, and file round trips
#'   through `terra`.
#' - the simulation and benchmarking cycle against analytical covariance,
#'   `fields::Matern()`, and `Kendall::MannKendall()` in 33 prespecified
#'   external controls. The full recorded experiment used 1,000 fields per
#'   spatial model and 500 paired replicates in each of eight scenarios;
#'   every control passed. Exact scope and limitations are retained under
#'   `inst/validation/`.
#'
#' TerrSet is retained as a complementary historical comparison rather
#' than ground truth because its contextual intermediate statistics are
#' not exposed. Reproducible external scripts and frozen results that
#' cannot run during `R CMD check` are stored under `inst/validation/`.
#'
#' Additional engineering controls include:
#' - **Runtime feedback**: iterative public computations use one
#'   dependency-free display reporting completed progress, elapsed duration
#'   and estimated time remaining. Indivisible operations report elapsed
#'   duration because a defensible remaining-time estimate is unavailable.
#'   Set `verbose = FALSE` to suppress runtime messages.
#' - **API stability**: the public function names, argument names, and
#'   argument order were audited for consistency across every exported
#'   function (e.g. `report` before `verbose` wherever both exist,
#'   `method` always immediately after the primary data argument,
#'   `seed` before `n_cores` wherever both exist) and are considered
#'   stable from this point forward -- breaking renames are not
#'   expected in future releases without a compelling reason.
#' - **Reproducibility**: every use of randomness in this package's own
#'   examples, tests, and simulations was audited to confirm an
#'   explicit `seed` is set wherever the result actually depends on it
#'   (calls that error out on invalid input before reaching any random
#'   step, or that use a `method` with no randomness involved, are the
#'   only exceptions, and were checked individually, not assumed safe).
#' - **Consistent `print()`/`summary()`/`plot()` styling**: every
#'   `"sptrends"` subclass's `print()` output opens with the same
#'   `<Title...>` bracketed convention, and `summary()` methods
#'   delegate to one shared underlying `_summary()` function per result
#'   type rather than duplicating formatting logic -- [compare_detections()]
#'   is the one deliberate exception, since it remains an ordinary data
#'   frame by design and is printed as one.
#' - **Fail-fast, specific error messages**: functions taking a
#'   `SpatRaster` time series check for zero complete-time-series cells
#'   immediately, before any computation, with a message naming what
#'   there was nothing to do -- rather than silently continuing through
#'   the full computation and eventually erroring (or worse, returning
#'   `NA`/`NaN`-filled output without erroring at all) somewhere
#'   downstream. Two real, previously-unguarded degenerate cases were
#'   found and fixed this way: [spatial_autocorrelation()] returning an
#'   uncaught, silent `NaN` for a perfectly constant raster, and
#'   [sim_trend_stack()] dividing by zero for a 1x1 grid under two of
#'   its three `trend_shape` options.
#' - **Performance**: reviewed for genuinely costly patterns (nested
#'   per-cell loops, repeated raster read/write round trips), not
#'   optimised indiscriminately at the expense of readability.
#'   [sim_trend_stack()] retains its batched legacy focal smoother for
#'   compatibility, while formal Gaussian, exponential and Matérn fields use
#'   circulant embedding and FFT rather than dense covariance matrices.
#' - **Memory model**: the principal analytical functions currently
#'   materialise the complete cell-by-time matrix with `terra::values()`.
#'   Vectorisation, sparse adjacency and optional PSOCK parallelism reduce
#'   runtime, but do not make the algorithms out-of-core. Peak memory grows
#'   with `ncell(x) * nlyr(x)` and can include several matrices of that
#'   size, especially for CMK and prewhitening. A reproducible scaling
#'   benchmark is installed at `inst/benchmarks/scalability.R`; users should
#'   test representative dimensions before processing very large rasters.
#' - **External validation**: this package's core statistics are
#'   checked, on real output rather than by inspection alone, against
#'   independent, long-established CRAN implementations of the same
#'   published methods -- the Mann-Kendall S statistic against
#'   `Kendall::MannKendall()` and `trend::mk.test()`, the Theil-Sen
#'   slope against `trend::sens.slope()`, and Yue-Pilon (2002)
#'   trend-free prewhitening's own pre-prewhitening slope against
#'   `modifiedmk::tfpwmk()` (confirmed, by reading its source directly,
#'   to implement the same published algorithm this package's own
#'   `prewhiten(method = "TFPW_Y")` does). Deterministic quantities
#'   with only one correct value (S, the Theil-Sen slope) matched
#'   these external packages exactly, down to a difference of `0`, not
#'   merely approximately; the p-value/z-statistic is checked only
#'   approximately, since different packages make different,
#'   individually legitimate choices about continuity correction and
#'   the variance-of-S formula under ties. See
#'   `test-external-validation.R` in the package's own test suite for
#'   the full detail.
#' - **Consistent visual identity across every map this package
#'   draws**: diverging palettes (Theil-Sen slope, the S/Sm trend
#'   statistic, direction-of-change results) always use
#'   `hcl.colors(50, "Blue-Red 3")`, deliberately colourblind-safer
#'   and domain-neutral -- not `"Green-Brown"`, which reads naturally
#'   as vegetation greening/browning but would be a strange default
#'   for, say, a temperature or precipitation series; that palette
#'   remains available and worth using explicitly for NDVI-like data
#'   specifically (see `?trend_maps`'s and `?workflow_rta`'s own
#'   `@examples` for exactly that). Every diverging map's own colour
#'   range is capped at 2 standard deviations rather than the single
#'   most extreme cell (`.robust_diverging_range()`, internal), so one
#'   outlier cannot wash out everyone else's colour into a narrow,
#'   barely-distinguishable band near the palette's midpoint. Purely
#'   binary significant/non-significant results (`fdr_significance_maps()`,
#'   `fdr_threshold_plot()`) use this package's own official brand
#'   navy (`"#071D4B"`, one of a fixed 5-colour palette --
#'   `.sptrends_brand`, internal -- taken from the exact colours used
#'   across every version of `man/figures/logo.png`, not re-derived
#'   from the image itself) for "significant" against a neutral grey
#'   for "not significant" -- the one colour a reader learns to
#'   recognise as this package's consistent "yes" signal, wherever it
#'   appears.
#'
#' @references
#' - Gutiérrez-Hernández, O. and García, L.V. (2025) Uncovering true
#'   significant trends in global greening. Remote Sensing Applications:
#'   Society and Environment, 37, 101377.
#'   \doi{10.1016/j.rsase.2024.101377}
#' - Gutiérrez-Hernández, O. and García, L.V. (2024) Robust Trend Analysis
#'   in Environmental Remote Sensing: A Case Study of Cork Oak Forest
#'   Decline. Remote Sensing, 16(20), 3886. \doi{10.3390/rs16203886}
#'
#' @keywords internal
"_PACKAGE"
