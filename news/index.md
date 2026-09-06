# Changelog

## sptrends 1.6.1

Full audit: functions, docs, website config, NEWS style, cran-comments.

- Removed 4 dead duplicate validation checks in `benchmark.R`/
  `validation.R` left over from the 1.6 audit integration.
- Added `_pkgdown.yml` entries for 5 functions exported in 1.5.3.
- Condensed NEWS from 1.5.4 onward to match the project’s usual
  telegraphic style.
- Updated `cran-comments.md` to reflect confirmed 1.6 results.
- Cross-checked `MK`/Theil-Sen against 5 external packages (Kendall,
  trend, modifiedmk, robslopes, zyp); no real discrepancies found.

## sptrends 1.6

### Bug fixes

- Fixed
  [`compute_anomalies()`](https://olive-r.github.io/sptrends/reference/compute_anomalies.md)
  climatology/standardisation alignment when `start_position` is not 1.
- Fixed
  [`read_ordered_stack()`](https://olive-r.github.io/sptrends/reference/read_ordered_stack.md)
  annual dates to use 1 January.
- Preserved constant residuals in Yue-Pilon prewhitening; skipped
  unavailable DW plots for Zhang’s method.
- Propagated worker counts in TST/RTA; honoured supplied clusters.
- Matched original-BKY plotting thresholds to its own decisions.
- Fixed zero F1 scores; benchmarks now apply `evaluation_mask`.
- Handled one-permutation summaries; rejected undefined General G inputs
  with fewer than two positive cells.
- Fixed missing-neighbour panels, single-cell histograms, and TST
  direction plots with BY or no slope estimate.

### Documentation

- Clarified cycle ordering, daily-data limits, calendar anchors.
- Updated website sources; synchronised exports with 1.5.3.
- Corrected estimator labels and MMK identification.

## sptrends 1.5.9

- Confirmed 100% coverage holds with `SPTRENDS_TEST_PARALLEL=true` set.
  No code changes.

## sptrends 1.5.8

- Fixed
  [`workflow_trends()`](https://olive-r.github.io/sptrends/reference/workflow_trends.md)’s
  example again: finds a valid cell programmatically, keeping real (not
  synthetic) data.

## sptrends 1.5.7

- Replaced the 1.5.5 example crop (had no valid cells) with synthetic
  data.

## sptrends 1.5.6

- Skipped 19 parallel-cluster tests by default
  (`SPTRENDS_TEST_PARALLEL=true` to run) to cut check time.

## sptrends 1.5.5

- Cropped
  [`workflow_trends()`](https://olive-r.github.io/sptrends/reference/workflow_trends.md)’s
  example to a small region to further cut check time.

## sptrends 1.5.4

- Used a null device for three `report = TRUE` tests to cut check time.

## sptrends 1.5.3

### Bug fixes

CRAN-requested fixes (submission feedback on 1.5.2):

- No longer modifies `.GlobalEnv`; uses `withr` for seed handling.
- Exported
  [`fdr_bh()`](https://olive-r.github.io/sptrends/reference/fdr_bh.md),
  [`fdr_by()`](https://olive-r.github.io/sptrends/reference/fdr_by.md),
  [`fdr_bky()`](https://olive-r.github.io/sptrends/reference/fdr_bky.md),
  [`prepare_cmk_neighbourhood()`](https://olive-r.github.io/sptrends/reference/prepare_cmk_neighbourhood.md),
  and
  [`plot_detection_comparison()`](https://olive-r.github.io/sptrends/reference/plot_detection_comparison.md).
- Removed internal (`:::`) calls from all examples.
- Replaced an interactive example’s `\dontrun{}` with
  `if (interactive())`.
- Added method references to the DESCRIPTION field.

Found while fixing the above:

- Fixed a conflicting `@inheritParams`/duplicate `@return` in
  [`fdr_bky()`](https://olive-r.github.io/sptrends/reference/fdr_bky.md)
  that corrupted the installed help database.
- Updated 4 stale tests asserting the newly-exported functions were
  still internal.

## sptrends 1.5.2

### Bug fixes

- Fixed all `path=`-based PNG-saving functions
  ([`trend_histograms()`](https://olive-r.github.io/sptrends/reference/trend_histograms.md),
  [`trend_maps()`](https://olive-r.github.io/sptrends/reference/trend_maps.md),
  [`fdr_pvalue_histogram()`](https://olive-r.github.io/sptrends/reference/fdr_pvalue_histogram.md),
  [`fdr_significance_maps()`](https://olive-r.github.io/sptrends/reference/fdr_significance_maps.md),
  [`fdr_comparison_barplot()`](https://olive-r.github.io/sptrends/reference/fdr_comparison_barplot.md),
  [`fdr_threshold_plot()`](https://olive-r.github.io/sptrends/reference/fdr_threshold_plot.md),
  [`fdr_direction_plot()`](https://olive-r.github.io/sptrends/reference/fdr_direction_plot.md),
  [`prewhiten_histograms()`](https://olive-r.github.io/sptrends/reference/prewhiten_histograms.md),
  [`prewhiten_maps()`](https://olive-r.github.io/sptrends/reference/prewhiten_maps.md),
  [`slope_map()`](https://olive-r.github.io/sptrends/reference/slope_map.md),
  [`spatial_autocorrelation()`](https://olive-r.github.io/sptrends/reference/spatial_autocorrelation.md)’s
  own plotting helpers,
  [`plot_detection_comparison()`](https://olive-r.github.io/sptrends/reference/plot_detection_comparison.md),
  and the `sptrends_simulation`/ `sptrends_benchmark`
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods),
  which silently produced no file at all on Debian/Linux (though not on
  Windows) – the internal helper relied on
  [`grDevices::dev.copy()`](https://rdrr.io/r/grDevices/dev2.html),
  which snapshots whatever device happens to already be active, and
  produces nothing in headless/batch environments with no active device.
  Now opens its own PNG device first and draws directly into it, which
  works identically regardless of platform or interactivity.

## sptrends 1.5.1

Resubmission of 1.5 to CRAN, with no code changes – version bump
required to avoid a filename collision on CRAN’s own incoming server
from an earlier submission attempt.

## sptrends 1.5

### New features

- [`read_ordered_stack()`](https://olive-r.github.io/sptrends/reference/read_ordered_stack.md)
  now supports explicit declaration of a series’ temporal structure
  (`files`, `time`, `cycle_type`, `start`, `end`, `time_anchor`), as an
  alternative to automatic detection from file names – recommended
  whenever the series is not simply annual. Eight unambiguous calendar
  conventions are supported: `"annual"`, `"monthly"`, `"16-day"`,
  `"semimonthly"`, `"10-day"`, `"8-day"`, `"weekly"`, and `"daily"`, all
  built with real calendar arithmetic (leap years included).
- [`compute_anomalies()`](https://olive-r.github.io/sptrends/reference/compute_anomalies.md)’s
  own `cycle_type` now uses the same vocabulary as
  [`read_ordered_stack()`](https://olive-r.github.io/sptrends/reference/read_ordered_stack.md)’s
  (`"monthly"`, `"16-day"`, `"semimonthly"`, `"10-day"`, `"8-day"`,
  `"weekly"`), so a value already used to read a stack can be reused
  directly here. `"annual"` and `"daily"` are deliberately not included
  – see
  [`?compute_anomalies`](https://olive-r.github.io/sptrends/reference/compute_anomalies.md)
  for why.

### Improvements

- Automatic mode’s gap detection is now scoped to simple annual
  sequences, avoiding false positives on combined numeric codes.
- Automatic mode now recommends explicit declaration for series that may
  not be annual.
- [`compute_anomalies()`](https://olive-r.github.io/sptrends/reference/compute_anomalies.md)
  now validates `cycle` fully (must be one finite integer `>= 2`),
  instead of only checking `cycle < 2`.

## sptrends 1.4.4

### Improvements

- `trend_test(method = "CMK")` now supports `connectivity = "rook"` as a
  genuine, working option, alongside the default `"queen"`.

## sptrends 1.4.3

### Documentation

- Added the CMK/`ConMK` external validation note to `README.md`.

## sptrends 1.4.2

### Documentation

- Documented a fresh, independent confirmation of CMK’s external
  validation against `ConMK`: the base statistic matches to
  floating-point precision, and the `continuity = TRUE` option
  reproduces `ConMK`’s own p-values exactly at the specific edge case
  (`Sm == 0`) where the two implementations would otherwise be expected
  to diverge.

## sptrends 1.4.1

### Improvements

- [`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md)
  now supports scenarios with different argument sets.
- Added validation for several
  [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)
  arguments.
- Clearer messages for
  [`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md)’s
  slope and prewhitening stages.

## sptrends 1.4

Passed a comprehensive external code audit; several improvements
resulted.

### Improvements

- Corrected method-comparison summaries to consistently treat error
  rates as “lower is better”.
- Added validation for
  [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)’s
  `ar1` and `noise_sd` arguments.
- Added the software’s own DOI to citation metadata.
- Unified the package website URL across `README.md` and `DESCRIPTION`.
- [`fdr_direction_summary()`](https://olive-r.github.io/sptrends/reference/fdr_direction_summary.md)
  now supports the `"BY"` method.
- Refined `NaN`/`NA` consistency in replicate-aggregated summaries.
- Corrected documentation for
  [`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)/[`workflow_rta()`](https://olive-r.github.io/sptrends/reference/workflow_rta.md)’s
  return value and documented a known limitation with irregular time
  spacing.

### Documentation

- Refreshed `cran-comments.md`.
- Softened prescriptive wording around `BY` in the FDR vignette.

## sptrends 1.3.4

### Improvements

- Improved reliability of categorical map rendering.

## sptrends 1.3.3

### Documentation

- Streamlined the beta-era `NEWS.md` history (0.96.4-0.97.1) too:
  shorter entries.

## sptrends 1.3.2

### Documentation

- Streamlined `NEWS.md` from 1.0 onward: shorter, improvement-focused
  entries.

## sptrends 1.3.1

### Improvements

- Improved robustness of categorical map rendering in edge cases.
- Improved tolerance of a regression test near signal boundaries.
- Direction maps: legend now lists “Increase” before “Decrease”.
- TST vignette example now calls
  [`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)
  directly.

## sptrends 1.3.0

### Improvements

- [`read_ordered_stack()`](https://olive-r.github.io/sptrends/reference/read_ordered_stack.md)
  now shows live reading progress.
- [`summary()`](https://rdrr.io/r/base/summary.html) on workflow results
  now leads with the FDR-corrected result, with uncorrected stats shown
  as diagnostic context.
- Documented reporting valid cell counts as good practice.
- Standardised valid-cell and raw/percentage reporting across every
  [`summary()`](https://rdrr.io/r/base/summary.html) table in the
  package.

## sptrends 1.2.5

### Improvements

- Improved “not significant” display consistency in smoothed direction
  plots.

## sptrends 1.2.4

### Improvements

- Refined category-colour assignment in categorical maps.

## sptrends 1.2.3

### Improvements

- Categorical maps now always render as discrete colour blocks.
- Unified category order and grey shade across every map, plot and bar
  chart. Full symbology audit completed.

## sptrends 1.2.2

### Improvements

- [`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)/[`workflow_rta()`](https://olive-r.github.io/sptrends/reference/workflow_rta.md)
  plots now show their own title.

## sptrends 1.2.1

### Improvements

- Refined the colour scheme into three distinct identities: diverging
  maps (ColorBrewer RdBu), categorical direction (Tableau red/blue), and
  binary significance (its own dedicated colour).

## sptrends 1.2.0

### Improvements

- Unified the colour scheme across the package’s plots.

## sptrends 1.1.5

### Quality assurance

- Improved test coverage in `tst-methods.R`, `rta.R` and
  `workflow_trends-methods.R`.

## sptrends 1.1.4

### Quality assurance

- Refined test coverage.

## sptrends 1.1.3

### Improvements

- Improved consistency of the direction panel across workflows.

## sptrends 1.1.2

### Documentation

- Added a “Quality assurance” section to `README.md`.

## sptrends 1.1.1

### Documentation

- Completed citations in the FDR-correction vignette.
- Audited citations and references across all six vignettes.

## sptrends 1.1

### Quality assurance

- Passed a comprehensive integral audit of the whole package.
- Confirmed via two independent rounds of external testing.

## sptrends 1.0.8

### Quality assurance

- Closed the last coverage gap; confirmed 100% test coverage and a clean
  `R CMD check`.

## sptrends 1.0.7

### Improvements

- Improved test reliability for direction-map testing.
- Added a round-trip fidelity test for
  [`read_netcdf_stack()`](https://olive-r.github.io/sptrends/reference/read_netcdf_stack.md).

## sptrends 1.0.6

### Improvements

- Improved test reliability with a stronger reference signal.

## sptrends 1.0.5

### Improvements

- Improved test compatibility with
  [`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md)
  and refined message formatting.

## sptrends 1.0.4

### Package identity

- Added README badges: license, lifecycle, Zenodo DOI, pkgdown site.

## sptrends 1.0.3

### Feature changes

- [`direction_map()`](https://olive-r.github.io/sptrends/reference/direction_map.md)
  gained an optional `slope` argument: trend direction can now be
  derived from the slope estimator’s own sign.
  [`workflow_trends()`](https://olive-r.github.io/sptrends/reference/workflow_trends.md),
  [`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)
  and
  [`workflow_rta()`](https://olive-r.github.io/sptrends/reference/workflow_rta.md)
  use this automatically when a slope result is available.

## sptrends 1.0.2

### Documentation

- Vignette and documentation improved.

## sptrends 1.0.1

### Documentation

- Polished `README.md`: official citation, author affiliations,
  reordered examples and modules.
- Polished `vignettes/a-getting-started.Rmd`: clarified TST origin,
  added references, reworked “Common mistakes”.

## sptrends 1.0

### Package identity

- Added `URL`/`BugReports` to `DESCRIPTION`.
- First stable release, confirmed clean on all win-builder platforms.

## sptrends 0.97.1

### Improvements

- Fixed four broken DOI links in vignette References sections (URL
  encoding).

### Documentation

- Added 11 words to `inst/WORDLIST`.

## sptrends 0.97

Last pre-release (beta) version. sptrends 1.0 is the first stable
release.

### Quality assurance

- Passed an extensive external test battery covering all 18 exported
  functions.

## sptrends 0.96.7

### Improvements

- Fixed `sim_trend_stack(trend_shape = "block")` ignoring
  `signal_location`/`signal_size`/`signal_angle`.

## sptrends 0.96.6

### Simulation and validation

- Corrected a validation-metric edge case.
- Confirmed simulator checks across spatial-dependence and
  multiple-testing scenarios.

## sptrends 0.96.5

### Simulation and validation

- Added field-significance power metric.
- Removed an R CMD check NOTE from console flushing.

## sptrends 0.96.4

### Simulation and validation

- Added progress/timing reporting to simulation, design, validation and
  benchmarking functions.
- Recorded a 33-control external validation of the simulation cycle.

## sptrends 0.96.3

### Simulation and validation

- Completed regression coverage for shared replicate-level validation
  inputs.

## sptrends 0.96.2

### Simulation and validation

- Expanded regression tests for benchmark scoring, covariance fallback,
  uncertainty graphics and replicate-specific validation inputs.

## sptrends 0.96.1

### Simulation and validation

- Corrected the real-data spatial validation test and limited benchmark
  stages to currently implemented method families.

## sptrends 0.96

### Simulation and validation

- Redesigned
  [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)
  and
  [`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md)
  for calibrated spatiotemporal simulation and known-truth evaluation.
  Added
  [`simulation_design()`](https://olive-r.github.io/sptrends/reference/simulation_design.md),
  [`benchmark_methods()`](https://olive-r.github.io/sptrends/reference/benchmark_methods.md),
  and
  [`benchmark_summary()`](https://olive-r.github.io/sptrends/reference/benchmark_summary.md)
  with unified S3 reporting and scenario-performance graphics.
- Standardised public examples and the documentation architecture of the
  simulation and benchmarking interfaces.

## sptrends 0.95.2

### Documentation

- Vignette and documentation improved.

## sptrends 0.95.1

### Documentation

- Vignette and documentation improved.

## sptrends 0.95

### Documentation

- Vignette and documentation improved.

## sptrends 0.94.10

### Vignettes

- Improved the multiple-testing and workflow vignettes and clarified the
  binarised workflow map.

## sptrends 0.94.9

### Vignettes

- Refined the multiple-testing vignette.

## sptrends 0.94.8

### Vignettes

- Improved the multiple-testing vignette.

## sptrends 0.94.7

### Vignettes

- Improved the package vignettes.

## sptrends 0.94.6

### Vignettes

- Improved the trend-magnitude vignette.

## sptrends 0.94.5

### Vignettes

- Improved the trend-inference vignette.

## sptrends 0.94.4

### Vignettes

- Improved the serial-correlation vignette.

## sptrends 0.94.3

### Vignettes

- Improved the data-loading and exploration vignette.

## sptrends 0.94.2

### Vignettes

- Improved the introductory vignettes.

## sptrends 0.94.1

### Documentation

- Simplified the README introduction and revised all introductory
  vignettes to use the bundled NDVI series, including a complete
  temporal mosaic, interactive animation code and linked navigation.
- Expanded the
  [`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md)
  example to show BH, BKY and BY explicitly, and listed additional CRAN
  packages used as external validation references.

## sptrends 0.94

### Package identity

- Renamed the package subtitle to *Statistical Inference for
  Spatiotemporal Trends in Gridded Data*, reflecting a general inference
  framework without restricting every component to robust statistical
  methods.

## sptrends 0.93.3

### Quality assurance

- Removed an invalid assumption about roxygen2’s internal placement of
  descriptive paragraphs from the documentation contract test.

## sptrends 0.93.2

### Quality assurance

- Corrected the documentation contract test to inspect `Function type`
  within the parsed description instead of relying on internal Rd
  ordering.

## sptrends 0.93.1

### Documentation

- Standardised `Function type` placement across the remaining documented
  helpers.

## sptrends 0.93

### Documentation

- Rebuilt six concise introductory vignettes around a common teaching
  structure, integrating TST and RTA into the trend-workflows guide.

## sptrends 0.92.2

### Documentation

- Completed the cross-references in four core help pages.

### Quality assurance

- Cleaned the remaining spelling, URL and style diagnostics.

## sptrends 0.92.1

### Documentation

- Clarified that the default CMK region follows Neeti and Eastman
  (2011), as implemented in TerrSet’s Kendall module.

### Quality assurance

- Made the documentation-architecture test compatible across supported
  `testthat` versions.

## sptrends 0.92

### New features

- Added configurable odd CMK neighbourhoods while retaining the 3 by 3
  region described by Neeti and Eastman (2011), as implemented in
  TerrSet’s Kendall module, as the default.

## sptrends 0.91.9

### Documentation

- Completed and regression-tested the common help-page architecture
  across the public API.

## sptrends 0.91.8

### Documentation

- Standardised the detailed help across the remaining public functions,
  workflows and shared S3 methods without removing methodological
  content.

## sptrends 0.91.7

### Documentation

- Standardised the detailed help for prewhitening, trend testing, slope
  estimation and FDR correction while preserving its methodological
  content.

## sptrends 0.91.6

### Documentation

- Confirmed a clean local `devtools::check()` for 0.91.5 itself: 0
  errors, 0 warnings, 0 notes.

## sptrends 0.91.5

### Documentation

- Confirmed win-builder R-oldrelease also clean: 0 errors, 0 warnings, 1
  (expected) note. All three win-builder platforms now confirmed clean.

## sptrends 0.91.4

### Documentation

- Confirmed win-builder R-release also clean: 0 errors, 0 warnings, 1
  (expected) note.

## sptrends 0.91.3

### Documentation

- Confirmed the LaTeX/PDF manual fix against a fresh win-builder R-devel
  run: 0 errors, 0 warnings, 1 (expected) note.

## sptrends 0.91.2

### Bug fixes

- Fixed a real win-builder R-devel error: a literal Unicode alpha
  character in `fdr.R`/`workflow_trends.R` broke PDF manual generation.
  Replaced with `\eqn{\alpha}`/plain text.

## sptrends 0.91.1

### Quality assurance

- Covered local verbose, parallel-cleanup and Gi\* printing branches.

## sptrends 0.91

### Quality assurance

- Hardened documentation tests and completed a static release audit.

### Documentation

- Standardised spatial-correlation terminology and citation metadata.

## sptrends 0.90.4

### Quality assurance

- Replaced formatting-sensitive help-text tests with API checks.

## sptrends 0.90.3

### Quality assurance

- Restored the consecutive release history.

## sptrends 0.90.2

### Quality assurance

- Made small-raster tests explicitly planar.

## sptrends 0.90.1

### Spatial autocorrelation

- Corrected sparse-matrix row sums in local calculations.

## sptrends 0.90

### Spatial autocorrelation

- Added permutation-based local Moran’s I and Getis-Ord Gi\*.
- Added local maps and FDR forwarding through
  [`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md).

## sptrends 0.89.1

### Documentation

- Reorganised the introductory and workflow vignettes around loading,
  exploration, configurable analysis and the published TST/RTA
  workflows.
- Presented all prewhitening, slope and FDR alternatives with clearer
  guidance on trend preservation, robustness, multiplicity and
  computational performance.
- Standardised the pedagogical distinction between raw *p*-values,
  per-test *α* and target FDR *q*.

### Quality assurance

- Added documentation-contract tests and resolved the remaining source
  line-length findings.

## sptrends 0.89

### Simulation

- Separated spatial scale (`smooth_radius`) from spatial intensity
  (`spatial_rho`) in
  [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md),
  preserving the established default output.

### API

- Exposed every
  [`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md)
  metric in Usage while keeping FWER opt-in for replicated validation.
- Added dependency-free elapsed-time messages to verbose analytical
  functions; workflow objects continue to retain per-stage timings.

### Documentation

- Added concise `Typical use` paths to every exported function and
  distinguished single-run from replicated benchmarking.
- Clarified that the simulator supports controlled method evaluation,
  not calibrated environmental-process simulation.

## sptrends 0.88.4

### Quality assurance

- Corrected the DCF structure of `.lintr` so `lintr::lint_package()` and
  the lintr stage of `goodpractice::gp()` can read the project
  configuration.

## sptrends 0.88.3

### API

- Displayed every supported method in public function usage while
  preserving established defaults; workflow FDR arguments now validate
  BH, BKY and BY before computation.
- Extended the default ordered-stack file pattern to common raster
  formats supported by `terra`, while retaining custom regular
  expressions for other GDAL-readable formats.

### Documentation

- Shortened the README to a general package overview and added Luis V.
  García’s ORCID to package metadata.
- Distinguished local FDR control at `q` (BH, BKY and BY) from local
  FWER control at `alpha` (permutation maxT) in the reserved
  spatial-autocorrelation interface.
- Documented OLS as fastest, repeated median as substantially slower,
  and Theil-Sen as the recommended general balance of robustness and
  computational cost.

### Quality assurance

- Excluded Windows shortcut files from source builds and added
  regression tests for API choices, metadata, raster-format defaults and
  global/local multiple-testing validation.

## sptrends 0.88

### API

- Added an explicit global/local scope and scope-specific S3 hierarchy
  to
  [`spatial_autocorrelation()`](https://olive-r.github.io/sptrends/reference/spatial_autocorrelation.md),
  preserving the global result while reserving a stable interface for
  independently validated local statistics.
- Reserved local multiple-testing choices (`none`, BH, BKY, BY and
  permutation-based maxT); global analyses reject non-`none` adjustments
  because they contain only one test.

### Documentation

- Reframed
  [`spatial_autocorrelation()`](https://olive-r.github.io/sptrends/reference/spatial_autocorrelation.md)
  as a general spatial diagnostic for environmental variables,
  residuals, coefficients and inferential fields; FDR-assumption
  assessment is now presented only as one qualified application.
- Documented that local Moran and Getis-Ord Gi\* remain unavailable
  until their permutation inference and multiple-testing behaviour have
  been validated.

## sptrends 0.87.2

### API

- Exposed configurable trend workflows through the single public
  interface
  [`workflow_trends()`](https://olive-r.github.io/sptrends/reference/workflow_trends.md);
  the temporary pre-release alias and compatibility class were removed.

### Bug fixes

- Removed a hardcoded exact-version check from an unrelated ORCID test;
  it broke on every routine version bump.
- Kept local quality-assurance logs out of source builds and made the
  FDR comparison-plot example use a minimal direct input.
- Added targeted regression tests for empty FDR method sets, malformed
  or incompatible spatial neighbourhoods, and inconsistent
  configurable-workflow FDR objects; removed an unreachable OLS
  time-axis branch already enforced by the common validator.
- Completed BY support in configurable-workflow printing, direction maps
  and slope-significance plots.
- Defined the default CMK result for completely tied, zero-variance
  neighbourhoods as `Z = 0`, `p = 1`.
- Precomputed spatial neighbourhoods now carry and validate raster
  geometry, connectivity and the valid-cell pattern before reuse, with
  type-compatible dimension checks.
- Added consistent validation for cores, probabilities, seeds, pair
  limits, convergence controls and Durbin-Watson bounds.
- Prevented stage argument lists from replacing workflow-managed inputs,
  methods, reporting or shared-cluster controls.
- Preserved an initially absent global RNG state during internal
  Theil-Sen pair sampling and restored graphical parameters after
  null-distribution plots.
- Restored the internal sequential fallback for `n_cores <= 1` while
  retaining strict positive-integer validation in exported workflows.
- Changed the optional Moran diagnostic from an automatic BH/BKY
  selector to a qualified assessment: positive association is compatible
  with BH, whereas other outcomes are inconclusive; BKY is no longer
  presented as a safeguard against arbitrary dependence.
- Added a shared time-axis validator to reject mismatched, non-finite,
  duplicated or unordered dates before slope estimation, prewhitening,
  OLS or MMK.
- Added explicit validation for FDR levels and p-value inputs, including
  empty and all-`NA` cases.
- Added explicit validation for permutation counts, cores, seeds and
  precomputed spatial-neighbourhood dimensions.

### Documentation

- Documented that current analytical functions are vectorised and
  optionally parallel but operate in memory, rather than claiming
  out-of-core scalability.
- Added `inst/benchmarks/scalability.R` for reproducible timing and
  memory-proxy measurements.
- Corrected the documented number of prewhitening procedures.

## sptrends 0.86

### Bug fixes

- Corrected CMK’s regional variance when `ties = TRUE`: each cell now
  contributes its own tie-corrected `VarS`, and every cross-covariance
  uses the corresponding pair of variances. The default no-ties
  calculation is unchanged.

### Quality assurance

- Added an analytical 3 x 3 RAMK regression test with heterogeneous tie
  patterns and independently calculated `Sm`, `VarSm`, and `p`.
- Added legitimate spelling terms reported by
  `spelling::spell_check_package()` to `inst/WORDLIST`.
- The final checks for this patched source are pending; results obtained
  for 0.85.10 before this correction are not carried forward.

## sptrends 0.85

### Documentation

- Documented the substantially faster CMK execution observed against
  TerrSet on the same validation raster, without presenting the
  platform-dependent result as a universal benchmark.
- Updated documentation for slope estimators, FDR methods and
  prewhitening methods.
- Updated vignettes, navigation, author metadata and spelling
  exceptions.

## sptrends 0.84

### Documentation

- Added reproducible ConMK and TerrSet validation materials under
  `inst/validation/`.
- Documented CMK’s continuity convention and the limits of the TerrSet
  comparison.
- Added an optional `rkt` comparison for independent validation of CMK’s
  regional score aggregation.

## sptrends 0.83

### New features

- Added `continuity` to CMK: `FALSE` follows the published equations;
  `TRUE` provides compatibility with ConMK’s convention.

### Documentation

- Documented CMK as the moving-window raster form of analytical RAMK,
  including their shared cross-correlation correction.
- Clarified that prewhitening is an upstream operation, not part of the
  CMK statistic itself.
- Added method-specific quality-assurance sections across the principal
  trend workflow functions.

## sptrends 0.82

### Bug fixes

- Corrected `VCTFPW` to use the published variance-ratio correction and
  95% lag-1-autocorrelation gate. Cells below the gate remain unchanged,
  the first time step is retained, and a `Modified` diagnostic is
  returned.
- Updated VCTFPW workflow reporting to show the number and percentage of
  modified cells.
- Masked complete VCTFPW series when their rescaling factor is invalid.

## sptrends 0.81

### Bug fixes

- Corrected repeated-median inspection and slope-comparison intercepts.
- Returned the repeated-median intercept as a `SpatRaster`.
- Corrected the repeated-median intercept calculation.

## sptrends 0.80

### New features

- Added repeated-median slopes and three-estimator comparison to
  [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md).

### Breaking changes

- [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md)
  no longer draws a shaded confidence band; the interval is now
  text-only in the legend.

## sptrends 0.79

### New features

- Added Siegel’s repeated median (`"RM"`) to
  [`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md).

### Bug fixes

- Added `"RM"` validation to
  [`workflow_trends()`](https://olive-r.github.io/sptrends/reference/workflow_trends.md).

## sptrends 0.78

### New features

- Added
  [`workflow_trends()`](https://olive-r.github.io/sptrends/reference/workflow_trends.md)
  for user-defined analytical workflows.
- Added a dedicated vignette for configurable trend workflows.

## sptrends 0.77

### New features

- Added Benjamini-Yekutieli (`"BY"`) to
  [`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md).

### Bug fixes

- Added BY support to workflow printing and plotting.

## sptrends 0.76

### New features

- Added `MMK` to
  [`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md).

### Bug fixes

- Handled non-positive MMK variance-correction factors.

## sptrends 0.75

### New features

- Added `TFPW_Z` and `VCTFPW` to
  [`prewhiten()`](https://olive-r.github.io/sptrends/reference/prewhiten.md).

### Bug fixes

- Corrected AR(1) estimation in `prewhiten(method = "TFPW_WS")`.

## sptrends 0.74

### Breaking changes

- Removed `prewhiten(method = "PW")`.

### Bug fixes

- Corrected method citations, links and Moran’s-I guidance.

## sptrends 0.73

### Breaking changes

- Changed all method identifiers to uppercase.

### New features

- Standardised package-wide `method` values to uppercase.

## sptrends 0.72

### Bug fixes

- Removed publisher `url` fields from `inst/CITATION`/`CITATION.cff`,
  keeping only `doi`. Both publisher URLs returned 403 to CRAN’s own
  automated checker, a known behaviour toward bots, not a dead link.
- Fixed a test using `expect_s4_class()` on an S3 object.
- Closed the package’s last remaining zero-coverage line (a
  cluster-sharing seed-handling branch).
- Investigated a suspected missing `\donttest{}` wrap; confirmed a false
  alarm (an earlier grep had missed the closing brace).

## sptrends 0.71.0

### New features

- Shared spatial-adjacency caching between
  [`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md)
  and
  [`spatial_autocorrelation()`](https://olive-r.github.io/sptrends/reference/spatial_autocorrelation.md).
  Both previously recomputed the identical structure independently; now
  computed once and reusable via a new `precomputed_neighbourhood`
  argument.
- [`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)/[`workflow_rta()`](https://olive-r.github.io/sptrends/reference/workflow_rta.md)
  gain a shared `n_cores`, reusing one cluster across their own internal
  parallel steps instead of building and tearing down a separate one for
  each.

## sptrends 0.70

### New features

- [`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md)
  gains `"fwer"` (family-wise error rate) as a requestable metric, only
  meaningful with `replicates = TRUE`.

## sptrends 0.69.0

### Bug fixes

- Fixed `summary.compare_detections()`/`..._replicates()` picking the
  wrong “best” method for `FPR`/`FDR`: both used
  [`which.max()`](https://rdrr.io/r/base/which.min.html)
  unconditionally, but lower is better for these two error-rate metrics.

### Documentation

- Clarified the distinction between `FDR` (a single run’s realised
  proportion) and `FDR_mean` (the actual rate estimate, averaged across
  replicates), and added the Benjamini & Hochberg (1995) reference to
  [`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md)’s
  own citation list.

## sptrends 0.68

### New features

- Added `method = "ols"` to
  [`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md).
  Provides the classical parametric alternative to CMK/MK, verified
  against [`stats::lm()`](https://rdrr.io/r/stats/lm.html).

## sptrends 0.67

### New features

- Exposed
  [`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md)’s
  own direction map, histogram, and bar chart via `report = TRUE` and
  `plot(which = ...)`. These reports previously existed only inside the
  full
  [`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)/
  [`workflow_rta()`](https://olive-r.github.io/sptrends/reference/workflow_rta.md)
  pipeline, not from
  [`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md)
  directly.

## sptrends 0.66

### New features

- Added a dedicated slope-estimation vignette.

## sptrends 0.65

### Breaking changes

- Renamed
  [`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md)’s
  and
  [`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md)’s
  own returned S3 classes (`"cmk"` to `"trend_test"`, `"theil"` to
  `"slope"`). Reflects that both functions now cover more than one
  method each.

## sptrends 0.64

### API changes

- Renamed
  [`fdr_bky()`](https://olive-r.github.io/sptrends/reference/fdr_bky.md)’s
  `variant` argument to `implementation`, and its `"definition6"` value
  to `"original"`. The old names did not say what they were a
  variant/implementation of.

## sptrends 0.63

### API changes

- Renamed
  [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md)’s
  `type` argument to `selection_type`. Avoided ambiguity with
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html)’s
  own unrelated `type` argument used nearby.

## sptrends 0.62

### API changes

- Renamed 7 internal reporting functions and 5 source files (e.g.
  `contextual_mk.R` to `trend_test.R`, `theilsen.R` to
  `slope_estimator.R`) to match the generic, multi-method scope their
  own functions and classes already had after the renames above.

## sptrends 0.61

### Bug fixes

- Fixed a real bug affecting every diverging map the package draws:
  `fill_range = TRUE` was missing from 5
  [`terra::plot()`](https://rspatial.github.io/terra/reference/plot.html)
  calls, so the highest-magnitude cells rendered blank instead of
  strongly coloured.

## sptrends 0.60

### Bug fixes

- Fixed
  [`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)/[`workflow_rta()`](https://olive-r.github.io/sptrends/reference/workflow_rta.md)’s
  own progress messages always saying “Contextual Mann-Kendall”
  regardless of the actual `method` used.

## sptrends 0.59

### Bug fixes

- Fixed inconsistent citation formatting (author name, DOI, title) for
  one paper cited from three separate files.

## sptrends 0.58

### Bug fixes

- Closed all remaining test-coverage gaps found by `covr`, including 4
  regressions introduced by adding `method = "ols"`; the package reached
  100% coverage except one line-group left deliberately untested
  (unreachable without unsafely mocking a base R function).

## sptrends 0.57

### Documentation

- Ran `lintr::lint_package()` for the first time. Of ~985 findings,
  fixed the 2 genuine issues found and recorded the rest as deliberate
  package conventions (short variable names mirroring published
  notation, a fixed-width indentation style) in a new `.lintr` file, so
  future runs do not re-surface the same false positives.
- Reorganised the documentation of most core functions
  ([`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md),
  [`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md),
  [`prewhiten()`](https://olive-r.github.io/sptrends/reference/prewhiten.md),
  [`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md),
  [`spatial_autocorrelation()`](https://olive-r.github.io/sptrends/reference/spatial_autocorrelation.md),
  [`read_ordered_stack()`](https://olive-r.github.io/sptrends/reference/read_ordered_stack.md),
  [`compute_anomalies()`](https://olive-r.github.io/sptrends/reference/compute_anomalies.md),
  [`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md),
  [`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md),
  [`workflow_rta()`](https://olive-r.github.io/sptrends/reference/workflow_rta.md),
  [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md),
  and the three shared S3 generics) for a consistent structure and to
  lead with general purpose rather than the specific methods available
  today.

## sptrends 0.56

### Documentation

- Corrected
  [`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md)’s
  own title (still describing only two of its three methods) and added a
  missing citation (Douglas, Vogel & Kroll, 2000) for the RAMK
  connection behind CMK’s own variance correction.
- Audited the whole package for vignettes/docs describing only one
  method where more than one now exists; added missing coverage in three
  vignettes and
  [`?sptrends`](https://olive-r.github.io/sptrends/reference/sptrends-package.md)’s
  own overview.
- Finalised the diverging-colour convention (navy, not cyan, for
  “significant”) across all affected plots.

## sptrends 0.55.9

### Bug fixes

- Fixed 2 real `R CMD check` NOTEs from `CITATION.cff` living at the
  package root: excluded it via `.Rbuildignore` (it is meant for
  GitHub/Zenodo, not R’s own
  [`citation()`](https://rdrr.io/r/utils/citation.html) mechanism, and
  does not need to be part of the built package).
- Wrapped
  [`compute_anomalies()`](https://olive-r.github.io/sptrends/reference/compute_anomalies.md)’s
  own example in `\donttest{}` (6.17s, over the 5-second threshold).

## sptrends 0.55.8

### Documentation

- Added `CITATION.cff` at the package root, with both authors, both
  published references, and the licence.

### Bug fixes

- Fixed `inst/CITATION` still referring to the old `tst()`/`rta()`
  function names after they were renamed.

## sptrends 0.55.7

### New features

- Added `verbose = TRUE` messaging to
  [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md),
  matching every other function’s own convention.
- Folded `summarise_replicates()` into
  [`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md)’s
  own `replicates = TRUE` argument. One function instead of two for the
  same task.
- `tst()` and `rta()` renamed to
  [`workflow_tst()`](https://olive-r.github.io/sptrends/reference/workflow_tst.md)
  and
  [`workflow_rta()`](https://olive-r.github.io/sptrends/reference/workflow_rta.md).
  Clarifies these are complete, opinionated workflows, not the only way
  to run each analytical step.
- [`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md),
  [`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md),
  and
  [`prewhiten()`](https://olive-r.github.io/sptrends/reference/prewhiten.md)
  now fail fast with an informative error on invalid input shapes,
  instead of a cryptic error partway through computation.
- [`spatial_autocorrelation()`](https://olive-r.github.io/sptrends/reference/spatial_autocorrelation.md)
  now errors informatively on two previously-confusing invalid input
  cases.
- Established a consistent visual identity (colours, logo, hex sticker)
  across the package, README, and pkgdown site.

### Bug fixes

- Fixed two real bugs in
  [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md),
  including a division-by-zero, found by new edge-case tests written for
  it,
  [`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md),
  and
  [`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md).
- Fixed a wrong DOI, caught by `devtools::check(remote = TRUE)`.
- Fixed a real `R CMD check` WARNING and a separate ERROR, both in
  [`compute_anomalies()`](https://olive-r.github.io/sptrends/reference/compute_anomalies.md)’s
  own documentation.
- Fixed 2 real test failures and a bug from the previous release’s
  `verbose =` addition, found running `test-external-validation.R` (a
  new suite comparing sptrends’ own core statistics against independent
  reference implementations).
- Removed a broken test added in the previous release, and fixed a
  broken `@examples` block and a malformed `@section` title.
- Fixed a visible numbering gap in the vignette titles, and used the
  wrong logo file from the previous release.

### Documentation

- Package title changed to remove “for”, and `DESCRIPTION`’s own
  description no longer repeats the package name.
- Added a runnable simulated-data example alongside the existing
  real-data one for 5 core functions, and wrapped several long-running
  examples in `\donttest{}`.
- New conceptual pipeline figure and a revised 5-step workflow diagram
  (from 6).
- New vignette covering the full workflow end-to-end, filling a real
  content gap; a separate RTA-only vignette folded into it.
- Reviewed the package in full against `goodpractice::gp()` and
  `lintr::lint_package()`; addressed genuine findings.

### Internal changes

- A large `NEWS.md` formatting inconsistency found and corrected across
  many earlier entries.
- Established a consistent brand colour palette, replacing an
  approximate logo blue used in the previous release.
- Closed several coverage gaps identified individually via
  `covr::zero_coverage()`, and removed a small amount of dead code in
  `prewhiten.R` found during that same review.
- 15 compound-semicolon statements split into separate lines for
  readability.

## sptrends 0.54

### Internal changes

- Minor internal maintenance and test upkeep.

## sptrends 0.53

### Internal changes

- Minor internal maintenance and test upkeep.

## sptrends 0.52

### Bug fixes

- Fixed 4 real test failures (in test code, not package behaviour):
  incorrect
  [`terra::values()`](https://rspatial.github.io/terra/reference/values.html)
  matrix-mode usage in two tests, a genuine bug in a hand-computed
  Theil-Sen reference implementation, and a named-vs-unnamed value
  mismatch in `expect_equal()`.

## sptrends 0.51

### Bug fixes

- Fixed, while completing a rename, a genuinely broken previous state:
  the earlier `moran_permutation_test()` to
  [`spatial_autocorrelation()`](https://olive-r.github.io/sptrends/reference/spatial_autocorrelation.md)
  rename had been left half-done – neither the old nor the new name was
  actually callable from outside the package, and
  [`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md)’s
  own `moran_check = TRUE` silently called the missing old name
  internally. Fixed throughout.

## sptrends 0.50.0

### New features

- `moran_permutation_test()` generalised into
  [`spatial_autocorrelation()`](https://olive-r.github.io/sptrends/reference/spatial_autocorrelation.md),
  gaining a second method (`method = "getis_ord"`, global Getis-Ord
  General G) alongside the existing Moran’s I. Follows the same
  `method =` pattern already used elsewhere in the package; verified
  numerically against the direct definition before being added. Returned
  field names changed to be generic (`$statistic`/ `$null_dist`, were
  `$I`/`$I_null`) – **breaking change for code reading the old field
  names.**

## sptrends 0.49.0

### New features

- [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)
  gains structural-break simulation (`break_type = "mean"`/`"slope"`),
  composing independently from the existing monotonic-trend simulation,
  with its own ground truth (`true_break`/`break_time`) for validating a
  future change-point method – none is implemented yet.
- [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md)
  gains `slope_method = "ols"`, matching
  [`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md)’s
  own second method; the confidence band’s “bowtie” shape was checked
  and confirmed correct (standard linear-trend theory), not a bug.
- [`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md)
  gains `method = "ols"`, a fast closed-form alternative to the default
  Theil-Sen, not a replacement for it (not robust to outliers the way
  Theil-Sen is).
- [`prewhiten()`](https://olive-r.github.io/sptrends/reference/prewhiten.md)
  gains `method = "yue_pilon"` (trend-free prewhitening) alongside the
  default `wang_swail`. Loses one time step (classic lag-1 differencing)
  and has a different diagnostics structure. Fixed three real
  integration bugs found while adding it –
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)/[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  `tst()`’s own [`print()`](https://rdrr.io/r/base/print.html), and
  [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md)
  all hard-assumed the other method’s own structure and would have
  errored on a `yue_pilon` result.

### Improvements

- [`read_ordered_stack()`](https://olive-r.github.io/sptrends/reference/read_ordered_stack.md)
  now stores detected years as real
  [`terra::time()`](https://rspatial.github.io/terra/reference/time.html)
  metadata;
  [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md)’s
  x-axis uses them automatically. Found and fixed a related bug while
  doing this: the prewhitened panel misaligned `t` by one step for
  `yue_pilon`’s own shorter output, silently, with no visible symptom.
- Fixed 3 real test failures from the previous release’s rename, missed
  because the old parameter name was passed through a forwarding
  argument list rather than called directly.

### Bug fixes

- Fixed a critical bug: `terra::values(x, mat = FALSE)[, 1]` errors on a
  single-layer raster (already a flat vector, not a matrix) – affected
  `prewhiten(method = "yue_pilon")` entirely (every call errored) and
  two `"ols"` tests.

### Internal changes

- [`direction_map()`](https://olive-r.github.io/sptrends/reference/direction_map.md)
  made internal – the FDR-corrected and uncorrected direction are the
  same computation with an optional filter, not two things needing
  separate public entry points. **Breaking: no longer exported**; use
  `plot(x, which = "direction")` instead.

## sptrends 0.48

### Internal changes

- Reviewed a `goodpractice::gp()` report; one genuine fix (a missing
  `add = TRUE` on an [`on.exit()`](https://rdrr.io/r/base/on.exit.html)
  call), the rest confirmed as tool limitations.

## sptrends 0.47.0

### Internal changes

- API redesign stage 1: `contextual_mann_kendall()`,
  `theil_sen_slope()`, and `wang_swail_prewhiten()` renamed to
  [`trend_test()`](https://olive-r.github.io/sptrends/reference/trend_test.md),
  [`slope_estimator()`](https://olive-r.github.io/sptrends/reference/slope_estimator.md),
  and
  [`prewhiten()`](https://olive-r.github.io/sptrends/reference/prewhiten.md),
  each gaining a `method =` argument ahead of a planned second method
  per function. **Breaking: update any direct calls to the old names**,
  and `neighbourhood = TRUE/FALSE` to `method = "cmk"`/`"mk"`. Fixed a
  real naming collision caught during the rename (`tst()`’s own
  `prewhiten` flag argument shadowed the new function of the same name).

## sptrends 0.46.1

### Internal changes

- Closed one more coverage gap; remaining `read_stack.R` gaps confirmed
  as unreachable defensive branches or an `ncdf4`-only test.

## sptrends 0.46.0

### New features

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) on
  `"tst"`/`"rta"` results gains 8 new `which` views
  (`slope_map`/`slope_direction`/`slope_hist`/`slope_bar`,
  `pvalue_map`/`pvalue_significance`/`pvalue_hist`/`pvalue_bar`), all
  uncorrected diagnostics, not a final result. Considered and rejected a
  bigger redesign (separate generics) first, since it would not
  generalise cleanly across every classed object type this package
  returns.

## sptrends 0.45.1

### Improvements

- `workflow_summary()`/`summary(x, which = "workflow")`, added in the
  same release, removed again immediately on reflection – it competed
  with [`summary()`](https://rdrr.io/r/base/summary.html)’s own simpler
  default rather than sitting alongside it. **No replacement**;
  `x$timing` has the same data.
- [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md)
  gains `show_neighbours`: a small-multiples grid of each aggregated
  neighbour cell’s own fit, drawn in a separate window, answering
  whether the clicked cell’s trend is representative of its
  neighbourhood.

### Internal changes

- Closed several coverage gaps found via `covr::zero_coverage()`;
  documented two `read_stack.R` branches as genuinely unreachable rather
  than chased with contrived tests.
- `method_citation()` removed entirely – redundant with each function’s
  own `@references`; fixed a real bug this surfaced along the way
  (duplicated, driftable citation text in `workflow_summary()`). **No
  replacement**; see each function’s own help page.
- [`classify_moran()`](https://olive-r.github.io/sptrends/reference/classify_moran.md)
  folded into
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)
  of a `"moran"` object; no longer exported separately.

## sptrends 0.44.0

### New features

- `workflow_summary()` folded into `summary(x, which = "workflow")` for
  `tst()`/`rta()` results (later removed again in 0.45.1).

### Internal changes

- `fdr_direction_map()` renamed to
  [`direction_map()`](https://olive-r.github.io/sptrends/reference/direction_map.md)
  and moved out of `fdr.R` into its own file – it isn’t itself an FDR
  method, it combines a trend test’s direction with an FDR result.
  **Breaking: rename any direct calls.**
- [`prepare_cmk_neighbourhood()`](https://olive-r.github.io/sptrends/reference/prepare_cmk_neighbourhood.md)
  no longer exported – a real, acknowledged capability loss for
  batch-processing many rasters sharing one grid geometry (no longer
  precomputable without `:::`).
- [`classify_moran()`](https://olive-r.github.io/sptrends/reference/classify_moran.md)
  no longer exported; its category now shown automatically in
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)
  of a `"moran"` object.
- Exported function count: 23 to 20.

### Documentation

- [`prepare_cmk_neighbourhood()`](https://olive-r.github.io/sptrends/reference/prepare_cmk_neighbourhood.md)’s
  batch-processing use case now demonstrated in examples and a vignette,
  not only described.

## sptrends 0.43.1

### Bug fixes

- Fixed a real test failure: `test-sptrends-methods.R` checked S3
  methods via
  [`getNamespaceExports()`](https://rdrr.io/r/base/ns-reflect.html),
  which never reflects S3 methods registered the standard way – fixed to
  use [`getS3method()`](https://rdrr.io/r/utils/getS3method.html), the
  actually-correct check.

## sptrends 0.43.0

### Improvements

- [`print.sptrends()`](https://olive-r.github.io/sptrends/reference/print.sptrends.md)/[`summary.sptrends()`](https://olive-r.github.io/sptrends/reference/summary.sptrends.md)/[`plot.sptrends()`](https://olive-r.github.io/sptrends/reference/plot.sptrends.md)
  now dispatch by naming convention (`.print_<class>()` lookup) instead
  of a hardcoded [`switch()`](https://rdrr.io/r/base/switch.html) –
  adding a future class no longer requires editing all three shared
  generics.
- API redesign: 18 separate
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)/[`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  methods (one per class) consolidated into the 3 shared generics above,
  mirroring `terra`’s own single-entry-point convention. **Breaking:
  code calling one of the 18 old method names directly** (not via
  generic dispatch) breaks; ordinary `print(x)`/`summary(x)`/`plot(x)`
  usage is unaffected. Exported function count: 43 to 28.

### Documentation

- [`fdr_bh()`](https://olive-r.github.io/sptrends/reference/fdr_bh.md)/[`fdr_bky()`](https://olive-r.github.io/sptrends/reference/fdr_bky.md)
  no longer exported (use `fdr_correction(method = ...)`);
  `moran_permutation_test()` now returns a classed object with its own
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)/[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  folding in two previously separate reporting functions;
  [`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md)
  likewise gains a class and a genuinely new
  [`summary()`](https://rdrr.io/r/base/summary.html). Exported function
  count: 28 to 23.

## sptrends 0.42.1

### Bug fixes

- Fixed 2 real `R CMD check` example errors: 15 internal functions
  called themselves by bare name in their own `@examples`, invisible in
  `devtools::test()` but not under `R CMD check`’s external-user
  context. Prefixed with `sptrends:::`.

## sptrends 0.42.0

### Improvements

- API redesign stages 4-5 of 5: 15 standalone reporting functions
  ([`trend_maps()`](https://olive-r.github.io/sptrends/reference/trend_maps.md),
  [`fdr_summary()`](https://olive-r.github.io/sptrends/reference/fdr_summary.md),
  and others) no longer exported – fully covered by the
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)/[`plot()`](https://rdrr.io/r/graphics/plot.default.html)
  methods added in earlier stages. **Breaking: switch to the equivalent
  S3 method** (see each function’s own docs for the exact mapping).
  Exported function count: 58 to 43.

## sptrends 0.41.1

### Bug fixes

- Fixed 4 real test failures missed in the previous stage’s own
  verification: a second test file also passed a classed result directly
  where the raw statistics raster was expected.

## sptrends 0.41.0

### Improvements

- API redesign stage 3 of 5: `contextual_mann_kendall()` now returns a
  classed object (`list(stats = <raster>, neighbourhood = <logical>)`)
  instead of a bare raster, with its own
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)/[`plot()`](https://rdrr.io/r/graphics/plot.default.html).
  **Breaking: update `trend$p` to `trend$stats$p`.**

## sptrends 0.40.0

### Improvements

- API redesign stages 1-2 of 5: `wang_swail_prewhiten()`,
  [`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md),
  and `theil_sen_slope()` now return classed objects with their own
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)/[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  matching `tst()`/`rta()`’s own pattern. `theil_sen_slope()`’s own
  change is breaking (was a bare raster, now
  `list(slope = ..., smoothed = ...)`). A deliberate 5-stage redesign,
  done before the associated papers are published rather than after.

## sptrends 0.39.2

### Documentation

- `R/` and all vignettes brought under the 80-character line-length
  limit (59 lines fixed), verified line by line, not just visually
  diffed – caught a real
  [`sprintf()`](https://rdrr.io/r/base/sprintf.html) mistake introduced
  during the cleanup itself before it reached the codebase. Test files
  deliberately left alone (developer-facing only, high volume).

## sptrends 0.39.1

### Internal changes

- `goodpractice::gp()` cleanup, low-risk items only: `stop(paste0())`
  replaced with [`stop()`](https://rdrr.io/r/base/stop.html)’s own
  concatenation, `fixed = TRUE` added to plain-string
  [`grepl()`](https://rdrr.io/r/base/grep.html) calls, several test
  assertions replaced with more specific equivalents, one real broken
  URL fixed in `README.md`.

## sptrends 0.39.0

### New features

- New “Preprocessing” pkgdown category, separating
  `wang_swail_prewhiten()`/[`compute_anomalies()`](https://olive-r.github.io/sptrends/reference/compute_anomalies.md)
  from “Core”/“Support” – documentation/organisation only, no
  behavioural change.

### Improvements

- **Breaking**: `tst()` no longer smooths its Theil-Sen slope by
  default, matching `rta()`’s own default and removing an unprincipled
  asymmetry (neither published method includes this smoothing). Pass
  `theil_sen_args = list(smooth_neighbourhood = TRUE)` to keep the old
  behaviour.

### Bug fixes

- Fixed a real, visually misleading bug found by inspecting an actual
  plot: queen-3x3 smoothing (`smooth_neighbourhood = TRUE`, and both
  `plot(which = "slope")` methods) used
  [`terra::focal()`](https://rspatial.github.io/terra/reference/focal.html)’s
  own default `na.policy = "all"`, letting a non-significant or no-data
  cell get painted with a neighbour’s colour/value. Fixed to `"omit"` in
  all three call sites, factored into one shared helper.

### Documentation

- Fixed 13 occurrences of a stale “one of the five pillars” phrase
  across 9 files, now inaccurate after prewhitening’s own
  reclassification.

## sptrends 0.38.1

### Bug fixes

- Fixed a test-suite bug: `expect_s3_class()` does not accept an `info`
  argument (unlike `expect_error()`); removed it from two new sweep
  tests.

## sptrends 0.38.0

### Improvements

- `rta()` gained `summary.rta()`/`plot.rta()`, closing an asymmetry with
  `tst()`.
- Added two systematic combinatorial sweep tests (`tst()` across 12
  argument combinations, `rta()` across 4).

### Internal changes

- Closed remaining coverage gaps found via `covr::zero_coverage()`
  across several reporting functions.

## sptrends 0.37.2

### Internal changes

- Covered
  [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md)’s
  interactive dispatch, reversing an earlier “this is the coverage
  ceiling” conclusion, by mocking
  [`terra::click()`](https://rspatial.github.io/terra/reference/click.html)/[`terra::draw()`](https://rspatial.github.io/terra/reference/draw.html)
  the same way `detectCores()` was already mocked elsewhere.

## sptrends 0.37.1

### Improvements

- Fixed non-ASCII characters in `R/method_citation.R`’s own string
  literals (accented names), caught by `R CMD check`’s portability check
  – replaced with Unicode escapes.

### Bug fixes

- Fixed a real bug: `tst()`/`rta()` called
  [`trend_summary()`](https://olive-r.github.io/sptrends/reference/trend_summary.md)
  a second, unconditional time with no way to suppress its own printed
  output, so `report = FALSE` never achieved full silence and
  `report = TRUE` printed the summary twice. Added `verbose =` to
  [`trend_summary()`](https://olive-r.github.io/sptrends/reference/trend_summary.md).

## sptrends 0.37.0

No new statistical methodology in this release – architecture, API and
documentation work only.

### New features

- `tst()`/`rta()` now record per-step timing in a new `$timing` list.
- New `workflow_summary()`: a Step/Method/Reference/Time table for a
  `tst()`/`rta()` result (later folded into
  [`summary()`](https://rdrr.io/r/base/summary.html), then removed – see
  later entries).
- New `method_citation()`: looks up a citation by short method name
  (later removed as redundant – see 0.45.1).
- New `theil_sen_summary()`/`theil_sen_map()`, and `theil_sen_slope()`
  gained `report = TRUE`, closing the one core function still missing
  companion reporting functions. Found and fixed a real bug while wiring
  this in:
  [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md)’s
  own internal call did not set `report = FALSE`, so every call would
  have started auto-printing a spurious one-cell summary/map.

### Improvements

- Added a “Methodological background” section to the top of each core
  method’s documentation, for a scannable orientation before the
  argument-level detail.
- `moran_permutation_test()` gained `report = TRUE` (closing the one
  remaining function requiring a separate manual call to see any report)
  and proper [`match.arg()`](https://rdrr.io/r/base/match.arg.html)
  validation on `connectivity` (previously failed silently deep inside
  [`terra::adjacent()`](https://rspatial.github.io/terra/reference/adjacent.html)
  on a typo).
- [`prewhiten_summary()`](https://olive-r.github.io/sptrends/reference/prewhiten_summary.md)
  now prints an informative narrative summary, matching the other
  `*_summary()` functions.

### Internal changes

- `tst()`/`rta()` results gain a shared `"sptrends"` superclass,
  non-breaking, for future common methods to dispatch on.

### Documentation

- Reorganised the pkgdown reference index into four blocks by user
  intent (“Published workflows”, “Core methods”, “Diagnostics”,
  “Utilities”) rather than internal function type.
- [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md)
  added to `README.md` for the first time.
- Rebalanced TST/RTA presentation further, per explicit feedback that an
  earlier pass still read as TST-first: rewrote the README opening,
  diagram, and citation section to present both as equal, parallel
  workflows; rewrote `rta()`’s own documentation to mirror `tst()`’s in
  full; added a new vignette walking through the complete RTA workflow,
  closing the one remaining asymmetry.

## sptrends 0.36.1

### New features

- Documentation review following `rta()`’s introduction:
  [`?sptrends`](https://olive-r.github.io/sptrends/reference/sptrends-package.md),
  `README.md`, `inst/CITATION`, and four vignettes updated to present
  TST and RTA as two workflows, not TST alone (`rta()` had been missing
  from several places despite already being complete and tested).
  `theil_sen_slope()`/[`fdr_bh()`](https://olive-r.github.io/sptrends/reference/fdr_bh.md)
  gained a “used by” note that both are shared infrastructure, not
  TST-specific.

## sptrends 0.36.0

### Improvements

- Completed `tst()`’s and `rta()`’s own `@references`: both previously
  cited only their own workflow paper, not the individual published
  methods (Wang & Swail, Neeti & Eastman, Theil, Sen, Benjamini and
  co-authors) each step builds on.

### Bug fixes

- Fixed a real, visually misleading bug:
  [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md)’s
  confidence band was anchored at `t = 0` (typically outside the
  observed range), widening in only one direction instead of the
  standard symmetric “bowtie” shape. Anchored at the centre of the
  observed range instead.

## sptrends 0.35.0

### New features

- New function `rta()`: a second, complete, published workflow (Robust
  Trend Analysis, Gutiérrez-Hernández & García 2024) alongside `tst()`,
  not replacing it – Theil-Sen, Contextual Mann-Kendall, and standard
  FDR-BH (not adaptive BKY), deliberately without prewhitening. Both
  differences documented at length rather than left implicit, since
  prewhitening’s own value before a Mann-Kendall test is a genuine,
  unresolved debate in the published literature (cited both sides). No
  Moran’s I step and no categorical slope classes, unlike the original
  published RTA method – kept consistent with how this package treats
  magnitude and diagnostics elsewhere. New `print.rta()`. Fixed an
  inaccurate existing cross-reference: `tst()`’s own docs had listed the
  RTA paper as “the same workflow”, no longer true now that RTA is
  separately implemented.
- Added a mosaic view and an optional year-by-year animation to the
  getting-started vignette’s data-loading section.

### Bug fixes

- Fixed a real test-suite bug: a shared test helper lived inside one
  test file, invisible to another that sorts earlier alphabetically.
  Moved to a `helper-sptrends.R`, always sourced first.

## sptrends 0.34.0

### New features

- Redesigned
  [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md):
  a single `neighbourhood` argument (default `TRUE`) replaces always
  overlaying two Theil-Sen lines. All aggregation modes now take the
  per-time-step median of raw values first, then fit one Theil-Sen slope
  to that aggregated series – deliberately different from the previous
  “estimate per cell, then median the slopes” order, since only the
  “aggregate first” order has a standard confidence-interval formula. As
  a consequence, `inspect_ts_cell(neighbourhood = TRUE)` no longer
  numerically matches
  `theil_sen_slope(..., smooth_neighbourhood = TRUE)` at the same cell –
  documented as a genuine, intended difference.
- A confidence interval is now always shown, using the standard
  Sen/Gilbert rank-based method; not offered for the old
  “median-of-slopes” mechanism, since no standard formula exists for it.

### Improvements

- [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md)
  gained an optional `prewhitened` argument, drawing a second panel
  using `wang_swail_prewhiten()`’s own output side by side with the raw
  one. Internally refactored into small, reusable helpers shared by both
  panels.
- Completed
  [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md)’s
  own references, and found and fixed a real gap while doing so:
  `contextual_mann_kendall()` itself had never explicitly cited Mann
  (1945)/Kendall (1975), the foundational statistic its own `S`/`VarS`
  computation rests on.

## sptrends 0.33.0

### New features

- New function
  [`inspect_ts_cell()`](https://olive-r.github.io/sptrends/reference/inspect_ts_cell.md):
  interactive click-to-inspect time series viewer, overlaying the
  clicked cell’s own Theil-Sen fit against one that borrows its queen
  neighbourhood – whether the two agree is itself the diagnostic. Split
  into a public interactive wrapper and an internal, fully testable core
  (the interactive part cannot be exercised by automated tests).
  Verified by direct comparison to match `theil_sen_slope()`’s own
  output exactly for the same cell.

## sptrends 0.32.0

### Improvements

- Replaced the bundled example dataset: annual mean temperature (CRU TS)
  out, annual mean NDVI (NOAA STAR Blended-VHP, 1982-2023) in –
  vegetation greening/browning is the motivating application of the TST
  method itself, so the bundled example now matches the package’s own
  subject matter. Updated everywhere the old dataset was referenced.

## sptrends 0.31.0

### Improvements

- Covered the `n_cores`-capping message safely for `R CMD check`, by
  mocking
  [`parallel::detectCores()`](https://rdrr.io/r/parallel/detectCores.html)
  rather than actually requesting more cores than the machine has.

## sptrends 0.30.1

### Improvements

- Simplified `plot.tst(which = "slope")`’s map title, dropping
  smoothing/neighbourhood wording from the displayed label.

### Internal changes

- Final coverage pass after a large refactor, closing several gaps
  identified by exact `covr::zero_coverage()` line numbers; confirmed
  (not forced) that the `n_cores`-capping message and two `read_stack.R`
  lines remain deliberately/structurally unreachable.

## sptrends 0.30.0

### Internal changes

- Finished the line-length (\>80 char) cleanup started in 0.29.1; from
  235 offending lines across 14 files down to 35 legitimate exceptions
  (a few-character string overruns, published reference tables, the
  package’s own title).

## sptrends 0.29.1

### Internal changes

- Started a line-length (\>80 char) cleanup across the codebase.

## sptrends 0.29.0

### New features

- New `variant` argument for
  [`fdr_bky()`](https://olive-r.github.io/sptrends/reference/fdr_bky.md)
  (`"multtest"`/`"definition6"`, two distinct published implementations
  of the same BKY procedure), forwarded through
  [`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md)/`tst()`.
  Confirmed `"multtest"`’s threshold is never stricter than
  `"definition6"`’s.

## sptrends 0.28.0

### Improvements

- Simplified `d-fdr-correction`’s narrative: Moran’s I check reframed as
  an optional diagnostic, not part of the main flow; now shows both
  FDR-BH and FDR-BKY direction maps side by side.

### Documentation

- Changed the FDR-corrected significance colour from royal blue to navy,
  and the general diverging palette from “RdBu” to “Blue-Red 3” (both
  [`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html)),
  across every function and vignette using it. Uncorrected significance
  stayed red deliberately, to reinforce it is diagnostic-only.
- `tst()`’s own examples now also show
  `plot(which = "significance"/ "slope")`, not only the default view.

## sptrends 0.27.1

### Improvements

- Fixed a test comparing raster values directly: the manually-built
  comparison object had a different auto-generated layer name, failing
  `expect_equal()` on the name alone, not the values. Fixed with
  [`unname()`](https://rdrr.io/r/base/unname.html).

## sptrends 0.27.0

### Documentation

- Restructured the vignette narrative into 5 steps: the getting-started
  vignette is now a true “hello world” (load and visualise only); the
  full TST demonstration moved into its own new vignette, after
  prewhitening/trend-test/FDR are each covered alone.

## sptrends 0.26.1

### Internal changes

- Fixed 3 real issues surfaced by `devtools::check()`: 7 tests used bare
  `values()` instead of
  [`terra::values()`](https://rspatial.github.io/terra/reference/values.html);
  a real portability bug in `moran_null_plot()` (Unicode arrow
  characters failing to convert on some locales, replaced with plain
  ASCII); and a
  [`max()`](https://rdrr.io/r/base/Extremes.html)-with-no-arguments
  warning in two all-NA guards, suppressed.

## sptrends 0.26.0

### New features

- New `plot.tst(which = "slope")` view: the Theil-Sen slope, masked to
  FDR-significant cells, smoothed by default at plot time only.
  Motivated by real data: comparing CMK’s neighbourhood-averaged
  statistic sign against Theil-Sen’s own sign on the bundled dataset
  showed 338 disagreeing cells, always pulled towards the locally
  dominant trend – this view deliberately uses Theil-Sen’s own sign.
- New `smooth_neighbourhood` argument for `theil_sen_slope()` (default
  `FALSE`): queen-3x3 median smoothing of the estimated slope, with an
  explicit caveat that – unlike CMK’s own published, validated method –
  this has no equivalent literature backing, since blending neighbouring
  cells’ magnitudes assumes they share a similar true slope, which does
  not always hold.

### Improvements

- **`tst()`’s own Theil-Sen step now smooths by default**, a deliberate
  methodological choice despite the caveats above; calling
  `theil_sen_slope()` directly is unaffected. Override via
  `theil_sen_args = list(smooth_neighbourhood = FALSE)`.
- Changed the package title to emphasise formal statistical inference as
  the actual differentiator.

### Bug fixes

- Fixed a real, potentially misleading plotting bug:
  [`trend_maps()`](https://olive-r.github.io/sptrends/reference/trend_maps.md)/
  [`prewhiten_maps()`](https://olive-r.github.io/sptrends/reference/prewhiten_maps.md)
  did not fix their diverging colour range symmetrically around zero,
  visually misrepresenting the sign boundary whenever positive and
  negative values were not equally extreme.

### Documentation

- Added a “Roadmap” section to `README.md`.

## sptrends 0.25.0

### Internal changes

- Final targeted coverage pass across the whole package (46 lines),
  using exact `covr::zero_coverage()` results. The single biggest
  recurring pattern: the `path =` (write-to-disk) argument had never
  been tested for most reporting functions – added across 10 of them.
  Closed several other real per-function gaps. Deliberately left
  uncovered: the `n_cores`-capping message and two structurally
  unreachable `read_stack.R` checks.

## sptrends 0.24.0

### Internal changes

- Used exact `covr::zero_coverage()` line numbers to close genuine gaps
  in `prewhiten.R` (Durbin-Watson out-of-range messages, several verbose
  branches, `path =` untested for two reporting functions) and
  `read_stack.R`
  ([`read_netcdf_stack()`](https://olive-r.github.io/sptrends/reference/read_netcdf_stack.md)’s
  own messages, never exercised since every existing test used
  `verbose = FALSE`). Noted, not forced: two checks appear structurally
  unreachable given existing upstream guards.

## sptrends 0.23.1

### Improvements

- Fixed a test based on a misreading of `wang_swail_prewhiten()`’s own
  code: the clamp only catches `rho` at or beyond exactly `+-1`, not a
  general cap at 0.99 – a legitimate `0.995`-ish estimate correctly
  failed the test’s own wrong assertion. Rewritten to check
  `abs(rho) < 1`, what the code actually guarantees.

## sptrends 0.23.0

### Internal changes

- More coverage for the remaining ~86%-covered files (`tst-methods.R`,
  `read_stack.R`, `prewhiten.R`); flagged one added test’s own assertion
  as weak, honestly, since execution could not confirm it actually drove
  the estimate past the clamping boundary.
- More coverage for `fdr.R` (worst-covered file at 84.65%): out-of-
  range input error, multi-layer input message, full `report = TRUE`
  branch, Moran’s I recommendation message, CSV-writing `path`.

## sptrends 0.22.1

### Bug fixes

- Fixed a real bug in the parallel chunk-splitting logic of
  `theil_sen_slope()`/`contextual_mann_kendall()`:
  [`cut()`](https://rdrr.io/r/base/cut.html) errors with only a single
  item to split, which happens at exactly one valid cell or two time
  layers. Fixed by switching to
  [`parallel::splitIndices()`](https://rdrr.io/r/parallel/splitIndices.html).

## sptrends 0.22.0

### Internal changes

- Further coverage for `theilsen.R`/`contextual_mk.R`: verbose messages
  previously executed but never asserted, and `n_cores` exceeding the
  number of valid cells without erroring.
- Found and fixed why `utils-internal.R`’s coverage was stuck at 58.5%:
  `covr::package_coverage()` does not set `NOT_CRAN`, so every
  `skip_on_cran()`-guarded test was silently skipped during coverage
  measurement. Nearly every parallel test used `n_cores` of 1-2 (already
  CRAN-compliant), so the guard was unnecessary caution – removed from
  all 9 affected tests.

## sptrends 0.21.1

### New features

- Fixed another type-mismatch test failure:
  [`stats::median()`](https://rdrr.io/r/stats/median.html) of an
  odd-length integer vector returns an integer, not a double, unlike the
  even-length case.

## sptrends 0.21.0

### Internal changes

- Continued the coverage sweep: `contextual_mk.R`’s full `report = TRUE`
  branch (previously only its own constituent functions were tested in
  isolation); `theilsen.R`’s subsampling branch
  (`n*(n-1)/2 > max_pairs`), previously completely untested;
  `utils-internal.R`’s remaining parallel-lapply arguments (`seed`,
  `packages`); `read_stack.R`’s NetCDF error path and diagnostic plot
  path.

## sptrends 0.20.1

### Improvements

- Removed a test requesting more cores than available to verify a
  capping message: `R CMD check` enforces its own separate hard limit on
  parallel workers, so the test errored before the function’s own logic
  ever ran.

## sptrends 0.20.0

### Internal changes

- Continued the coverage sweep on the worst-covered files:
  `utils-internal.R` (31.7%, every internal helper tested directly),
  `pipeline.R`/`tst()` (71.8%, the flagship function, largely untouched
  by earlier rounds), `theilsen.R` (76.1%, including a hand-worked exact
  reference value, not just a sanity check).

## sptrends 0.19.1

### Bug fixes

- Fixed 13 test failures from the previous round’s new tests: several
  genuine type mismatches
  ([`terra::values()`](https://rspatial.github.io/terra/reference/values.html)
  on a single layer returns a vector not a matrix;
  [`terra::nlyr()`](https://rspatial.github.io/terra/reference/dimensions.html)
  returns a double, not integer), two malformed test fixtures (a
  file-name collision, an NetCDF write missing `overwrite = TRUE`), and
  one over-strict `expect_silent()` where warnings are legitimately
  expected.

## sptrends 0.19.0

### Internal changes

- Added substantial test suites for the three remaining lowest-coverage
  core files: `prewhiten.R` (gating behaviour, `dw_method`, error
  conditions), `fdr.R` (raster input, all reporting functions), and
  `contextual_mk.R` (hand-worked exact `S`/variance values, isolated-
  cell fallback, connectivity comparison).
- Added a test suite for
  [`read_ordered_stack()`](https://olive-r.github.io/sptrends/reference/read_ordered_stack.md)
  (previously zero tests, unlike
  [`read_netcdf_stack()`](https://olive-r.github.io/sptrends/reference/read_netcdf_stack.md)),
  plus further NetCDF tests guarded by `skip_if_not_installed("ncdf4")`.
- Added a substantial test suite for `moran.R` (previously 40.8%
  coverage): hand-worked exact Moran’s I values under both connectivity
  types on a checkerboard pattern, confirming they genuinely differ and
  are each individually correct.

## sptrends 0.18.1

### Internal changes

- Fixed a leftover [`unname()`](https://rdrr.io/r/base/unname.html)
  omission in two tests comparing an extracted raster value against a
  literal number.

## sptrends 0.18.0

### Internal changes

- Added a full test suite for
  [`compute_anomalies()`](https://olive-r.github.io/sptrends/reference/compute_anomalies.md)
  (previously 0% coverage), verified against hand-worked expected
  values, not just “does it run”.

## sptrends 0.17.1

### Bug fixes

- Fixed a real bug from the previous round’s `testthat` idiom cleanup:
  [`terra::SpatRaster`](https://rspatial.github.io/terra/reference/SpatRaster-class.html)
  is S4, not S3 – 4 `expect_s3_class()` calls fixed to
  `expect_s4_class()`.

## sptrends 0.17.0

### Internal changes

- Addressed `goodpractice::gp()`’s findings: defensive
  `on.exit(add = TRUE)` fixes, removed placeholder `URL`/`BugReports`
  fields pointing to a repository that doesn’t exist yet, fixed a real
  NA-safety bug introduced while flattening a nested
  [`ifelse()`](https://rdrr.io/r/base/ifelse.html) in
  `fdr_direction_map()`, several idiomatic cleanups
  ([`anyDuplicated()`](https://rdrr.io/r/base/duplicated.html), inverted
  negated conditions), and confirmed two flagged findings as false
  positives.
- Added `goodpractice` and `lintr` to `Suggests`, formalising the
  package’s quality-evaluation toolkit.

## sptrends 0.16.0

### Improvements

- Clarified `alpha`’s own multi-value documentation: `0.05` is the
  conventional standard actually used for the map/reference message,
  `0.1`/`0.01` shown alongside for context, not as interchangeable
  defaults.

### Bug fixes

- Fixed a real bug affecting default output throughout the package: with
  the default `alpha = c(0.1, 0.05, 0.01)`, reporting functions used
  `min(alpha) = 0.01` (the strictest) instead of the conventional
  `0.05`, silently under-reporting significance in every default run.
  Fixed via a new `.reference_alpha()` helper.

### Internal changes

- Added `covr` to `Suggests`, formalising test coverage as a quality
  standard.

### Documentation

- Pedagogical sweep: every `@examples` block across the package now has
  at least one explanatory comment, verified programmatically.

## sptrends 0.15.0

### Internal changes

- Renamed `sptrends_example()` to
  [`example_data()`](https://olive-r.github.io/sptrends/reference/example_data.md)
  (inconsistent prefix), and two identifiers to British spelling
  (`summarise_replicates()`, `standardise`) – done now, before the
  package reaches CRAN.

### Documentation

- Made several functions’ and vignettes’ examples more pedagogical, and
  switched to the bundled real dataset where a synthetic one was less
  intuitive.
- Final exhaustive British-English sweep across the whole package;
  confirmed clean beyond the two renames above and identifiers that must
  keep American spelling (not ours to rename).

## sptrends 0.14.0

### Documentation

- `tst()`’s own examples and all four vignettes switched to the bundled
  real temperature dataset instead of
  [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md),
  a more recognisable example for the package’s flagship function and
  walkthroughs specifically.

## sptrends 0.13.1

### Improvements

- Fixed a stale test asserting `"accuracy"` should error, left over from
  before `Accuracy` became a valid metric.

## sptrends 0.13.0

### New features

- New bundled real-world example dataset (annual mean temperature, CRU
  TS v4.10, 1976-2025) and new function
  [`example_data()`](https://olive-r.github.io/sptrends/reference/example_data.md)
  to access it, doubling as a realistic
  [`read_ordered_stack()`](https://olive-r.github.io/sptrends/reference/read_ordered_stack.md)
  example.

### Improvements

- [`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md)
  gained `Accuracy` and `MCC` (Matthews correlation coefficient,
  generally preferred under class imbalance) as two more metrics, both
  included by default.

### Internal changes

- [`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md)’s
  `truth` argument renamed to `ground_truth`.

### Documentation

- New “Simulation and benchmarking” README section showcasing
  [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md) +
  [`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md)
  together.

## sptrends 0.12.0

### New features

- New functions `summarise_replicates()` and
  [`plot_detection_comparison()`](https://olive-r.github.io/sptrends/reference/plot_detection_comparison.md),
  completing the validation toolkit.
- New function
  [`compare_detections()`](https://olive-r.github.io/sptrends/reference/compare_detections.md):
  a confusion-matrix comparison of one or more trend-detection results
  against a known ground truth, deliberately agnostic to where either
  side comes from.

### Improvements

- [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)
  gained `noise_dist = "gaussian"/"t"`, generating heavy-tailed noise –
  the condition under which rank-based methods’ robustness argument over
  a parametric one like OLS actually shows up.

### Bug fixes

- Fixed a real bug in
  [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)’s
  own trend-assignment logic: sign flips were applied cell-by-cell at
  random (“salt-and-pepper”), silently defeating CMK’s whole rationale
  (borrowing strength from a plausibly-trending neighbourhood). Switched
  to coarse, spatially contiguous blocks.
- Fixed a real documentation bug: an internal helper had ended up inside
  [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)’s
  own roxygen block, documenting the wrong function as exported.

## sptrends 0.11.0

### Documentation

- Fixed the package logo’s background: solid white instead of
  transparent, showing a visible rectangle on any non-white background.

## sptrends 0.10.0

### New features

- [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)
  gained `noise_sd`, `trend_fraction` (proportion of cells keeping a
  non-zero true slope; `0` gives a complete null field for
  false-positive-rate checks), and `trend_shape`
  (“radial”/“gradient”/“block”).

### Improvements

- **Breaking**:
  [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)
  now returns `list(series = ..., true_slope = ...)` instead of the
  raster directly – what makes it a simulator with known ground truth,
  comparable directly against a fitted slope or a method’s own
  significance calls.

### Documentation

- Added a package logo.

## sptrends 0.9.0

### Documentation

- [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)’s
  examples now demonstrate `smooth_radius`’s own effect directly via
  Moran’s I (~0.07 vs ~0.65).

## sptrends 0.8.0

### New features

- [`sim_trend_stack()`](https://olive-r.github.io/sptrends/reference/sim_trend_stack.md)
  now genuinely generates spatial autocorrelation between neighbouring
  cells, via a new `smooth_radius` argument – previously documented but
  not actually implemented (noise was independent per cell).

### Improvements

- `moran_permutation_test()` documentation now explains why the
  permutation approach matters for this package’s own use case (bounded
  and/or skewed trend statistics, not normally distributed).
- Removed `sptrends_pipeline()` entirely – a deprecated alias for
  `tst()`, never part of a public release.
- [`read_netcdf_stack()`](https://olive-r.github.io/sptrends/reference/read_netcdf_stack.md)
  now issues a real [`warning()`](https://rdrr.io/r/base/warning.html),
  not just a documentation note, when the detected time step looks
  sub-annual.
- **Behaviour change**: `tst()`’s default `fdr_method` changed from
  `c("BH", "BKY")` to `"BKY"` only, the method actually used in the
  published TST methodology.

## sptrends 0.7.0

### Improvements

- Benjamini & Yekutieli (2001) now cited explicitly in all 11 `fdr_*`
  functions, not only
  [`fdr_bh()`](https://olive-r.github.io/sptrends/reference/fdr_bh.md).
- `tst()`’s own argument order regrouped for logical proximity (`alpha`
  with `moran_check`, `fdr_method` with `q`); named-argument code is
  unaffected.
- Strengthened the `alpha`/`q` note: applying a per-test `alpha` cell by
  cell across a raster is the multiple-testing error FDR exists to fix,
  not a harmless simplification.

### Internal changes

- Prose switched to British English throughout the package. Left
  unchanged, deliberately: identifiers that must keep American spelling
  (not ours to rename) and every quoted paper title.
- Renamed “multiple comparisons” to “multiple testing” throughout.

### Documentation

- Corrected a co-author’s affiliation, mistakenly duplicated from the
  other author’s.
- `Description` rewritten twice more for clarity, naming each workflow
  component explicitly before citing the TST paper.
- Fixed vignette ordering in
  [`browseVignettes()`](https://rdrr.io/r/utils/browseVignettes.html)
  (sorts by title, not file name) by prefixing the getting-started
  vignette’s own title.

## sptrends 0.6.0

### Internal changes

- Vignettes renamed with an `a-`/`b-`/`c-`/`d-` prefix so
  [`browseVignettes()`](https://rdrr.io/r/utils/browseVignettes.html)
  lists them in reading order.

### Documentation

- Added a “Citation” section to `README.md`.

## sptrends 0.5.0

### Improvements

- Every exported function’s help page now has a “Function type” section,
  and a runnable `@examples` block (previously missing on 21 of 34
  functions).

### Documentation

- Strengthened the package’s identity messaging (official R
  implementation of the published TST framework) across `README.md` and
  package docs; finalised `Title`/`Description` around a
  reproducible-framework framing; added `_pkgdown.yml` for a grouped
  future reference index; fixed a duplicated header in `README.md`.

## sptrends 0.4.0

### Improvements

- `Title` finalised, deliberately not naming TST or “monotonic” in the
  title itself, since TST is meant to be the flagship workflow, not
  necessarily the only one going forward.

### Documentation

- `Description` rewritten to lead with the reproducible-framework
  framing before naming TST.
- Corrected imprecise use of “adaptive” – only BKY is adaptive in the
  technical sense, BH is not.
- New “Getting started” vignette.

## sptrends 0.3.0

### New features

- `tst()` now includes Theil-Sen as a fourth, optional step, completing
  the four-component TST workflow as published.

### Improvements

- Every exported function’s help page now states its Function type
  (Core/Support/Reporting-derived).
- `tst()` now returns a classed object with
  [`print()`](https://rdrr.io/r/base/print.html)/[`summary()`](https://rdrr.io/r/base/summary.html)/
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

### Internal changes

- New main function `tst()`, renamed from `sptrends_pipeline()`
  (matching the published workflow’s own name); the old name kept as a
  deprecated alias.

### Documentation

- `DESCRIPTION` and `README.md` rewritten to lead with the general
  monotonic-trend framing before naming TST by name.

## sptrends 0.2.0

### New features

- New `theil_sen_slope()`, `fdr_direction_map()`/`_plot()`/`_summary()`,
  [`compute_anomalies()`](https://olive-r.github.io/sptrends/reference/compute_anomalies.md),
  and NetCDF support for
  [`read_ordered_stack()`](https://olive-r.github.io/sptrends/reference/read_ordered_stack.md)
  plus a new
  [`read_netcdf_stack()`](https://olive-r.github.io/sptrends/reference/read_netcdf_stack.md).

### Improvements

- Added second author; added `inst/CITATION`; added the package-level
  help page; `report = TRUE` added as the default across the core
  functions; `moran_check` added to
  [`fdr_correction()`](https://olive-r.github.io/sptrends/reference/fdr_correction.md)/the
  pipeline; parallel computation (`n_cores`) added to the CMK and
  Moran’s I steps; documentation overhaul with verified/corrected DOIs
  throughout.

### Internal changes

- `DESCRIPTION` rewritten around the three-step workflow.

### Documentation

- Vignettes restructured from one combined vignette into three
  step-by-step ones.

## sptrends 0.1.0

Initial release, built from four standalone scripts to reproduce the TST
workflow described at the top of this file: selective AR(1) prewhitening
(Wang & Swail, 2001), Contextual Mann-Kendall (Neeti & Eastman, 2011),
pixel-wise FDR correction (Benjamini-Hochberg 1995;
Benjamini-Krieger-Yekutieli 2006), and Moran’s I permutation testing,
plus
[`compute_anomalies()`](https://olive-r.github.io/sptrends/reference/compute_anomalies.md),
raster/NetCDF reading, a convenience pipeline wrapper, and a synthetic
data generator for examples and tests.
