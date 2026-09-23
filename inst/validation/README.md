# External validation of `trend_test(method = "CMK")`

This folder documents the external comparisons discussed in
`?trend_test`'s own "External validation" section, in a form future
maintainers or reviewers can inspect or re-run without depending on
this session's own conversation history.

## Files

- `VALIDATION_MATRIX.md` -- concise map from each methodological family to
  its independent-package comparison or analytical reference test.
- `SIMULATION_CYCLE_0.96.3.md` -- retained protocol, numerical results,
  corrected graphical interpretation and limitations for the 33-control
  external validation of simulation and benchmarking.
- `mk-type-i-by-ar1.png` and `mk-power-by-ar1.png` -- corrected external
  validation figures separating null Type I error from signal-region power.
- `benchmark_mk_cmk_spatial_dependence.R` -- reproducible Monte Carlo
  comparison in which MK and CMK receive the same realisations, known truth
  and evaluation cells under several levels of spatial dependence.

- `compare_conmk.R` -- the reproducible script for comparing
  `sptrends`'s `S`/`Sm`/`VarSm`/`p` against
  `ConMK::contextual_mann_kendall()`
  (<https://github.com/antiphon/ConMK>, install via
  `devtools::install_github("antiphon/ConMK")`; not on CRAN, so this
  script is not run automatically by `R CMD check` and requires that
  package installed manually first).
- `conmk_comparison_100cells.csv` -- the frozen numeric result retained
  from the original run, for all 100 cells of a
  `sim_trend_stack(nrow = 10, ncol = 10, n_time = 15, seed = 1)` raster:
  it records `Sm`/`theirs_S` (which match to floating-point precision
  throughout) and `ours_p`/`theirs_p`. It predates the addition of the
  two variance columns to `compare_conmk.R`; therefore, this frozen CSV
  does not itself substantiate a numerical `VarSm`/`s2` match. Re-running
  the current script refreshes the file with both variance estimates.
  For the recorded p-values,
  `ours_p`/`theirs_p` (`theirs_p` is larger for 97 of the 100 cells --
  the continuity-correction effect described below; the 3 exceptions,
  cells 1/12/59, all have `Sm` near zero, where subtracting 1 crosses
  zero and increases rather than decreases `|Z|`, the same mechanism
  producing the opposite-signed effect at that specific edge -- see
  `?trend_test`).
- `trend_test()`'s own `continuity = TRUE` argument reproduces
  `theirs_p` in `conmk_comparison_100cells.csv` for 91 of the
  100 cells -- verified by an automated test (`test-cmk-extra.R`), not
  only checked by hand once when this comparison was first performed.
  The remaining 9 cells (1, 2, 3, 11, 12, 13, 21, 22, 23) all sit
  within one queen-neighbourhood step of a constant, zero-variance
  cell (this specific simulated raster happens to produce a small
  constant-valued cluster at cells 1/2/11/12) -- the cross-correlation
  term either implementation needs is undefined (0/0) for those, and
  this package's own convention for it is not asserted to match
  `ConMK`'s own, unverified, handling of the same undefined case.
  When `Sm = 0`, the option deliberately follows `ConMK`'s `+1` branch;
  it is a compatibility convention rather than the package default.
- `terrset_comparison_10cells.csv` -- the same comparison against
  TerrSet's own Kendall module (both classic `MK`, no spatial pooling,
  and the contextual `CMK` case, via its
  `KENDALL_Crosscorrelation` module), for the first 10 cells of the
  same raster. Limited to 10 cells (not the full 100 `ConMK`
  comparison covers) because TerrSet's own `.rst` output was
  transcribed manually from its GUI, not scripted -- see the "External
  validation" section in `?trend_test` for why the `MK` case matched
  exactly while the `CMK` case is reported as inconclusive, not a
  confirmed discrepancy, given TerrSet's own closed source.
  In the like-for-like run on the same raster, `sptrends` completed CMK
  substantially faster than TerrSet. No universal speed ratio is
  claimed because timings depend on hardware, storage, raster size,
  software versions and parallel settings.
- The optional automated test in
  `tests/testthat/test-external-validation.R` also treats a 3 x 3
  neighbourhood as nine regional blocks and verifies that
  `rkt::rkt()` returns the same aggregated score, `S / 9 = Sm`.
  It deliberately does not require equal corrected variances or
  p-values: `rkt` implements the Hirsch--Slack inter-block correction,
  whereas CMK uses the Neeti--Eastman contextual variance formula.

## Why this lives in `inst/`, not `tests/`

The frozen ConMK and TerrSet comparisons do not run during
`R CMD check`: `ConMK` is not installable from CRAN, and TerrSet is
proprietary desktop software with no scriptable interface used here.
The `rkt` comparison does run when that suggested package is available
and otherwise skips cleanly. The files in this directory are a
permanent, inspectable record of validation exercises already
performed -- if `sptrends`'s own `CMK` formula ever changes,
re-running `compare_conmk.R` (with `ConMK` installed) is how to check
whether it still agrees with this comparison, not something this
package's own test suite verifies on every run.
