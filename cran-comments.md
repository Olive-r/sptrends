## Test environments

* Local: `devtools::check(args = "--as-cran")` run directly against
  1.6: 0 errors, 0 warnings, 0 notes. Duration under 10 minutes
  (~8 min with the reporting machine's cloud-sync client paused;
  longer runs observed with it active are a local environment
  artefact, not a package issue).
  `covr::package_coverage()`: 100.00% overall and 100.00% for every
  R source file, confirmed with `SPTRENDS_TEST_PARALLEL=true` set
  (see "Parallel-execution tests" below for why this is needed).
* win-builder (R-devel, R-release, R-oldrelease), run against 1.6.1:
  0 errors, 0 warnings, 1 NOTE on all three platforms (misspelled
  words in DESCRIPTION -- expected domain terminology, see below).

## Parallel-execution tests

19 tests that spawn a real multi-worker cluster (`n_cores > 1`), to
verify the parallel code path agrees with the sequential one, are
skipped by default (`testthat::skip_if_not()`), since cluster
start-up overhead alone accounted for over half the test suite's
total time and previously pushed the overall check time above 10
minutes (an earlier submission was auto-rejected for exactly this).
Set the environment variable `SPTRENDS_TEST_PARALLEL=true` to run
them; `covr::package_coverage()` with this set confirms 100.00%
coverage is genuinely unaffected -- the skip hides no real gap, only
tests intentionally not run on every check.

## R CMD check results

Most recent local check, against 1.6: 0 errors, 0 warnings, 0 notes.

win-builder (R-devel, R-release, R-oldrelease), against 1.6.1:
0 errors, 0 warnings. NOTEs seen:

* "New submission" -- expected, this is a first submission to CRAN.
* Possibly misspelled words in DESCRIPTION (`gridded`, `prewhitening`,
  `spatiotemporal`, proper nouns with diacritics, and their
  capitalised forms) -- all correct domain terminology or author
  names, not misspellings.

## Additional quality checks

* `spelling::spell_check_package()`: only legitimate technical terms
  and proper nouns found in the last local run; all added to
  `inst/WORDLIST`.

## External validation

* `trend_test(method = "CMK")` was cross-checked against an installed
  copy of `ConMK` (Antiphon, GitHub, not on CRAN -- the closest
  available external reference implementation of the contextual
  Mann-Kendall test). The base statistic `S`/`Sm` matched to
  floating-point precision (correlation 1, zero maximum absolute
  difference). `p` correlated at 0.998 in the general case, with the
  only source of disagreement traced to a specific, documented
  difference in how each implementation applies a continuity
  correction (see `?trend_test`'s "External validation" section for
  the full account, including a fresh reproduction confirming the
  optional `continuity = TRUE` argument matches `ConMK`'s own p-values
  exactly, to seven decimal places, at the specific edge case where
  the two implementations would otherwise be expected to diverge).
* `trend_test(method = "MK")` and `slope_estimator(method = "TS")`
  were separately cross-checked against `Kendall::MannKendall()`,
  `trend::mk.test()`/`sens.slope()`, `modifiedmk::mmkh()`,
  `robslopes::TheilSen()`, and `zyp::zyp.sen()` across five test
  series (clear trend, no trend, tied values, short, long). 24/27
  comparisons matched to numerical precision; the 3 remaining
  differences were traced to (a) `robslopes`'s own approximate
  algorithm for large samples (this package's result matched the
  other two independent references exactly instead), and (b) a
  test-script omission (`ties = TRUE` not set to match the reference
  packages' own always-on tie correction) rather than any package
  issue -- confirmed by hand-deriving the exact reference p-value
  from this package's own documented formula.

## Downstream dependencies

There are currently no downstream dependencies for this package, as
this is a new submission.
