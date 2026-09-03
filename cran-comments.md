## Test environments

* Local: `devtools::check(args = "--as-cran")` run directly against
  1.5.9: 0 errors, 0 warnings, 0 notes. Duration under 10 minutes.
  `covr::package_coverage()`: 100.00% overall and 100.00% for every
  R source file, confirmed with `SPTRENDS_TEST_PARALLEL=true` set
  (see "Parallel-execution tests" below for why this is needed).
* win-builder (R-devel, R-release, R-oldrelease), most recently run
  against 1.5.8: 0 errors, 0 warnings. Remaining NOTEs (misspelled
  words in DESCRIPTION) are expected, domain terminology -- see
  below.

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

Most recent local check, against 1.5.9: 0 errors, 0 warnings,
0 notes.

win-builder (R-devel, R-release, R-oldrelease), against 1.5.8:
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

## Downstream dependencies

There are currently no downstream dependencies for this package, as
this is a new submission.
