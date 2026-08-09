## Test environments

* Local: `devtools::check()` run directly against 1.3.4: 0 errors,
  0 warnings, 0 notes. `covr::package_coverage()`: 100.00% overall
  and 100.00% for every R source file.
* win-builder (R-devel, R-release, R-oldrelease) was last run against
  0.97.1, a much earlier version -- **needs to be re-run against the
  current version before submission.** That earlier run confirmed
  0 errors, 0 warnings, and only the expected new-submission note
  across all three platforms, including a since-fixed URL issue;
  those specific results should not be quoted as current.

## R CMD check results

Most recent completed local check, against 1.3.4: 0 errors,
0 warnings, 0 notes.

win-builder has not yet been re-run against 1.3.4 -- do this before
submission and update this file with the actual result.

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
