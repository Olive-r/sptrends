## Test environments

* Local: `devtools::check()` run directly against 1.5: 0 errors,
  0 warnings, 0 notes. `covr::package_coverage()`: 100.00% overall
  and 100.00% for every R source file.
* win-builder (R-devel, R-release, R-oldrelease), run against 1.5:
  0 errors, 0 warnings, 3 NOTEs on R-devel and R-release, 2 NOTEs
  on R-oldrelease (see below). A `LICENSE` file present but not
  referenced in `DESCRIPTION` was among these NOTEs; it has since
  been excluded from the build via `.Rbuildignore` (CRAN does not
  permit bundling a copy of a standard license, and `GPL (>= 3)`
  alone is already a complete, standard licence specification) --
  confirmed resolved in a subsequent local `--as-cran` check
  (0 errors, 0 warnings, 0 notes), though win-builder itself has
  not been re-run after this specific fix.

## R CMD check results

Most recent completed local check, against 1.5 (after the
`.Rbuildignore` fix described above): 0 errors, 0 warnings, 0 notes.

win-builder (R-devel, R-release, R-oldrelease), against 1.5, run
before the `.Rbuildignore` fix: 0 errors, 0 warnings. NOTEs seen:

* "New submission" -- expected, this is a first submission to CRAN.
* Possibly misspelled words in DESCRIPTION (`gridded`, `prewhitening`,
  `spatiotemporal`, and their capitalised forms) -- all correct
  domain terminology for this package's subject matter, not
  misspellings.
* `LICENSE` file not mentioned in DESCRIPTION -- since fixed, see
  above; not expected to reappear.
* A timeout connecting to <https://www.gnu.org/licenses/gpl-3.0> (an
  existing, working URL) appeared on R-devel and R-release but not
  on R-oldrelease at a different check time, consistent with a
  temporary network issue on the checking infrastructure rather than
  the URL itself.

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
