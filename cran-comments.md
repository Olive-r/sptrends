## Test environments

* Local: `devtools::check()` run directly against 1.0.8: 0 errors,
  0 warnings, 0 notes (20 minutes 31.4 seconds).
  `covr::package_coverage()`: 100.00% overall and 100.00% for every R
  source file.
* win-builder R-devel -- run against 0.97.1: **1 NOTE only** (the
  expected new-submission note and known false-positive spelling
  flags -- `Gridded`, `Spatiotemporal`, `gridded`, `prewhitening`,
  `spatiotemporal`). The URL fix is confirmed: the `Found the
  following (possibly) invalid URLs` section that flagged the
  Wang & Swail (2001) DOI as a 404 against 0.97 no longer appears.
* win-builder R-release (R 4.6.1) -- run against 0.97.1: **1 NOTE
  only** (same expected new-submission note and spelling flags as
  R-devel). URL fix confirmed here too: no invalid-URL section.
* win-builder R-oldrelease (R 4.5.3) -- run against 0.97.1: **1 NOTE
  only** (same expected new-submission note and spelling flags as
  R-devel/R-release). URL fix confirmed here too: no invalid-URL
  section.

**All three win-builder platforms (R-devel, R-release, R-oldrelease)
are confirmed clean against 0.97.1, including the URL fix.**

## R CMD check results

Most recent completed local check, against 1.0.8: 0 errors, 0 warnings,
0 notes.

The win-builder checks listed above were run against 0.97.1 and each
returned 0 errors, 0 warnings, and 1 expected CRAN incoming note,
unavoidable for a new submission:

* New submission.

## Additional quality checks

* `spelling::spell_check_package()`: only legitimate technical terms
  and proper nouns found in the last local run; all added to
  `inst/WORDLIST`.

## Downstream dependencies

There are currently no downstream dependencies for this package, as
this is a new submission.
