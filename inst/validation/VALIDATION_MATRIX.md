# Numerical validation matrix

This file maps the main methodological families to independent or analytical
controls. Automated comparisons are skipped cleanly when an optional reference
package is unavailable.

| Family | Control | Location |
|---|---|---|
| Mann-Kendall | `Kendall::MannKendall()` and `trend::mk.test()` | `test-external-validation.R` |
| Contextual Mann-Kendall | `rkt::rkt()`, ConMK and TerrSet records | `test-external-validation.R`, `inst/validation/` |
| Theil-Sen slope | `trend::sens.slope()` | `test-external-validation.R` |
| OLS test and slope | `stats::lm()` analytical contracts | trend and slope test files |
| Repeated-median slope | direct Siegel repeated-median formula | `test-rm-slope.R` |
| TFPW-Y | `modifiedmk::tfpwmk()` | `test-external-validation.R` |
| MMK | equations and source behaviour of `modifiedmk` | `test-mmk.R` |
| BH and BY | `stats::p.adjust()` | FDR test files |
| BKY | published two-stage step-up equations | FDR test files |
| Spatial covariance simulation | analytical covariance plus `fields::Matern()` | `test-simulation-formal.R` |
| Detection, FDR and FWER metrics | hand-calculated confusion matrices and replicate experiments | `test-validation.R` |
| Complete simulation cycle | 33-control external run; 1,000 fields per spatial model and paired `Kendall` MK benchmarks | `SIMULATION_CYCLE_0.96.3.md` |

Agreement with another implementation is not substituted for a methodological
definition. When conventions legitimately differ, the tests document the
quantity expected to agree and avoid asserting equality for a different
estimand.
