# External validation of trend_test(method = "CMK") against ConMK
# (Antiphon, GitHub: https://github.com/antiphon/ConMK, not on CRAN).
#
# Requires ConMK installed manually first:
#   devtools::install_github("antiphon/ConMK", build_vignettes = TRUE)
#
# Not run automatically by R CMD check -- see README.md in this same
# folder for why, and for what this comparison found.

library(sptrends)
library(ConMK)
library(terra)
library(raster)

# Same simulated raster used throughout this comparison.
r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 15, seed = 1)$series

# sptrends's own CMK.
ours <- trend_test(r, method = "CMK", report = FALSE, verbose = FALSE)

# ConMK's own contextual Mann-Kendall, on the same raster converted to
# the older raster::RasterStack representation ConMK itself expects.
r_stack <- raster::stack(r)
theirs <- ConMK::contextual_mann_kendall(r_stack)

result <- data.frame(
  cell     = seq_len(terra::ncell(r)),
  ours_Sm  = terra::values(ours$stats$Sm),
  theirs_S = raster::values(theirs$S),
  ours_VarSm  = terra::values(ours$stats$VarSm),
  theirs_s2   = raster::values(theirs$s2),
  ours_p   = terra::values(ours$stats$p),
  theirs_p = raster::values(theirs$p)
)

# Sm/theirs_S and VarSm/theirs_s2 are expected to match to
# floating-point precision; ours_p/theirs_p are expected to differ
# systematically -- see ?trend_test's own "External validation"
# section for why (ConMK applies a continuity correction this
# package's own reading of Eqs. 7-10, Neeti & Eastman 2011, does not
# support for the neighbourhood-averaged Sm).
write.csv(result, "conmk_comparison_100cells.csv", row.names = FALSE)
