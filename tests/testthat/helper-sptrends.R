# Shared test helpers, in a helper-*.R file on purpose: testthat sources
# every helper-*.R before ANY test-*.R file, regardless of alphabetical
# order -- unlike a plain top-level function defined inside a test-*.R
# file, which is only available to test files that happen to run AFTER
# it. .make_ar1_raster() is used from more than one test file, so it
# belongs here, not inside any single one of them (a real bug this fixed:
# it used to live in test-prewhiten-extra.R, invisible to
# test-inspect-ts-cell.R, which sorts alphabetically before it).

.make_ar1_raster <- function(rho, n_time = 40, nrow = 4, ncol = 4, sd_e = 1, seed = 1) {
  set.seed(seed)
  ncell_ <- nrow * ncol
  noise <- matrix(0, ncell_, n_time)
  noise[, 1] <- stats::rnorm(ncell_, sd = sd_e)
  for (k in 2:n_time) {
    noise[, k] <- rho * noise[, k - 1] + stats::rnorm(ncell_, sd = sd_e * sqrt(1 - rho^2))
  }
  r <- terra::rast(nrows = nrow, ncols = ncol, xmin = 0, xmax = ncol, ymin = 0, ymax = nrow)
  layers <- lapply(seq_len(n_time), function(k) terra::setValues(r, noise[, k] + 10))
  do.call(c, layers)
}
