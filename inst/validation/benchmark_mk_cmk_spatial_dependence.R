# Reproducible MK-versus-CMK experiment under calibrated spatial dependence.
# This script deliberately evaluates raw trend-test behaviour; it does not
# apply FDR or FWER correction.

rho_levels <- c(0, 0.3, 0.6, 0.9)
signal_sizes <- c(5, 20, 50)

scenarios <- list()
for (rho in rho_levels) {
  for (size in signal_sizes) {
    scenario_name <- sprintf("rho_%0.1f_size_%d", rho, size)
    scenarios[[scenario_name]] <- list(
      nrow = 100, ncol = 100, n_time = 34,
      trend_shape = "square", signal_size = size,
      signal_location = "centre", trend_strength = 0.03,
      trend_fraction = 1, constant_block = FALSE,
      ar1 = 0, noise_sd = 1,
      spatial_model = if (rho == 0) "independent" else "gaussian",
      spatial_rho = if (rho == 0) 1 else rho)
  }
}

methods <- list(
  MK = function(series, simulation) {
    result <- trend_test(
      series, method = "MK", report = FALSE, verbose = FALSE)
    list(significant = result$stats$p <= 0.05,
         direction = result$stats$S)
  },
  CMK = function(series, simulation) {
    result <- trend_test(
      series, method = "CMK", window_size = 3,
      report = FALSE, verbose = FALSE)
    list(significant = result$stats$p <= 0.05,
         direction = result$stats$Sm)
  })

# Increase n_replicates to 1000 for the final Monte Carlo experiment.
comparison <- benchmark_methods(
  scenarios = scenarios, methods = methods, n_replicates = 100,
  seed = 2026,
  metrics = c("type_i", "type_ii", "type_iii", "global_power",
              "within_image_power", "directional_power"))

comparison

# Performance profiles of the kind used in simulation studies: one curve per
# method, spatial dependence on the x axis, and signal size in separate panels.
plot(
  comparison,
  metric = c("GlobalPower", "WithinImagePower", "TypeI",
             "TypeII", "TypeIII"),
  scenario = "spatial_rho", facet = "signal_size",
  type = "profile", interval = "ci")

benchmark_summary(comparison)
