# Reproducible time-and-memory scaling benchmark for sptrends.
#
# This script is installed for manual quality assurance; it is not run by
# R CMD check. The reported matrix size is a transparent lower-bound proxy,
# not peak resident memory: trend_test(method = "CMK") and prewhiten() create
# several working matrices of comparable dimensions.

benchmark_sptrends_scalability <- function(
    dimensions = c(50L, 100L, 200L),
    n_time = 20L,
    methods = c("MK", "CMK", "TFPW_WS"),
    n_cores = 1L,
    seed = 1L) {
  stopifnot(
    all(dimensions > 0),
    n_time >= 4L,
    n_cores >= 1L,
    all(methods %in% c("MK", "CMK", "TFPW_WS"))
  )

  rows <- list()
  k <- 0L
  for (dimension in dimensions) {
    simulated <- sptrends::sim_trend_stack(
      nrow = dimension,
      ncol = dimension,
      n_time = n_time,
      smooth_radius = 0,
      constant_block = FALSE,
      seed = seed
    )$series

    n_cells <- terra::ncell(simulated)
    matrix_mib <- n_cells * n_time * 8 / 1024^2

    for (method in methods) {
      gc()
      elapsed <- system.time({
        result <- if (method %in% c("MK", "CMK")) {
          sptrends::trend_test(
            simulated,
            method = method,
            n_cores = n_cores,
            report = FALSE,
            verbose = FALSE
          )
        } else {
          sptrends::prewhiten(
            simulated,
            method = method,
            report = FALSE,
            verbose = FALSE
          )
        }
      })[["elapsed"]]

      k <- k + 1L
      rows[[k]] <- data.frame(
        dimension = dimension,
        cells = n_cells,
        time_steps = n_time,
        method = method,
        n_cores = n_cores,
        input_matrix_MiB = matrix_mib,
        result_object_MiB = as.numeric(object.size(result)) / 1024^2,
        elapsed_seconds = elapsed
      )
      rm(result)
    }
    rm(simulated)
  }

  do.call(rbind, rows)
}

# Example:
# timings <- benchmark_sptrends_scalability(
#   dimensions = c(50L, 100L, 200L),
#   n_time = 20L,
#   n_cores = 1L
# )
# print(timings)
