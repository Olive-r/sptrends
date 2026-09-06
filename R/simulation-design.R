#' Build a factorial design of simulation scenarios
#'
#' Creates named argument lists for [sim_trend_stack()] by crossing temporal,
#' spatial, signal and noise conditions. It separates experimental design from
#' data generation so the complete scenario grid can be inspected and retained.
#'
#' **Function type:** **Benchmarking function** -- defines simulation
#' scenarios; it does not generate data or perform inference.
#'
#' @section Typical use:
#' Define the factors that should vary, add shared settings through `constants`,
#' and pass the returned list to [benchmark_methods()].
#'
#' @section Methodological details:
#' **Experimental design**
#'
#' Every combination is retained. This makes comparisons across spatial and
#' temporal dependence explicit and prevents methods from being evaluated on
#' accidentally different scenario sets. Values that must remain grouped,
#' such as a two-number `signal_size`, should be wrapped in a list.
#'
#' **Computational considerations**
#'
#' This function only constructs argument lists. Memory use grows with the
#' product of the numbers of supplied factor levels; data generation remains
#' deferred to [benchmark_methods()].
#'
#' **Limitations**
#'
#' A full factorial design can become unnecessarily large. Users should vary
#' scientifically relevant factors and keep fixed settings in `constants`.
#'
#' **Quality assurance**
#'
#' Tests verify factorial completeness, deterministic ordering, grouped
#' values, validation failures, metadata retention and S3 presentation. The
#' complete external simulation-cycle validation passed all 33 prespecified
#' controls; see `inst/validation/` for the retained protocol and results.
#'
#' @param ... Named vectors or lists of factor levels to cross.
#' @param constants Named list of arguments shared by every scenario.
#' @param prefix Character prefix used for generated scenario names.
#' @param verbose Logical. If `TRUE`, reports progress, elapsed time and the
#'   estimated time remaining while scenarios are assembled.
#'
#' @return A named list of argument lists suitable for the `scenarios`
#'   argument of [benchmark_methods()], with classes
#'   `"sptrends_simulation_design"` and `"sptrends"` for unified
#'   printing, summaries, and plotting.
#'
#' @examples
#' design <- simulation_design(
#'   spatial_model = c("independent", "exponential"),
#'   spatial_rho = c(0.3, 0.7), ar1 = c(0, 0.5),
#'   trend_strength = c(0, 0.05),
#'   constants = list(nrow = 20, ncol = 20, n_time = 20,
#'                    constant_block = FALSE), verbose = FALSE)
#' length(design)
#'
#' @family validation functions
#' @seealso [sim_trend_stack()], [benchmark_methods()]
#' @export
simulation_design <- function(..., constants = list(), prefix = "scenario",
                              verbose = TRUE) {
  finish_timer <- .sptrends_elapsed_timer("simulation_design()", verbose)
  on.exit(finish_timer(), add = TRUE)
  factors <- list(...)
  if (length(factors) == 0L || is.null(names(factors)) ||
      any(names(factors) == "")) {
    stop("Simulation factors supplied through '...' must be named.")
  }
  if (!is.list(constants) ||
      (length(constants) > 0L &&
       (is.null(names(constants)) || any(names(constants) == "")))) {
    stop("'constants' must be a named list.")
  }
  overlap <- intersect(names(factors), names(constants))
  if (length(overlap) > 0L) {
    stop("Factors and constants duplicate: ", paste(overlap, collapse = ", "))
  }
  if (!is.character(prefix) || length(prefix) != 1L || is.na(prefix) ||
      !nzchar(prefix)) {
    stop("'prefix' must be one non-empty character value.")
  }
  empty <- vapply(factors, length, integer(1)) == 0L
  if (any(empty)) {
    stop("Every simulation factor must contain at least one level.")
  }

  indices <- expand.grid(
    lapply(factors, seq_along), KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE)
  progress <- .sptrends_progress(
    nrow(indices), "Building simulation scenarios", verbose)
  on.exit(.sptrends_progress_close(progress), add = TRUE)
  scenarios <- lapply(seq_len(nrow(indices)), function(row_index) {
    selected <- lapply(seq_along(factors), function(column_index) {
      factors[[column_index]][[indices[row_index, column_index]]]
    })
    names(selected) <- names(factors)
    scenario <- c(selected, constants)
    .sptrends_progress_step(progress, row_index)
    scenario
  })
  .sptrends_progress_close(progress)
  width <- nchar(length(scenarios))
  names(scenarios) <- sprintf(
    paste0(prefix, "_%0", width, "d"), seq_along(scenarios))
  attr(scenarios, "factors") <- names(factors)
  attr(scenarios, "constants") <- names(constants)
  class(scenarios) <- c(
    "sptrends_simulation_design", "sptrends", "list")
  scenarios
}
