#' Benchmark statistical methods across known-truth simulation scenarios
#'
#' Coordinates reproducible Monte Carlo experiments without tying the
#' benchmark to sptrends implementations. Each method is an ordinary function,
#' so methods from other packages can be evaluated under exactly the same
#' simulated realisations and truth fields.
#'
#' **Function type:** **Benchmarking function** -- coordinates simulation,
#' method execution, scoring and timing; it is not an inferential workflow.
#'
#' @section Typical use:
#' Define named scenarios, define package or external method wrappers, and run
#' the same methods on every generated realisation. Built-in scoring is
#' available for prewhitening, trend tests, slopes, FDR and FWER. A custom
#' evaluator supports arbitrary result structures without changing the Monte
#' Carlo engine.
#'
#' @section Methodological details:
#' **Paired method comparison**
#'
#' A replicate seed is generated once and shared by every method within that
#' replicate. Consequently, method differences are paired within identical
#' data, known truth, prepared input and evaluation domain rather than
#' confounded with different random fields or cell subsets. The returned table
#' retains replicate-level results; uncertainty summaries must be calculated
#' across those independent replicates, not across raster cells.
#'
#' **Statistical assumptions**
#'
#' Performance estimates are conditional on the simulated scenarios and the
#' common evaluation domain. Monte Carlo replicates, not raster cells, are
#' the independent units used to estimate repeated-sampling behaviour.
#'
#' **Computational considerations**
#'
#' Runtime grows with scenarios, replicates and methods. Scenario results are
#' retained at replicate level so expensive experiments can be summarised or
#' plotted without rerunning the methods.
#'
#' **Limitations**
#'
#' Built-in scorers require the documented truth components. Method-specific
#' objects or cluster-level targets require a custom `evaluator`; incompatible
#' method-specific evaluation domains are deliberately rejected.
#'
#' **Quality assurance**
#'
#' Tests cover every supported stage, paired inputs, external simulators,
#' common masks, reproducible seeds, timing, failures, summaries and plots.
#' In an independent full run, 500 replicates were evaluated in each of eight
#' spatiotemporal scenarios. Cell-level MK decisions and directions agreed
#' exactly between sptrends and `Kendall::MannKendall()` for every retained
#' performance metric, while paired seeds, summaries and graphics passed all
#' recorded controls. This validates the benchmark orchestration and scoring;
#' it does not imply that every method compared by future users is equivalent.
#'
#' @param scenarios Named list of argument lists passed to `simulator`.
#' @param methods Named list of functions. Each receives `(input, simulation)`
#'   and returns a transformed raster series for `stage = "prewhitening"`, a
#'   detection object for trend/FDR/FWER stages, a slope
#'   raster/vector for `stage = "slope"`, or an object accepted by `evaluator`.
#' @param n_replicates Positive integer number of realisations per scenario.
#' @param stage One of `"prewhitening"`, `"trend_test"`, `"slope"`, `"fdr"`,
#'   `"fwer"`, or `"custom"`.
#'   `"detection"` is retained as an alias for `"trend_test"`.
#' @param simulator Function used to generate one known-truth realisation. It
#'   must accept the scenario arguments plus `seed`, and return at least
#'   `series`, `true_signal`, `true_direction`, and `true_slope` components
#'   when the corresponding built-in scorer is used.
#' @param prepare Optional function called once per replicate as
#'   `prepare(series, simulation)`. Its result becomes the identical `input`
#'   passed to every method. Use it to compute one common p-value or statistic
#'   map before comparing FDR or FWER methods.
#' @param evaluator Optional scoring function called as
#'   `evaluator(outputs, simulation)`. Required for `stage = "custom"`.
#' @param seed Integer seed that deterministically generates replicate seeds.
#' @param metrics Metrics passed to [compare_detections()] for detection
#'   benchmarks.
#' @param evaluation_mask Optional common mask or a function of `simulation`
#'   returning one. The same mask is applied to every method so all methods
#'   are compared on exactly the same cells.
#' @param verbose Logical. If `TRUE`, report method-level progress, elapsed
#'   duration and estimated time remaining across the complete Monte Carlo
#'   experiment. Simulator-level verbosity is disabled automatically when the
#'   simulator exposes a `verbose` argument and the scenario does not override
#'   it.
#'
#' @return A data frame with one row per scenario, replicate and method,
#'   including explicit scenario factors, known-truth composition, elapsed
#'   time and stage-appropriate accuracy metrics. It has classes
#'   `"sptrends_benchmark"` and `"sptrends"`, providing unified
#'   `print()`, `summary()`, and scenario-performance `plot()` methods.
#'
#' @examples
#' methods <- list(
#'   MK = function(series, simulation) {
#'     fit <- trend_test(series, method = "MK", report = FALSE,
#'                       verbose = FALSE)
#'     list(significant = fit$stats$p <= 0.05,
#'          direction = fit$stats$S)
#'   }
#' )
#' scenarios <- list(null = list(nrow = 6, ncol = 6, n_time = 8,
#'                               trend_fraction = 0))
#' result <- benchmark_methods(scenarios, methods, n_replicates = 2,
#'                             seed = 1, verbose = FALSE)
#' result
#'
#' @family validation functions
#' @seealso [simulation_design()], [sim_trend_stack()],
#'   [compare_detections()], [benchmark_summary()]
#' @export
benchmark_methods <- function(scenarios, methods, n_replicates = 100L,
                              stage = c("trend_test", "prewhitening", "slope",
                                        "fdr", "fwer", "custom", "detection"),
                              simulator = sim_trend_stack, prepare = NULL,
                              evaluator = NULL,
                              seed = 1L,
                              metrics = c("type_i", "type_ii", "type_iii",
                                          "field_power",
                                          "global_power",
                                          "within_image_power",
                                          "directional_power", "fdr"),
                              evaluation_mask = NULL, verbose = TRUE) {
  finish_timer <- .sptrends_elapsed_timer("benchmark_methods()", verbose)
  on.exit(finish_timer(), add = TRUE)
  stage <- match.arg(stage)
  if (stage == "detection") stage <- "trend_test"
  if (!is.list(scenarios) || length(scenarios) == 0L) {
    stop("'scenarios' must be a non-empty list of argument lists.")
  }
  if (is.null(names(scenarios)) || any(names(scenarios) == "")) {
    stop("'scenarios' must be named.")
  }
  if (!is.list(methods) || length(methods) == 0L ||
      is.null(names(methods)) || any(names(methods) == "") ||
      !all(vapply(methods, is.function, logical(1)))) {
    stop("'methods' must be a non-empty named list of functions.")
  }
  if (length(n_replicates) != 1L || !is.numeric(n_replicates) ||
      is.na(n_replicates) || n_replicates < 1 ||
      n_replicates != floor(n_replicates)) {
    stop("'n_replicates' must be one positive integer.")
  }
  if (stage == "custom" && !is.function(evaluator)) {
    stop("'evaluator' must be supplied for stage = 'custom'.")
  }
  if (!is.null(prepare) && !is.function(prepare)) {
    stop("'prepare' must be NULL or a function.")
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed)) {
    stop("'seed' must be one finite numeric value.")
  }

  # Validate every scenario before starting potentially expensive
  # simulations. Scenarios are not required to share the same argument
  # names (a null scenario need not set spatial_rho; different
  # spatial_model choices take different parameters) -- without this
  # upfront check, a mismatch between scenarios' own argument names
  # only surfaced as a do.call(rbind, ...) column-mismatch error at
  # the very end, after every scenario/replicate/method combination
  # had already run.
  reserved_scenario_names <- c(
    "Scenario", "Replicate", "Seed", "Elapsed", "Method",
    "DomainCells", "TrueNulls", "Pi0", "SignalProportion"
  )
  for (scenario_name in names(scenarios)) {
    scenario_arguments <- scenarios[[scenario_name]]
    if (!is.list(scenario_arguments)) {
      stop(sprintf("Scenario '%s' is not an argument list.",
                   scenario_name))
    }
    if (length(scenario_arguments) > 0L &&
        (is.null(names(scenario_arguments)) ||
         any(names(scenario_arguments) == "") ||
         anyDuplicated(names(scenario_arguments)))) {
      stop(sprintf("Scenario '%s' must contain uniquely named arguments.",
                   scenario_name))
    }
    overlap_reserved <- intersect(names(scenario_arguments),
                                  reserved_scenario_names)
    if (length(overlap_reserved) > 0L) {
      stop(sprintf("Scenario '%s' uses reserved argument name(s): %s.",
                   scenario_name, paste(overlap_reserved, collapse = ", ")))
    }
  }
  # Common schema for heterogeneous scenario argument lists -- every
  # scenario's own column set gets padded to this shared union later,
  # rather than requiring identical argument names across scenarios.
  scenario_argument_names <- unique(unlist(
    lapply(scenarios, names), use.names = FALSE
  ))

  set.seed(seed)
  replicate_seeds <- sample.int(.Machine$integer.max,
                                length(scenarios) * n_replicates)
  rows <- list()
  row_index <- 0L
  seed_index <- 0L
  progress_index <- 0L
  progress <- .sptrends_progress(
    length(scenarios) * n_replicates * length(methods),
    "Benchmarking methods", verbose
  )
  on.exit(.sptrends_progress_close(progress), add = TRUE)
  for (scenario_name in names(scenarios)) {
    arguments <- scenarios[[scenario_name]]
    for (replicate_index in seq_len(n_replicates)) {
      seed_index <- seed_index + 1L
      simulation_arguments <- utils::modifyList(
        arguments, list(seed = replicate_seeds[[seed_index]])
      )
      simulator_formals <- names(formals(simulator))
      if ("verbose" %in% simulator_formals &&
          !"verbose" %in% names(simulation_arguments)) {
        simulation_arguments$verbose <- FALSE
      }
      simulation <- do.call(simulator, simulation_arguments)
      shared_input <- if (is.function(prepare)) {
        prepare(simulation$series, simulation)
      } else {
        simulation$series
      }
      outputs <- list()
      timings <- numeric(length(methods))
      names(timings) <- names(methods)
      for (method_name in names(methods)) {
        started <- proc.time()[["elapsed"]]
        outputs[[method_name]] <- methods[[method_name]](
          shared_input, simulation)
        timings[[method_name]] <- proc.time()[["elapsed"]] - started
        progress_index <- progress_index + 1L
        .sptrends_progress_step(progress, progress_index)
      }
      current_mask <- if (is.function(evaluation_mask)) {
        evaluation_mask(simulation)
      } else {
        evaluation_mask
      }
      scores <- .benchmark_score(outputs, simulation, stage, evaluator,
                                 metrics, current_mask)
      .validate_benchmark_scores(scores, names(methods))
      scenario_data <- .benchmark_scenario_data(arguments, simulation)
      missing_arguments <- setdiff(scenario_argument_names,
                                   names(scenario_data))
      for (missing_name in missing_arguments) {
        scenario_data[[missing_name]] <- NA
      }
      truth_fields_order <- c("DomainCells", "TrueNulls", "Pi0",
                              "SignalProportion")
      scenario_data <- scenario_data[
        c(scenario_argument_names,
          intersect(truth_fields_order, names(scenario_data)))
      ]
      overlap <- intersect(names(scenario_data), names(scores))
      if (length(overlap) > 0L) {
        stop("Scenario and score columns overlap: ",
             paste(overlap, collapse = ", "), ".")
      }
      for (method_name in names(methods)) {
        row_index <- row_index + 1L
        method_score <- scores[scores$Method == method_name, , drop = FALSE]
        rows[[row_index]] <- cbind(
          data.frame(Scenario = scenario_name,
                     Replicate = replicate_index,
                     Seed = replicate_seeds[[seed_index]],
                     Elapsed = timings[[method_name]]),
          scenario_data,
          method_score)
      }
    }
  }
  .sptrends_progress_close(progress)
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  truth_fields <- c(
    "DomainCells", "TrueNulls", "Pi0", "SignalProportion")
  scenario_fields <- scenario_argument_names
  attr(result, "scenario_fields") <- scenario_fields
  attr(result, "truth_fields") <- truth_fields
  attr(result, "stage") <- stage
  class(result) <- c("sptrends_benchmark", "sptrends", "data.frame")
  result
}

#' Retain scenario factors and known-truth composition as columns
#' @noRd
.benchmark_scenario_data <- function(arguments, simulation) {
  reserved <- c(
    "Scenario", "Replicate", "Seed", "Elapsed", "Method",
    "DomainCells", "TrueNulls", "Pi0", "SignalProportion")
  if (length(intersect(names(arguments), reserved)) > 0L) {
    stop("Scenario arguments use a reserved benchmark column name.")
  }
  values <- lapply(arguments, function(value) {
    if (inherits(value, "SpatRaster")) {
      return(sprintf("<SpatRaster[%d cells]>", terra::ncell(value)))
    }
    if (is.matrix(value)) {
      return(sprintf("<matrix[%dx%d]>", nrow(value), ncol(value)))
    }
    if (length(value) == 1L && !is.list(value)) return(value)
    if (is.atomic(value) && length(value) <= 4L) {
      return(paste(value, collapse = ","))
    }
    sprintf("<%s[%d]>", class(value)[1L], length(value))
  })
  truth <- if (!is.null(simulation$true_signal)) {
    if (inherits(simulation$true_signal, "SpatRaster")) {
      terra::values(simulation$true_signal, mat = FALSE)
    } else {
      as.vector(simulation$true_signal)
    }
  } else if (!is.null(simulation$true_slope)) {
    if (inherits(simulation$true_slope, "SpatRaster")) {
      terra::values(simulation$true_slope, mat = FALSE) != 0
    } else {
      as.vector(simulation$true_slope) != 0
    }
  } else {
    NULL
  }
  if (is.null(truth)) {
    m <- if (inherits(simulation$series, "SpatRaster")) {
      terra::ncell(simulation$series)
    } else {
      NA_integer_
    }
    m1 <- NA_integer_
  } else {
    valid <- !is.na(truth)
    m <- sum(valid)
    m1 <- sum(truth[valid] != 0)
  }
  derived <- list(
    DomainCells = m,
    TrueNulls = if (is.na(m1)) NA_integer_ else m - m1,
    Pi0 = if (is.na(m1) || m == 0L) NA_real_ else (m - m1) / m,
    SignalProportion = if (is.na(m1) || m == 0L) {
      NA_real_
    } else {
      m1 / m
    })
  as.data.frame(c(values, derived), stringsAsFactors = FALSE)
}

#' Validate the method rows returned by a benchmark scorer
#' @noRd
.validate_benchmark_scores <- function(scores, method_names) {
  if (!is.data.frame(scores) || !"Method" %in% names(scores)) {
    stop("A benchmark evaluator must return a data frame with Method.")
  }
  if (anyDuplicated(scores$Method) ||
      !setequal(as.character(scores$Method), method_names)) {
    stop("Benchmark scores must contain exactly one row per method.")
  }
  invisible(TRUE)
}

#' @noRd
.benchmark_score <- function(outputs, simulation, stage, evaluator, metrics,
                             evaluation_mask) {
  if (is.function(evaluator)) return(evaluator(outputs, simulation))
  detection_stages <- c("trend_test", "fdr", "fwer")
  if (stage %in% detection_stages) {
    result <- compare_detections(
      outputs, simulation$true_signal,
      metrics = setdiff(metrics, "fwer"),
      truth_direction = simulation$true_direction,
      evaluation_mask = evaluation_mask)
    discoveries <- result$TP + result$FP
    result$FalseDiscoveryProportion <- ifelse(
      discoveries == 0, 0, result$FP / discoveries)
    result$AnyFalsePositive <- as.numeric(result$FP > 0)
    return(result)
  }
  if (stage == "prewhitening") {
    return(.benchmark_prewhitening(outputs, simulation))
  }
  if (is.null(simulation$true_slope)) {
    stop("stage = 'slope' requires the simulator's output to include ",
         "'true_slope' (missing here -- check your 'simulator' function).")
  }
  truth <- terra::values(simulation$true_slope, mat = FALSE)
  rows <- lapply(names(outputs), function(method_name) {
    estimate <- outputs[[method_name]]
    if (inherits(estimate, "SpatRaster")) {
      estimate <- terra::values(estimate, mat = FALSE)
    }
    estimate <- as.vector(estimate)
    if (length(estimate) != length(truth)) {
      stop(sprintf("Slope output '%s' and true_slope differ in length.",
                   method_name))
    }
    ok <- is.finite(estimate) & is.finite(truth)
    error <- estimate[ok] - truth[ok]
    data.frame(
      Method = method_name,
      Bias = mean(error), MAE = mean(abs(error)),
      RMSE = sqrt(mean(error^2)),
      DirectionError = mean(sign(estimate[ok]) != sign(truth[ok])))
  })
  do.call(rbind, rows)
}

#' Score transformed series against known temporal structure
#' @noRd
.benchmark_prewhitening <- function(outputs, simulation) {
  if (is.null(simulation$true_slope)) {
    stop("stage = 'prewhitening' requires the simulator's output to ",
         "include 'true_slope' (missing here -- check your 'simulator' ",
         "function).")
  }
  true_slope <- terra::values(simulation$true_slope, mat = FALSE)
  rows <- lapply(names(outputs), function(method_name) {
    output <- outputs[[method_name]]
    transformed <- if (is.list(output) && !is.null(output$series)) {
      output$series
    } else {
      output
    }
    if (!inherits(transformed, "SpatRaster")) {
      stop(sprintf(
        "Prewhitening output '%s' must contain a SpatRaster series.",
        method_name))
    }
    values <- terra::values(transformed, mat = TRUE)
    if (nrow(values) != length(true_slope)) {
      stop(sprintf(
        "Prewhitening output '%s' and truth differ in cell count.",
        method_name))
    }
    if (ncol(values) < 3L) {
      stop(sprintf(
        "Prewhitening output '%s' needs at least three time steps.",
        method_name))
    }
    time <- seq_len(ncol(values))
    centred_time <- time - mean(time)
    denominator <- sum(centred_time^2)
    estimated_slope <- as.vector(values %*% centred_time / denominator)
    residuals <- values - outer(estimated_slope, time)
    lag1 <- vapply(seq_len(nrow(residuals)), function(index) {
      series <- residuals[index, ]
      if (any(!is.finite(series)) || stats::sd(series) == 0) {
        return(NA_real_)
      }
      stats::cor(series[-length(series)], series[-1L])
    }, numeric(1))
    slope_error <- estimated_slope - true_slope
    modified_fraction <- NA_real_
    if (is.list(output) && inherits(output$diagnostics, "SpatRaster") &&
        "Modified" %in% names(output$diagnostics)) {
      modified <- terra::values(
        output$diagnostics$Modified, mat = FALSE)
      modified_fraction <- mean(modified != 0, na.rm = TRUE)
    }
    data.frame(
      Method = method_name,
      ResidualACF1 = .benchmark_mean(lag1),
      AbsResidualACF1 = .benchmark_mean(abs(lag1)),
      SlopeBias = .benchmark_mean(slope_error),
      SlopeMAE = .benchmark_mean(abs(slope_error)),
      SlopeRMSE = sqrt(.benchmark_mean(slope_error^2)),
      ModifiedFraction = modified_fraction,
      OutputLength = ncol(values))
  })
  do.call(rbind, rows)
}

#' Summarise a method benchmark across Monte Carlo replicates
#'
#' Aggregates the replicate-level output of [benchmark_methods()] while
#' preserving scenarios and methods. For detection-like stages, empirical FDR
#' is the mean false-discovery proportion and empirical FWER is the proportion
#' of replicates containing at least one false positive.
#'
#' **Function type:** **Benchmarking function** -- summarises known-truth
#' experiments; it does not perform statistical inference on user data.
#'
#' @section Typical use:
#' Run [benchmark_methods()] and pass its result directly to this function.
#'
#' @section Methodological details:
#' **Monte Carlo aggregation**
#'
#' Replicates, rather than raster cells, are the independent Monte Carlo units.
#' Therefore FDR and FWER are aggregated across replicate-level false
#' discovery proportions and false-positive indicators, respectively.
#'
#' **Computational considerations**
#'
#' Aggregation operates on the retained benchmark table and does not rerun
#' simulations or methods.
#'
#' **Limitations**
#'
#' Summary precision depends on the number of independent replicates and the
#' range of scenarios evaluated. It does not generalise beyond those designs.
#'
#' **Quality assurance**
#'
#' Tests verify grouping, scenario retention, means, standard deviations,
#' empirical FDR and FWER, CSV output and invalid input handling. The retained
#' external validation also verifies the summaries against 4,000 paired MK
#' replicates and 33 independent simulation-cycle controls.
#'
#' @param x Replicate-level result returned by [benchmark_methods()].
#' @param path Character or `NULL`. If supplied, the summary is written as a
#'   CSV file at this path.
#' @param verbose Logical. If `TRUE`, reports progress, elapsed time and the
#'   estimated time remaining while scenario-method groups are summarised.
#'
#' @return A data frame with one row per scenario and method, retained
#'   scenario factors, the number of
#'   replicates, means and standard deviations of numerical metrics, and
#'   `EmpiricalFDR`/`EmpiricalFWER` when detection metrics are present.
#'
#' @examples
#' x <- structure(
#'   data.frame(Scenario = rep("null", 2), Replicate = 1:2,
#'              Method = rep("method", 2), FP = c(1, 0),
#'              FalseDiscoveryProportion = c(1, 0),
#'              AnyFalsePositive = c(1, 0)),
#'   class = c("sptrends_benchmark", "data.frame"))
#' benchmark_summary(x, verbose = FALSE)
#'
#' @family validation functions
#' @seealso [benchmark_methods()], [compare_detections()]
#' @export
benchmark_summary <- function(x, path = NULL, verbose = TRUE) {
  finish_timer <- .sptrends_elapsed_timer("benchmark_summary()", verbose)
  on.exit(finish_timer(), add = TRUE)
  if (!inherits(x, "sptrends_benchmark") || !is.data.frame(x)) {
    stop("'x' must be a result returned by benchmark_methods().")
  }
  required <- c("Scenario", "Method", "Replicate")
  if (!all(required %in% names(x))) {
    stop("Benchmark result lacks Scenario, Method, or Replicate columns.")
  }
  scenario_fields <- attr(x, "scenario_fields")
  scenario_fields <- intersect(scenario_fields, names(x))
  group_names <- c("Scenario", scenario_fields, "Method")
  groups <- unique(x[, group_names, drop = FALSE])
  progress <- .sptrends_progress(
    nrow(groups), "Summarising benchmark groups", verbose)
  on.exit(.sptrends_progress_close(progress), add = TRUE)
  rows <- lapply(seq_len(nrow(groups)), function(index) {
    selected <- x$Scenario == groups$Scenario[index] &
      x$Method == groups$Method[index]
    subset <- x[selected, , drop = FALSE]
    numeric_names <- names(subset)[vapply(subset, is.numeric, logical(1))]
    numeric_names <- setdiff(
      numeric_names, c("Replicate", "Seed", scenario_fields))
    row <- groups[index, , drop = FALSE]
    row$n_replicates <- nrow(subset)
    for (name in numeric_names) {
      row[[paste0(name, "_mean")]] <- .benchmark_mean(subset[[name]])
      row[[paste0(name, "_sd")]] <- .benchmark_sd(subset[[name]])
    }
    if ("FalseDiscoveryProportion" %in% names(subset)) {
      row$EmpiricalFDR <- mean(subset$FalseDiscoveryProportion)
    }
    if ("AnyFalsePositive" %in% names(subset)) {
      row$EmpiricalFWER <- mean(subset$AnyFalsePositive)
    }
    .sptrends_progress_step(progress, index)
    row
  })
  .sptrends_progress_close(progress)
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  class(result) <- c("sptrends_benchmark_summary", "data.frame")
  if (!is.null(path)) utils::write.csv(result, path, row.names = FALSE)
  result
}

#' @noRd
.benchmark_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

#' @noRd
.benchmark_sd <- function(x) {
  if (sum(!is.na(x)) < 2L) return(NA_real_)
  stats::sd(x, na.rm = TRUE)
}
