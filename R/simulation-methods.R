# S3 presentation methods for simulation and benchmarking -----------------

#' @noRd
.print_sptrends_simulation <- function(x, ...) {
  p <- x$parameters
  truth <- terra::values(x$true_signal, mat = FALSE)
  valid <- !is.na(truth)
  cat("<sptrends simulation>\n")
  cat(sprintf("Grid: %d x %d | time steps: %d | valid cells: %d\n",
              p$nrow, p$ncol, p$n_time, sum(valid)))
  cat(sprintf(
    "Temporal AR(1): %.3g | spatial model: %s | signal: %.1f%%\n",
    p$ar1, p$spatial_model,
    100 * mean(truth[valid] != 0)))
  invisible(x)
}

#' @noRd
.summary_sptrends_simulation <- function(object, path = NULL, ...) {
  slope <- terra::values(object$true_slope, mat = FALSE)
  signal <- terra::values(object$true_signal, mat = FALSE)
  valid <- !is.na(signal)
  m <- sum(valid)
  m1 <- sum(signal[valid] != 0)
  result <- data.frame(
    DomainCells = m,
    TrueSignals = m1,
    PctTrueSignals = if (m == 0L) NA_real_ else round(100 * m1 / m, 2),
    TrueNulls = m - m1,
    Pi0 = if (m == 0L) NA_real_ else (m - m1) / m,
    MinimumSlope = min(slope, na.rm = TRUE),
    MaximumSlope = max(slope, na.rm = TRUE),
    AR1 = object$parameters$ar1,
    SpatialModel = object$parameters$spatial_model)
  print.data.frame(result, row.names = FALSE)
  if (!is.null(path)) utils::write.csv(result, path, row.names = FALSE)
  invisible(result)
}

#' @noRd
.plot_sptrends_simulation <- function(
    x, which = c("truth", "slope", "direction", "breaks", "series"),
    path = NULL, ...) {
  which <- match.arg(which)
  plot_chosen <- function() {
    if (which == "truth") {
      .safe_categorical_plot(
        x$true_signal, values = c(0, 1),
        colours = c("grey90", "#0072B2"),
        labels = c("No trend", "True trend"),
        main = "Known trend region", ...)
    } else if (which == "slope") {
      terra::plot(x$true_slope, main = "Known slope", ...)
    } else if (which == "direction") {
      .safe_categorical_plot(
        x$true_direction, values = c(-1, 0, 1),
        colours = c("#D55E00", "grey90", "#009E73"),
        labels = c("Decrease", "No trend", "Increase"),
        main = "Known trend direction", ...)
    } else if (which == "breaks") {
      .safe_categorical_plot(
        x$true_break, values = c(0, 1),
        colours = c("grey90", "#CC79A7"),
        labels = c("No break", "True break"),
        main = "Known break region", ...)
    } else {
      terra::plot(x$series, main = names(x$series), ...)
    }
  }
  plot_chosen()
  if (!is.null(path)) .save_current_plot(path, plot_chosen)
  invisible(x)
}

#' @noRd
.print_sptrends_simulation_design <- function(x, ...) {
  cat("<sptrends simulation design>\n")
  cat(sprintf("Scenarios: %d | varied factors: %s\n",
              length(x), paste(attr(x, "factors"), collapse = ", ")))
  constants <- attr(x, "constants")
  if (length(constants) > 0L) {
    cat("Shared settings:", paste(constants, collapse = ", "), "\n")
  }
  invisible(x)
}

#' @noRd
.summary_sptrends_simulation_design <- function(object, path = NULL, ...) {
  factors <- attr(object, "factors")
  result <- data.frame(
    Factor = factors,
    Levels = vapply(factors, function(name) {
      length(unique(vapply(object, function(scenario) {
        paste(scenario[[name]], collapse = ",")
      }, character(1))))
    }, integer(1)))
  print.data.frame(result, row.names = FALSE)
  if (!is.null(path)) utils::write.csv(result, path, row.names = FALSE)
  invisible(result)
}

#' @noRd
.plot_sptrends_simulation_design <- function(x, path = NULL, ...) {
  factors <- attr(x, "factors")
  result <- data.frame(
    Factor = factors,
    Levels = vapply(factors, function(name) {
      length(unique(vapply(x, function(scenario) {
        paste(scenario[[name]], collapse = ",")
      }, character(1))))
    }, integer(1)))
  plot_bars <- function() {
    graphics::barplot(
      result$Levels, names.arg = result$Factor,
      ylab = "Number of levels", main = "Simulation design", ...)
  }
  plot_bars()
  if (!is.null(path)) .save_current_plot(path, plot_bars)
  invisible(x)
}

#' @noRd
.print_sptrends_benchmark <- function(x, ...) {
  cat("<sptrends benchmark>\n")
  cat(sprintf(
    "Stage: %s | scenarios: %d | methods: %d | replicates: %d\n",
    attr(x, "stage"), length(unique(x$Scenario)),
    length(unique(x$Method)), length(unique(x$Replicate))))
  cat(sprintf("Rows: %d | elapsed method time: %.3f s\n",
              nrow(x), sum(x$Elapsed, na.rm = TRUE)))
  invisible(x)
}

#' @noRd
.summary_sptrends_benchmark <- function(object, ...) {
  result <- benchmark_summary(object, ...)
  print.data.frame(result, row.names = FALSE)
  invisible(result)
}

#' Aggregate one benchmark metric for plotting
#' @noRd
.benchmark_plot_summary <- function(x, metric, scenario, group, facet) {
  columns <- unique(c(scenario, group, facet))
  keys <- interaction(x[, columns, drop = FALSE], drop = TRUE, lex.order = TRUE)
  indices <- split(seq_len(nrow(x)), keys)
  rows <- lapply(indices, function(index) {
    values <- as.numeric(x[[metric]][index])
    row <- x[index[1L], columns, drop = FALSE]
    row$Mean <- .benchmark_mean(values)
    row$SD <- .benchmark_sd(values)
    row$N <- sum(!is.na(values))
    row
  })
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Resolve aliases used by Monte Carlo summaries
#' @noRd
.benchmark_plot_metric <- function(x, metric) {
  aliases <- c(
    EmpiricalFDR = "FalseDiscoveryProportion",
    EmpiricalFWER = "AnyFalsePositive")
  resolved <- if (metric %in% names(aliases)) aliases[[metric]] else metric
  if (!resolved %in% names(x) || !is.numeric(x[[resolved]])) {
    stop("Unknown numerical benchmark metric: '", metric, "'.")
  }
  resolved
}

#' Draw benchmark performance curves or grouped bars
#' @noRd
.plot_benchmark_panel <- function(
    x, metric, metric_label, scenario, group, facet_value,
    type, interval, level, ...) {
  facet <- attr(x, "plot_facet")
  data <- if (is.null(facet)) {
    x
  } else {
    x[x[[facet]] == facet_value, , drop = FALSE]
  }
  summary <- .benchmark_plot_summary(
    data, metric, scenario, group, facet = NULL)
  x_values <- unique(data[[scenario]])
  if (is.numeric(x_values)) x_values <- sort(x_values)
  groups <- unique(data[[group]])
  colours <- grDevices::hcl.colors(length(groups), "Dark 3")
  standard_error <- summary$SD / sqrt(pmax(summary$N, 1))
  spread <- if (interval == "sd") summary$SD else standard_error
  multiplier <- if (interval == "ci") {
    value <- stats::qt(1 - (1 - level) / 2, df = pmax(summary$N - 1, 1))
    value[summary$N < 2L] <- NA_real_
    value
  } else if (interval == "none") {
    0
  } else {
    1
  }
  if (interval == "none") {
    summary$Lower <- summary$Mean
    summary$Upper <- summary$Mean
  } else {
    summary$Lower <- summary$Mean - multiplier * spread
    summary$Upper <- summary$Mean + multiplier * spread
  }
  rate_metrics <- c(
    "Sensitivity", "Specificity", "Precision", "Accuracy", "F1",
    "FPR", "FDR", "TypeI", "TypeII", "TypeIII", "FieldPower",
    "GlobalPower", "WithinImagePower", "DirectionalPower",
    "FalseDiscoveryProportion", "AnyFalsePositive")
  if (metric %in% rate_metrics) {
    summary$Lower <- pmax(0, summary$Lower)
    summary$Upper <- pmin(1, summary$Upper)
  }
  finite_bounds <- c(summary$Lower, summary$Upper)
  finite_bounds <- finite_bounds[is.finite(finite_bounds)]
  y_range <- if (length(finite_bounds) == 0L) {
    c(0, 1)
  } else {
    range(finite_bounds)
  }
  if (diff(y_range) == 0) {
    margin <- if (y_range[1L] == 0) 0.05 else abs(y_range[1L]) * 0.05
    y_range <- y_range + c(-margin, margin)
  }
  title <- metric_label
  if (!is.null(facet)) title <- paste0(title, " | ", facet_value)

  if (type == "bar") {
    y_range <- range(c(0, y_range))
    means <- matrix(
      NA_real_, nrow = length(groups), ncol = length(x_values),
      dimnames = list(groups, as.character(x_values)))
    lower <- upper <- means
    for (i in seq_len(nrow(summary))) {
      row <- match(summary[[group]][i], groups)
      column <- match(summary[[scenario]][i], x_values)
      means[row, column] <- summary$Mean[i]
      lower[row, column] <- summary$Lower[i]
      upper[row, column] <- summary$Upper[i]
    }
    positions <- graphics::barplot(
      means, beside = TRUE, col = colours, border = NA,
      names.arg = x_values, ylim = y_range,
      xlab = scenario, ylab = metric_label, main = title, ...)
    if (interval != "none") {
      drawable <- is.finite(lower) & is.finite(upper) & lower != upper
      if (any(drawable)) {
        graphics::arrows(
          positions[drawable], lower[drawable],
          positions[drawable], upper[drawable],
          angle = 90, code = 3, length = 0.04)
      }
    }
  } else {
    positions <- seq_along(x_values)
    graphics::plot(
      positions, rep(NA_real_, length(positions)), type = "n",
      xaxt = "n", xlab = scenario, ylab = metric_label,
      ylim = y_range, main = title, ...)
    graphics::axis(1, at = positions, labels = x_values)
    for (i in seq_along(groups)) {
      selected <- summary[[group]] == groups[i]
      current <- summary[selected, , drop = FALSE]
      ordering <- match(current[[scenario]], x_values)
      current <- current[base::order(ordering), , drop = FALSE]
      current_x <- match(current[[scenario]], x_values)
      bounds <- c(current$Lower, current$Upper)
      if (interval != "none" && nrow(current) > 1L &&
          all(is.finite(bounds))) {
        graphics::polygon(
          c(current_x, rev(current_x)),
          c(current$Lower, rev(current$Upper)),
          col = grDevices::adjustcolor(colours[i], alpha.f = 0.18),
          border = NA)
      }
      graphics::lines(
        current_x, current$Mean, type = "b", pch = i,
        col = colours[i], lwd = 2)
    }
  }
  graphics::legend(
    "topright", legend = groups, col = colours,
    lty = if (type == "bar") 0 else 1,
    pch = if (type == "bar") 15 else seq_along(groups), bty = "n")
}

#' Draw distributions of replicate-level performance
#' @noRd
.plot_benchmark_boxplot <- function(
    x, metric, metric_label, scenario, group, facet_value, ...) {
  facet <- attr(x, "plot_facet")
  data <- if (is.null(facet)) x else x[x[[facet]] == facet_value, ]
  labels <- interaction(data[[scenario]], data[[group]], sep = "\n")
  title <- metric_label
  if (!is.null(facet)) title <- paste0(title, " | ", facet_value)
  graphics::boxplot(
    split(data[[metric]], labels), las = 2,
    ylab = metric_label, main = title, ...)
}

#' Draw a two-factor performance heatmap
#' @noRd
.plot_benchmark_heatmap <- function(
    x, metric, metric_label, scenario, facet, group_value, ...) {
  data <- x[x[[attr(x, "plot_group")]] == group_value, , drop = FALSE]
  summary <- .benchmark_plot_summary(
    data, metric, scenario, attr(x, "plot_group"), facet)
  x_values <- unique(data[[scenario]])
  y_values <- unique(data[[facet]])
  matrix_values <- matrix(
    NA_real_, nrow = length(x_values), ncol = length(y_values))
  for (i in seq_len(nrow(summary))) {
    matrix_values[
      match(summary[[scenario]][i], x_values),
      match(summary[[facet]][i], y_values)] <- summary$Mean[i]
  }
  graphics::image(
    seq_along(x_values), seq_along(y_values), matrix_values,
    axes = FALSE, xlab = scenario, ylab = facet,
    main = paste(metric_label, "|", group_value),
    col = grDevices::hcl.colors(30, "YlOrRd", rev = TRUE), ...)
  graphics::axis(1, at = seq_along(x_values), labels = x_values)
  graphics::axis(2, at = seq_along(y_values), labels = y_values)
  graphics::box()
}

#' @noRd
.plot_sptrends_benchmark <- function(
    x, metric = NULL, scenario = NULL, group = "Method", facet = NULL,
    type = c("line", "bar", "boxplot", "heatmap", "profile"),
    interval = c("ci", "se", "sd", "none"), level = 0.95,
    path = NULL, ...) {
  type <- match.arg(type)
  interval <- match.arg(interval)
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("'level' must be one number strictly between zero and one.")
  }
  scenario_fields <- attr(x, "scenario_fields")
  varying <- scenario_fields[vapply(
    x[scenario_fields], function(value) length(unique(value)) > 1L,
    logical(1))]
  if (is.null(scenario)) {
    scenario <- if (length(varying) > 0L) varying[1L] else "Scenario"
  }
  for (name in c(scenario, group, facet)) {
    if (!is.null(name) && !name %in% names(x)) {
      stop("Unknown benchmark grouping column: '", name, "'.")
    }
  }
  excluded <- c(
    "Replicate", "Seed", "Elapsed", "DomainCells", "TrueNulls",
    "Pi0", "SignalProportion", "TP", "FP", "TN", "FN",
    scenario_fields)
  candidates <- names(x)[vapply(x, is.numeric, logical(1))]
  candidates <- setdiff(candidates, excluded)
  if (is.null(metric)) {
    if (length(candidates) == 0L) {
      stop("No numerical performance metric is available to plot.")
    }
    metric <- if ("WithinImagePower" %in% candidates) {
      "WithinImagePower"
    } else {
      candidates[1L]
    }
  }
  if (length(metric) > 1L && type != "profile") {
    stop("Multiple metrics require type = 'profile'.")
  }
  if (type == "profile" && length(metric) == 1L) {
    type <- "line"
  }
  resolved <- vapply(metric, function(value) {
    .benchmark_plot_metric(x, value)
  }, character(1))
  facets <- if (is.null(facet)) NA else unique(x[[facet]])
  panels <- if (type == "heatmap") {
    if (is.null(facet)) stop("A heatmap requires a second factor in 'facet'.")
    unique(x[[group]])
  } else {
    rep(facets, times = length(resolved))
  }
  panel_count <- if (type == "heatmap") {
    length(panels)
  } else {
    length(resolved) * length(facets)
  }
  plot_panels <- function() {
    old <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old), add = TRUE)
    if (panel_count > 1L) {
      graphics::par(mfrow = grDevices::n2mfrow(panel_count))
    }
    attr(x, "plot_facet") <- facet
    attr(x, "plot_group") <- group
    if (type == "heatmap") {
      for (group_value in panels) {
        .plot_benchmark_heatmap(
          x, resolved[1L], metric[1L], scenario, facet, group_value, ...)
      }
    } else {
      for (metric_index in seq_along(resolved)) {
        for (facet_value in facets) {
          if (type == "boxplot") {
            .plot_benchmark_boxplot(
              x, resolved[metric_index], metric[metric_index], scenario,
              group, facet_value, ...)
          } else {
            panel_type <- if (type == "profile") "line" else type
            .plot_benchmark_panel(
              x, resolved[metric_index], metric[metric_index], scenario,
              group, facet_value, panel_type, interval, level, ...)
          }
        }
      }
    }
  }
  plot_panels()
  if (!is.null(path)) .save_current_plot(path, plot_panels)
  invisible(x)
}
