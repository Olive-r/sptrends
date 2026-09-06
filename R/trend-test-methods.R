#' @noRd
.print_trend_test <- function(x, ...) {
  variant <- if (isTRUE(x$neighbourhood)) {
    size <- if (is.null(x$window_size)) 3L else x$window_size
    sprintf("Contextual Mann-Kendall (%dx%d)", size, size)
  } else if (identical(x$method, "MMK")) {
    "modified Mann-Kendall"
  } else if ("beta" %in% names(x$stats)) {
    "OLS"
  } else {
    "classic Mann-Kendall"
  }
  cat(sprintf("<%s result>\n", variant))
  p_vals <- terra::values(x$stats$p, mat = FALSE)
  n_valid <- sum(!is.na(p_vals))
  n_sig <- sum(p_vals <= 0.05, na.rm = TRUE)
  cat(sprintf(
    paste0("Cells tested: %d | significant at alpha=0.05 (uncorrected): ",
           "%d (%.1f%%)\n"),
    n_valid, n_sig, 100 * n_sig / n_valid))
  invisible(x)
}

#' @noRd
.summary_trend_test <- function(object, ...) {
  trend_summary(object$stats, ...)
}

#' @noRd
.plot_trend_test <- function(x, which = c("maps", "histograms"),
                              alpha = 0.05, ...) {
  which <- match.arg(which)
  if (which == "maps") {
    trend_maps(x$stats, alpha = alpha, ...)
  } else {
    trend_histograms(x$stats, ...)
  }
  invisible(x)
}
