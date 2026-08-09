#' @noRd
.print_fdr <- function(x, ...) {
  cat(sprintf("<FDR correction result (%s)>\n",
              paste(x$method, collapse = ", ")))
  n_valid <- sum(!is.na(x$p))
  cat(sprintf("Valid cells: %d | raw: %d significant\n",
              n_valid, sum(x$reject_raw, na.rm = TRUE)))
  if (!is.null(x$reject_BH)) {
    cat(sprintf("FDR-BH: %d significant\n", sum(x$reject_BH, na.rm = TRUE)))
  }
  if (!is.null(x$reject_BKY)) {
    cat(sprintf("FDR-BKY: %d significant\n", sum(x$reject_BKY, na.rm = TRUE)))
  }
  if (!is.null(x$reject_BY)) {
    cat(sprintf("FDR-BY: %d significant\n", sum(x$reject_BY, na.rm = TRUE)))
  }
  invisible(x)
}

#' @noRd
.summary_fdr <- function(object, ...) {
  fdr_summary(object, ...)
}

#' @noRd
.plot_fdr <- function(x, which = c("significance", "pvalue_histogram",
                                   "comparison", "threshold"), ...) {
  which <- match.arg(which)
  if (which == "significance") {
    if (is.null(x$rasters)) {
      stop("'x' has no rasters -- fdr_correction() must have been run on ",
           "a SpatRaster to draw significance maps.")
    }
    fdr_significance_maps(x)
  } else if (which == "pvalue_histogram") {
    fdr_pvalue_histogram(x$p)
  } else if (which == "comparison") {
    fdr_comparison_barplot(x)
  } else {
    if (!"BKY" %in% x$method) {
      stop("'x' has no BKY results -- fdr_correction() must have been ",
           "run with method including \"BKY\" to draw the threshold plot.")
    }
    fdr_threshold_plot(x)
  }
  invisible(x)
}
