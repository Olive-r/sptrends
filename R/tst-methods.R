#' @noRd
.print_tst <- function(x, ...) {
  cat("<True Significant Trends (TST) result>\n")

  if (!is.null(x$prewhiten)) {
    if (identical(x$prewhiten$method, "TFPW_Y")) {
      rho_vals <- terra::values(x$prewhiten$diagnostics$Rho, mat = FALSE)
      n_valid <- sum(!is.na(rho_vals))
      cat(sprintf("Prewhitening: Yue-Pilon, all %d valid cells\n", n_valid))
    } else {
      mod <- terra::values(x$prewhiten$diagnostics$Modified, mat = FALSE)
      n_valid <- sum(!is.na(mod))
      n_mod <- sum(mod == 1, na.rm = TRUE)
      cat(sprintf("Prewhitening: %d of %d cells modified (%.1f%%)\n",
                  n_mod, n_valid, 100 * n_mod / n_valid))
    }
  } else {
    cat("Prewhitening: skipped\n")
  }

  s_name <- if ("Sm" %in% names(x$trend)) {
    "Sm"
  } else if ("S" %in% names(x$trend)) {
    "S"
  } else {
    "beta"
  }
  p_vals <- terra::values(x$trend$p, mat = FALSE)
  n_cells <- sum(!is.na(p_vals))
  cat(sprintf("Trend test: %d cells (%s statistic)\n", n_cells, s_name))

  if (!is.null(x$theil_sen)) {
    slope_vals <- terra::values(x$theil_sen, mat = FALSE)
    cat(sprintf("Theil-Sen slope: median %.4g (range %.4g to %.4g)\n",
                stats::median(slope_vals, na.rm = TRUE),
                min(slope_vals, na.rm = TRUE), max(slope_vals, na.rm = TRUE)))
  }

  if (!is.null(x$fdr)) {
    if (!is.null(x$fdr$reject_BH)) {
      n_sig <- sum(x$fdr$reject_BH, na.rm = TRUE)
      cat(sprintf("Significant after FDR-BH: %d (%.1f%%)\n",
                  n_sig, 100 * n_sig / n_cells))
    }
    if (!is.null(x$fdr$reject_BKY)) {
      n_sig <- sum(x$fdr$reject_BKY, na.rm = TRUE)
      cat(sprintf("Significant after FDR-BKY: %d (%.1f%%)\n",
                  n_sig, 100 * n_sig / n_cells))
    }
    if (!is.null(x$fdr$moran_assessment)) {
      cat(sprintf("Moran's I assessment: %s\n",
                  x$fdr$moran_assessment))
    }
  } else {
    cat("FDR correction: skipped\n")
  }

  cat("Use summary() for details, plot() for a map.\n")
  invisible(x)
}

#' @noRd
.summary_tst <- function(object, ...) {
  cat("=== Trend test (uncorrected) ===\n")
  print(object$trend_summary_table)

  if (!is.null(object$theil_sen)) {
    cat("\n=== Theil-Sen slope ===\n")
    print(summary(terra::values(object$theil_sen, mat = FALSE)))
  }

  fdr_tab <- NULL
  if (!is.null(object$fdr)) {
    cat("\n=== FDR correction ===\n")
    fdr_tab <- fdr_summary(object$fdr)
  }

  invisible(list(trend = object$trend_summary_table, fdr = fdr_tab))
}

#' @noRd
.plot_tst <- function(x, which = c("direction", "significance", "trend",
                                   "slope", "slope_map", "slope_direction",
                                   "slope_hist", "slope_bar", "pvalue_map",
                                   "pvalue_significance", "pvalue_hist",
                                   "pvalue_bar"),
                      method = c("BKY", "BH"), smooth = TRUE, ...) {
  which <- match.arg(which)
  method <- match.arg(method)

  if (which == "direction") {
    if (is.null(x$fdr)) {
      stop("No FDR correction in this workflow_tst object -- rerun ",
           "workflow_tst() with fdr_method set, or use which = 'trend'.")
    }
    direction <- direction_map(x$trend, x$fdr, method = method,
                                    verbose = FALSE)
    fdr_direction_plot(direction)
  } else if (which == "significance") {
    if (is.null(x$fdr)) {
      stop("No FDR correction in this workflow_tst object -- rerun ",
           "workflow_tst() with fdr_method set, or use which = 'trend'.")
    }
    fdr_significance_maps(x$fdr)
  } else if (which == "slope") {
    if (is.null(x$theil_sen)) {
      stop("No Theil-Sen slope in this workflow_tst object -- rerun ",
           "workflow_tst() with theil_sen = TRUE.")
    }
    if (is.null(x$fdr)) {
      stop("No FDR correction in this workflow_tst object -- rerun ",
           "workflow_tst() with fdr_method set.")
    }
    reject <- if (method == "BKY") x$fdr$reject_BKY else x$fdr$reject_BH
    if (is.null(reject)) {
      stop(sprintf(
        paste0("FDR-%s is not available in this workflow_tst object -- ",
               "rerun workflow_tst() with that method included, or set ",
               "'method' to match."),
        method))
    }

    # x$theil_sen may have been computed with smoothing applied at
    # source, if theil_sen_args = list(smooth_neighbourhood = TRUE) was
    # explicitly requested (not workflow_tst()'s own default -- see
    # ?workflow_tst) -- do
    # not smooth it a second time here, and be honest in the subtitle
    # about which stage did it (or, if smooth = FALSE was requested but
    # source smoothing was already applied, that raw values were never
    # kept to go back to).
    already_smoothed <- isTRUE(x$theil_sen_smoothed)
    apply_smooth_now <- isTRUE(smooth) && !already_smoothed

    slope_sig <- .mask_and_smooth_slope(x$theil_sen, reject, apply_smooth_now)

    if (already_smoothed && !isTRUE(smooth)) {
      message(
        "Note: smooth = FALSE was requested, but this slope was already ",
        "computed with smoothing applied at source (theil_sen_args = ",
        "list(smooth_neighbourhood = TRUE) was requested) -- the unsmoothed ",
        "values were never kept to show instead. Rerun with theil_sen_args ",
        "= list(smooth_neighbourhood = FALSE) for a genuinely unsmoothed ",
        "result.")
    }
    subtitle <- sprintf("Theil-Sen slope, FDR-%s significant cells", method)

    range_lim <- .robust_diverging_range(slope_sig)
    # fill_range = TRUE: without it, cells beyond the robust range would
    # render as blank (NA) instead of saturating to the extreme colour --
    # see slope_map()'s own comment on this for the full reasoning.
    terra::plot(slope_sig, col = grDevices::hcl.colors(50, "Blue-Red 3"),
                range = range_lim, fill_range = TRUE, main = subtitle)
  } else if (which %in% c("slope_map", "slope_direction", "slope_hist",
                           "slope_bar")) {
    if (is.null(x$theil_sen)) {
      stop("No Theil-Sen slope in this workflow_tst object -- rerun ",
           "workflow_tst() with theil_sen = TRUE.")
    }
    switch(which,
      slope_map       = slope_map(x$theil_sen),
      slope_direction = slope_direction_map(x$theil_sen, ...),
      slope_hist      = slope_histogram(x$theil_sen, ...),
      slope_bar       = slope_direction_barplot(x$theil_sen, ...)
    )
  } else if (which %in% c("pvalue_map", "pvalue_significance", "pvalue_hist",
                           "pvalue_bar")) {
    switch(which,
      pvalue_map           = .plot_pvalue_map(x$trend, ...),
      pvalue_significance  = .plot_pvalue_significance(x$trend, ...),
      pvalue_hist          = .plot_pvalue_hist(x$trend, ...),
      pvalue_bar           = .plot_pvalue_bar(x$trend, ...)
    )
  } else {
    trend_maps(x$trend)
  }

  invisible(x)
}
