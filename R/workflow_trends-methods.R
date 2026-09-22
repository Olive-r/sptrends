# Internal methods for workflow_trends objects.
#' @noRd
.print_workflow_trends <- function(x, ...) {
  cat("<workflow_trends result>\n")

  if (!is.null(x$prewhiten)) {
    if (identical(x$prewhiten$method, "TFPW_Y")) {
      rho_vals <- terra::values(x$prewhiten$diagnostics$Rho, mat = FALSE)
      n_valid <- sum(!is.na(rho_vals))
      cat(sprintf("Prewhitening (%s): all %d valid cells\n",
                  x$prewhiten$method, n_valid))
    } else if ("Modified" %in% names(x$prewhiten$diagnostics)) {
      mod <- terra::values(x$prewhiten$diagnostics$Modified, mat = FALSE)
      n_valid <- sum(!is.na(mod))
      n_mod <- sum(mod == 1, na.rm = TRUE)
      cat(sprintf("Prewhitening (%s): %d of %d cells modified (%.1f%%)\n",
                  x$prewhiten$method, n_mod, n_valid, 100 * n_mod / n_valid))
    } else {
      # Unreachable with the 4 prewhiten() methods this package
      # currently offers (TFPW_WS, TFPW_Z, TFPW_Y, VCTFPW): TFPW_Y is
      # caught by the first branch above, and TFPW_WS/TFPW_Z/VCTFPW
      # all now carry their own "Modified" diagnostic field, caught by
      # the second branch. Kept as a defensive fallback for any future
      # prewhiten() method that has neither shape, not dead code to
      # delete -- marked nocov rather than chased with a contrived
      # test that would need to fabricate a diagnostics shape no real
      # prewhiten() call can currently produce.
      cat(sprintf("Prewhitening (%s): done\n", x$prewhiten$method)) # nocov
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

  if (!is.null(x$slope)) {
    slope_vals <- terra::values(x$slope, mat = FALSE)
    cat(sprintf("Slope: median %.4g (range %.4g to %.4g)\n",
                stats::median(slope_vals, na.rm = TRUE),
                min(slope_vals, na.rm = TRUE), max(slope_vals, na.rm = TRUE)))
  }

  if (!is.null(x$fdr)) {
    for (nm in c("reject_BH", "reject_BKY", "reject_BY")) {
      if (!is.null(x$fdr[[nm]])) {
        n_sig <- sum(x$fdr[[nm]], na.rm = TRUE)
        cat(sprintf("Significant after FDR-%s: %d (%.1f%%)\n",
                    sub("^reject_", "", nm), n_sig, 100 * n_sig / n_cells))
      }
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
.summary_workflow_trends <- function(object, ...) {
  fdr_tab <- NULL
  if (!is.null(object$fdr)) {
    cat("=== FDR correction (the actual result) ===\n")
    fdr_tab <- fdr_summary(object$fdr)
  }

  if (!is.null(object$slope)) {
    cat("\n=== Slope ===\n")
    print(summary(terra::values(object$slope, mat = FALSE)))
  }

  if (!is.null(object$fdr)) {
    cat("\n=== Trend test, uncorrected (diagnostic only -- not the final result; see FDR correction above) ===\n")
  } else {
    cat("=== Trend test (uncorrected -- no FDR correction was run for this workflow) ===\n")
  }
  print(object$trend_summary_table)

  invisible(list(trend = object$trend_summary_table, fdr = fdr_tab))
}

#' @noRd
.plot_workflow_trends <- function(x, which = c("direction", "significance",
                                                 "trend", "slope",
                                                 "slope_map",
                                                 "slope_direction",
                                                 "slope_hist", "slope_bar",
                                                 "pvalue_map",
                                                 "pvalue_significance",
                                                 "pvalue_hist", "pvalue_bar"),
                                    smooth = TRUE, ...) {
  which <- match.arg(which)

  # workflow_trends(), unlike workflow_tst(), only ever computes ONE
  # FDR method per object (whichever fdr_method was requested at
  # construction) -- so, unlike .plot_tst()'s own "method" argument
  # (a real choice, since workflow_tst() can hold both BH and BKY at
  # once), there is only one valid method to use here, detected from
  # whichever rejection vector is actually present rather than
  # asked for as a separate argument.
  .fdr_method_used <- function(fdr_result) {
    present <- c(
      BH = !is.null(fdr_result$reject_BH),
      BKY = !is.null(fdr_result$reject_BKY),
      BY = !is.null(fdr_result$reject_BY)
    )
    methods <- names(present)[present]
    if (length(methods) != 1L) {
      stop("The workflow FDR result must contain exactly one of ",
           "reject_BH, reject_BKY, or reject_BY.")
    }
    methods
  }

  if (which == "direction") {
    if (is.null(x$fdr)) {
      stop("No FDR correction in this workflow_trends object -- rerun ",
           "workflow_trends() with fdr_method set, or use which = ",
           "'trend'.")
    }
    method_used_dir <- .fdr_method_used(x$fdr)
    slope_for_direction <- if (is.null(x$slope)) {
      NULL
    } else {
      suppressMessages(
        .smooth_slope_for_direction(x$slope, isTRUE(smooth))
      )
    }
    direction <- direction_map(x$trend, x$fdr, slope = slope_for_direction,
                                method = method_used_dir,
                                verbose = FALSE)
    fdr_direction_plot(direction)
  } else if (which == "significance") {
    if (is.null(x$fdr)) {
      stop("No FDR correction in this workflow_trends object -- rerun ",
           "workflow_trends() with fdr_method set, or use which = ",
           "'trend'.")
    }
    fdr_significance_maps(x$fdr)
  } else if (which == "slope") {
    if (is.null(x$slope)) {
      stop("No slope in this workflow_trends object -- rerun ",
           "workflow_trends() with slope_method set.")
    }
    if (is.null(x$fdr)) {
      stop("No FDR correction in this workflow_trends object -- rerun ",
           "workflow_trends() with fdr_method set.")
    }
    method_used <- .fdr_method_used(x$fdr)
    reject <- x$fdr[[paste0("reject_", method_used)]]
    slope_sig <- .mask_and_smooth_slope(x$slope, reject, isTRUE(smooth))
    subtitle <- sprintf("Slope, FDR-%s significant cells", method_used)
    range_lim <- .robust_diverging_range(slope_sig)
    terra::plot(slope_sig, col = .sptrends_diverging_palette(50),
                range = range_lim, fill_range = TRUE, main = subtitle)
  } else if (which %in% c("slope_map", "slope_direction", "slope_hist",
                           "slope_bar")) {
    if (is.null(x$slope)) {
      stop("No slope in this workflow_trends object -- rerun ",
           "workflow_trends() with slope_method set.")
    }
    switch(which,
      slope_map       = slope_map(x$slope),
      slope_direction = slope_direction_map(x$slope, ...),
      slope_hist      = slope_histogram(x$slope, ...),
      slope_bar       = slope_direction_barplot(x$slope, ...)
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
    trend_maps(x$trend, alpha = .reference_alpha(c(0.1, 0.05, 0.01)))
  }

  invisible(x)
}
