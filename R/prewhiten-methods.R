#' @noRd
.print_prewhiten <- function(x, ...) {
  if (identical(x$method, "TFPW_Y")) {
    cat("<Yue, Pilon & Cavadias (2002) trend-free prewhitening result>\n")
    rho <- terra::values(x$diagnostics$Rho, mat = FALSE)
    n_total <- sum(!is.na(rho))
    cat(sprintf("Prewhitened: all %d valid cells (no DW gate)\n", n_total))
  } else if (x$method %in% c("TFPW_WS", "TFPW_Z")) {
    label <- if (identical(x$method, "TFPW_Z")) {
      "<Zhang, Vincent, Hogg & Niitsoo (2000) prewhitening result>\n"
    } else {
      "<Wang & Swail (2001) prewhitening result>\n"
    }
    cat(label)
    mod <- terra::values(x$diagnostics$Modified, mat = FALSE)
    n_total <- sum(!is.na(mod))
    n_mod <- sum(mod == 1, na.rm = TRUE)
    cat(sprintf("Prewhitened: %d of %d valid cells (%.1f%%)\n",
                n_mod, n_total, 100 * n_mod / n_total))
  } else {
    # VCTFPW is the only remaining value match.arg() allows in
    # prewhiten() itself -- no further fallback branch is reachable
    # here, since that match.arg() call already restricts method to
    # exactly TFPW_Y/TFPW_WS/TFPW_Z/VCTFPW before this function is
    # ever invoked with a genuine prewhiten() result.
    cat(paste0(
      "<Wang, Chen, Becker & Liu (2015) variance-corrected ",
      "prewhitening result>\n"))
    mod <- terra::values(x$diagnostics$Modified, mat = FALSE)
    n_total <- sum(!is.na(mod))
    n_mod <- sum(mod == 1, na.rm = TRUE)
    cat(sprintf(
      "Prewhitened: %d of %d valid cells (95%% lag-1 ACF gate; %.1f%%)\n",
      n_mod, n_total, 100 * n_mod / n_total))
  }
  cat("Use summary() for diagnostic detail, or inspect $diagnostics ",
      "directly.\n", sep = "")
  invisible(x)
}

#' @noRd
.summary_prewhiten <- function(object, ...) {
  if (identical(object$method, "TFPW_Y")) {
    .TFPW_Y_summary(object$diagnostics, ...)
  } else if (object$method %in% c("TFPW_WS", "TFPW_Z")) {
    prewhiten_summary(object$diagnostics, ...)
  } else {
    # VCTFPW: no dedicated summary function -- its own diagnostics
    # (Rho + Beta_corrected + Modified) do not match the shape
    # prewhiten_summary() expects for the four-layer
    # (DW_initial/Rho/Modified/Clamped)
    # methods. Message this clearly, matching the same honesty
    # prewhiten()'s own report = TRUE already gives at construction
    # time, rather than erroring here when the object is inspected
    # afterwards instead.
    message("No dedicated summary() for method = \"", object$method,
            "\" yet -- its own diagnostics do not match the shape this ",
            "function expects. Try summary(terra::values(object$",
            "diagnostics)) directly, or terra::plot(object$diagnostics).")
    invisible(NULL)
  }
}

#' @noRd
.plot_prewhiten <- function(x, which = c("maps", "histograms"), ...) {
  which <- match.arg(which)
  if (identical(x$method, "TFPW_Y")) {
    if (which == "maps") {
      .TFPW_Y_maps(x$diagnostics)
    } else {
      .TFPW_Y_histograms(x$diagnostics)
    }
    invisible(x)
  } else if (x$method %in% c("TFPW_WS", "TFPW_Z")) {
    if (which == "maps") {
      prewhiten_maps(x$diagnostics)
    } else {
      prewhiten_histograms(x$diagnostics)
    }
    invisible(x)
  } else {
    # VCTFPW: see .summary_prewhiten()'s own comment above for why it
    # has no dedicated reporting function yet.
    message("No dedicated plot() for method = \"", x$method, "\" yet -- ",
            "its own diagnostics do not match the shape this function ",
            "expects. Try terra::plot(x$diagnostics) directly.")
    invisible(x)
  }
}
