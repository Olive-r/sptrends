#' @noRd
.print_slope <- function(x, ...) {
  label <- switch(x$method,
    TS  = "Theil-Sen",
    OLS = "Ordinary least squares",
    RM  = "Repeated median (Siegel)",
    x$method
  )
  cat(sprintf("<%s slope result%s>\n", label,
              if (isTRUE(x$smoothed)) " (queen-3x3 smoothed)" else ""))
  vals <- terra::values(x$slope, mat = FALSE)
  ok <- !is.na(vals)
  cat(sprintf("Valid cells: %d, range [%.4g, %.4g]\n",
              sum(ok), min(vals[ok]), max(vals[ok])))
  invisible(x)
}

#' @noRd
.summary_slope <- function(object, ...) {
  slope_summary(object$slope, ...)
}

#' @noRd
.plot_slope <- function(x, which = c("map", "direction", "histogram", "bar"),
                         ...) {
  which <- match.arg(which)
  switch(which,
    map       = slope_map(x$slope, ...),
    direction = slope_direction_map(x$slope, ...),
    histogram = slope_histogram(x$slope, ...),
    bar       = slope_direction_barplot(x$slope, ...)
  )
  invisible(x)
}
