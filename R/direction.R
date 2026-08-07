#' Binarised trend map after multiple-testing correction
#'
#' Combines the sign of the trend statistic (`Sm`/`S`) with a chosen
#' FDR-corrected (or raw) rejection vector to classify each cell as an
#' increase, a decrease, or a non-significant result -- a binarised map
#' for reporting *after* multiple-comparison correction, as
#' opposed to [trend_maps()], which uses the uncorrected p-value.
#'
#' **Function type:** **Support function** -- computes something real
#' (a genuinely new
#' raster, combining trend direction with significance), but is not one
#' of the core building blocks of TST or RTA itself; it is a
#' post-processing step that consumes the output of two of them
#' ([trend_test()] and [fdr_correction()]) together. Not exported --
#' the binarised trend direction is the same underlying
#' computation as the uncorrected direction already reachable via
#' `plot(x, which = "trend")`, just with a significance filter applied
#' on top; reachable directly via `plot(x, which = "direction")` for
#' `workflow_tst()`/`workflow_rta()` results, or with `:::` for
#' programmatic use.
#'
#' @param trend The `$stats` field of [trend_test()]'s output.
#' @param fdr_result Output of [fdr_correction()], run on `trend$stats$p`.
#' @param slope Optional single-layer `SpatRaster` (e.g.
#'   [slope_estimator()]'s own `$slope`) whose sign determines direction
#'   instead of the trend test's own statistic. `NULL` (default): use
#'   `trend`'s own `Sm`/`S`/`beta`, as before. Direction from the slope and
#'   direction from the trend statistic usually agree but are not
#'   guaranteed to: a cell with a near-flat slope of its own can still
#'   inherit its neighbours' sign in `Sm` under CMK's neighbourhood
#'   averaging. The significance mask (which cells are shown as
#'   increasing/decreasing at all, from `fdr_result`) is identical either
#'   way -- only the source of the *sign* for already-significant cells
#'   changes.
#' @param method `"BH"` (default), `"BKY"`, `"BY"`, or `"raw"` (uncorrected,
#'   included for comparison) -- which rejection vector in `fdr_result` to
#'   use as the significance mask.
#' @param verbose Logical. Print a one-line count summary.
#'
#' @return A single-layer `terra::SpatRaster`, values `-1` (significant
#'   decrease), `0` (not significant), `1` (significant increase), named
#'   `"binarised_trend_map"`.
#'
#' @examples
#' \donttest{
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' trend <- trend_test(r, report = FALSE, verbose = FALSE)
#' fdr_result <- fdr_correction(trend$stats$p, report = FALSE, verbose = FALSE)
#'
#' # Combines "is it significant" (from fdr_result) with "which way is
#' # it going" (from trend) into a single binarised trend map.
#' direction <- sptrends:::direction_map(trend$stats, fdr_result, method = "BH")
#' }
#'
#' @references
#' This combination of FDR-corrected significance with trend direction is
#' this package's own contribution, not from an external method paper;
#' cited here as the source of the overall TST workflow it belongs to:
#' - Gutiérrez-Hernández, O. and García, L.V. (2025) Uncovering true
#'   significant trends in global greening. Remote Sensing Applications:
#'   Society and Environment, 37, 101377. \doi{10.1016/j.rsase.2024.101377}
#'
#' Underlying theoretical justification for the FDR-BH assumption this
#' significance is based on:
#' - Benjamini, Y., & Yekutieli, D. (2001) The control of the false
#'   discovery rate in multiple testing under dependency. Annals of
#'   Statistics, 29(4), 1165-1188. \doi{10.1214/aos/1013699998}
#' @keywords internal
direction_map <- function(trend, fdr_result, slope = NULL,
                               method = c("BH", "BKY", "BY", "raw"),
                               verbose = TRUE) {
  method <- match.arg(method)

  if (!is.null(slope)) {
    if (!inherits(slope, "SpatRaster")) {
      stop("'slope' must be a terra SpatRaster (e.g. slope_estimator()'s ",
           "own $slope), or NULL.")
    }
    Sv <- terra::values(slope, mat = FALSE)
    source_used <- "slope"
  } else {
    s_name <- if ("Sm" %in% names(trend)) {
      "Sm"
    } else if ("S" %in% names(trend)) {
      "S"
    } else {
      "beta"
    }
    Sv <- terra::values(trend[[s_name]], mat = FALSE)
    source_used <- s_name
  }

  reject <- switch(method,
    BH  = fdr_result$reject_BH,
    BKY = fdr_result$reject_BKY,
    BY  = fdr_result$reject_BY,
    raw = fdr_result$reject_raw
  )
  if (is.null(reject)) {
    stop(sprintf(
      paste0("fdr_result has no rejection vector for method='%s' -- ",
             "rerun fdr_correction() with that method included."),
      method))
  }

  na_mask <- is.na(reject) | is.na(Sv)
  cond_inc <- reject & Sv > 0
  cond_dec <- reject & Sv < 0
  cond_inc[is.na(cond_inc)] <- FALSE
  cond_dec[is.na(cond_dec)] <- FALSE

  cls_vals <- rep(0, length(Sv))
  cls_vals[cond_inc] <- 1
  cls_vals[cond_dec] <- -1
  cls_vals[na_mask] <- NA

  r1 <- trend[[1]]
  out <- terra::setValues(r1, cls_vals)
  names(out) <- "binarised_trend_map"

  if (verbose) {
    tab <- table(factor(cls_vals, levels = c(-1, 0, 1),
                         labels = c("decrease", "no_change", "increase")))
    message(sprintf(
      paste0("Binarised trend map (%s) [direction from %s] -- increase: %d | ",
             "decrease: %d | not significant: %d"),
      method, source_used, tab[["increase"]], tab[["decrease"]],
      tab[["no_change"]]))
  }

  out
}
