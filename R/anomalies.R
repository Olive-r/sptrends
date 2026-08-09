#' Remove the seasonal cycle from raster time series
#'
#' Monotonic trend tests assume that the input series does not contain a
#' regular periodic component -- the ranking of observations should
#' primarily reflect long-term change, not where in the year a value
#' falls. When a seasonal cycle is present (e.g. monthly means over many
#' years), that cycle itself dominates the ranking and can obscure or
#' distort the monotonic signal a trend test is looking for.
#'
#' This function removes that cycle by computing the mean value at each
#' position within it (e.g. the mean of all Januaries, all Februaries,
#' ...) -- the **climatology** -- and subtracting it from every layer at
#' that position, leaving an **anomaly series**: successive observations
#' that are directly comparable to one another, with the seasonal
#' pattern removed. This is appropriate input for [trend_test()] or
#' [prewhiten()], both of which assume no periodic component -- in a
#' typical workflow, this function is applied first, before
#' [prewhiten()], since prewhitening's own AR(1) model assumes the
#' series it receives has no remaining seasonal structure of its own.
#'
#' **Function type:** **Preprocessing function** -- prepares seasonal
#' raster time series before prewhitening or trend analysis. It is not
#' itself a trend test or slope estimator.
#'
#' @section Typical use:
#' ```
#' seasonal raster time series
#'     |
#' compute_anomalies()
#'     |
#' anomaly time series (`result$anomalies`)
#'     |
#' prewhiten(), trend_test(), or workflow_trends()
#' ```
#' Apply this step only when `x` contains a recurring seasonal cycle.
#' [workflow_tst()] and [workflow_trends()] start from prewhitening or
#' trend analysis, so pass `result$anomalies`, rather than the original
#' seasonal series, when deseasonalisation is required.
#'
#' @section Methodological details:
#' **How it works**
#'
#' A simple mean per cycle position, over the whole series -- **not** a
#' moving climatology, a LOESS-smoothed seasonal curve, or a spline fit.
#' This is deliberately the simplest well-defined choice, not a
#' limitation to work around: for the purpose of removing a fixed
#' seasonal pattern before a monotonic-trend test, a fixed per-position
#' mean is sufficient, and a more elaborate seasonal model would add
#' complexity without changing what this function is for.
#'
#' **Statistical assumptions**
#'
#' The seasonal cycle has a known, fixed length and corresponding cycle
#' positions are comparable across repetitions. The estimated fixed
#' climatology is assumed to represent the recurring component that
#' should be removed before later analysis.
#'
#' **Limitations**
#'
#' This function does not detrend the series in any other sense: the
#' climatology is computed directly over the raw values at each cycle
#' position, so a strong underlying trend is already partly absorbed
#' into the climatology mean itself (e.g. if summers have been getting
#' warmer throughout the record, the "mean summer" climatology reflects
#' an average across that warming, not any single year's summer). This
#' is expected, not a defect -- it is the trend test applied afterwards
#' that isolates the trend itself; this function's only job is removing
#' the periodic component. Likewise, a genuine regime shift (an abrupt,
#' one-time change in the seasonal pattern partway through the series,
#' as opposed to a gradual trend) would show up mixed into the
#' climatology rather than being detected as such -- this function has
#' no mechanism to distinguish a regime shift from ordinary seasonal
#' variability; that is a different kind of question from the one it
#' answers.
#'
#' **Quality assurance**
#'
#' Tests verify monthly and arbitrary-cycle climatologies against direct
#' calculations, centred and standardised anomalies, layer names,
#' retained geometry, missing values, zero-variance cycles, and invalid
#' inputs. Workflow tests confirm that anomaly outputs remain compatible
#' with later preprocessing and trend stages. See `?sptrends` for the
#' common release-check protocol.
#'
#' @param x A `terra::SpatRaster` with `nlyr(x)` a multiple of `cycle` (or
#'   not -- a partial final cycle is allowed, but its climatology mean will
#'   be based on fewer years for those positions).
#' @param cycle Integer. Length of the seasonal cycle in layers -- e.g.
#'   `12` for monthly data with an annual cycle (the default), `4` for
#'   quarterly/seasonal data, `365` for daily data with an annual cycle.
#' @param standardise Logical. If `FALSE` (default), returns raw anomalies
#'   (`x - climatology_mean`), in the original units -- preserves the
#'   variable's own physical units, and is the more common choice when
#'   that is what downstream reporting needs. If `TRUE`, also
#'   divides by the per-cycle-position standard deviation (z-scores:
#'   `(x - climatology_mean) / climatology_sd`), making the anomaly
#'   magnitude comparable across cycle positions that have very different
#'   natural variability (e.g. winter vs. summer temperature variance) --
#'   useful specifically when comparing or combining anomalies across
#'   cycle positions whose natural spread differs.
#'   Cycle positions with zero standard deviation (constant across all
#'   years) get `NA` anomalies rather than a division by zero.
#' @param verbose Logical. Print progress messages and elapsed time.
#'
#' @return Returns a list with:
#'   \item{anomalies}{A `SpatRaster`, same number of layers as `x` -- raw
#'     or standardised depending on `standardise`.}
#'   \item{climatology}{A `SpatRaster` with `cycle` layers (the mean field
#'     for each position in the cycle).}
#'   \item{climatology_sd}{Only if `standardise = TRUE`: a `SpatRaster`
#'     with `cycle` layers (the standard deviation field for each position
#'     in the cycle).}
#'   A plain list, not a classed `"sptrends"` object -- unlike
#'   [prewhiten()] or [trend_test()], this function's output is
#'   typically fed straight into the next preprocessing or inferential
#'   step rather than inspected on its own via `print()`/`summary()`/
#'   `plot()`.
#'
#' @seealso [prewhiten()] for temporal dependence treatment after anomaly
#'   construction; [trend_test()] and [workflow_trends()] for subsequent
#'   trend inference; [sim_trend_stack()] for controlled example data.
#'
#' @examples
#' # The bundled NDVI series is annual and therefore has no subannual
#' # climatological cycle to remove.
#' r <- read_ordered_stack(example_data("vhp_ndvi"))
#' terra::nlyr(r)
#' # Apply compute_anomalies() only to real observations with a genuine
#' # cycle, for example cycle = 12 for monthly data.
#'
#' @family preprocessing functions
#' @references
#' General references for the anomaly/standardisation concept (removing a
#' periodic mean, optionally scaling by its standard deviation) -- not a
#' single named method with one original paper, but a standard technique
#' in climatology and atmospheric science:
#' - Wilks, D.S. (2019) Statistical Methods in the Atmospheric Sciences
#'   (4th edn). Elsevier/Academic Press. No DOI available (book).
#' - Mather, P.M. (1999) Computer Processing of Remotely-Sensed Images.
#'   John Wiley and Sons.
#' @export
compute_anomalies <- function(x, cycle = 12, standardise = FALSE,
                               verbose = TRUE) {
  finish_timer <- .sptrends_elapsed_timer("compute_anomalies()", verbose)
  on.exit(finish_timer(), add = TRUE)
  if (!inherits(x, "SpatRaster")) stop("'x' must be a terra SpatRaster.")
  n <- terra::nlyr(x)
  if (cycle < 2) stop("'cycle' must be >= 2.")
  if (n < cycle) {
    stop(sprintf(
      paste0("Fewer layers (%d) than cycle length (%d) -- cannot estimate ",
             "a climatology."),
      n, cycle))
  }
  if (n %% cycle != 0 && verbose) {
    message(sprintf(
      paste0("Note: %d layers is not a multiple of cycle=%d -- the last ",
             "partial cycle uses fewer years for its climatology mean."),
      n, cycle))
  }

  position <- ((seq_len(n) - 1) %% cycle) + 1

  if (verbose) {
    message(sprintf("Computing climatology (%d cycle positions)...", cycle))
  }
  clim <- terra::tapp(x, index = position, fun = "mean", na.rm = TRUE)
  names(clim) <- paste0("clim_", seq_len(cycle))

  if (verbose) message("Subtracting climatology from each layer...")
  clim_per_layer <- clim[[position]]
  anomalies <- x - clim_per_layer
  names(anomalies) <- paste0(names(x), "_anom")

  out <- list(anomalies = anomalies, climatology = clim)

  if (isTRUE(standardise)) {
    if (verbose) {
      message("Standardising by per-cycle-position standard deviation...")
    }
    clim_sd <- terra::tapp(x, index = position, fun = "sd", na.rm = TRUE)
    names(clim_sd) <- paste0("sd_", seq_len(cycle))
    clim_sd_per_layer <- clim_sd[[position]]
    clim_sd_safe <- terra::ifel(clim_sd_per_layer == 0, NA, clim_sd_per_layer)
    anomalies <- anomalies / clim_sd_safe
    names(anomalies) <- paste0(names(x), "_zanom")
    out$anomalies <- anomalies
    out$climatology_sd <- clim_sd
  }

  if (verbose) message("Done.")
  out
}
