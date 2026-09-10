#' Path to sptrends' bundled example dataset
#'
#' sptrends is built around one idea: spatiotemporal trend analysis of
#' gridded (raster) data -- a stack of layers over the same spatial grid,
#' one layer per time step. `example_data()` gives you a real dataset in
#' exactly that shape, bundled with the package, so every function's
#' `@examples`, every vignette, and your own first attempts at the
#' package have something real to run on without downloading anything.
#'
#' The bundled dataset is annual mean NDVI (Normalized Difference
#' Vegetation Index) -- derived from the NOAA STAR Blended Vegetation
#' Health Product, 1982-2023, global land 100 km Eckert IV equal-area
#' grid, ~3.5 MB -- a folder of one GeoTIFF per year, which is exactly
#' the layout [read_ordered_stack()] expects. So `example_data()`
#' doubles as a realistic example of *that* function's intended use (a
#' folder of yearly rasters), not just a shortcut to a pre-loaded R
#' object. NDVI trend analysis (vegetation "greening" and "browning")
#' is also the motivating application of the True Significant Trends
#' workflow itself -- see the primary reference in `?workflow_tst`.
#'
#' **Function type:** **Support function** -- locates package example
#' files; it performs no statistical analysis.
#'
#' @section Typical use:
#' ```
#' example_data("vhp_ndvi")
#'     |
#' read_ordered_stack()
#'     |
#' bundled annual raster time series
#'     |
#' workflow_tst(), workflow_rta(), or workflow_trends()
#' ```
#' Call `example_data()` without an argument first to list all bundled
#' paths.
#'
#' @section Methodological details:
#' **Data source and licence**
#'
#' Derived from the Blended Vegetation Health Product (Blended-VHP),
#' provided by NOAA's Center for Satellite Applications and Research
#' (STAR) -- see `example_data("vhp_ndvi/Readme.txt")` for full
#' provenance, processing steps, and the required acknowledgement. This
#' derived subset is for demonstration purposes only; it is not a
#' substitute for the original product, which has global land coverage
#' at the original 4 km resolution, weekly observations, and additional
#' variables (brightness temperature, Vegetation Condition Index,
#' Vegetation Health Index).
#'
#' **Limitations**
#'
#' The bundled raster is a small, coarsened demonstration dataset and
#' must not be treated as a replacement for the original NOAA product.
#'
#' **Quality assurance**
#'
#' Tests verify file listing, path resolution, informative failure for
#' absent paths, and end-to-end compatibility of the bundled dataset
#' with [read_ordered_stack()]. See `?sptrends` for the package-wide
#' release-check protocol.
#'
#' @param path Character or `NULL`. If `NULL` (default), lists every
#'   bundled example file (as paths relative to the package's `extdata`
#'   directory) -- use this first, to see what is available, before
#'   deciding what to ask for. If a specific relative path from that
#'   listing (e.g. `"vhp_ndvi"`, the whole dataset folder, or one
#'   particular `.tif` file inside it), the full absolute path to it on
#'   this machine.
#'
#' @return If `path = NULL`, a character vector of relative paths.
#'   Otherwise, a single absolute file path (as from `system.file()`).
#'   Errors if the requested `path` does not exist.
#'
#' @examples
#' # example_data() with no arguments lists every file bundled with
#' # the package, so you can see what is available before using any of
#' # it -- you do not need to know the file names in advance.
#' example_data()
#'
#' \donttest{
#' # Passing "vhp_ndvi" (the folder name from the listing above) gives
#' # the full path to that folder on your machine. read_ordered_stack()
#' # reads every GeoTIFF inside it and stacks them into one multi-layer
#' # raster, one layer per year, already sorted chronologically -- this
#' # is the "gridded time series" shape sptrends is built around.
#' # report = FALSE just turns off the automatic year-order check plot,
#' # to keep this example quiet; leave it on (the default) when
#' # exploring interactively, as a sanity check on the file order.
#' # Wrapped in \donttest{} rather than run unconditionally: reading
#' # all 42 real GeoTIFFs takes several seconds, unlike every other
#' # example in this package, which uses small synthetic rasters --
#' # still runs when a user tries this example directly, just not
#' # timed as part of R CMD check's own example suite.
#' r <- read_ordered_stack(example_data("vhp_ndvi"), report = FALSE)
#'
#' # nlyr() ("number of layers") confirms how many years came through --
#' # 42, one per year from 1982 to 2023.
#' terra::nlyr(r)
#'
#' # r[[1]] is the first layer (year 1982) on its own. terra's default
#' # colour scheme is not designed for a vegetation index, so we ask for
#' # a better one: grDevices::hcl.colors(50, "Greens 3") gives 50 shades
#' # from pale to deep green, matching how NDVI maps are normally read
#' # (higher NDVI, more/denser vegetation).
#' terra::plot(r[[1]], col = rev(grDevices::hcl.colors(50, "Greens 3")))
#' }
#'
#' @references
#' NOAA Center for Satellite Applications and Research (STAR). Blended
#' Vegetation Health Product (Blended-VHP).
#' \url{https://www.star.nesdis.noaa.gov/smcd/emb/vci/VH/vh_ftp.php}
#' @family example data functions
#' @export
example_data <- function(path = NULL) {
  if (is.null(path)) {
    return(list.files(system.file("extdata", package = "sptrends"),
                       recursive = TRUE))
  }
  system.file("extdata", path, package = "sptrends", mustWork = TRUE)
}
