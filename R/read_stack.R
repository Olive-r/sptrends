#' @noRd
.order_for_stacking <- function(x) {
  order(x)
}

#' Read and chronologically order a folder of raster files
#'
#' `list.files()` sorts alphabetically, not numerically: with file names
#' that lack leading zeros (`"image1"`, ..., `"image10"`), `"image10"`
#' sorts before `"image2"` -- the same problem regardless of whether the
#' numbering is a year, a time step, or any other sequential index.
#' Using that order silently invalidates a trend analysis with no
#' visible error anywhere in the code. This function instead extracts
#' an explicit ordering number from each file name and sorts
#' numerically by it, then prints and (optionally) plots a verification
#' of the detected order.
#'
#' Works with raster formats readable by [terra::rast()], subject to the
#' GDAL drivers available in the user's `terra` installation. The default
#' pattern covers common GeoTIFF, NetCDF, native raster, ERDAS Imagine,
#' virtual raster and ASCII-grid extensions; other supported formats can
#' be selected through `pattern`. Each file must represent one time step.
#' For a single multi-temporal NetCDF file instead, use
#' [read_netcdf_stack()].
#'
#' **Why this matters more than it looks**: a trend analysis assumes
#' that successive raster layers represent the true chronological
#' sequence. If the temporal order is wrong, every subsequent
#' statistical result becomes invalid even though the analysis
#' completes without any error -- there is nothing in a Mann-Kendall
#' test's own output that could reveal a shuffled input series. Unlike
#' generic file readers, this function deliberately refuses to proceed
#' when the temporal order cannot be established unambiguously, rather
#' than silently reverting to alphabetical file names. This same
#' philosophy -- surface a silent risk rather than let it pass
#' unnoticed -- recurs throughout this package: multiple-testing
#' correction ([fdr_correction()]), the spatial-autocorrelation
#' diagnostic behind it ([spatial_autocorrelation()]), and the
#' monotonic-trend-only assumption checked informally in "Monotonic
#' trends only" sections elsewhere ([trend_test()],
#' [slope_estimator()]) are the same instinct applied to different
#' risks.
#'
#' This function does not use `list.files()`'s own ordering, does not
#' attempt to parse arbitrary or ambiguous date formats, and does not
#' guess when no candidate pattern extracts a unique number from every
#' file -- it stops instead, on the view that it is better to stop than
#' to silently analyse the wrong chronology.
#'
#' **Function type:** **Data import function** -- builds the ordered
#' `SpatRaster` expected by the analytical functions. It performs no
#' statistical inference.
#'
#' @section Typical use:
#' ```
#' folder with one raster file per time step
#'     |
#' read_ordered_stack()
#'     |
#' chronologically ordered raster time series
#'     |
#' compute_anomalies() if seasonal, then a trend workflow
#' ```
#'
#' @section Methodological details:
#' **Temporal ordering**
#'
#' Ordering is derived from explicit numeric labels in file names, never
#' from alphabetical order. Ambiguous, duplicated, or incomplete labels
#' stop the import rather than trigger an undocumented fallback.
#'
#' **Monthly or seasonal data**
#'
#' This function only orders layers chronologically -- it does not know
#' or care whether the data has a seasonal cycle (e.g. monthly values with
#' an annual signal). If it does, deseasonalise with [compute_anomalies()]
#' *before* passing the result to [trend_test()],
#' [slope_estimator()], or [workflow_tst()], all of which assume a
#' monotonic trend,
#' not a periodic one. Ordering and deseasonalisation solve different
#' problems.
#'
#' **Limitations**
#'
#' Files must share compatible raster geometry and each file must
#' represent one time step. The function deliberately does not attempt
#' to infer arbitrary calendar formats that are not captured by the
#' supplied or candidate regular expressions.
#'
#' **Quality assurance**
#'
#' Tests cover year and year-month filename parsing, chronological
#' ordering, duplicate and ambiguous labels, mixed geometries, variable
#' selection, preserved `terra` geometry/time metadata, and informative
#' failures for empty or invalid inputs. See `?sptrends` for the
#' package-wide release-check protocol.
#'
#' @param dir Character. Path to a folder containing one raster file per
#'   time step, all with the same extent, resolution, and CRS.
#' @param pattern Character. Regular expression used to list candidate
#'   files. The default matches common raster extensions: `.tif`, `.tiff`,
#'   `.nc`, `.grd`, `.img`, `.vrt`, and `.asc`. Supply another expression
#'   for any additional format supported by `terra::rast()`/GDAL. The
#'   extension identifies candidate files; actual readability and
#'   compatible geometry are still validated when the stack is opened.
#' @param order_regex Character or `NULL`. A regular expression with a
#'   single capture group, `"(...)"`, that extracts the ordering number
#'   from each file name -- e.g. `"year([0-9]+)"` extracts `25` from
#'   `"year25.tif"`. If `NULL` (default), several common patterns are
#'   tried automatically (year-like `"year25"`, MODIS-style `"_A2000065"`,
#'   a bare 4-digit year, or any run of digits) and the first one that
#'   extracts a **unique** number from every file is used. If none does,
#'   the function stops rather than silently falling back to alphabetical
#'   order.
#' @param candidate_regex Character vector of patterns tried automatically
#'   when `order_regex` is `NULL`. Most users should never need to modify
#'   this.
#' @param var Character or `NULL`. For NetCDF input only: variable name to
#'   read from each file if each file has more than one variable. Leave
#'   `NULL` for ordinary single-variable raster formats.
#' @param report Logical. If `TRUE` (default), draw a diagnostic plot of
#'   stack position vs. detected order number (a perfect diagonal line
#'   confirms correct ordering). Deviations from the diagonal immediately
#'   reveal files that would otherwise have been read in the wrong
#'   temporal order.
#' @param verbose Logical. Print progress messages and elapsed time.
#'
#' @return A `terra::SpatRaster` with one layer per input file, ordered
#'   chronologically, with layer names taken from the (de-duplicated) file
#'   names, and the detected order numbers stored as proper time metadata
#'   (`terra::time(result, "years")`) rather than discarded after the
#'   verification step above -- used, for instance, as the default `t`
#'   in [inspect_ts_cell()], and by any future function in this package
#'   requiring explicit time coordinates.
#'
#' @examples
#' \donttest{
#' # The bundled dataset contains one real NDVI GeoTIFF per year.
#' s <- read_ordered_stack(example_data("vhp_ndvi"))
#' terra::nlyr(s)
#' terra::time(s, "years")
#' terra::plot(s[[1]], main = "NDVI, first year")
#' }
#'
#' @family Data import functions
#' @export
read_ordered_stack <- function(
    dir, pattern = "\\.(tif|tiff|nc|grd|img|vrt|asc)$",
    order_regex = NULL,
    candidate_regex = c(
      "year([0-9]+)", "_A([0-9]{4})",
      "(19[0-9]{2}|20[0-9]{2})",
      "([0-9]{4})", "([0-9]+)"),
    var = NULL, report = TRUE, verbose = TRUE) {
  finish_timer <- .sptrends_elapsed_timer("read_ordered_stack()", verbose)
  on.exit(finish_timer(), add = TRUE)
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No files matching '", pattern, "' in: ", dir)

  file_names <- basename(files)

  try_regex <- function(regex) {
    numbers <- suppressWarnings(as.numeric(
      sub(paste0(".*", regex, ".*"), "\\1", file_names)
    ))
    if (length(numbers) != length(file_names) || anyNA(numbers)) return(NULL)
    if (anyDuplicated(numbers) > 0) return(NULL)
    numbers
  }

  if (!is.null(order_regex)) {
    order_numbers <- try_regex(order_regex)
    if (is.null(order_numbers)) {
      stop(sprintf(paste(
        "TEMPORAL ORDER NOT VERIFIABLE -- stopping on purpose.",
        "The pattern '%s' you supplied did not extract a UNIQUE number from",
        "every file in '%s'. Check your file names or try another pattern.",
        sep = "\n"
      ), order_regex, dir))
    }
  } else {
    order_numbers <- NULL
    pattern_used <- NULL
    for (candidate in candidate_regex) {
      result <- try_regex(candidate)
      if (!is.null(result)) {
        order_numbers <- result
        pattern_used <- candidate
        break
      }
    }
    if (is.null(order_numbers)) {
      stop(sprintf(paste(
        "TEMPORAL ORDER NOT VERIFIABLE -- stopping on purpose.",
        "These patterns were tried automatically without success: %s.",
        "None extracted a UNIQUE number from every file name in '%s'.",
        "Alphabetical order is NEVER used as a fallback.",
        "Supply 'order_regex' manually with a pattern that captures the",
        "real ordering number of your files (one capture group, \"(...)\").",
        sep = "\n"
      ), paste(candidate_regex, collapse = ", "), dir))
    }
    if (verbose) {
      message(sprintf("Temporal order auto-detected with pattern '%s'.",
                       pattern_used))
    }
  }

  # A thin, purely internal wrapper around order() -- exists so this
  # exact permutation step can be safely mocked in a test (mocking this
  # package's own internal function) without touching base::order()
  # itself (which every other part of R also depends on, and mocking
  # broadly would risk destabilising the test run in ways unrelated to
  # this package's own logic).
  ord <- .order_for_stacking(order_numbers)
  files <- files[ord]
  ordered_numbers <- order_numbers[ord]

  valid_sequence <- anyDuplicated(ordered_numbers) == 0 &&
    !any(diff(ordered_numbers) <= 0)

  if (verbose) {
    order_table <- data.frame(
      stack_position = seq_along(files),
      detected_number = ordered_numbers,
      file = basename(files)
    )
    message("Temporal order verification (mandatory, cannot be skipped):")
    message(paste(utils::capture.output(print(order_table, row.names = FALSE)),
                  collapse = "\n"))
  }

  # Defensive: try_regex() above already rejects any candidate with
  # duplicate or non-finite numbers (returns NULL, caught earlier as
  # "NOT VERIFIABLE"), so order_numbers is already unique by the time
  # it reaches .order_for_stacking() -- under a correct permutation,
  # sorting unique numbers always gives strictly positive diffs. This
  # branch exists as an explicit safety net for the case where
  # .order_for_stacking() itself misbehaves (e.g. a future refactor
  # introduces a bug in it) -- silent data corruption (files stacked in
  # a subtly wrong temporal order) would be far worse than trusting an
  # unchecked permutation. Exercised directly in the test suite by
  # mocking .order_for_stacking() to return a deliberately wrong
  # permutation, not left untested.
  if (!valid_sequence) {
    stop(paste(
      "INVALID TEMPORAL ORDER -- stopping on purpose.",
      "The sequence of detected numbers has repeats or is not strictly",
      "increasing after sorting. This points to a real problem with the",
      "file names -- check it before continuing.",
      sep = "\n"
    ))
  } else if (verbose && any(diff(ordered_numbers) != 1)) {
    message("Note: there are gaps in the numbering ",
            "(not all years/steps are consecutive).")
    message("Check whether this is intentional (missing years) or a ",
            "misnamed file.")
  }

  if (isTRUE(report)) {
    graphics::plot(seq_along(ordered_numbers), ordered_numbers,
         type = "b", pch = 19, col = "steelblue",
         xlab = "Position in the stack",
         ylab = "Detected order number",
         main = "Temporal order check\n(should be a perfect diagonal)")
    graphics::abline(a = ordered_numbers[1] - 1, b = 1, col = "red", lty = 2)
    graphics::legend("topleft",
           legend = c("Detected order", "Expected diagonal"),
           col = c("steelblue", "red"), lty = c(1, 2), pch = c(19, NA),
           bty = "n")
  }

  if (!is.null(var)) {
    layers <- lapply(files, terra::rast, subds = var)
    s <- do.call(c, layers)
  } else {
    s <- terra::rast(files)
  }
  if (!all(terra::inMemory(s))) s <- s + 0

  names_ <- tools::file_path_sans_ext(basename(files))
  # Also hard to trigger honestly through normal usage: for two files'
  # stripped names to collide, their un-stripped names would need to
  # collide everywhere except inside the final extension itself (the
  # only part tools::file_path_sans_ext() removes) -- but that same
  # near-identical stem is exactly what would make any reasonable
  # order_regex/candidate_regex extract the same number from both,
  # which is already rejected earlier (try_regex()'s own uniqueness
  # check) before layer names are ever built. Kept as a genuine safety
  # net regardless -- a silently-overwritten layer name is a worse
  # failure mode than an untested branch, and this doesn't depend on
  # anything about how the numbers were extracted remaining true.
  if (anyDuplicated(names_) > 0) {
    if (verbose) {
      message("Note: duplicate layer names -- made unique automatically.")
    }
    names_ <- make.unique(names_, sep = "_")
  }
  names(s) <- names_

  # The detected temporal order numbers (already verified above, table
  # and diagonal plot included) are otherwise computed and then thrown
  # away -- stored here as proper time metadata so downstream functions
  # (inspect_ts_cell() in particular) can use the real detected years
  # instead of a generic 1:n time-step index by default. tstep="years"
  # is the closest fit for what candidate_regex/order_regex extract in
  # practice (a bare year, or a "year N" sequential index) -- if that
  # index isn't literally a calendar year for your data, terra will
  # still display the values correctly, just not as calendar dates.
  terra::time(s, tstep = "years") <- ordered_numbers

  if (verbose) {
    message(sprintf("Stack built: %d layers, %d x %d cells.",
                     terra::nlyr(s), terra::nrow(s), terra::ncol(s)))
  }
  s
}

#' Read and chronologically order a single multi-temporal NetCDF file
#'
#' Wraps `terra::rast()` for the common case of one NetCDF file holding an
#' entire time series (e.g. reanalysis or climate model output), verifying
#' that the layers come out in chronological order using the file's own
#' time dimension (via `terra::time()`) rather than assuming the on-disk
#' layer order is already correct.
#'
#' **Why this matters more than it looks**: a trend analysis assumes
#' that successive raster layers represent the true chronological
#' sequence. A NetCDF file's own internal layer order and its time
#' dimension are two separate pieces of metadata, written independently
#' -- nothing in the file format itself guarantees they agree, and if
#' they do not, every subsequent statistical result becomes invalid
#' even though the analysis completes without any error. This function
#' reorders by the time dimension explicitly rather than trusting the
#' on-disk order, the same "surface a silent risk rather than let it
#' pass unnoticed" instinct behind [read_ordered_stack()]'s own,
#' stricter, filename-based check (see its own documentation for the
#' fuller reasoning, which applies here too).
#'
#' **Function type:** **Data import function** -- builds an ordered
#' `SpatRaster` from one multi-temporal NetCDF file. It performs no
#' statistical inference.
#'
#' @section Typical use:
#' ```
#' NetCDF file with a time dimension
#'     |
#' read_netcdf_stack()
#'     |
#' chronologically ordered raster time series
#'     |
#' compute_anomalies() if seasonal, then a trend workflow
#' ```
#'
#' @section Methodological details:
#' **Temporal ordering**
#'
#' Layers are reordered explicitly from the NetCDF time coordinate; the
#' function does not assume that physical layer order and temporal
#' metadata already agree.
#'
#' **Monthly or seasonal data**
#'
#' This function only orders layers chronologically -- it does not remove
#' or otherwise account for a seasonal cycle (e.g. monthly reanalysis
#' values with an annual signal, very common in NetCDF climate data). If
#' the detected time step looks sub-annual (based on `terra::timeInfo()`),
#' a `warning()` is issued reminding you to deseasonalise with
#' [compute_anomalies()] *before* passing the result to
#' [trend_test()], [slope_estimator()], or [workflow_tst()], all of
#' which assume a monotonic trend, not a periodic one. Ordering and
#' deseasonalisation solve different problems.
#'
#' **Limitations**
#'
#' A usable, unambiguous time coordinate is required. With multiple
#' variables, `var` must identify the intended field; this function does
#' not guess among scientifically different variables.
#'
#' **Quality assurance**
#'
#' Tests verify time-coordinate extraction and ordering, variable
#' selection, layer naming, preserved `terra` geometry/time metadata,
#' and failures for absent or ambiguous NetCDF variables. See
#' `?sptrends` for the package-wide release-check protocol.
#'
#' @param path Character. Path to the `.nc` file.
#' @param var Character or `NULL`. Variable name to read, if the file has
#'   more than one. If `NULL` and the file has a single variable, that one
#'   is used; if it has several, the function stops and lists them.
#' @param report Logical. If `TRUE` (default), draw a diagnostic plot of
#'   stack position vs. time value (a perfect diagonal line confirms
#'   correct ordering; deviations reveal layers that would otherwise
#'   have been read in the wrong temporal order).
#' @param verbose Logical. Print progress messages and elapsed time.
#'
#' @return A `terra::SpatRaster`, ordered chronologically, with layer names
#'   taken from the time values, and the time dimension itself preserved
#'   as proper time metadata (readable via `terra::time()`) rather than
#'   discarded after the verification step above.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("ncdf4", quietly = TRUE)) {
#'   # Convert the bundled environmental series to one temporary NetCDF.
#'   r <- read_ordered_stack(example_data("vhp_ndvi"))
#'   path <- tempfile(fileext = ".nc")
#'   terra::writeCDF(r, path, varname = "ndvi", overwrite = TRUE)
#'   s <- read_netcdf_stack(path)
#'   terra::nlyr(s)
#'   terra::time(s)
#'   unlink(path)
#' }
#' }
#'
#' @family Data import functions
#' @export
read_netcdf_stack <- function(path, var = NULL, report = TRUE, verbose = TRUE) {
  finish_timer <- .sptrends_elapsed_timer("read_netcdf_stack()", verbose)
  on.exit(finish_timer(), add = TRUE)
  if (!file.exists(path)) stop("File not found: ", path)

  if (is.null(var)) {
    var_names <- tryCatch(names(terra::sds(path)), error = function(e) NULL)
    if (!is.null(var_names) && length(var_names) > 1) {
      stop(sprintf(
        paste0("This NetCDF file has %d variables and 'var' is NULL -- ",
               "specify which one to read. Available: %s"),
        length(var_names), paste(var_names, collapse = ", ")
      ))
    }
    r <- terra::rast(path)
  } else {
    r <- terra::rast(path, subds = var)
  }

  time_vals <- terra::time(r)
  if (all(is.na(time_vals))) {
    stop(paste(
      "TEMPORAL ORDER NOT VERIFIABLE -- stopping on purpose.",
      paste0("This NetCDF file has no usable time dimension ",
             "(terra::time() is all NA)."),
      "Fix the file's time metadata, or use read_ordered_stack() on a folder",
      "of single-layer rasters with parseable file names instead.",
      sep = "\n"
    ))
  }

  ord <- order(time_vals)
  r <- r[[ord]]
  ordered_time <- time_vals[ord]

  valid_sequence <- anyDuplicated(ordered_time) == 0 &&
    !any(diff(as.numeric(ordered_time)) <= 0)

  if (verbose) {
    order_table <- data.frame(stack_position = seq_along(ordered_time),
                               time = as.character(ordered_time))
    message("Temporal order verification (mandatory, cannot be skipped):")
    message(paste(utils::capture.output(print(order_table, row.names = FALSE)),
                  collapse = "\n"))
  }

  if (!valid_sequence) {
    stop(paste(
      "INVALID TEMPORAL ORDER -- stopping on purpose.",
      "The time values in this file have repeats or are not strictly",
      "increasing after sorting -- check the file's time dimension.",
      sep = "\n"
    ))
  }

  time_step_unit <- tryCatch(terra::timeInfo(r)$step,
                              error = function(e) NA_character_)
  time_gaps <- diff(ordered_time)
  gap_days <- if (inherits(time_gaps, "difftime")) {
    suppressWarnings(
      stats::median(as.numeric(time_gaps, units = "days"), na.rm = TRUE)
    )
  } else {
    suppressWarnings(stats::median(as.numeric(time_gaps), na.rm = TRUE))
  }
  # Only trust the measured gap-in-days when the underlying time unit is
  # itself days/months (for "years" or unknown units, a numeric diff of 1
  # means "1 year", not "1 day" -- checking gap_days there would misfire).
  sub_annual <- isTRUE(!is.na(time_step_unit) && time_step_unit == "months") ||
    isTRUE(!is.na(time_step_unit) && time_step_unit %in% c("days", "seconds") &&
             !is.na(gap_days) && gap_days > 0 && gap_days < 300)
  if (isTRUE(sub_annual)) {
    gap_text <- if (!is.na(gap_days)) {
      sprintf("~%.0f days apart", gap_days)
    } else {
      "sub-annual"
    }
    warning(sprintf(paste(
      "Time steps are %s -- this looks like sub-annual",
      "(e.g. monthly) data, which often has a seasonal cycle. If it does,",
      "deseasonalise with compute_anomalies() before passing this stack to",
      "trend_test(), slope_estimator(), or workflow_tst(): all three",
      "assume a monotonic trend, not a periodic one."
    ), gap_text), call. = FALSE)
  }

  if (isTRUE(report)) {
    graphics::plot(seq_along(ordered_time), as.numeric(ordered_time),
         type = "b", pch = 19, col = "steelblue",
         xlab = "Position in the stack",
         ylab = "Time value (numeric)",
         main = "Temporal order check (NetCDF)\n(should be a perfect diagonal)")
  }

  names(r) <- as.character(ordered_time)
  if (verbose) {
    message(sprintf("Stack built from NetCDF: %d layers, %d x %d cells.",
                     terra::nlyr(r), terra::nrow(r), terra::ncol(r)))
  }
  r
}
