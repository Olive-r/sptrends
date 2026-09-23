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
#' **Explicit declaration (`files`, `time`, `cycle_type`)**
#'
#' Automatic detection above only extracts one ordering number per file
#' name, which is reliable for genuinely annual series but not for finer
#' cadences: the same numeric shape in a file name can mean genuinely
#' different things across real datasets -- `"19820102"` is 2 January in
#' some daily products, but the second half of January under PKU-GIMMS
#' NDVI's own semimonthly convention. `files` (an explicit, already
#' correctly ordered vector) sidesteps this by never interpreting file
#' names at all; combined with `time` (fully explicit dates) or
#' `cycle_type` (a named calendar convention) it is recommended whenever
#' the series is not simply annual.
#'
#' **Supported cycle types**
#'
#' Eight unambiguous calendar conventions, all built with genuine calendar
#' arithmetic (leap years included) rather than naive interval division.
#' Fixed compositing intervals that do not follow one of these exact
#' conventions (e.g. a genuinely continuous 8-day or 16-day interval) are
#' out of scope for `cycle_type` -- supply `time` explicitly instead.
#'
#' \describe{
#'   \item{`"annual"`}{One date per year, 1 January. 1 sample/year.
#'     E.g. the bundled `example_data("vhp_ndvi")` dataset itself.}
#'   \item{`"monthly"`}{One period per calendar month, beginning on the 1st.
#'     12 samples/year. E.g. CRU TS, TerraClimate, ERA5 monthly means.}
#'   \item{`"16-day"`}{23 fixed periods per year, resetting every
#'     1 January (never continuing across a year boundary), starting
#'     on year-day 1, 17, 33, ..., 353 -- the last period of the year
#'     is shorter (13 or 14 days) to fit within the calendar year.
#'     Matches MODIS' own 16-day compositing convention (e.g.
#'     MOD13Q1), verified directly against real product file names.}
#'   \item{`"semimonthly"`}{Two periods per calendar month, beginning on the 1st and
#'     the 16th. 24 samples/year. Matches PKU-GIMMS NDVI's own
#'     "half-month" convention -- not the same cadence as a continuous
#'     14-day interval.}
#'   \item{`"10-day"`}{Three calendar periods per month, starting on
#'     the 1st, 11th and 21st (the last one running to the end of the
#'     month, so its own length varies: 8 to 11 days). 36 samples/year.
#'     The standard "dekad" convention, e.g. SPOT-VEGETATION and
#'     several FEWS NET agricultural products -- not a continuous
#'     10-day interval, which would not align with month boundaries
#'     and would give a different total (37, not 36).}
#'   \item{`"8-day"`}{46 fixed periods per year, resetting every
#'     1 January, starting on year-day 1, 9, 17, ..., 361 -- the last
#'     period of the year is shorter, capturing the remainder.
#'     Matches MODIS' own 8-day compositing convention (e.g.
#'     MCD15A2H, the LAI/FPAR product).}
#'   \item{`"weekly"`}{52 fixed, complete 7-day periods per year
#'     (Earth Trends Modeler's own number), starting on year-day 1, 8,
#'     ..., 358 -- unlike `"8-day"`/`"16-day"`, the final one or two days
#'     of the year belong to no period and are skipped straight into
#'     next year's first period, rather than forming a shorter final
#'     period.}
#'   \item{`"daily"`}{One date per real calendar day, including 29
#'     February in leap years. 365 or 366 samples/year -- the only one
#'     of these eight where the yearly total itself changes with leap
#'     years (the other seven keep a fixed count per year; only the
#'     exact date of late-period boundaries shifts). E.g. ERA5,
#'     CHIRPS daily precipitation.}
#' }
#'
#' `"weekly"`'s own remainder-discarding convention (unlike the other
#' fixed-interval types here) is not this function's own inconsistency
#' -- it follows Earth Trends Modeler's own table exactly for each
#' named type individually, and that table is not internally
#' consistent on this specific point across its own listed
#' conventions.
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
#' represent one time step. Automatic detection deliberately does not
#' attempt to infer arbitrary calendar formats that are not captured by
#' the supplied or candidate regular expressions -- use `files` with
#' `time` or `cycle_type` for those cases instead of expecting automatic
#' detection to guess correctly.
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
#' @param files Character vector of file paths, in the exact order you
#'   want them read as layers -- no sorting of any kind is ever applied.
#'   Use instead of `dir` whenever automatic detection from file names is
#'   not reliable enough, which in practice means: whenever the series is
#'   not simply annual. The same numeric shape in a file name can mean
#'   genuinely different things across real datasets (e.g. `"19820102"`
#'   is 2 January in some daily products, but the second half of January
#'   in PKU-GIMMS' own semimonthly convention) -- automatic detection
#'   cannot resolve that safely, only declaring it explicitly can.
#'   Exactly one of `dir` or `files` must be supplied.
#' @param time Optional. A `Date`, `POSIXct`, or numeric vector, one value
#'   per layer, in the same order as `files`. Must be strictly increasing
#'   with no missing or duplicated values. The supplied file order and
#'   time vector are treated as authoritative -- checked for internal
#'   consistency (right length, valid dates, correct order), not that you
#'   paired the right date with the right file; that correspondence
#'   remains your responsibility. Mutually exclusive with `cycle_type`.
#'   Requires `files`.
#' @param cycle_type Optional. One of `"annual"`, `"monthly"`,
#'   `"16-day"`, `"semimonthly"`, `"10-day"`, `"8-day"`, `"weekly"`,
#'   or `"daily"` -- generates real calendar dates instead of
#'   you supplying `time` yourself, for series that genuinely follow one
#'   of these conventions (see "Supported cycle types" above for the
#'   generative rule and real-product examples behind each). If your
#'   series does not follow one of these exact conventions, supply
#'   `time` explicitly instead. Requires `files` and `start`.
#' @param start Required with `cycle_type`. A single `Date` marking the
#'   beginning of the first period. Must fall on a boundary valid for the
#'   declared `cycle_type` (day 1 for `"monthly"`; day 1 or 16 for
#'   `"semimonthly"`; day 1, 11 or 21 for `"10-day"`; year-day 1, 17, 33,
#'   ..., 353 for `"16-day"`; year-day 1, 9, 17, ..., 361 for `"8-day"`;
#'   year-day 1, 8, ..., 358 for `"weekly"`; 1 January for `"annual"`;
#'   any date for `"daily"`) --
#'   this checks that your own declaration is internally consistent, not
#'   that it matches your real product. If your product genuinely starts
#'   elsewhere, it does not follow this `cycle_type` convention --
#'   supply `files` + `time` instead.
#' @param end Optional, only with `cycle_type`. A single `Date`, the final
#'   inclusive end of the last declared period. If omitted, the number of
#'   layers found in `files` determines how many dates are generated from
#'   `start` onward. If supplied, the calendar expected between `start`
#'   and `end` is built first, and the number of layers found must match
#'   that expected count exactly -- this is what detects an incomplete
#'   series (e.g. a declared full year missing one month's file), which
#'   omitting `end` cannot catch on its own.
#' @param time_anchor One of `"start"`, `"centre"` (default), or `"end"`.
#'   Only applicable with `cycle_type = "monthly"`, `"semimonthly"`,
#'   `"10-day"`, `"16-day"`, `"8-day"` or `"weekly"`.
#'   Annual layers are always dated 1 January; daily layers use their
#'   own date. Do not supply `time_anchor` for these two conventions.
#'   Controls which point within
#'   each period's own date
#'   range is assigned as that layer's date -- irrelevant for MK/CMK
#'   (rank-based), but affects the numeric time values OLS/MMK use
#'   directly. For a period with an even number of days, `"centre"`
#'   deterministically picks the earlier of the two central days.
#'
#' @return A `terra::SpatRaster` with one layer per time step, ordered
#'   chronologically -- one layer per input file for ordinary single-layer
#'   formats (the typical case), but a single input file can contribute
#'   more than one layer for multi-layer formats (e.g. a NetCDF file with
#'   several time steps of its own; `files` and the explicit modes fully
#'   support this, tracking which layers came from which file for the
#'   verbose order-check table). Layer names are taken from the
#'   (de-duplicated) file names, and temporal order is stored as proper
#'   time metadata (`terra::time(result)`) rather than discarded after
#'   the verification step above -- used, for instance, as the default
#'   `t` in [inspect_ts_cell()], and by any future function in this
#'   package requiring explicit time coordinates.
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
    dir = NULL, pattern = "\\.(tif|tiff|nc|grd|img|vrt|asc)$",
    order_regex = NULL,
    candidate_regex = c(
      "year([0-9]+)", "_A([0-9]{4})",
      "(19[0-9]{2}|20[0-9]{2})",
      "([0-9]{4})", "([0-9]+)"),
    var = NULL, report = TRUE, verbose = TRUE,
    files = NULL, time = NULL, cycle_type = NULL,
    start = NULL, end = NULL, time_anchor = "centre") {
  finish_timer <- .sptrends_elapsed_timer("read_ordered_stack()", verbose)
  on.exit(finish_timer(), add = TRUE)

  # ------------------------------------------------------------
  # Argument-combination validation -- before touching any file.
  # See ?read_ordered_stack, "Automatic vs. explicit mode", for why
  # explicit declaration (files/time/cycle_type) exists alongside
  # automatic detection from file names: the same numeric shape in a
  # file name can mean genuinely different things across real
  # datasets (e.g. "19820102" is 2 January in some daily products,
  # but the second half of January in PKU-GIMMS' own semimonthly
  # convention) -- automatic detection cannot resolve that safely,
  # only the user declaring it explicitly can.
  # ------------------------------------------------------------
  if (is.null(dir) == is.null(files)) {
    stop("Supply exactly one of 'dir' or 'files'.")
  }
  if (!is.null(time) && !is.null(cycle_type)) {
    stop("'time' and 'cycle_type' are mutually exclusive -- supply ",
         "one or the other, never both.")
  }
  if (!is.null(cycle_type) && is.null(files)) {
    stop("'cycle_type' requires 'files' as an explicit, already ",
         "correctly ordered vector -- automatic sorting is never ",
         "assumed for this path.")
  }
  if (!is.null(cycle_type) && is.null(start)) {
    stop("'cycle_type' requires 'start' (a single Date) -- there is ",
         "no way to generate calendar dates without knowing where ",
         "the declared series begins.")
  }
  if (!is.null(time) && is.null(files)) {
    stop("'time' requires 'files' too (explicit order) -- no file ",
         "order is ever assumed implicitly.")
  }
  if (is.null(cycle_type) && (!is.null(start) || !is.null(end))) {
    stop("'start' and 'end' are only applicable with 'cycle_type'.")
  }
  if (!is.null(cycle_type)) {
    cycle_type <- match.arg(cycle_type,
      c("annual", "monthly", "16-day", "semimonthly", "10-day",
        "8-day", "weekly", "daily"))
  }
  time_anchor_given <- !missing(time_anchor)
  time_anchor <- match.arg(time_anchor, c("start", "centre", "end"))
  if (is.null(cycle_type) && time_anchor_given) {
    stop("'time_anchor' is only applicable with 'cycle_type'.")
  }
  if (!is.null(cycle_type) &&
      !cycle_type %in% c("monthly", "semimonthly", "10-day", "16-day",
                         "8-day", "weekly") &&
      time_anchor_given) {
    stop("'time_anchor' only applies to 'monthly', 'semimonthly', ",
         "'10-day', '16-day', '8-day' and 'weekly' -- ",
         "'annual' uses 1 January and 'daily' uses each day's date.")
  }
  if (!is.null(start)) {
    if (!inherits(start, "Date") || length(start) != 1L || is.na(start)) {
      stop("'start' must be one non-missing Date.")
    }
  }
  if (!is.null(end)) {
    if (!inherits(end, "Date") || length(end) != 1L || is.na(end)) {
      stop("'end' must be one non-missing Date.")
    }
    if (end < start) stop("'end' must not be before 'start'.")
  }
  if (!is.null(cycle_type)) {
    day_ <- as.integer(format(start, "%d"))
    yday_ <- as.integer(format(start, "%j"))
    start_ok <- switch(cycle_type,
      annual      = format(start, "%m-%d") == "01-01",
      monthly     = day_ == 1L,
      semimonthly = day_ %in% c(1L, 16L),
      "10-day"     = day_ %in% c(1L, 11L, 21L),
      "16-day"     = yday_ %in% seq(1L, 353L, by = 16L),
      "8-day"      = yday_ %in% seq(1L, 361L, by = 8L),
      "weekly"     = yday_ %in% seq(1L, by = 7L, length.out = 52L),
      daily       = TRUE
    )
    if (!isTRUE(start_ok)) {
      stop(sprintf(paste(
        "'start' (%s) is not a valid period boundary for",
        "cycle_type='%s'. If your product genuinely starts on a",
        "different day, supply 'files' + 'time' explicitly instead --",
        "it does not follow the calendar convention this function",
        "generates."
      ), start, cycle_type))
    }
  }
  if (!is.null(files)) {
    if (!is.character(files) || length(files) == 0L || anyNA(files)) {
      stop("'files' must be a non-empty character vector.")
    }
    if (anyDuplicated(files)) {
      stop("'files' must not contain duplicate paths.")
    }
    if (any(!file.exists(files))) {
      stop("Some paths supplied in 'files' do not exist.")
    }
  }

  if (!is.null(files)) {
    return(.read_ordered_stack_explicit(
      files = files, time = time, cycle_type = cycle_type,
      start = start, end = end, time_anchor = time_anchor,
      var = var, report = report, verbose = verbose
    ))
  }

  detected_files <- list.files(dir, pattern = pattern, full.names = TRUE)
  if (length(detected_files) == 0) {
    stop("No files matching '", pattern, "' in: ", dir)
  }

  file_names <- basename(detected_files)

  try_regex <- function(regex) {
    numbers <- suppressWarnings(as.numeric(
      sub(paste0(".*", regex, ".*"), "\\1", file_names)
    ))
    if (length(numbers) != length(file_names) || anyNA(numbers)) return(NULL)
    if (anyDuplicated(numbers) > 0) return(NULL)
    numbers
  }

  pattern_used <- NULL
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
      message(
        "Automatic mode: order detected from file names. For higher ",
        "reliability -- especially if the series is not annual -- ",
        "supplying 'files' explicitly (with 'time' or 'cycle_type') ",
        "is recommended. See ?read_ordered_stack."
      )
    }
  }

  # A thin, purely internal wrapper around order() -- exists so this
  # exact permutation step can be safely mocked in a test (mocking this
  # package's own internal function) without touching base::order()
  # itself (which every other part of R also depends on, and mocking
  # broadly would risk destabilising the test run in ways unrelated to
  # this package's own logic).
  ord <- .order_for_stacking(order_numbers)
  detected_files <- detected_files[ord]
  ordered_numbers <- order_numbers[ord]

  valid_sequence <- anyDuplicated(ordered_numbers) == 0 &&
    !any(diff(ordered_numbers) <= 0)

  if (verbose) {
    order_table <- data.frame(
      stack_position = seq_along(detected_files),
      detected_number = ordered_numbers,
      file = basename(detected_files)
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
    # Gap detection is only meaningful for a simple, unambiguously
    # annual sequence of consecutive integers (e.g. years): a
    # combined code (like a YYYYMM-style "202312", "202401" for
    # consecutive December/January) would differ by far more than 1
    # between genuinely consecutive entries, and flagging that as a
    # "gap" would be a false alarm, not a real one. Scoped to the two
    # candidate patterns that extract a genuinely simple, single-field
    # sequential index ("year([0-9]+)" and a bare 4-digit year) --
    # not the looser fallback patterns, which are more likely to be
    # part of a larger combined code.
    looks_like_plain_years <- !is.null(pattern_used) &&
      (pattern_used == "year([0-9]+)" ||
       (pattern_used == "(19[0-9]{2}|20[0-9]{2})" &&
          all(ordered_numbers >= 1900 & ordered_numbers <= 2099)))
    if (looks_like_plain_years) {
      message("Note: there are gaps in the numbering ",
              "(not all years/steps are consecutive).")
      message("Check whether this is intentional (missing years) or a ",
              "misnamed file.")
    } else {
      message("Note: detected numbers are not evenly spaced by 1 -- this ",
              "is expected for non-simple-year patterns (e.g. combined ",
              "codes) and is not flagged as a gap. Inspect the order ",
              "table above manually if in doubt.")
    }
  }

  if (isTRUE(report)) {
    graphics::plot(seq_along(ordered_numbers), ordered_numbers,
         type = "b", pch = 19, col = "steelblue",
         xlab = "Position in the stack",
         ylab = "Detected order number",
         main = "Temporal order check\n(must increase chronologically)")
    graphics::abline(a = ordered_numbers[1] - 1, b = 1, col = "red", lty = 2)
    graphics::legend("topleft",
           legend = c("Detected order", "Expected diagonal"),
           col = c("steelblue", "red"), lty = c(1, 2), pch = c(19, NA),
           bty = "n")
  }

  pb <- .sptrends_progress(length(detected_files), "Reading files", verbose)
  if (!is.null(var)) {
    layers <- vector("list", length(detected_files))
    for (i in seq_along(detected_files)) {
      layers[[i]] <- terra::rast(detected_files[i], subds = var)
      .sptrends_progress_step(pb, i)
    }
    .sptrends_progress_close(pb)
    s <- do.call(c, layers)
  } else {
    layers <- vector("list", length(detected_files))
    for (i in seq_along(detected_files)) {
      layers[[i]] <- terra::rast(detected_files[i])
      .sptrends_progress_step(pb, i)
    }
    .sptrends_progress_close(pb)
    s <- do.call(c, layers)
  }
  if (!all(terra::inMemory(s))) s <- s + 0

  names_ <- tools::file_path_sans_ext(basename(detected_files))
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
         main = "Temporal order check (NetCDF)\n(must increase chronologically)")
  }

  names(r) <- as.character(ordered_time)
  if (verbose) {
    message(sprintf("Stack built from NetCDF: %d layers, %d x %d cells.",
                     terra::nlyr(r), terra::nrow(r), terra::ncol(r)))
  }
  r
}

# Generates real calendar dates for one of the eight unambiguous
# cycle_type conventions, following genuine calendar arithmetic
# (month lengths, leap years) throughout -- never naive interval
# division. Two modes: give 'n' (generate exactly that many periods
# from 'start' onward), or give 'end' (generate the full calendar
# between 'start' and 'end', however many periods that produces).
#' @noRd
.generate_calendar_dates <- function(cycle_type, start, n = NULL,
                                      end = NULL, time_anchor = "centre") {
  last_day_of_month <- function(d) {
    y <- as.integer(format(d, "%Y"))
    m <- as.integer(format(d, "%m"))
    if (m == 12L) {
      as.Date(sprintf("%d-12-31", y))
    } else {
      as.Date(sprintf("%d-%02d-01", y, m + 1L)) - 1
    }
  }

  period_starts <- vector("list", 0)
  period_ends <- vector("list", 0)
  cursor <- start

  repeat {
    if (cycle_type == "annual") {
      p_end <- as.Date(sprintf("%d-12-31", as.integer(format(cursor, "%Y"))))
      next_cursor <- as.Date(sprintf(
        "%d-01-01", as.integer(format(cursor, "%Y")) + 1L))
    } else if (cycle_type == "monthly") {
      p_end <- last_day_of_month(cursor)
      next_cursor <- p_end + 1
    } else if (cycle_type == "semimonthly") {
      day_ <- as.integer(format(cursor, "%d"))
      if (day_ == 1L) {
        p_end <- as.Date(paste0(format(cursor, "%Y-%m-"), "15"))
        next_cursor <- p_end + 1
      } else {
        p_end <- last_day_of_month(cursor)
        next_cursor <- as.Date(sprintf(
          "%d-%02d-01", as.integer(format(p_end + 1, "%Y")),
          as.integer(format(p_end + 1, "%m"))))
      }
    } else if (cycle_type == "10-day") {
      day_ <- as.integer(format(cursor, "%d"))
      if (day_ == 1L) {
        p_end <- as.Date(paste0(format(cursor, "%Y-%m-"), "10"))
        next_cursor <- p_end + 1
      } else if (day_ == 11L) {
        p_end <- as.Date(paste0(format(cursor, "%Y-%m-"), "20"))
        next_cursor <- p_end + 1
      } else {
        p_end <- last_day_of_month(cursor)
        next_cursor <- as.Date(sprintf(
          "%d-%02d-01", as.integer(format(p_end + 1, "%Y")),
          as.integer(format(p_end + 1, "%m"))))
      }
    } else if (cycle_type == "16-day") {
      # 23 periods/year, capturing the remainder as a shorter final
      # period -- matches ETM's own table (23) and real MOD13Q1 file
      # names (year-day 353 is the last period's start).
      year_ <- as.integer(format(cursor, "%Y"))
      yday_ <- as.integer(format(cursor, "%j"))
      starts_this_year <- seq(1L, 353L, by = 16L)
      idx <- which(starts_this_year == yday_)
      if (idx == length(starts_this_year)) {
        p_end <- as.Date(sprintf("%d-12-31", year_))
        next_cursor <- as.Date(sprintf("%d-01-01", year_ + 1L))
      } else {
        next_yday <- starts_this_year[idx + 1L]
        jan1 <- as.Date(sprintf("%d-01-01", year_))
        p_end <- jan1 + (next_yday - 1L) - 1L
        next_cursor <- jan1 + (next_yday - 1L)
      }
    } else if (cycle_type == "8-day") {
      # 46 periods/year, capturing the remainder as a shorter final
      # period -- matches ETM's own table exactly.
      year_ <- as.integer(format(cursor, "%Y"))
      yday_ <- as.integer(format(cursor, "%j"))
      starts_this_year <- seq(1L, 361L, by = 8L)
      idx <- which(starts_this_year == yday_)
      if (idx == length(starts_this_year)) {
        p_end <- as.Date(sprintf("%d-12-31", year_))
        next_cursor <- as.Date(sprintf("%d-01-01", year_ + 1L))
      } else {
        next_yday <- starts_this_year[idx + 1L]
        jan1 <- as.Date(sprintf("%d-01-01", year_))
        p_end <- jan1 + (next_yday - 1L) - 1L
        next_cursor <- jan1 + (next_yday - 1L)
      }
    } else if (cycle_type == "weekly") {
      # 52 periods/year, ETM's own number -- only complete 7-day
      # periods count (1, 8, ..., 358); the last day of the year
      # (365, or 365-366 in a leap year) belongs to no period under
      # this convention, unlike "8-day"/"16-day", which
      # capture the remainder instead. Matches the same asymmetry
      # already found and accepted for ETM's own table (confirmed
      # not internally consistent across its own listed conventions).
      year_ <- as.integer(format(cursor, "%Y"))
      yday_ <- as.integer(format(cursor, "%j"))
      starts_this_year <- seq(1L, by = 7L, length.out = 52L)
      idx <- which(starts_this_year == yday_)
      jan1 <- as.Date(sprintf("%d-01-01", year_))
      if (idx == length(starts_this_year)) {
        p_end <- jan1 + (starts_this_year[idx] - 1L) + 6L
        next_cursor <- as.Date(sprintf("%d-01-01", year_ + 1L))
      } else {
        next_yday <- starts_this_year[idx + 1L]
        p_end <- jan1 + (next_yday - 1L) - 1L
        next_cursor <- jan1 + (next_yday - 1L)
      }
    } else {
      p_end <- cursor
      next_cursor <- cursor + 1
    }

    period_starts[[length(period_starts) + 1L]] <- cursor
    period_ends[[length(period_ends) + 1L]] <- p_end

    if (!is.null(n) && length(period_starts) >= n) break
    if (!is.null(end) && p_end >= end) break
    cursor <- next_cursor
  }

  period_starts <- do.call(c, period_starts)
  period_ends <- do.call(c, period_ends)

  result <- if (cycle_type %in% c("annual", "daily")) {
    period_starts
  } else switch(time_anchor,
    start  = period_starts,
    end    = period_ends,
    centre = period_starts +
      (as.numeric(period_ends - period_starts)) %/% 2L
  )
  list(dates = result, period_start = period_starts,
       period_end = period_ends)
}

# Orchestrates the three explicit sub-paths of read_ordered_stack():
# files alone (ordinal time), files + time (fully explicit), and
# files + cycle_type + start [+ end] (generated calendar dates).
# Argument-combination validation already happened in the caller --
# this function assumes 'files' is a valid, already-ordered vector.
#' @noRd
.read_ordered_stack_explicit <- function(files, time, cycle_type, start,
                                          end, time_anchor, var, report,
                                          verbose) {
  pb <- .sptrends_progress(length(files), "Reading files", verbose)
  layers <- vector("list", length(files))
  file_of_layer <- character(0)
  for (i in seq_along(files)) {
    if (!is.null(var)) {
      layers[[i]] <- terra::rast(files[i], subds = var)
    } else {
      layers[[i]] <- terra::rast(files[i])
    }
    file_of_layer <- c(file_of_layer,
                       rep(files[i], terra::nlyr(layers[[i]])))
    .sptrends_progress_step(pb, i)
  }
  .sptrends_progress_close(pb)
  s <- do.call(c, layers)
  if (!all(terra::inMemory(s))) s <- s + 0
  n_capas <- terra::nlyr(s)

  if (!is.null(time)) {
    # -------------------- files + time --------------------
    if (!(inherits(time, "Date") || inherits(time, "POSIXct") ||
          is.numeric(time))) {
      stop("'time' must be Date, POSIXct or numeric.")
    }
    if (anyNA(time)) stop("'time' must not contain NA.")
    if (is.numeric(time) && any(!is.finite(time))) {
      stop("'time' must be finite if numeric.")
    }
    if (anyDuplicated(time)) stop("'time' must not contain duplicates.")
    if (any(diff(time) <= 0)) stop("'time' must be strictly increasing.")
    if (length(time) != n_capas) {
      stop(sprintf(paste(
        "'time' has %d values but the files contain %d layers in",
        "total (a single file can hold more than one layer)."
      ), length(time), n_capas))
    }
    terra::time(s) <- time
    if (verbose) {
      message(
        "Explicit mode: the supplied file order and 'time' vector ",
        "are treated as authoritative."
      )
    }

  } else if (!is.null(cycle_type)) {
    # -------------------- files + cycle_type + start [+ end] ----
    if (is.null(end)) {
      generated <- .generate_calendar_dates(
        cycle_type, start, n = n_capas, time_anchor = time_anchor
      )
      fechas <- generated$dates
    } else {
      generated <- .generate_calendar_dates(
        cycle_type, start, end = end, time_anchor = time_anchor
      )
      fechas <- generated$dates
      ultimo_period_end <- generated$period_end[length(generated$period_end)]
      if (ultimo_period_end != end) {
        stop(sprintf(paste(
          "'end' (%s) does not fall exactly on the end of a period for",
          "cycle_type='%s' -- the last generated period actually ends",
          "%s. 'end' must be the exact inclusive end of the last",
          "declared period, not just an upper bound within it."
        ), end, cycle_type, ultimo_period_end))
      }
      if (length(fechas) != n_capas) {
        stop(sprintf(paste(
          "Declared range %s to %s (cycle_type='%s') expects %d",
          "layers; found %d in 'files'."
        ), start, end, cycle_type, length(fechas), n_capas))
      }
    }
    terra::time(s) <- fechas
    if (verbose) {
      message(sprintf(
        "Declared series: %s, starting %s. Expected layers: %d. Found: %d -- matches.",
        cycle_type, start, length(fechas), n_capas
      ))
    }

  } else {
    # -------------------- files alone -- ordinal time ----------
    terra::time(s, tstep = "years") <- seq_len(n_capas)
  }

  if (verbose) {
    order_table <- data.frame(
      Layer = seq_len(n_capas),
      File = basename(file_of_layer),
      Time = as.character(terra::time(s))
    )
    message("Temporal order verification (mandatory, cannot be skipped):")
    message(paste(utils::capture.output(print(order_table, row.names = FALSE)),
                  collapse = "\n"))
  }

  if (isTRUE(report)) {
    time_numeric <- suppressWarnings(as.numeric(terra::time(s)))
    graphics::plot(seq_len(n_capas), time_numeric,
         type = "b", pch = 19, col = "steelblue",
         xlab = "Position in the stack", ylab = "Time value",
         main = "Temporal order check\n(must increase chronologically)")
  }

  names_ <- tools::file_path_sans_ext(basename(file_of_layer))
  if (anyDuplicated(names_) > 0) {
    if (verbose) {
      message("Note: duplicate layer names -- made unique automatically.")
    }
    names_ <- make.unique(names_, sep = "_")
  }
  names(s) <- names_

  if (verbose) {
    message(sprintf("Stack built: %d layers, %d x %d cells.",
                     terra::nlyr(s), terra::nrow(s), terra::ncol(s)))
  }
  s
}
