# Internal helpers shared across sptrends functions. None of these are
# exported; they replace the duplicated boilerplate (timers, progress bars,
# safe categorical plotting, PNG export) that appeared independently in each
# of the four original scripts this package was built from.

# This package's own fixed 5-colour brand palette (man/figures/logo.png),
# used wherever a plot needs this package's own signature colour rather
# than a generic one -- currently just the "significant" colour in
# fdr_significance_maps() and fdr_threshold_plot(). Deliberately a short,
# named, exact list rather than re-deriving shades from the logo image
# each time a new use for it comes up.
.sptrends_brand <- list(
  navy = "#bf3688",
  royal_blue = "#1565D8",
  azure = "#179BEA",
  cyan = "#e2a0c9",
  white = "#FFFFFF"
)

# This package's own diverging colour ramp for continuous statistic maps
# (Sm/S/beta, rho, Theil-Sen slope, etc.): genuine ColorBrewer RdBu
# endpoints (#2166ac blue, white midpoint, #b2182b red), not an
# hcl.colors() approximation -- a single source of truth so every
# diverging map in the package uses exactly the same three colours.
#' @noRd
.sptrends_diverging_palette <- function(n = 50) {
  grDevices::colorRampPalette(c("#2166ac", "white", "#b2182b"))(n)
}

#' @noRd
.reference_alpha <- function(alpha) {
  # The conventional 0.05 significance level, if present in the vector of
  # thresholds being reported -- not the strictest (minimum) one, which
  # would silently understate significance relative to what "alpha=0.05"
  # normally means. Falls back to the strictest value only if 0.05 was
  # not among the thresholds requested.
  if (0.05 %in% alpha) 0.05 else min(alpha)
}

# Validate the time coordinate once and use the same contract everywhere
# a user can supply it. Pairwise-slope and prewhitening formulae require
# unique, chronologically increasing values; accepting duplicates or a
# recycled vector would otherwise produce Inf/NaN values or, worse, a
# numerically plausible result attached to the wrong dates.
#' @noRd
.validate_time_axis <- function(t, n, arg = "t") {
  if (!is.numeric(t)) {
    stop("'", arg, "' must be a numeric vector.")
  }
  if (length(t) != n) {
    stop("'", arg, "' must have length ", n,
         " (one value per raster layer).")
  }
  if (anyNA(t) || any(!is.finite(t))) {
    stop("'", arg, "' must contain only finite, non-missing values.")
  }
  if (anyDuplicated(t)) {
    stop("'", arg, "' must not contain duplicate time values.")
  }
  if (any(diff(t) <= 0)) {
    stop("'", arg, "' must be strictly increasing in chronological order.")
  }
  as.numeric(t)
}

#' @noRd
.validate_positive_integer <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1 || x > .Machine$integer.max || x != floor(x)) {
    stop("'", arg, "' must be one positive integer.")
  }
  as.integer(x)
}

#' @noRd
.validate_odd_window_size <- function(x, arg = "window_size") {
  x <- .validate_positive_integer(x, arg)
  if (x < 3L || x %% 2L == 0L) {
    stop("'", arg, "' must be one odd integer greater than or equal to 3.")
  }
  x
}

#' @noRd
.validate_probability <- function(x, arg, vector = FALSE) {
  length_ok <- if (vector) length(x) >= 1L else length(x) == 1L
  if (!is.numeric(x) || !length_ok || anyNA(x) || any(!is.finite(x)) ||
      any(x <= 0 | x >= 1)) {
    qualifier <- if (vector) "one or more" else "one"
    stop("'", arg, "' must contain ", qualifier,
         " finite numeric value", if (vector) "s" else "",
         " strictly between 0 and 1.")
  }
  as.numeric(x)
}

#' @noRd
.validate_seed <- function(seed, arg = "seed") {
  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
       !is.finite(seed))) {
    stop("'", arg, "' must be NULL or one finite numeric value.")
  }
  seed
}

#' @noRd
.validate_max_pairs <- function(x, arg = "max_pairs") {
  if (length(x) != 1L || !is.numeric(x) || is.na(x) ||
      (!is.infinite(x) &&
       (!is.finite(x) || x < 1 || x > .Machine$integer.max ||
        x != floor(x))) ||
      (is.infinite(x) && x < 0)) {
    stop("'", arg, "' must be one positive integer or Inf.")
  }
  if (is.infinite(x)) Inf else as.integer(x)
}

#' @noRd
.validate_positive_numeric <- function(x, arg) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x <= 0) {
    stop("'", arg, "' must be one finite positive numeric value.")
  }
  as.numeric(x)
}

#' @noRd
.check_reserved_args <- function(args, reserved, arg) {
  if (!is.list(args) ||
      (length(args) > 0L &&
       (is.null(names(args)) || anyNA(names(args)) ||
        any(names(args) == "")))) {
    stop("'", arg, "' must be a fully named list.")
  }
  conflicts <- intersect(names(args), reserved)
  if (length(conflicts)) {
    stop("'", arg, "' cannot override workflow-managed argument",
         if (length(conflicts) > 1L) "s" else "", ": ",
         paste(conflicts, collapse = ", "), ".")
  }
  invisible(args)
}

#' @noRd
.with_timer <- function(label, expr, verbose = TRUE) {
  t0 <- Sys.time()
  result <- force(expr)
  if (verbose) {
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    message(sprintf(">> [%s] completed in %.2f s", label, dt))
  }
  result
}

# Return a closure suitable for on.exit(). It reports time spent even when
# execution stops with an error, so the message deliberately says "elapsed"
# rather than implying successful completion.
#' @noRd
.sptrends_elapsed_timer <- function(label, verbose = TRUE) {
  start <- proc.time()[["elapsed"]]
  force(label)
  force(verbose)
  function() {
    if (isTRUE(verbose)) {
      elapsed <- proc.time()[["elapsed"]] - start
      message(sprintf(">> [%s] elapsed: %.2f s", label, elapsed))
    }
    invisible(NULL)
  }
}

#' @noRd
.sptrends_progress <- function(total, label, verbose = TRUE) {
  if (!verbose || total <= 0) return(NULL)
  progress <- new.env(parent = emptyenv())
  progress$total <- as.numeric(total)
  progress$label <- as.character(label)
  progress$start <- proc.time()[["elapsed"]]
  progress$last_update <- -Inf
  progress$last_step <- 0
  progress$closed <- FALSE
  class(progress) <- "sptrends_progress"
  cat(sprintf(
    "%s | progress: 0.0%% | elapsed: 00:00 | remaining: estimating",
    progress$label
  ), file = stderr())
  utils::flush.console()
  progress
}

#' @noRd
.sptrends_progress_step <- function(pb, i) {
  if (is.null(pb)) return(invisible(NULL))
  if (!inherits(pb, "sptrends_progress") || isTRUE(pb$closed)) {
    return(invisible(NULL))
  }
  i <- min(max(as.numeric(i), 0), pb$total)
  now <- proc.time()[["elapsed"]]
  elapsed <- max(now - pb$start, 0)
  should_update <- pb$last_step == 0 || i >= pb$total ||
    now - pb$last_update >= 0.25
  if (!should_update) return(invisible(NULL))
  remaining <- if (i > 0 && elapsed > 0) {
    elapsed * (pb$total - i) / i
  } else {
    NA_real_
  }
  remaining_text <- if (is.finite(remaining)) {
    .sptrends_format_duration(remaining)
  } else {
    "estimating"
  }
  cat(sprintf(
    "\r%s | progress: %5.1f%% (%d/%d) | elapsed: %s | remaining: %s",
    pb$label, 100 * i / pb$total, as.integer(i), as.integer(pb$total),
    .sptrends_format_duration(elapsed), remaining_text
  ), file = stderr())
  utils::flush.console()
  pb$last_update <- now
  pb$last_step <- i
  invisible(NULL)
}

#' @noRd
.sptrends_progress_close <- function(pb) {
  if (is.null(pb)) return(invisible(NULL))
  if (!inherits(pb, "sptrends_progress") || isTRUE(pb$closed)) {
    return(invisible(NULL))
  }
  cat("\n", file = stderr())
  utils::flush.console()
  pb$closed <- TRUE
  invisible(NULL)
}

#' @noRd
.sptrends_format_duration <- function(seconds) {
  if (!is.numeric(seconds) || length(seconds) != 1L ||
      !is.finite(seconds) || seconds < 0) {
    return("unknown")
  }
  seconds <- as.integer(round(seconds))
  hours <- seconds %/% 3600L
  minutes <- (seconds %% 3600L) %/% 60L
  remainder <- seconds %% 60L
  if (hours > 0L) {
    sprintf("%02d:%02d:%02d", hours, minutes, remainder)
  } else {
    sprintf("%02d:%02d", minutes, remainder)
  }
}

# Safe categorical raster plotting: terra::plot() infers the number of
# categories actually present in the data. If a fixed-length legend is
# supplied but a category is completely absent in this particular run
# (e.g. zero significant cells), terra errors out instead of just omitting
# that legend entry. This builds the legend only from values that are
# actually present.
#
# IMPORTANT: 'values' must be given in ASCENDING numeric order, with
# 'colours'/'labels' aligned to that same order (colours[i]/labels[i]
# describes values[i]). 'idx <- values %in% present' is a logical FILTER,
# not a reordering -- colours[idx]/labels[idx] preserve whatever order
# they were written in, dropping only the entries for absent categories.
# terra::plot(type = "classes") itself assigns colours to the raster's
# own ascending sorted unique values, one-to-one. Writing 'values' in
# any other order (e.g. to try to control which category appears first
# in the legend) silently misassigns colours to the wrong cells -- this
# happened once in early 1.2.x development and was caught and fixed;
# do not repeat that mistake. Legend *display* order is a separate,
# unrelated concern (controlled by terra::plot()'s own 'sort'/'reverse'
# arguments, not by this function's own vectors).
#' @noRd
.safe_categorical_plot <- function(x, values, colours, labels, main = NULL,
                                    ...) {
  present <- sort(unique(stats::na.omit(
    as.numeric(terra::values(x, mat = FALSE))
  )))
  idx <- values %in% present
  if (!any(idx)) {
    message(sprintf(
      "Note: no valid cells to plot for '%s' -- map skipped.",
      if (!is.null(main)) main else "this raster"
    ))
    return(invisible(NULL))
  }
  terra::plot(x, col = colours[idx], main = main,
              plg = list(legend = labels[idx]), ...)
}

# Save the most recently drawn plot as a PNG at an explicit path. Never
# called automatically by any exported function -- only when the caller
# supplies a path.
#' @noRd
.save_current_plot <- function(path, width = 1000, height = 800, res = 120) {
  grDevices::dev.copy(grDevices::png, filename = path, width = width,
                       height = height, res = res)
  invisible(grDevices::dev.off())
}

# Cross-platform parallel lapply: uses a PSOCK cluster (works on Windows,
# macOS, and Linux alike -- unlike parallel::mclapply, which is fork-based
# and unavailable on Windows). Falls back to plain lapply when n_cores <= 1.
# 'export_vars'/'export_env' let the caller ship the objects each worker
# needs (PSOCK workers start with an empty environment, nothing is shared
# by default). 'seed' uses the L'Ecuyer-CMRG stream for reproducible
# parallel RNG, distinct per worker.
#' @noRd
.sptrends_parallel_lapply <- function(X, FUN, n_cores = 1,
                                       export_vars = character(0),
                                       export_env = parent.frame(), seed = NULL,
                                       packages = character(0),
                                       shared_cluster = NULL) {
  n_cores <- .validate_positive_integer(n_cores, "n_cores")
  # A pre-made, already-running cluster (typically built once by a
  # calling workflow function -- see .sptrends_shared_cluster() below --
  # and reused across several of this package's own parallel steps in
  # the same call) skips this function's own create/tear-down entirely.
  # Its own n_cores is then irrelevant here: the shared cluster's size
  # was already decided by whoever built it.
  if (!is.null(shared_cluster)) {
    if (length(packages) > 0) {
      parallel::clusterCall(shared_cluster, function(pkgs) {
        for (p in pkgs) requireNamespace(p, quietly = TRUE)
        NULL
      }, packages)
    }
    if (length(export_vars) > 0) {
      parallel::clusterExport(shared_cluster, export_vars, envir = export_env)
    }
    if (!is.null(seed)) {
      parallel::clusterSetRNGStream(shared_cluster, iseed = seed)
    }
    return(parallel::parLapply(shared_cluster, X, FUN))
  }

  if (n_cores <= 1L) return(lapply(X, FUN))

  max_cores <- parallel::detectCores(logical = TRUE)
  if (!is.na(max_cores) && n_cores > max_cores) {
    message(sprintf(
      paste0("Note: n_cores=%d requested but only %d logical cores ",
             "detected -- using %d."),
      n_cores, max_cores, max_cores))
    n_cores <- max_cores
  }

  cl <- parallel::makeCluster(n_cores, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)

  if (length(packages) > 0) {
    parallel::clusterCall(cl, function(pkgs) {
      for (p in pkgs) requireNamespace(p, quietly = TRUE)
      NULL
    }, packages)
  }
  if (length(export_vars) > 0) {
    parallel::clusterExport(cl, export_vars, envir = export_env)
  }
  if (!is.null(seed)) {
    parallel::clusterSetRNGStream(cl, iseed = seed)
  }

  parallel::parLapply(cl, X, FUN)
}

# Builds one PSOCK cluster meant to be shared across several of this
# package's own parallel steps within a single workflow_tst()/
# workflow_rta() call, instead of each step creating and tearing down
# its own -- avoiding the repeated process-spawn overhead of doing that
# 2-3 times in a row for what is, from the user's own perspective, one
# parallel request (workflow_tst(n_cores = 4), not "parallelise the CMK
# step, then separately parallelise the Theil-Sen step"). Returns NULL
# (not a cluster) when n_cores <= 1, so callers can pass its result
# straight through to .sptrends_parallel_lapply()'s own shared_cluster
# argument unconditionally, with sequential execution the natural
# result when there is no cluster to share.
#' @noRd
.sptrends_shared_cluster <- function(n_cores) {
  # Internal compatibility contract: callers that use this helper directly
  # may use zero (or another finite value <= 1) to request no cluster.
  # Exported workflows validate n_cores as a positive integer before they
  # reach this point.
  if (is.numeric(n_cores) && length(n_cores) == 1L && !is.na(n_cores) &&
      is.finite(n_cores) && n_cores <= 1) {
    return(NULL)
  }
  n_cores <- .validate_positive_integer(n_cores, "n_cores")

  max_cores <- parallel::detectCores(logical = TRUE)
  if (!is.na(max_cores) && n_cores > max_cores) {
    message(sprintf(
      paste0("Note: n_cores=%d requested but only %d logical cores ",
             "detected -- using %d."),
      n_cores, max_cores, max_cores))
    n_cores <- max_cores
  }
  parallel::makeCluster(n_cores, type = "PSOCK")
}

# Shared spatial adjacency builder behind both trend_test()'s own
# prepare_cmk_neighbourhood() and spatial_autocorrelation()'s own
# internal neighbourhood setup -- previously two independent, near-
# identical implementations (same terra::adjacent() call, same sparse
# Matrix construction), diverging only in which one or two fields of
# the same underlying adjacency matrix each caller actually needed
# (nb_count -- valid-neighbour count per cell -- for CMK's own
# variance adjustment; S0 -- the total sum of W, i.e. sum(nb_count) --
# for Moran's I/Getis-Ord's own normalisation). Unified here so a
# raster's spatial structure is computed once, not twice, when both
# functions are run on the same geometry with the default 3 by 3 structure
# (exactly analyses that use both diagnostics and canonical CMK). Both
# callers' own precomputed arguments share that structure. Broader CMK
# windows deliberately carry a different signature and cannot be reused by
# a diagnostic expecting the default adjacency.
#' @noRd
.prepare_spatial_neighbourhood <- function(x, ok, connectivity = "queen",
                                            window_size = 3L) {
  connectivity <- match.arg(connectivity, c("queen", "rook"))
  window_size <- .validate_odd_window_size(window_size)
  ncell_ <- terra::ncell(x)
  directions <- connectivity
  if (window_size != 3L) {
    directions <- matrix(0L, nrow = window_size, ncol = window_size)
    centre <- (window_size + 1L) %/% 2L
    if (connectivity == "queen") {
      directions[] <- 1L
    } else {
      directions[centre, ] <- 1L
      directions[, centre] <- 1L
    }
    directions[centre, centre] <- 0L
  }
  adj <- terra::adjacent(x, cells = 1:ncell_, directions = directions,
                          pairs = TRUE)
  from <- adj[, 1]
  to <- adj[, 2]
  valid <- ok[from] & ok[to]
  from <- from[valid]
  to <- to[valid]

  W <- Matrix::sparseMatrix(i = from, j = to, x = 1, dims = c(ncell_, ncell_))
  nb_count <- as.numeric(Matrix::rowSums(W))
  signature <- list(
    nrow = terra::nrow(x),
    ncol = terra::ncol(x),
    ncell = ncell_,
    extent = as.vector(terra::ext(x)),
    resolution = terra::res(x),
    crs = terra::crs(x, proj = TRUE),
    connectivity = connectivity,
    window_size = window_size,
    ok = as.logical(ok)
  )
  list(W = W, nb_count = nb_count, S0 = sum(nb_count),
       signature = signature)
}

#' @noRd
.validate_spatial_neighbourhood <- function(neighbourhood, x, ok,
                                             connectivity = "queen",
                                             window_size = 3L) {
  window_size <- .validate_odd_window_size(window_size)
  ncell_ <- terra::ncell(x)
  if (!is.list(neighbourhood) ||
      !all(c("W", "nb_count", "signature") %in% names(neighbourhood))) {
    stop("'precomputed_neighbourhood' must be an object returned by ",
         "prepare_cmk_neighbourhood() or the package's spatial ",
         "neighbourhood builder.")
  }
  w_dim <- dim(neighbourhood$W)
  if (!inherits(neighbourhood$W, "Matrix") ||
      length(w_dim) != 2L ||
      !all(as.numeric(w_dim) == c(ncell_, ncell_))) {
    stop("'precomputed_neighbourhood$W' must be a ", ncell_, " x ",
         ncell_, " Matrix.")
  }
  if (!is.numeric(neighbourhood$nb_count) ||
      length(neighbourhood$nb_count) != ncell_ ||
      any(!is.finite(neighbourhood$nb_count)) ||
      any(neighbourhood$nb_count < 0)) {
    stop("'precomputed_neighbourhood$nb_count' must contain one finite, ",
         "non-negative value per raster cell.")
  }

  sig <- neighbourhood$signature
  expected <- list(
    nrow = terra::nrow(x),
    ncol = terra::ncol(x),
    ncell = ncell_,
    extent = as.vector(terra::ext(x)),
    resolution = terra::res(x),
    crs = terra::crs(x, proj = TRUE),
    connectivity = connectivity,
    window_size = window_size,
    ok = as.logical(ok)
  )
  geometry_ok <- isTRUE(all.equal(as.numeric(sig$nrow),
                                  as.numeric(expected$nrow), tolerance = 0)) &&
    isTRUE(all.equal(as.numeric(sig$ncol),
                     as.numeric(expected$ncol), tolerance = 0)) &&
    isTRUE(all.equal(as.numeric(sig$ncell),
                     as.numeric(expected$ncell), tolerance = 0)) &&
    isTRUE(all.equal(sig$extent, expected$extent, tolerance = 0)) &&
    isTRUE(all.equal(sig$resolution, expected$resolution, tolerance = 0)) &&
    identical(sig$crs, expected$crs)
  if (!geometry_ok) {
    stop("'precomputed_neighbourhood' was built for different raster geometry.")
  }
  if (!identical(sig$connectivity, connectivity)) {
    stop("'precomputed_neighbourhood' uses connectivity='",
         sig$connectivity, "', but connectivity='", connectivity,
         "' is required.")
  }
  if (!identical(sig$window_size, window_size)) {
    found <- if (is.null(sig$window_size)) "unknown" else sig$window_size
    stop("'precomputed_neighbourhood' uses window_size=", found,
         ", but window_size=", window_size, " is required.")
  }
  if (!identical(sig$ok, expected$ok)) {
    stop("'precomputed_neighbourhood' was built for a different valid-cell ",
         "pattern.")
  }
  invisible(neighbourhood)
}

#' @noRd
.mask_and_smooth_slope <- function(theil_sen, reject, smooth) {
  reject_r <- terra::setValues(theil_sen, as.numeric(reject))
  slope_sig <- terra::mask(theil_sen, reject_r, maskvalues = c(0, NA))
  if (isTRUE(smooth)) {
    w <- matrix(1, nrow = 3, ncol = 3)
    # na.policy = "omit": a cell that was masked out (not significant)
    # must stay NA -- terra::focal()'s own default ("all") would
    # otherwise paint a non-significant cell with a colour just because
    # it happens to sit next to a significant one, visually
    # misrepresenting which cells actually passed the significance
    # test. Shared by the "tst" and "rta" cases of plot.sptrends()'s
    # which = "slope" view, so the fix (and its test) lives once, not
    # duplicated in both files.
    slope_sig <- terra::focal(slope_sig, w = w, fun = "median",
                               na.rm = TRUE, na.policy = "omit")
  }
  slope_sig
}

# Smooths the slope WITHOUT pre-masking by significance first --
# deliberately distinct from .mask_and_smooth_slope() above. That
# function is correct for the "slope" panel (showing a raw slope
# VALUE at a non-significant cell would be misleading), but is wrong
# as input to direction_map(): that function does its OWN
# significance-based masking internally (reject & Sv > 0/< 0), and
# treats "Sv is NA" as "no valid underlying data at all", not "this
# cell was masked out for being non-significant". Feeding it an
# already-significance-masked slope collapses every non-significant
# cell to NA instead of the intended 0 ("not significant") category,
# so the grey "not significant" class silently vanishes from the
# rendered map -- found and fixed as a real bug, not a hypothetical
# concern.
#' @noRd
.smooth_slope_for_direction <- function(theil_sen, smooth) {
  if (!isTRUE(smooth)) return(theil_sen)
  w <- matrix(1, nrow = 3, ncol = 3)
  terra::focal(theil_sen, w = w, fun = "median",
               na.rm = TRUE, na.policy = "omit")
}

#' @noRd
.sptrends_dispatch <- function(generic, x, ...) {
  # Naming convention, not an enumerated list: adding a new "sptrends"
  # subclass (e.g. a future published workflow or core method) only
  # ever requires writing .print_<newclass>(), .summary_<newclass>(),
  # and .plot_<newclass>() -- print.sptrends()/summary.sptrends()/
  # plot.sptrends() themselves never need editing again, since they
  # look the function up by this convention rather than enumerating
  # every known class in a switch(). See sptrends-methods.R.
  cls <- class(x)[1]
  fn_name <- paste0(".", generic, "_", cls)
  fn <- get0(fn_name, envir = asNamespace("sptrends"), inherits = FALSE)
  if (is.null(fn)) {
    stop("Unknown 'sptrends' object class: '", cls, "'.")
  }
  fn(x, ...)
}

#' @noRd
.robust_diverging_range <- function(x, k = 2) {
  # A single extreme cell can otherwise dictate the colour range for
  # every diverging map this package draws (Theil-Sen slope, S/Sm
  # statistic, direction-of-change results, Rho): terra::plot()'s
  # own range = c(-max_abs, max_abs) stretches from the most extreme
  # value seen, so one outlier compresses everyone else's colour into
  # a narrow, washed-out band near the palette's midpoint. Using
  # min(max_abs, k standard deviations) instead keeps the range no
  # WIDER than the actual data (so it never shows a colour more
  # extreme than what is really there), but shrinks it when the data's
  # typical spread is much tighter than its single most extreme cell.
  #
  # IMPORTANT for every caller: terra::plot()'s own DEFAULT for values
  # outside 'range' is to colour them NA (blank/invisible), not
  # saturate them to the palette's extreme colour -- the opposite of
  # what this range is meant to achieve. Every terra::plot() call using
  # this range MUST also pass fill_range = TRUE, or genuinely extreme
  # cells (the exact ones this range narrows around) will silently
  # vanish from the map instead of showing up strongly coloured. This
  # was found and fixed as a real bug in this package's own maps, not
  # a hypothetical concern.
  rng <- terra::global(x, "range", na.rm = TRUE)
  max_abs <- suppressWarnings(max(abs(rng$min), abs(rng$max), na.rm = TRUE))
  if (!is.finite(max_abs) || max_abs == 0) return(c(-1, 1))

  sd_val <- suppressWarnings(terra::global(x, "sd", na.rm = TRUE)$sd)
  if (!is.finite(sd_val) || sd_val == 0) return(c(-max_abs, max_abs))

  limit <- min(max_abs, k * sd_val)
  c(-limit, limit)
}

# Vectorised autocorrelation function across every row of Y at once, for
# every lag from 1 to lag.max, matching R's own acf(x, lag.max)$acf[-1]
# definition (mean-centred per row, normalised by that row's own lag-0
# variance -- not renormalised separately for each lag's own shrinking
# sample). Returns a matrix, rows = cells, columns = lag 1..lag.max.
# Used by trend_test(method = "MMK")'s own rank-autocorrelation step
# (Hamed and Rao, 1998), which needs every lag, not just lag-1 the way
# prewhiten()'s own .lag1_acf_vectorised() does.
#' @noRd
.acf_multilag_vectorised <- function(Y, lag.max) {
  n <- ncol(Y)
  Yc <- Y - rowMeans(Y)
  den <- rowSums(Yc^2)
  # A row with zero variance has no autocorrelation to speak of at any
  # lag -- same convention as .lag1_acf_vectorised()'s own den == 0
  # guard, to avoid NaN propagating into later steps rather than a
  # single, sensible value.
  den[den == 0] <- NA_real_
  out <- matrix(NA_real_, nrow(Y), lag.max)
  for (k in seq_len(lag.max)) {
    lhs <- Yc[, (k + 1):n, drop = FALSE]
    rhs <- Yc[, 1:(n - k), drop = FALSE]
    num <- rowSums(lhs * rhs)
    out[, k] <- num / den
  }
  out[is.na(out)] <- 0
  out
}
