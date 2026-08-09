#' Compare detection methods against a known ground truth
#'
#' Evaluates how well one or more trend-detection methods recover a
#' known truth -- the scientific question this function exists to
#' answer is simple: *which method actually recovers the true trends
#' better?* A confusion-matrix comparison is how it answers that
#' question, not the point of the function itself. Built for
#' [sim_trend_stack()]'s `true_slope`, but deliberately agnostic to
#' where either side comes from. `detections` and `ground_truth` are
#' just logical vectors (or objects that reduce to one): this works
#' equally well comparing `workflow_tst()`
#' against a classic Mann-Kendall, a linear model, a method from another
#' package, or an algorithm you implemented yourself -- as long as each is
#' reduced to a "significant / not significant" call per cell.
#'
#' **Function type:** **Validation function** -- scores detection methods
#' against known truth. It is not part of the inferential pipeline.
#'
#' @section Typical use:
#' Single-run validation:
#' ```
#' sim_trend_stack()
#'     |
#' run one or more detection methods on sim$series
#'     |
#' compare_detections()
#'     |
#' one table of cell-wise detection metrics
#' ```
#' Replicated validation is an alternative route, not a step after the
#' single-run call:
#' ```
#' repeat simulation and detection across seeds
#'     |
#' collect detections and ground truths in lists
#'     |
#' compare_detections(replicates = TRUE)
#'     |
#' mean and standard deviation of metrics + optional empirical FWER
#' ```
#' See the examples below for both routes worked through in full.
#'
#' @section Methodological details:
#' **Why use simulated data?**
#'
#' Real environmental datasets have no known "correct" answer, so
#' detection performance cannot be measured directly -- there is no
#' ground truth to score against for a real NDVI or temperature raster.
#' Simulation studies ([sim_trend_stack()]) provide a known ground truth
#' instead, allowing sensitivity, specificity, precision, and the other
#' metrics below to be computed objectively, something no amount of
#' visual inspection of a real map can substitute for.
#'
#' **A platform for comparing methods, not just running them**
#'
#' Together, [sim_trend_stack()], this function, and the two published
#' workflows this package offers ([workflow_tst()] and
#' [workflow_rta()]) mean sptrends is not only a tool for running a
#' trend analysis, but also a platform for validating and comparing
#' *how well* different methods do it -- including methods this package
#' did not itself implement, since both `detections` and `ground_truth`
#' are accepted in plain, method-agnostic forms (see the `detections`
#' argument below).
#'
#' **Methods and metric selection**
#'
#' All requested metrics derive from the same cell-wise confusion matrix.
#' Sensitivity and specificity describe conditional detection behaviour;
#' precision and FDR describe the rejection set; MCC provides a balanced
#' summary under class imbalance. Empirical FDR and FWER require repeated
#' simulations and therefore are only interpretable in replicated mode.
#'
#' **Statistical assumptions and limitations**
#'
#' Detection and truth must refer to the same cells in the same order and
#' reduce to binary calls. Results describe performance under the supplied
#' simulated scenarios; they do not establish performance for every real
#' environmental process. A single replicate estimates realised
#' proportions, not repeated-sampling error rates.
#'
#' **Computational considerations**
#'
#' Scoring is vectorised over the common evaluation domain and is generally
#' inexpensive relative to generating data or fitting the compared methods.
#' Replicated mode retains method-level summaries rather than raster outputs.
#'
#' **Quality assurance**
#'
#' Tests compare confusion counts and every derived metric with
#' analytically tractable cases, including undefined denominators.
#' Additional tests cover vector/raster inputs, replicated aggregation,
#' empirical FDR and FWER, validation failures, method ordering and S3
#' printing, summaries and plots. In the retained external validation,
#' `sptrends` and `Kendall` produced identical MK decisions, directions and
#' derived performance metrics across 4,000 paired Monte Carlo replicates.
#' See `inst/validation/` and `?sptrends` for scope and limitations.
#'
#' @param detections When `replicates = FALSE` (default): a named list,
#'   one element per method being compared. Each element is the method's
#'   significance call: a logical vector (`TRUE` = called significant), a
#'   single-layer `terra::SpatRaster` of 0/1 or `TRUE`/`FALSE` (coerced
#'   with `!= 0`), or a numeric/integer vector coerced the same way. The
#'   list names become the `Method` column. When `replicates = TRUE`: a
#'   list of such named lists, one per replicate (e.g. one per random
#'   seed) -- every inner list must use the same method names.
#' @param ground_truth When `replicates = FALSE`: the ground truth, in
#'   any of the same forms as one element of `detections` -- most
#'   commonly [sim_trend_stack()]'s `true_slope` (a `SpatRaster`,
#'   non-zero where a true trend exists) or its raw values. Coerced to
#'   logical the same way as `detections`. When `replicates = TRUE`:
#'   either a single ground truth reused for every replicate (the common
#'   case, when the same underlying truth is being detected repeatedly
#'   under different noise draws), or a list the same length as
#'   `detections` with one ground truth per replicate (for the rarer
#'   case where the truth itself is re-simulated every time too, e.g.
#'   [sim_trend_stack()]'s own `true_slope` with a different `seed` on
#'   each call).
#' @param metrics Character vector of derived metrics to include, from
#'   `c("sensitivity", "specificity", "precision", "accuracy", "f1",
#'   "mcc", "fpr", "fdr", "fwer", "type_i", "type_ii", "type_iii",
#'   "field_power", "global_power", "within_image_power",
#'   "directional_power")`
#'   (case-insensitive). The omitted default remains the original eight
#'   confusion-matrix metrics. The
#'   confusion matrix counts (`TP`, `FP`, `TN`, `FN`) are always included
#'   regardless of `metrics` -- every requested metric is computed from
#'   that same confusion matrix, not recomputed independently. `"fwer"`
#'   is structurally different from the other eight -- see its own
#'   entry in "Value" below -- and can only be requested when
#'   `replicates = TRUE`; it errors otherwise, since a single run's own
#'   false-positive count is either zero or not, with no per-replicate
#'   ratio to compute in the first place.
#'   `"type_i"` is the background false-positive proportion; `"type_ii"`
#'   is the missed-signal proportion; `"type_iii"` is the wrong-direction
#'   proportion among detected true signals. `"field_power"` records whether
#'   a replicate produces at least one rejection anywhere in the evaluation
#'   domain. It can equal one because of a false positive alone.
#'   `"global_power"` is stricter: it records whether a replicate detects at
#'   least one true signal, whereas
#'   `"within_image_power"` is the proportion of signal cells detected.
#' @param replicates Logical. `FALSE` (default): score a single run, as
#'   described above. `TRUE`: score every replicate the same way, then
#'   aggregate to the mean and standard deviation of every column, per
#'   method -- a single seed can favour or disfavour a method by chance,
#'   so this is what actually supports a claim like "CMK detects more
#'   true trends" (see `@examples` below). Generating the data and
#'   computing each method's detections for every replicate still stays
#'   as an ordinary loop you write, since that step depends entirely on
#'   which methods you are comparing and cannot be generalised without
#'   either guessing at your workflow or adding a callback-style
#'   interface -- `replicates = TRUE` only does the scoring and
#'   aggregation, for the per-replicate detections/ground-truth you
#'   already computed.
#' @param directions Optional named list of estimated direction vectors, one
#'   per method. Alternatively, each detection may be a list containing
#'   `significant` and `direction`.
#' @param truth_direction Optional known direction vector, normally
#'   [sim_trend_stack()]'s `true_direction`. Required for Type III error and
#'   directional power.
#' @param evaluation_mask Optional logical vector/raster selecting one common
#'   evaluation domain for every method. Missing and `FALSE` cells are not
#'   scored. Method-specific masks are rejected because methods must be
#'   compared against the same cells and the same ground truth.
#' @param verbose Logical. If `TRUE`, reports progress, elapsed time and the
#'   estimated time remaining while runs are scored.
#'
#' @return **When `replicates = FALSE`**: a data frame with one row per
#'   method: `Method`, `TP`, `FP`,
#'   `TN`, `FN`, and the requested metric columns among `Sensitivity`
#'   (a.k.a. recall or power: `TP / (TP + FN)`), `Specificity`
#'   (`TN / (TN + FP)`), `Precision` (`TP / (TP + FP)`), `Accuracy`
#'   (`(TP + TN) / (TP + FP + TN + FN)` -- included because it is widely
#'   expected, though it can be misleading whenever true trends are rare
#'   relative to the whole raster, which is common in these simulations),
#'   `F1` (harmonic mean of precision and sensitivity), `MCC` (Matthews
#'   correlation coefficient -- a single balanced summary of all four
#'   confusion-matrix cells at once, generally preferred over `Accuracy`
#'   or `F1` under the class imbalance a low `trend_fraction` produces),
#'   and `FPR` (false positive rate, `FP / (FP + TN)`). `FDR`, in this
#'   single-run table, is the realised false discovery *proportion*
#'   (Benjamini & Hochberg, 1995) in this one result
#'   (`FP / (FP + TP)`) -- not the false discovery *rate* itself, which
#'   is formally defined as the expectation of that proportion over many
#'   repetitions. `q` in [fdr_correction()] bounds that expectation, not
#'   the value in any single realisation -- see `replicates = TRUE`
#'   below for what actually estimates the expectation this proportion is
#'   a single draw from. By standard convention, this realised proportion is
#'   zero when no discoveries are made; retaining those zeroes is essential
#'   when estimating FDR across replicates. Other ratios with an undefined
#'   denominator (e.g. `Precision` when a method calls nothing significant,
#'   or `MCC` when any confusion-matrix margin is zero) are returned as `NA`,
#'   not `NaN` or `0`.
#'   `FieldPower`, `GlobalPower`, and `WithinImagePower` distinguish,
#'   respectively, whether any cell was rejected, whether at least one true
#'   signal cell was detected, and the fraction of true signal cells detected.
#'
#'   **When `replicates = TRUE`**: a data frame with one row per method:
#'   `Method`, `n_replicates`, and `<column>_mean`/`<column>_sd` for
#'   every column above (`TP_mean`, `TP_sd`, `Sensitivity_mean`,
#'   `Sensitivity_sd`, and so on). `FDR_mean` here is the actual estimate
#'   of the false discovery rate itself (the average of the per-replicate
#'   proportion described above, across replicates) -- the quantity `q`
#'   in [fdr_correction()] is meant to bound, unlike the single-run
#'   `FDR` column above. If `"fwer"` was requested, also a plain `FWER`
#'   column (no `_mean`/`_sd` pair): the proportion of replicates with at
#'   least one false positive (`FP > 0`) -- the family-wise error rate
#'   (Benjamini & Hochberg, 1995's own point of contrast for FDR), a
#'   single already-aggregated rate rather than a per-replicate value
#'   with its own mean and spread.
#'
#' @examples
#' # Simulated data only, deliberately -- compare_detections() needs a
#' # KNOWN ground truth to score against (sim$true_slope below), which
#' # by definition only sim_trend_stack() can provide; real environmental
#' # data has no equivalent "true" answer to compare a detection against.
#' \donttest{
#' sim <- sim_trend_stack(nrow = 15, ncol = 15, n_time = 15, seed = 1)
#'
#' # Two variants of the same test: classic Mann-Kendall (no
#' # neighbourhood averaging) vs. the Contextual version.
#' trend_mk  <- trend_test(sim$series, method = "MK",
#'                          report = FALSE, verbose = FALSE)
#' trend_cmk <- trend_test(sim$series, method = "CMK",
#'                          report = FALSE, verbose = FALSE)
#'
#' # Since sim$true_slope is known exactly (this is simulated data), we
#' # can score each method's raw significance calls against the truth.
#' compare_detections(
#'   detections = list(MK = trend_mk$stats$p <= 0.05,
#'                      CMK = trend_cmk$stats$p <= 0.05),
#'   ground_truth = sim$true_slope, verbose = FALSE
#' )
#'
#' # replicates = TRUE: repeat the comparison across several random
#' # seeds and aggregate -- a single run can favour a method just by
#' # chance. Both detections and the true_slope ground truth change
#' # every replicate here, so ground_truth is a list too, not a single
#' # shared vector.
#' detections_list <- list()
#' truths_list <- list()
#' for (s in 1:10) {
#'   sim_s <- sim_trend_stack(nrow = 12, ncol = 12, n_time = 12, seed = s)
#'   mk_s  <- trend_test(sim_s$series, method = "MK",
#'                        report = FALSE, verbose = FALSE)
#'   cmk_s <- trend_test(sim_s$series, method = "CMK",
#'                        report = FALSE, verbose = FALSE)
#'   detections_list[[s]] <- list(MK = mk_s$stats$p <= 0.05,
#'                                 CMK = cmk_s$stats$p <= 0.05)
#'   truths_list[[s]] <- sim_s$true_slope
#' }
#' compare_detections(detections_list, truths_list, replicates = TRUE,
#'                    verbose = FALSE)
#'
#' # FWER ("did at least one false positive happen in this run at all?")
#' # only means anything across many replicates -- request it alongside
#' # the usual metrics with replicates = TRUE; it errors with
#' # replicates = FALSE, since a single run cannot estimate a rate.
#' compare_detections(detections_list, truths_list, replicates = TRUE,
#'                     metrics = c("sensitivity", "fdr", "fwer"),
#'                     verbose = FALSE)
#' }
#'
#' @references
#' Fawcett, T. (2006) An introduction to ROC analysis. Pattern
#' Recognition Letters, 27(8), 861-874. \doi{10.1016/j.patrec.2005.10.010}
#'
#' Source of the Matthews correlation coefficient (`MCC`):
#' - Matthews, B.W. (1975) Comparison of the predicted and observed
#'   secondary structure of T4 phage lysozyme. Biochimica et Biophysica
#'   Acta (BBA) - Protein Structure, 405(2), 442-451.
#'   \doi{10.1016/0005-2795(75)90109-9}
#'
#' Source of the false discovery rate concept `FDR`/`FDR_mean` estimate,
#' and of the family-wise error rate `FWER` is contrasted against (see
#' "Value" above for both):
#' - Benjamini, Y., & Hochberg, Y. (1995) Controlling the False Discovery
#'   Rate: A Practical and Powerful Approach to Multiple Testing. Journal
#'   of the Royal Statistical Society: Series B, 57, 289-300.
#'   \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' @family validation functions
#' @export
compare_detections <- function(detections, ground_truth,
                                metrics = c("sensitivity", "specificity",
                                            "precision", "accuracy",
                                            "f1", "mcc", "fpr", "fdr",
                                            "fwer", "type_i", "type_ii",
                                            "type_iii", "field_power",
                                            "global_power",
                                            "within_image_power",
                                            "directional_power"),
                                replicates = FALSE, directions = NULL,
                                truth_direction = NULL,
                                evaluation_mask = NULL, verbose = TRUE) {
  finish_timer <- .sptrends_elapsed_timer("compare_detections()", verbose)
  on.exit(finish_timer(), add = TRUE)
  default_metrics <- c("sensitivity", "specificity", "precision",
                       "accuracy", "f1", "mcc", "fpr", "fdr")
  if (missing(metrics)) metrics <- default_metrics
  metrics <- tolower(metrics)
  extended_metrics <- c("type_i", "type_ii", "type_iii",
                        "field_power", "global_power",
                        "within_image_power",
                        "directional_power")
  valid_metrics <- c(default_metrics, "fwer", extended_metrics)
  bad <- setdiff(metrics, valid_metrics)
  if (length(bad) > 0) {
    stop(sprintf("Unknown metric(s): %s. Valid options: %s.",
                 paste(bad, collapse = ", "),
                 paste(valid_metrics, collapse = ", ")))
  }
  if ("fwer" %in% metrics && !isTRUE(replicates)) {
    stop("'fwer' only has a meaning across many replicates (the ",
         "proportion of replicates with at least one false positive) -- ",
         "it cannot be computed from a single run. Set replicates = TRUE, ",
         "or remove 'fwer' from 'metrics'.")
  }

  if (isTRUE(replicates)) {
    if (!is.list(detections) || length(detections) == 0 ||
          !all(vapply(detections, is.list, logical(1)))) {
      stop("When replicates = TRUE, 'detections' must be a list of ",
           "named lists (one inner list per replicate) -- not a single ",
           "named list of detection vectors, which is what replicates = ",
           "FALSE (the default) expects.")
    }
    n_rep <- length(detections)
    # A single shared ground truth (one vector/raster reused for every
    # replicate) is the common case; a list of length n_rep is for the
    # rarer case where the truth itself changes per replicate too (e.g.
    # sim_trend_stack()'s own true_slope, freshly simulated each time).
    gt_is_list <- is.list(ground_truth) && !inherits(ground_truth, "SpatRaster")
    if (gt_is_list && length(ground_truth) != n_rep) {
      stop("When 'ground_truth' is a list (one truth per replicate), it ",
           "must have the same length as 'detections' (", n_rep, ").")
    }
    # "fwer" is not a per-replicate column at all (a single run's FP
    # count is either 0 or not -- there is no per-replicate "FWER" to
    # compute) -- it is derived later, in .aggregate_replicate_tables(),
    # directly from the FP column every per-replicate table already has
    # regardless of 'metrics'. Excluded here so .score_one_run() (which
    # only knows about genuine per-replicate metrics) never sees it.
    per_replicate_metrics <- setdiff(metrics, "fwer")
    progress <- .sptrends_progress(
      n_rep, "Scoring detection replicates", verbose)
    on.exit(.sptrends_progress_close(progress), add = TRUE)
    tables <- lapply(seq_len(n_rep), function(i) {
      det_i <- detections[[i]]
      if (is.null(names(det_i)) || any(names(det_i) == "")) {
        stop(sprintf(
          "detections[[%d]] must be a named list (one name per method).",
          i))
      }
      gt_i <- if (gt_is_list) ground_truth[[i]] else ground_truth
      dir_i <- .replicate_component(directions, i, n_rep, names(det_i))
      truth_dir_i <- .replicate_component(truth_direction, i, n_rep)
      mask_i <- .replicate_component(evaluation_mask, i, n_rep)
      table <- .score_one_run(
        det_i, gt_i, per_replicate_metrics, directions = dir_i,
        truth_direction = truth_dir_i, evaluation_mask = mask_i)
      .sptrends_progress_step(progress, i)
      table
    })
    .sptrends_progress_close(progress)
    result <- .aggregate_replicate_tables(
      tables, compute_fwer = "fwer" %in% metrics)
    class(result) <- c("compare_detections_replicates", "compare_detections",
                        "sptrends", "data.frame")
    return(result)
  }

  if (!is.list(detections) || is.null(names(detections)) ||
      any(names(detections) == "")) {
    stop("'detections' must be a named list (one name per method).")
  }
  progress <- .sptrends_progress(1L, "Scoring detection methods", verbose)
  on.exit(.sptrends_progress_close(progress), add = TRUE)
  result <- .score_one_run(
    detections, ground_truth, metrics, directions = directions,
    truth_direction = truth_direction, evaluation_mask = evaluation_mask)
  .sptrends_progress_step(progress, 1L)
  .sptrends_progress_close(progress)
  class(result) <- c("compare_detections", "sptrends", "data.frame")
  result
}

#' @noRd
.score_one_run <- function(detections, ground_truth, metrics,
                           directions = NULL, truth_direction = NULL,
                           evaluation_mask = NULL) {
  to_vector <- function(x) {
    if (inherits(x, "SpatRaster")) x <- terra::values(x, mat = FALSE)
    as.vector(x)
  }
  to_logical <- function(x) {
    x <- to_vector(x)
    if (is.logical(x)) return(x)
    x != 0
  }

  truth_logical <- to_logical(ground_truth)
  truth_dir <- if (is.null(truth_direction)) {
    NULL
  } else {
    sign(to_vector(truth_direction))
  }
  if (is.list(evaluation_mask) &&
      !inherits(evaluation_mask, "SpatRaster")) {
    stop("For one run, 'evaluation_mask' must be one common vector or ",
         "raster, not a method-specific list.")
  }
  shared_mask <- if (is.null(evaluation_mask)) {
    rep(TRUE, length(truth_logical))
  } else {
    to_logical(evaluation_mask)
  }
  if (length(shared_mask) != length(truth_logical)) {
    stop("'evaluation_mask' and 'ground_truth' must have equal lengths.")
  }
  if (!is.null(truth_dir) && length(truth_dir) != length(truth_logical)) {
    stop("'truth_direction' and 'ground_truth' must have equal lengths.")
  }

  rows <- lapply(names(detections), function(nm) {
    eval_mask <- shared_mask
    detection <- detections[[nm]]
    if (is.list(detection) && !inherits(detection, "SpatRaster")) {
      if (is.null(detection$significant)) {
        stop(sprintf("'%s' is a list but has no 'significant' element.", nm))
      }
      sig <- to_logical(detection$significant)
      method_direction <- detection$direction
    } else {
      sig <- to_logical(detection)
      method_direction <- if (!is.null(directions)) directions[[nm]] else NULL
    }
    if (length(sig) != length(truth_logical)) {
      stop(sprintf(
        "'%s' has length %d, but 'ground_truth' has length %d.",
        nm, length(sig), length(truth_logical)))
    }
    ok <- !is.na(sig) & !is.na(truth_logical) & !is.na(eval_mask) & eval_mask
    sig <- sig[ok]
    truth_ok <- truth_logical[ok]
    truth_dir_ok <- if (is.null(truth_dir)) NULL else truth_dir[ok]
    direction_ok <- if (is.null(method_direction)) {
      NULL
    } else {
      method_direction <- sign(to_vector(method_direction))
      if (length(method_direction) != length(truth_logical)) {
        stop(sprintf(
          "Direction for '%s' and ground truth differ in length.", nm))
      }
      method_direction[ok]
    }

    TP <- sum(sig & truth_ok)
    FP <- sum(sig & !truth_ok)
    TN <- sum(!sig & !truth_ok)
    FN <- sum(!sig & truth_ok)

    safe_div <- function(num, den) if (den == 0) NA_real_ else num / den

    sensitivity <- safe_div(TP, TP + FN)
    specificity <- safe_div(TN, TN + FP)
    precision   <- safe_div(TP, TP + FP)
    accuracy    <- safe_div(TP + TN, TP + FP + TN + FN)
    f1 <- if (is.na(precision) || is.na(sensitivity) ||
                (precision + sensitivity) == 0) {
      NA_real_
    } else {
      2 * precision * sensitivity / (precision + sensitivity)
    }
    mcc_denom <- as.double(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)
    mcc <- if (mcc_denom == 0) {
      NA_real_
    } else {
      (TP * TN - FP * FN) / sqrt(mcc_denom)
    }
    fpr <- safe_div(FP, FP + TN)
    discoveries <- TP + FP
    fdr <- if (discoveries == 0L) 0 else FP / discoveries
    type_i <- fpr
    type_ii <- safe_div(FN, TP + FN)
    field_power <- as.numeric((TP + FP) > 0)
    global_power <- as.numeric(TP > 0)
    within_image_power <- sensitivity
    if (is.null(direction_ok) || is.null(truth_dir_ok)) {
      type_iii <- NA_real_
      directional_power <- NA_real_
    } else {
      true_detected <- sig & truth_ok
      wrong_direction <- true_detected & direction_ok != truth_dir_ok
      correct_direction <- true_detected & direction_ok == truth_dir_ok
      type_iii <- safe_div(sum(wrong_direction), sum(true_detected))
      directional_power <- safe_div(sum(correct_direction), sum(truth_ok))
    }

    all_cols <- data.frame(
      Method = nm, TP = TP, FP = FP, TN = TN, FN = FN,
      Sensitivity = sensitivity, Specificity = specificity,
      Precision = precision, Accuracy = accuracy, F1 = f1, MCC = mcc,
      FPR = fpr, FDR = fdr, TypeI = type_i, TypeII = type_ii,
      TypeIII = type_iii, FieldPower = field_power,
      GlobalPower = global_power,
      WithinImagePower = within_image_power,
      DirectionalPower = directional_power
    )
    metric_name_map <- c(sensitivity = "Sensitivity",
                          specificity = "Specificity",
                          precision = "Precision", accuracy = "Accuracy",
                          f1 = "F1", mcc = "MCC", fpr = "FPR", fdr = "FDR",
                          type_i = "TypeI", type_ii = "TypeII",
                          type_iii = "TypeIII",
                          field_power = "FieldPower",
                          global_power = "GlobalPower",
                          within_image_power = "WithinImagePower",
                          directional_power = "DirectionalPower")
    keep_cols <- c("Method", "TP", "FP", "TN", "FN",
                   unname(metric_name_map[metrics]))
    all_cols[, keep_cols, drop = FALSE]
  })

  do.call(rbind, rows)
}

#' Resolve an optional component that may be shared or vary by replicate
#' @noRd
.replicate_component <- function(x, i, n_rep, method_names = NULL) {
  if (is.null(x)) return(NULL)
  if (!is.null(method_names) && is.list(x) &&
      setequal(names(x), method_names)) {
    return(x)
  }
  if (is.list(x) && !inherits(x, "SpatRaster") && length(x) == n_rep) {
    return(x[[i]])
  }
  x
}


#' @noRd
.aggregate_replicate_tables <- function(tables, compute_fwer = FALSE) {
  combined <- do.call(rbind, tables)
  method_order <- unique(combined$Method)
  metric_cols <- setdiff(names(combined), "Method")

  rows <- lapply(method_order, function(m) {
    sub <- combined[combined$Method == m, metric_cols, drop = FALSE]
    row <- data.frame(Method = m, n_replicates = nrow(sub))
    for (col in metric_cols) {
      # mean(na.rm = TRUE) on an all-NA column returns NaN (0/0), while
      # sd(na.rm = TRUE) on the same all-NA column already returns a
      # clean NA -- an inconsistency between the two columns for the
      # exact same "no valid replicate values for this metric" case.
      # Converting NaN to NA_real_ here matches sd()'s own convention
      # and this package's established pattern elsewhere (see
      # .summary_compare_detections()'s own all-NA handling) of
      # surfacing a clean NA rather than a bare NaN.
      col_mean <- mean(sub[[col]], na.rm = TRUE)
      row[[paste0(col, "_mean")]] <- if (is.nan(col_mean)) {
        NA_real_
      } else {
        col_mean
      }
      row[[paste0(col, "_sd")]] <- stats::sd(sub[[col]], na.rm = TRUE)
    }
    if (isTRUE(compute_fwer)) {
      # FWER (family-wise error rate): the probability of at least one
      # false positive occurring in a single run -- structurally
      # different from every other metric here, which are each a
      # per-replicate *ratio* subsequently averaged across replicates.
      # A single replicate's own FP count is either zero or not; there
      # is no per-replicate "FWER" to compute in the first place. The
      # proportion of replicates with FP > 0 *is* the estimate of FWER
      # itself, computed directly here from the FP column every
      # per-replicate table already has -- not a mean/sd pair like the
      # other metrics, since FWER is already a single rate, not a
      # per-replicate value being averaged.
      row[["FWER"]] <- mean(sub$FP > 0, na.rm = TRUE)
    }
    row
  })

  do.call(rbind, rows)
}


#' @noRd
.print_compare_detections <- function(x, ...) {
  print.data.frame(x, ...)
  invisible(x)
}

#' @noRd
.summary_compare_detections <- function(object, ...) {
  numeric_cols <- names(object)[vapply(object, is.numeric, logical(1))]
  numeric_cols <- setdiff(numeric_cols, c("TP", "FP", "TN", "FN"))
  # For every metric here, a HIGHER value is the better outcome, except
  # FPR, FDR, TypeI, TypeII and TypeIII (all are error rates: a LOWER
  # value is better) -- picking the "winner" with which.max() for
  # these would name the method with the *most* errors as the best
  # one. TypeI is literally the same quantity as FPR under a different
  # name (see the "type_i <- fpr" assignment where these are
  # computed); TypeII is the missed-signal rate; TypeIII is the
  # wrong-direction rate among detections.
  lower_is_better <- c("FPR", "FDR", "TypeI", "TypeII", "TypeIII")
  winners <- vapply(numeric_cols, function(col) {
    vals <- object[[col]]
    if (all(is.na(vals))) return(NA_character_)
    pick <- if (col %in% lower_is_better) which.min else which.max
    object$Method[pick(vals)]
  }, character(1))
  tab <- data.frame(metric = numeric_cols, best_method = winners,
                     row.names = NULL)
  cat("Best-scoring method per metric:\n")
  print.data.frame(tab, row.names = FALSE)
  invisible(tab)
}

#' @noRd
.plot_compare_detections <- function(x, ...) {
  plot_detection_comparison(x, ...)
  invisible(x)
}

#' @noRd
.print_compare_detections_replicates <- function(x, ...) {
  print.data.frame(x, ...)
  invisible(x)
}

#' @noRd
.summary_compare_detections_replicates <- function(object, ...) {
  mean_cols <- grep("_mean$", names(object), value = TRUE)
  mean_cols <- setdiff(mean_cols, c("TP_mean", "FP_mean", "TN_mean",
                                     "FN_mean"))
  # Same "lower is better" exception as .summary_compare_detections()
  # above, checked against the "_mean"-stripped metric name.
  lower_is_better <- c("FPR", "FDR", "TypeI", "TypeII", "TypeIII")
  winners <- vapply(mean_cols, function(col) {
    vals <- object[[col]]
    if (all(is.na(vals))) return(NA_character_)
    metric_name <- sub("_mean$", "", col)
    pick <- if (metric_name %in% lower_is_better) which.min else which.max
    object$Method[pick(vals)]
  }, character(1))
  # Strip the "_mean" suffix for display -- "best method per metric",
  # matching .summary_compare_detections()'s own table, not "best
  # method per metric_mean" which would just be visual noise repeated
  # on every row.
  tab <- data.frame(metric = sub("_mean$", "", mean_cols),
                     best_method = winners, row.names = NULL)
  # FWER is not "_mean"-suffixed (see .aggregate_replicate_tables()'s own
  # comment on why), and unlike every metric above, a *lower* FWER is the
  # better outcome -- which.min(), not which.max(), or this would name
  # the method with the *most* false positives as the "winner".
  if ("FWER" %in% names(object)) {
    fwer_vals <- object[["FWER"]]
    fwer_winner <- if (all(is.na(fwer_vals))) {
      NA_character_
    } else {
      object$Method[which.min(fwer_vals)]
    }
    tab <- rbind(tab, data.frame(metric = "FWER", best_method = fwer_winner,
                                  row.names = NULL))
  }
  cat(sprintf("Best-scoring method per metric (mean across %d replicates):\n",
              object$n_replicates[1]))
  print.data.frame(tab, row.names = FALSE)
  invisible(tab)
}

#' @noRd
.plot_compare_detections_replicates <- function(x, ...) {
  # plot_detection_comparison() expects plain metric column names
  # (e.g. "Sensitivity"), not the "_mean"/"_sd"-suffixed ones this
  # aggregated table has -- build that plain-named view from the mean
  # columns (the point estimate across replicates) and reuse the same
  # plotting function, rather than duplicating its bar-chart logic
  # here for a second time.
  mean_cols <- grep("_mean$", names(x), value = TRUE)
  plain <- x[, c("Method", mean_cols), drop = FALSE]
  names(plain) <- sub("_mean$", "", names(plain))
  plot_detection_comparison(plain, ...)
  invisible(x)
}

#' Bar plot of a compare_detections() comparison
#'
#' A grouped bar plot of the metric columns from a [compare_detections()]
#' result (either a single run, or the mean columns of a `replicates =
#' TRUE` aggregated result), one group of bars per method.
#'
#' **Function type:** **Reporting/derived function** -- summarises or
#' plots the output of
#' another function; it does not compute any new statistic. Not
#' exported -- reachable from outside the package via `plot()`.
#'
#' @param comparison Output of [compare_detections()], or a data frame
#'   with the same plain (non-suffixed) metric column names -- e.g. the
#'   mean-columns view [compare_detections()]'s own `replicates = TRUE`
#'   plot dispatch constructs internally.
#' @param metrics Character vector of column names in `comparison` to
#'   plot. Default: the six proportion-based metrics from
#'   [compare_detections()] (`TP`/`FP`/`TN`/`FN` are counts on a different
#'   scale and are not plotted by default).
#' @param path Character or `NULL`. If supplied, a PNG is written there.
#'
#' @return `NULL`, invisibly.
#'
#' @examples
#' sim <- sim_trend_stack(nrow = 12, ncol = 12, n_time = 12, seed = 1)
#' trend_mk  <- trend_test(sim$series, method = "MK",
#'                          report = FALSE, verbose = FALSE)
#' trend_cmk <- trend_test(sim$series, method = "CMK",
#'                          report = FALSE, verbose = FALSE)
#' comparison <- compare_detections(
#'   detections = list(MK = trend_mk$stats$p <= 0.05,
#'                      CMK = trend_cmk$stats$p <= 0.05),
#'   ground_truth = sim$true_slope
#' )
#'
#' # A grouped bar chart of the table above -- one group of bars per
#' # method, one bar per metric.
#' sptrends:::plot_detection_comparison(comparison)
#'
#' @family validation functions
#' @keywords internal
plot_detection_comparison <- function(comparison,
                                       metrics = c("Sensitivity", "Specificity",
                                                   "Precision", "F1", "FPR",
                                                   "FDR"),
                                       path = NULL) {
  metrics <- metrics[metrics %in% names(comparison)]
  if (length(metrics) == 0) {
    stop("None of the requested 'metrics' are present in 'comparison'. ",
         "Did you mean to pass the output of compare_detections()?")
  }

  vals <- t(as.matrix(comparison[, metrics, drop = FALSE]))
  storage.mode(vals) <- "double"
  colnames(vals) <- comparison$Method

  cols <- grDevices::hcl.colors(length(metrics), "Set 2")
  op <- graphics::par(mar = c(5, 4, 4, 8) + 0.1, xpd = TRUE)
  on.exit(graphics::par(op), add = TRUE)
  bp <- graphics::barplot(vals, beside = TRUE, col = cols, border = "white",
                           ylab = "Value", main = "Detection comparison")
  legend_x <- max(bp) + 1
  graphics::legend(legend_x, max(vals, na.rm = TRUE), legend = metrics,
                    fill = cols, bty = "n", cex = 0.85)

  if (!is.null(path)) .save_current_plot(path)
  invisible(NULL)
}
