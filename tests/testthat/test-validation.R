test_that("compare_detections computes a correct confusion matrix on a known case", {
  truth <- c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE)
  perfect <- truth
  none <- rep(FALSE, 6)

  tab <- compare_detections(list(perfect = perfect, none = none), ground_truth = truth)

  expect_identical(tab$Method, c("perfect", "none"))

  row_perfect <- tab[tab$Method == "perfect", ]
  expect_identical(row_perfect$TP, 3L)
  expect_identical(row_perfect$FP, 0L)
  expect_identical(row_perfect$TN, 3L)
  expect_identical(row_perfect$FN, 0L)
  expect_equal(row_perfect$Sensitivity, 1)
  expect_equal(row_perfect$Specificity, 1)

  row_none <- tab[tab$Method == "none", ]
  expect_identical(row_none$TP, 0L)
  expect_identical(row_none$FN, 3L)
  expect_equal(row_none$Sensitivity, 0)
  expect_true(is.na(row_none$Precision))  # 0/0 -- no positives predicted at all
  expect_equal(row_none$FDR, 0)
})

test_that("FDP is zero when a replicate contains no discoveries", {
  truth <- c(TRUE, FALSE, FALSE, FALSE)
  detections <- list(
    list(method = c(FALSE, FALSE, FALSE, FALSE)),
    list(method = c(FALSE, TRUE, FALSE, FALSE)),
    list(method = c(TRUE, FALSE, FALSE, FALSE))
  )
  result <- compare_detections(
    detections, truth, replicates = TRUE,
    metrics = "fdr", verbose = FALSE
  )
  expect_equal(result$FDR_mean, 1 / 3)
  expect_equal(result$FDR_sd, stats::sd(c(0, 1, 0)))
})

test_that("compare_detections works with SpatRaster inputs (e.g. sim_trend_stack)", {
  sim <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 1)
  trend <- trend_test(sim$series, report = FALSE, verbose = FALSE)
  sig <- trend$stats$p <= 0.05  # SpatRaster of TRUE/FALSE

  tab <- compare_detections(list(raw = sig), ground_truth = sim$true_slope)
  expect_equal(nrow(tab), 1)
  expect_true(all(c("TP", "FP", "TN", "FN", "Sensitivity", "FDR") %in% names(tab)))
})

test_that("compare_detections respects the metrics argument", {
  truth <- c(TRUE, FALSE, TRUE, FALSE)
  sig <- c(TRUE, FALSE, FALSE, FALSE)
  tab <- compare_detections(list(m = sig), ground_truth = truth, metrics = c("sensitivity", "f1"))
  expect_true(all(c("Sensitivity", "F1") %in% names(tab)))
  expect_false(any(c("Specificity", "Precision", "FPR", "FDR") %in% names(tab)))
})

test_that("compare_detections errors on an unnamed list", {
  expect_error(compare_detections(list(TRUE, FALSE), ground_truth = c(TRUE, FALSE)), "named list")
})

test_that("compare_detections errors on an unknown metric", {
  expect_error(
    compare_detections(list(m = c(TRUE, FALSE)), ground_truth = c(TRUE, FALSE), metrics = "recall"),
    "Unknown metric"
  )
})

test_that("compare_detections computes Accuracy and MCC correctly on a known case", {
  # TP=3, FP=1, TN=2, FN=0 -> accuracy = 5/6; MCC has a known closed form here
  sig   <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
  truth <- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
  tab <- compare_detections(list(m = sig), ground_truth = truth)

  expect_identical(tab$TP, 3L)
  expect_identical(tab$FP, 1L)
  expect_identical(tab$TN, 2L)
  expect_identical(tab$FN, 0L)
  expect_equal(tab$Accuracy, 5 / 6)

  expected_mcc <- (3 * 2 - 1 * 0) / sqrt((3 + 1) * (3 + 0) * (2 + 1) * (2 + 0))
  expect_equal(tab$MCC, expected_mcc)
})

test_that("compare_detections returns NA for MCC when a confusion-matrix margin is zero", {
  # everything predicted positive, and truth is all positive -> TN=FP=0
  sig <- c(TRUE, TRUE, TRUE)
  truth <- c(TRUE, TRUE, TRUE)
  tab <- compare_detections(list(m = sig), ground_truth = truth)
  expect_true(is.na(tab$MCC))
})

test_that("compare_detections errors on mismatched lengths", {
  expect_error(
    compare_detections(list(m = c(TRUE, FALSE, TRUE)), ground_truth = c(TRUE, FALSE)),
    "length"
  )
})

test_that("compare_detections(replicates = TRUE) aggregates correctly across replicates", {
  truth <- c(TRUE, TRUE, FALSE, FALSE)
  detections_list <- list(
    list(m = c(TRUE, FALSE, FALSE, FALSE)),  # sens = 0.5
    list(m = c(TRUE, TRUE, FALSE, FALSE))    # sens = 1.0
  )

  summary_tab <- compare_detections(detections_list, ground_truth = truth,
                                     replicates = TRUE)

  expect_s3_class(summary_tab, "compare_detections_replicates")
  expect_identical(summary_tab$Method, "m")
  expect_equal(summary_tab$n_replicates, 2)
  expect_equal(summary_tab$Sensitivity_mean, 0.75)
  expect_gt(summary_tab$Sensitivity_sd, 0)
})

test_that("compare_detections(replicates = TRUE) preserves method order and handles several methods", {
  truth <- c(TRUE, FALSE, TRUE, FALSE)
  one_rep <- list(B = c(TRUE, FALSE, TRUE, FALSE), A = c(FALSE, FALSE, TRUE, FALSE))
  summary_tab <- compare_detections(list(one_rep, one_rep), ground_truth = truth,
                                     replicates = TRUE)
  expect_identical(summary_tab$Method, c("B", "A"))
})

test_that("compare_detections(replicates = TRUE) accepts a ground_truth list that varies per replicate", {
  detections_list <- list(
    list(m = c(TRUE, FALSE, FALSE, FALSE)),
    list(m = c(TRUE, TRUE, FALSE, FALSE))
  )
  truths_list <- list(c(TRUE, TRUE, FALSE, FALSE), c(TRUE, TRUE, FALSE, FALSE))
  summary_tab <- compare_detections(detections_list, ground_truth = truths_list,
                                     replicates = TRUE)
  expect_equal(summary_tab$n_replicates, 2)
})

test_that("compare_detections(replicates = TRUE) errors on invalid 'detections'", {
  truth <- c(TRUE, FALSE, TRUE, FALSE)
  # A single named list of vectors (the replicates = FALSE shape), not
  # a list of such lists -- the most likely real mistake to make here.
  expect_error(
    compare_detections(list(m = c(TRUE, FALSE, TRUE, FALSE)),
                        ground_truth = truth, replicates = TRUE),
    "list of.*named lists"
  )
  expect_error(
    compare_detections(list(), ground_truth = truth, replicates = TRUE),
    "list of.*named lists"
  )
})

test_that("compare_detections(replicates = TRUE) errors when ground_truth list length does not match detections", {
  detections_list <- list(list(m = c(TRUE, FALSE)), list(m = c(TRUE, FALSE)))
  truths_list <- list(c(TRUE, FALSE))  # only 1, but detections has 2
  expect_error(
    compare_detections(detections_list, ground_truth = truths_list,
                        replicates = TRUE),
    "same length"
  )
})

test_that("plot_detection_comparison runs without error on a compare_detections() result", {
  truth <- c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE)
  tab <- compare_detections(list(MK = truth, CMK = c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE)), ground_truth = truth)
  expect_silent(plot_detection_comparison(tab))
})

test_that("plot_detection_comparison errors when no requested metric is present", {
  tab <- data.frame(Method = "m", TP = 1, FP = 0, TN = 1, FN = 0)
  expect_error(plot_detection_comparison(tab), "None of the requested")
})

test_that("plot()/print()/summary() dispatch correctly for a 'compare_detections_replicates' result", {
  truth <- c(TRUE, TRUE, FALSE, FALSE)
  detections_list <- list(
    list(m = c(TRUE, FALSE, FALSE, FALSE)),
    list(m = c(TRUE, TRUE, FALSE, FALSE))
  )
  result <- compare_detections(detections_list, ground_truth = truth,
                                replicates = TRUE)
  expect_error(print(result), NA)
  expect_error(suppressMessages(summary(result)), NA)
  expect_error(plot(result), NA)
})

test_that("plot_detection_comparison can write a PNG via path", {
  truth <- c(TRUE, TRUE, FALSE, FALSE, TRUE, FALSE)
  tab <- compare_detections(list(m = truth), ground_truth = truth)
  path <- tempfile(fileext = ".png")
  plot_detection_comparison(tab, path = path)
  expect_true(file.exists(path))
  unlink(path)
})

test_that("compare_detections(replicates = TRUE) errors clearly when one specific replicate's own inner list is unnamed, not just when the whole detections list is malformed", {
  truth <- c(TRUE, FALSE, TRUE, FALSE)
  detections_list <- list(
    list(m = c(TRUE, FALSE, TRUE, FALSE)),  # replicate 1: correctly named
    list(c(TRUE, FALSE, TRUE, FALSE))       # replicate 2: unnamed -- the bug
  )
  expect_error(
    compare_detections(detections_list, ground_truth = truth,
                        replicates = TRUE),
    "detections\\[\\[2\\]\\].*named list"
  )
})

test_that("summary.compare_detections_replicates() returns the best-scoring method per metric, computed on the aggregated table's own _mean columns", {
  truth <- c(TRUE, TRUE, FALSE, FALSE)
  detections_list <- list(
    list(a = c(TRUE, FALSE, FALSE, FALSE),   # weaker: sens = 0.5
         b = c(TRUE, TRUE, FALSE, FALSE)),   # stronger: sens = 1.0
    list(a = c(TRUE, FALSE, FALSE, FALSE),
         b = c(TRUE, TRUE, FALSE, FALSE))
  )
  result <- compare_detections(detections_list, ground_truth = truth,
                                replicates = TRUE)

  winners <- summary(result)
  expect_true("best_method" %in% names(winners))
  # Method "b" should win on every metric where it strictly dominates "a".
  sens_row <- winners[winners$metric == "Sensitivity", ]
  expect_equal(sens_row$best_method, "b")
})

test_that("summary.compare_detections_replicates() returns NA_character_ for a metric that is NA across every replicate (a method that never predicts positive, so TP + FP = 0 makes Precision undefined every time)", {
  truth <- c(TRUE, TRUE, FALSE, FALSE)
  detections_list <- list(
    list(never_positive = c(FALSE, FALSE, FALSE, FALSE)),
    list(never_positive = c(FALSE, FALSE, FALSE, FALSE))
  )
  result <- compare_detections(detections_list, ground_truth = truth,
                                replicates = TRUE)

  winners <- summary(result)
  precision_row <- winners[winners$metric == "Precision", ]
  expect_true(is.na(precision_row$best_method))
})

test_that("the aggregated table itself reports a clean NA (not NaN) for '_mean' when a metric is NA across every replicate -- a real inconsistency an external audit found: sd(na.rm=TRUE) already gave clean NA for the same case, but mean(na.rm=TRUE) gave NaN", {
  truth <- c(TRUE, TRUE, FALSE, FALSE)
  detections_list <- list(
    list(never_positive = c(FALSE, FALSE, FALSE, FALSE)),
    list(never_positive = c(FALSE, FALSE, FALSE, FALSE))
  )
  result <- compare_detections(detections_list, ground_truth = truth,
                                replicates = TRUE)

  precision_mean <- result$Precision_mean[result$Method == "never_positive"]
  expect_true(is.na(precision_mean))
  expect_false(is.nan(precision_mean))
})

test_that("compare_detections(metrics = 'fwer') errors clearly when replicates = FALSE", {
  truth <- c(TRUE, FALSE)
  det <- list(a = c(TRUE, FALSE))
  expect_error(
    compare_detections(det, ground_truth = truth, metrics = "fwer"),
    "only has a meaning across many replicates"
  )
})

test_that("compare_detections(replicates = TRUE, metrics = 'fwer') computes the correct proportion of replicates with at least one false positive", {
  truth <- c(TRUE, TRUE, FALSE, FALSE)
  # Method "noisy": FP > 0 in 2 of 3 replicates -> FWER = 2/3.
  # Method "clean": FP == 0 in every replicate -> FWER = 0.
  detections_list <- list(
    list(noisy = c(TRUE, TRUE, TRUE, FALSE),   # 1 FP
         clean = c(TRUE, TRUE, FALSE, FALSE)),  # 0 FP
    list(noisy = c(TRUE, TRUE, TRUE, FALSE),   # 1 FP
         clean = c(TRUE, TRUE, FALSE, FALSE)),  # 0 FP
    list(noisy = c(TRUE, TRUE, FALSE, FALSE),  # 0 FP (this replicate)
         clean = c(TRUE, TRUE, FALSE, FALSE))   # 0 FP
  )
  result <- compare_detections(detections_list, ground_truth = truth,
                                replicates = TRUE,
                                metrics = c("sensitivity", "fwer"))
  expect_true("FWER" %in% names(result))
  noisy_row <- result[result$Method == "noisy", ]
  clean_row <- result[result$Method == "clean", ]
  expect_equal(noisy_row$FWER, 2 / 3)
  expect_equal(clean_row$FWER, 0)

  # summary() must pick the LOWER FWER as "best" (which.min, not
  # which.max, unlike every other metric here) -- "clean" should win.
  winners <- summary(result)
  fwer_row <- winners[winners$metric == "FWER", ]
  expect_equal(fwer_row$best_method, "clean")
})

test_that("summary.compare_detections() correctly picks the method with the LOWER FPR/FDR as best (which.min, not which.max) -- a real, previously-existing bug fixed here, not a hypothetical one", {
  truth <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
  # Method "sloppy": 2 false positives (FP = 2) among 3 true negatives.
  # Method "careful": 0 false positives.
  detections <- list(
    sloppy  = c(TRUE, TRUE, TRUE, TRUE, FALSE),   # 2 FP (positions 3,4)
    careful = c(TRUE, TRUE, FALSE, FALSE, FALSE)  # 0 FP
  )
  result <- compare_detections(detections, ground_truth = truth,
                                metrics = c("sensitivity", "fpr", "fdr"))

  winners <- summary(result)
  fpr_row <- winners[winners$metric == "FPR", ]
  fdr_row <- winners[winners$metric == "FDR", ]
  # "careful" has the lower (better) FPR and FDR -- before this fix,
  # which.max() would have wrongly named "sloppy" (the higher FPR/FDR)
  # as the "best" method for both.
  expect_equal(fpr_row$best_method, "careful")
  expect_equal(fdr_row$best_method, "careful")
  # Sensitivity is unaffected by this fix (higher is genuinely better,
  # which.max() was always correct for it) -- both methods tie here
  # (both catch both true positives), so either coming back is fine;
  # this just confirms the fix did not somehow break the unaffected case.
  sens_row <- winners[winners$metric == "Sensitivity", ]
  expect_true(sens_row$best_method %in% c("sloppy", "careful"))
})

test_that("summary.compare_detections() correctly picks the method with the LOWER TypeI/TypeII as best (which.min, not which.max) -- a real bug that shipped alongside the FPR/FDR fix above but was missed for these three error-rate metrics until an external audit caught it", {
  truth <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
  # Method "misses_more": only catches 1 of 4 true signals (TypeII =
  # 3/4 = 0.75, a high miss rate) and has 1 false positive (TypeI =
  # FPR = 1/2 = 0.5).
  # Method "misses_less": catches 3 of 4 true signals (TypeII = 1/4 =
  # 0.25) and has 0 false positives (TypeI = 0).
  detections <- list(
    misses_more = c(TRUE, FALSE, FALSE, FALSE, TRUE, FALSE),
    misses_less = c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
  )
  result <- compare_detections(detections, ground_truth = truth,
                                metrics = c("type_i", "type_ii"))

  winners <- summary(result)
  type_i_row <- winners[winners$metric == "TypeI", ]
  type_ii_row <- winners[winners$metric == "TypeII", ]
  # "misses_less" has the lower (better) TypeI and TypeII -- before
  # this fix, which.max() would have wrongly named "misses_more" (the
  # higher, worse error rate) as the "best" method for both.
  expect_equal(type_i_row$best_method, "misses_less")
  expect_equal(type_ii_row$best_method, "misses_less")
})

test_that("summary.compare_detections_replicates() correctly picks the method with the LOWER FPR_mean/FDR_mean as best (which.min, not which.max) -- same real bug as the single-run case above, fixed in both places", {
  truth <- c(TRUE, TRUE, FALSE, FALSE, FALSE)
  detections_list <- list(
    list(sloppy = c(TRUE, TRUE, TRUE, TRUE, FALSE),
         careful = c(TRUE, TRUE, FALSE, FALSE, FALSE)),
    list(sloppy = c(TRUE, TRUE, TRUE, TRUE, FALSE),
         careful = c(TRUE, TRUE, FALSE, FALSE, FALSE))
  )
  result <- compare_detections(detections_list, ground_truth = truth,
                                replicates = TRUE,
                                metrics = c("sensitivity", "fpr", "fdr"))

  winners <- summary(result)
  fpr_row <- winners[winners$metric == "FPR", ]
  fdr_row <- winners[winners$metric == "FDR", ]
  expect_equal(fpr_row$best_method, "careful")
  expect_equal(fdr_row$best_method, "careful")
})

test_that("compare_detections Usage exposes every metric without changing the omitted default", {
  displayed <- eval(formals(compare_detections)$metrics)
  expect_setequal(
    displayed,
    c("sensitivity", "specificity", "precision", "accuracy",
      "f1", "mcc", "fpr", "fdr", "fwer", "type_i", "type_ii",
      "type_iii", "field_power", "global_power", "within_image_power",
      "directional_power")
  )

  truth <- c(TRUE, FALSE, TRUE, FALSE)
  result <- compare_detections(list(method = truth), truth)
  expect_false("FWER" %in% names(result))
  expect_true(all(c("Sensitivity", "Specificity", "Precision",
                    "Accuracy", "F1", "MCC", "FPR", "FDR") %in%
                  names(result)))
})

test_that("extended trend metrics have analytically exact values", {
  truth <- c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE)
  truth_direction <- c(1, 1, -1, 0, 0, 0)
  detection <- list(
    significant = c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE),
    direction = c(1, -1, 0, 1, 0, 0))
  result <- compare_detections(
    list(method = detection), truth,
    metrics = c("type_i", "type_ii", "type_iii", "field_power",
                "global_power", "within_image_power",
                "directional_power"),
    truth_direction = truth_direction)
  expect_equal(result$TypeI, 1 / 3)
  expect_equal(result$TypeII, 1 / 3)
  expect_equal(result$TypeIII, 1 / 2)
  expect_equal(result$FieldPower, 1)
  expect_equal(result$GlobalPower, 1)
  expect_equal(result$WithinImagePower, 2 / 3)
  expect_equal(result$DirectionalPower, 1 / 3)
})

test_that("field and true-global power remain distinct", {
  truth <- c(TRUE, FALSE, FALSE, FALSE)
  false_only <- c(FALSE, TRUE, FALSE, FALSE)
  result <- compare_detections(
    list(method = false_only), truth,
    metrics = c("field_power", "global_power"),
    verbose = FALSE
  )
  expect_equal(result$FieldPower, 1)
  expect_equal(result$GlobalPower, 0)
})

test_that("field power is supported by detection comparison plots", {
  truth <- c(TRUE, FALSE, FALSE, FALSE)
  result <- compare_detections(
    list(false_only = c(FALSE, TRUE, FALSE, FALSE)), truth,
    metrics = c("field_power", "global_power"), verbose = FALSE
  )
  expect_error(
    sptrends:::plot_detection_comparison(
      result, metrics = c("FieldPower", "GlobalPower")
    ),
    NA
  )
})

test_that("one evaluation mask is applied identically to every method", {
  truth <- c(TRUE, FALSE, FALSE, FALSE)
  detections <- list(
    MK = c(TRUE, TRUE, FALSE, FALSE),
    CMK = c(TRUE, TRUE, FALSE, FALSE))
  result <- compare_detections(
    detections, truth, metrics = "type_i",
    evaluation_mask = c(TRUE, FALSE, TRUE, TRUE))
  expect_equal(result$TypeI, c(0, 0))
})

test_that("method-specific evaluation masks are rejected", {
  truth <- c(TRUE, FALSE, FALSE, FALSE)
  detections <- list(
    MK = c(TRUE, TRUE, FALSE, FALSE),
    CMK = c(TRUE, TRUE, FALSE, FALSE))
  masks <- list(
    MK = rep(TRUE, 4),
    CMK = c(TRUE, FALSE, TRUE, TRUE))
  expect_error(
    compare_detections(
      detections, truth, metrics = "type_i", evaluation_mask = masks),
    "one common")
})

test_that("replicated global and within-image power aggregate correctly", {
  truth <- c(TRUE, TRUE, FALSE, FALSE)
  detections <- list(
    list(MK = c(FALSE, FALSE, FALSE, FALSE)),
    list(MK = c(TRUE, FALSE, FALSE, FALSE)),
    list(MK = c(TRUE, TRUE, FALSE, FALSE)))
  result <- compare_detections(
    detections, truth, replicates = TRUE,
    metrics = c("global_power", "within_image_power"))
  expect_equal(result$GlobalPower_mean, 2 / 3)
  expect_equal(result$WithinImagePower_mean, 0.5)
})

test_that("replicated field power counts any rejection", {
  truth <- c(TRUE, TRUE, FALSE, FALSE)
  detections <- list(
    list(MK = c(FALSE, FALSE, FALSE, FALSE)),
    list(MK = c(FALSE, FALSE, TRUE, FALSE)),
    list(MK = c(TRUE, FALSE, FALSE, FALSE))
  )
  result <- compare_detections(
    detections, truth, replicates = TRUE,
    metrics = c("field_power", "global_power"), verbose = FALSE
  )
  expect_equal(result$FieldPower_mean, 2 / 3)
  expect_equal(result$GlobalPower_mean, 1 / 3)
})

test_that("benchmark_methods scores detections and records each replicate", {
  scenarios <- list(
    signal = list(nrow = 4, ncol = 4, n_time = 4,
                  trend_shape = "square", signal_size = 2,
                  trend_fraction = 1, constant_block = FALSE,
                  spatial_model = "independent"))
  methods <- list(
    oracle = function(series, simulation) {
      list(significant = simulation$true_signal,
           direction = simulation$true_direction)
    })
  result <- benchmark_methods(
    scenarios, methods, n_replicates = 3, seed = 4,
    metrics = c("type_i", "type_ii", "directional_power"))
  expect_s3_class(result, "sptrends_benchmark")
  expect_equal(nrow(result), 3)
  expect_true(all(result$TypeI == 0))
  expect_true(all(result$TypeII == 0))
  expect_true(all(result$DirectionalPower == 1))
  expect_true(all(result$Elapsed >= 0))
})

test_that("benchmark_methods scores slope estimators against true slopes", {
  scenarios <- list(
    signal = list(nrow = 4, ncol = 4, n_time = 4,
                  trend_fraction = 1, constant_block = FALSE,
                  spatial_model = "independent"))
  methods <- list(
    oracle = function(series, simulation) simulation$true_slope)
  result <- benchmark_methods(
    scenarios, methods, n_replicates = 2, stage = "slope", seed = 5)
  expect_true(all(result$Bias == 0))
  expect_true(all(result$MAE == 0))
  expect_true(all(result$RMSE == 0))
  expect_true(all(result$DirectionError == 0))
})

test_that("benchmark architecture exposes only implemented stages", {
  expect_identical(
    eval(formals(benchmark_methods)$stage),
    c("trend_test", "prewhitening", "slope", "fdr", "fwer", "custom",
      "detection"))
})

test_that("benchmark methods receive the identical simulated realisation", {
  observed <- list()
  methods <- list(
    first = function(series, simulation) {
      observed$first <<- terra::values(series, mat = TRUE)
      list(significant = simulation$true_signal,
           direction = simulation$true_direction)
    },
    second = function(series, simulation) {
      observed$second <<- terra::values(series, mat = TRUE)
      list(significant = simulation$true_signal,
           direction = simulation$true_direction)
    })
  scenarios <- list(
    one = list(nrow = 3, ncol = 3, n_time = 4,
               constant_block = FALSE, spatial_model = "independent"))
  benchmark_methods(
    scenarios, methods, n_replicates = 1, seed = 8,
    metrics = "type_i")
  expect_identical(observed$first, observed$second)
})

test_that("prepare runs once and shares one derived input across methods", {
  prepare_calls <- 0L
  observed <- list()
  prepare <- function(series, simulation) {
    prepare_calls <<- prepare_calls + 1L
    terra::values(series[[1]], mat = FALSE)
  }
  methods <- list(
    first = function(input, simulation) {
      observed$first <<- input
      list(significant = simulation$true_signal,
           direction = simulation$true_direction)
    },
    second = function(input, simulation) {
      observed$second <<- input
      list(significant = simulation$true_signal,
           direction = simulation$true_direction)
    })
  scenarios <- list(
    one = list(nrow = 3, ncol = 3, n_time = 4,
               constant_block = FALSE, spatial_model = "independent"))
  benchmark_methods(
    scenarios, methods, n_replicates = 2, prepare = prepare,
    seed = 8, metrics = "type_i")
  expect_equal(prepare_calls, 2L)
  expect_identical(observed$first, observed$second)
})

test_that("prewhitening stage accepts every transformed-series method", {
  scenarios <- list(
    ar1 = list(nrow = 3, ncol = 3, n_time = 8, ar1 = 0.5,
               constant_block = FALSE, spatial_model = "independent"))
  methods <- list(
    unchanged = function(series, simulation) series,
    packaged = function(series, simulation) list(series = series))
  result <- benchmark_methods(
    scenarios, methods, n_replicates = 1,
    stage = "prewhitening", seed = 9)
  expect_equal(result$Method, c("unchanged", "packaged"))
  expect_true(all(c("ResidualACF1", "SlopeRMSE", "OutputLength") %in%
                  names(result)))
})

test_that("FDR and FWER share detection scoring", {
  scenarios <- list(
    null = list(nrow = 3, ncol = 3, n_time = 4,
                trend_fraction = 0, constant_block = FALSE,
                spatial_model = "independent"))
  methods <- list(
    empty = function(series, simulation) {
      list(significant = simulation$true_signal,
           direction = simulation$true_direction)
    })
  for (stage in c("fdr", "fwer")) {
    result <- benchmark_methods(
      scenarios, methods, n_replicates = 1, stage = stage,
      metrics = c("fdr", "fwer"), seed = 10)
    expect_equal(result$FDR, 0)
    expect_equal(result$FDR, result$FalseDiscoveryProportion)
    expect_equal(result$FalseDiscoveryProportion, 0)
    expect_equal(result$AnyFalsePositive, 0)
  }
})

test_that("benchmark_summary estimates FDR and FWER across replicates", {
  x <- structure(
    data.frame(
      Scenario = rep("null", 3), Replicate = 1:3, Seed = 1:3,
      Elapsed = 0, Method = rep("method", 3), FP = c(1, 1, 0),
      FDR = c(1, 1, 0),
      FalseDiscoveryProportion = c(1, 1, 0),
      AnyFalsePositive = c(1, 1, 0)),
    class = c("sptrends_benchmark", "data.frame"))
  result <- benchmark_summary(x)
  expect_equal(result$FDR_mean, result$EmpiricalFDR)
  expect_equal(result$EmpiricalFDR, 2 / 3)
  expect_equal(result$EmpiricalFWER, 2 / 3)
  expect_equal(result$n_replicates, 3)
})

test_that("varying realised truth composition does not split replicates", {
  x <- structure(
    data.frame(
      Scenario = rep("variable_truth", 3), Replicate = 1:3,
      Seed = 1:3, Elapsed = 0, trend_fraction = 0.5,
      DomainCells = 100, TrueNulls = c(40, 50, 60),
      Pi0 = c(0.4, 0.5, 0.6), SignalProportion = c(0.6, 0.5, 0.4),
      Method = "method", TypeI = c(0.04, 0.05, 0.06)),
    scenario_fields = "trend_fraction",
    truth_fields = c(
      "DomainCells", "TrueNulls", "Pi0", "SignalProportion"),
    class = c("sptrends_benchmark", "sptrends", "data.frame"))

  result <- benchmark_summary(x)
  expect_equal(nrow(result), 1L)
  expect_equal(result$n_replicates, 3L)
  expect_equal(result$Pi0_mean, 0.5)
})

test_that("custom evaluators permit arbitrary future method outputs", {
  scenarios <- list(
    one = list(nrow = 2, ncol = 2, n_time = 3,
               constant_block = FALSE, spatial_model = "independent"))
  methods <- list(arbitrary = function(series, simulation) list(value = 7))
  evaluator <- function(outputs, simulation) {
    data.frame(Method = names(outputs), CustomScore = outputs[[1]]$value)
  }
  result <- benchmark_methods(
    scenarios, methods, n_replicates = 1, stage = "custom",
    evaluator = evaluator, seed = 11)
  expect_equal(result$CustomScore, 7)
})

test_that("benchmark evaluators must score every method exactly once", {
  scenarios <- list(
    one = list(nrow = 2, ncol = 2, n_time = 3,
               constant_block = FALSE, spatial_model = "independent"))
  methods <- list(
    first = function(series, simulation) 1,
    second = function(series, simulation) 2)
  incomplete <- function(outputs, simulation) {
    data.frame(Method = "first", Score = 1)
  }
  duplicated <- function(outputs, simulation) {
    data.frame(Method = c("first", "first"), Score = c(1, 2))
  }

  expect_error(
    benchmark_methods(
      scenarios, methods, n_replicates = 1, stage = "custom",
      evaluator = incomplete),
    "exactly one row")
  expect_error(
    benchmark_methods(
      scenarios, methods, n_replicates = 1, stage = "custom",
      evaluator = duplicated),
    "exactly one row")
})

test_that("benchmark_methods validates simulator seeds", {
  expect_error(
    benchmark_methods(
      list(one = list()), list(method = function(x, simulation) x),
      seed = NA_real_),
    "finite numeric")
})

test_that("scenario and score columns cannot silently collide", {
  scenarios <- list(
    one = list(nrow = 2, ncol = 2, n_time = 3, TypeI = 7,
               constant_block = FALSE, spatial_model = "independent"))
  simulator <- function(nrow, ncol, n_time, TypeI, constant_block,
                        spatial_model, seed) {
    sim_trend_stack(
      nrow = nrow, ncol = ncol, n_time = n_time,
      constant_block = constant_block, spatial_model = spatial_model,
      seed = seed)
  }
  methods <- list(oracle = function(series, simulation) {
    list(significant = simulation$true_signal,
         direction = simulation$true_direction)
  })

  expect_error(
    benchmark_methods(
      scenarios, methods, n_replicates = 1, simulator = simulator,
      metrics = "type_i"),
    "columns overlap")
})
test_that("detection comparisons report all three timing indicators", {
  output <- capture.output(
    compare_detections(
      detections = list(method = c(TRUE, FALSE, TRUE)),
      ground_truth = c(TRUE, FALSE, FALSE), verbose = TRUE
    ),
    type = "message"
  )
  text <- paste(output, collapse = " ")
  expect_match(text, "progress:", fixed = TRUE)
  expect_match(text, "remaining:", fixed = TRUE)
  expect_match(text, "elapsed:", fixed = TRUE)
})
