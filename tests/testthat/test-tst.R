test_that("workflow_tst() returns a classed object with the expected structure", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 15, ar1 = 0.4, seed = 1)$series
  result <- workflow_tst(r, report = FALSE, verbose = FALSE)

  expect_s3_class(result, "tst")
  expect_true(all(c("prewhiten", "trend", "trend_summary_table", "fdr") %in% names(result)))
})

test_that("print.tst runs without error and returns invisibly", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, report = FALSE, verbose = FALSE)
  expect_output(print(result), "True Significant Trends")
})

test_that("summary.tst runs without error", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, report = FALSE, verbose = FALSE)
  expect_output(summary(result), "Trend test")
})

test_that("plot.tst errors clearly when which='direction' but no FDR was run", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result <- workflow_tst(r, fdr_method = NULL, report = FALSE, verbose = FALSE)
  expect_error(plot(result), "No FDR correction")
})

test_that("workflow_tst() runs cleanly across every combination of prewhiten x theil_sen x fdr_method (a systematic sweep, not just individual cases)", {
  grid <- expand.grid(
    prewhiten = c(TRUE, FALSE),
    theil_sen = c(TRUE, FALSE),
    fdr_method = list(NULL, "BH", "BKY", "BY")
  )
  r <- sim_trend_stack(nrow = 6, ncol = 6, n_time = 10, seed = 100)$series

  for (i in seq_len(nrow(grid))) {
    prewhiten_i <- grid$prewhiten[[i]]
    theil_sen_i <- grid$theil_sen[[i]]
    fdr_method_i <- grid$fdr_method[[i]]

    result <- workflow_tst(r, prewhiten = prewhiten_i, theil_sen = theil_sen_i,
                  fdr_method = fdr_method_i, report = FALSE, verbose = FALSE)

    label <- sprintf("prewhiten=%s, theil_sen=%s, fdr_method=%s",
                      prewhiten_i, theil_sen_i,
                      if (is.null(fdr_method_i)) "NULL" else fdr_method_i)

    expect_s3_class(result, "tst")
    expect_identical(is.null(result$prewhiten), !prewhiten_i, info = label)
    expect_identical(is.null(result$theil_sen), !theil_sen_i, info = label)
    expect_identical(is.null(result$fdr), is.null(fdr_method_i), info = label)

    # print() must not error for any combination in the grid.
    expect_error(print(result), NA, info = label)
  }
})
