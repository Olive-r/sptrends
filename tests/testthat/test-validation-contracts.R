test_that("public numerical controls reject malformed scalar values", {
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 6, seed = 901)$series

  expect_error(trend_test(r, n_cores = 1.5, report = FALSE,
                          verbose = FALSE), "positive integer")
  expect_error(trend_test(r, alpha = c(0.05, NA), report = FALSE,
                          verbose = FALSE), "strictly between 0 and 1")
  expect_error(slope_estimator(r, max_pairs = 0, report = FALSE,
                               verbose = FALSE), "positive integer or Inf")
  expect_error(slope_estimator(r, seed = Inf, report = FALSE,
                               verbose = FALSE), "finite numeric")
  expect_error(prewhiten(r, eps = 0, report = FALSE, verbose = FALSE),
               "finite positive")
  expect_error(prewhiten(r, itmax = 2.5, report = FALSE, verbose = FALSE),
               "positive integer")
  expect_error(prewhiten(r, dw_low = 3, dw_high = 2,
                         report = FALSE, verbose = FALSE),
               "0 <= dw_low < dw_high <= 4")
  expect_error(inspect_ts_cell(r, conf_level = 1),
               "strictly between 0 and 1")
})

test_that("workflow stage lists cannot replace workflow-managed arguments", {
  r <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 6, seed = 902)$series

  expect_error(
    workflow_trends(r, trend_args = list(n_cores = 2),
                    report = FALSE, verbose = FALSE),
    "cannot override workflow-managed argument"
  )
  expect_error(
    workflow_tst(r, cmk_args = list(x = r),
                 report = FALSE, verbose = FALSE),
    "cannot override workflow-managed argument"
  )
  expect_error(
    workflow_rta(r, theil_sen_args = list(shared_cluster = NULL),
                 report = FALSE, verbose = FALSE),
    "cannot override workflow-managed argument"
  )
  expect_error(
    workflow_trends(r, slope_args = list(1),
                    report = FALSE, verbose = FALSE),
    "fully named list"
  )
})
