test_that("workflow_tst() and workflow_rta() results share the 'sptrends' superclass", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 1)$series
  result_tst <- workflow_tst(r, report = FALSE, verbose = FALSE)
  result_rta <- workflow_rta(r, report = FALSE, verbose = FALSE)
  expect_identical(class(result_tst), c("tst", "sptrends"))
  expect_identical(class(result_rta), c("rta", "sptrends"))
  expect_s3_class(result_tst, "sptrends")
  expect_s3_class(result_rta, "sptrends")
})

test_that("workflow_tst() and workflow_rta() record per-step timing", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 2)$series
  result_tst <- workflow_tst(r, report = FALSE, verbose = FALSE)
  result_rta <- workflow_rta(r, report = FALSE, verbose = FALSE)

  expect_true(all(c("prewhiten", "CMK", "theil_sen", "fdr") %in% names(result_tst$timing)))
  expect_true(all(vapply(result_tst$timing, function(t) is.numeric(t) && t >= 0, logical(1))))

  expect_true(all(c("theil_sen", "CMK", "fdr") %in% names(result_rta$timing)))
  expect_true(all(vapply(result_rta$timing, function(t) is.numeric(t) && t >= 0, logical(1))))
})

test_that("workflow_tst()'s timing omits steps that did not run", {
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 12, seed = 3)$series
  result <- workflow_tst(r, prewhiten = FALSE, theil_sen = FALSE, fdr_method = NULL,
                report = FALSE, verbose = FALSE)
  expect_null(result$timing$prewhiten)
  expect_null(result$timing$theil_sen)
  expect_null(result$timing$fdr)
  expect_false(is.null(result$timing$CMK))
})

test_that(".mask_and_smooth_slope() never paints a non-significant cell, even with significant neighbours (na.policy fix)", {
  # A 5x5 raster where only the centre cell (cell 13) is "significant";
  # all its 8 queen neighbours are not. Before the na.policy = "omit"
  # fix, terra::focal()'s own default ("all") would assign the centre's
  # neighbours a smoothed value too, since they each have at least one
  # valid (significant) neighbour -- painting cells that were never
  # actually significant.
  r <- terra::rast(nrows = 5, ncols = 5, vals = seq_len(25))
  reject <- rep(FALSE, 25)
  reject[13] <- TRUE  # centre cell only

  result <- sptrends:::.mask_and_smooth_slope(r, reject, smooth = TRUE)
  vals <- terra::values(result, mat = FALSE)

  # the centre cell (the only significant one) must have a real value
  expect_false(is.na(vals[13]))
  # every other cell -- including the centre's own queen neighbours,
  # which the bug would have "leaked" a value into -- must stay NA
  expect_true(all(is.na(vals[-13])))
})

test_that(".mask_and_smooth_slope() with smooth = FALSE just masks, no focal() involved", {
  r <- terra::rast(nrows = 3, ncols = 3, vals = 1:9)
  reject <- c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, TRUE)
  result <- sptrends:::.mask_and_smooth_slope(r, reject, smooth = FALSE)
  vals <- terra::values(result, mat = FALSE)
  expect_equal(is.na(vals), !reject)
  expect_equal(vals[reject], (1:9)[reject])
})

test_that("workflow_tst(n_cores = 2) builds a single shared cluster and gives identical results to n_cores = 1", {
  .skip_unless_parallel_tests()
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 70)$series

  result_seq <- workflow_tst(r, prewhiten = FALSE, n_cores = 1,
                              report = FALSE, verbose = FALSE)
  result_par <- workflow_tst(r, prewhiten = FALSE, n_cores = 2,
                              report = FALSE, verbose = FALSE)

  expect_equal(terra::values(result_seq$trend, mat = FALSE),
               terra::values(result_par$trend, mat = FALSE))
  expect_equal(terra::values(result_seq$theil_sen, mat = FALSE),
               terra::values(result_par$theil_sen, mat = FALSE))
})

test_that("workflow_tst(n_cores = 2)'s top-level cluster takes precedence over a conflicting n_cores set inside cmk_args/theil_sen_args, rather than each step building its own on top of the shared one", {
  .skip_unless_parallel_tests()
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 71)$series

  # cmk_args/theil_sen_args each ask for n_cores = 5 -- if the shared
  # top-level cluster were NOT taking precedence, each step would try
  # to build its own 5-worker cluster instead of reusing the 2-worker
  # shared one, which this test cannot observe directly without deeper
  # mocking -- but it can at least confirm the call succeeds and gives
  # a genuine result rather than erroring or silently returning nothing,
  # which a broken precedence (e.g. two conflicting clusters fighting
  # over the same connection) would be likely to do.
  result <- workflow_tst(r, prewhiten = FALSE, n_cores = 2,
                          cmk_args = list(n_cores = 5),
                          theil_sen_args = list(n_cores = 5),
                          report = FALSE, verbose = FALSE)
  expect_false(is.null(result$trend))
  expect_false(is.null(result$theil_sen))
})

test_that("workflow_rta(n_cores = 2) builds a single shared cluster and gives identical results to n_cores = 1", {
  .skip_unless_parallel_tests()
  r <- sim_trend_stack(nrow = 8, ncol = 8, n_time = 10, seed = 72)$series

  result_seq <- workflow_rta(r, n_cores = 1, report = FALSE, verbose = FALSE)
  result_par <- workflow_rta(r, n_cores = 2, report = FALSE, verbose = FALSE)

  expect_equal(terra::values(result_seq$trend, mat = FALSE),
               terra::values(result_par$trend, mat = FALSE))
  expect_equal(terra::values(result_seq$theil_sen, mat = FALSE),
               terra::values(result_par$theil_sen, mat = FALSE))
})

test_that("published workflows forward a configurable CMK window", {
  r <- sim_trend_stack(nrow = 7, ncol = 7, n_time = 8, seed = 908)$series
  direct <- trend_test(r, window_size = 5L, report = FALSE,
                       verbose = FALSE)$stats

  tst <- workflow_tst(
    r, prewhiten = FALSE, theil_sen = FALSE, fdr_method = NULL,
    cmk_args = list(window_size = 5L), report = FALSE, verbose = FALSE
  )
  rta <- workflow_rta(
    r, cmk_args = list(window_size = 5L),
    report = FALSE, verbose = FALSE
  )

  expect_equal(terra::values(tst$trend, mat = TRUE),
               terra::values(direct, mat = TRUE), tolerance = 0)
  expect_equal(terra::values(rta$trend, mat = TRUE),
               terra::values(direct, mat = TRUE), tolerance = 0)
})
