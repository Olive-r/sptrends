test_that("directional validation rejects inconsistent structures", {
  truth <- c(TRUE, FALSE, TRUE)
  expect_error(
    compare_detections(
      list(m = list(direction = c(1, 0, -1))), truth
    ),
    "no 'significant'"
  )
  expect_error(
    compare_detections(
      list(m = truth), truth, truth_direction = c(1, -1)
    ),
    "truth_direction"
  )
  expect_error(
    compare_detections(
      list(m = truth), truth,
      directions = list(m = c(1, -1))
    ),
    "Direction"
  )
  expect_error(
    compare_detections(
      list(m = truth), truth,
      evaluation_mask = c(TRUE, FALSE)
    ),
    "equal lengths"
  )
})

test_that("replicated validation accepts shared and varying directions", {
  truth <- c(TRUE, TRUE, FALSE)
  detections <- list(
    list(m = c(TRUE, TRUE, FALSE)),
    list(m = c(TRUE, FALSE, FALSE))
  )
  directions <- list(
    list(m = c(1, -1, 0)),
    list(m = c(1, 1, 0))
  )
  truth_directions <- list(c(1, 1, 0), c(1, 1, 0))
  result <- compare_detections(
    detections, truth, replicates = TRUE,
    metrics = c("type_iii", "directional_power"),
    directions = directions,
    truth_direction = truth_directions
  )
  expect_equal(result$n_replicates, 2L)
  expect_true("TypeIII_mean" %in% names(result))
})

test_that("replicate components distinguish methods from replicates", {
  by_method <- list(MK = 1:2, CMK = 3:4)
  expect_identical(
    sptrends:::.replicate_component(
      by_method, i = 1L, n_rep = 2L,
      method_names = c("MK", "CMK")
    ),
    by_method
  )

  by_replicate <- list(first = 1:2, second = 3:4)
  expect_identical(
    sptrends:::.replicate_component(
      by_replicate, i = 2L, n_rep = 2L,
      method_names = c("MK", "CMK")
    ),
    3:4
  )

  shared_vector <- c(1, 0, -1)
  expect_identical(
    sptrends:::.replicate_component(
      shared_vector, i = 2L, n_rep = 2L,
      method_names = c("MK", "CMK")
    ),
    shared_vector
  )
})
