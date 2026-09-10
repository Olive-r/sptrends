test_that("simulation_design creates the complete factorial grid", {
  design <- simulation_design(
    spatial_model = c("independent", "exponential"),
    spatial_rho = c(0.3, 0.7),
    ar1 = c(0, 0.5),
    trend_strength = c(0, 0.05),
    constants = list(nrow = 10, ncol = 12, n_time = 20))

  expect_length(design, 16L)
  expect_identical(names(design), sprintf("scenario_%02d", 1:16))
  expect_true(all(vapply(design, `[[`, numeric(1), "nrow") == 10))
  combinations <- vapply(design, function(scenario) {
    paste(scenario$spatial_model, scenario$spatial_rho,
          scenario$ar1, scenario$trend_strength, sep = "|")
  }, character(1))
  expect_length(unique(combinations), 16L)
})

test_that("simulation_design preserves grouped values and constants", {
  design <- simulation_design(
    signal_size = list(c(4, 8), c(8, 12)),
    signal_location = c("centre", "random"),
    constants = list(signal_axis_ratio = 0.5), prefix = "shape")

  expect_length(design, 4L)
  expect_identical(design[[1]]$signal_size, c(4, 8))
  expect_true(all(vapply(
    design, `[[`, numeric(1), "signal_axis_ratio") == 0.5))
  expect_identical(names(design), paste0("shape_", 1:4))
})

test_that("simulation_design rejects ambiguous specifications", {
  expect_error(simulation_design(c(1, 2)), "must be named")
  expect_error(simulation_design(ar1 = numeric()), "at least one")
  expect_error(
    simulation_design(ar1 = c(0, 0.5), constants = list(ar1 = 0)),
    "duplicate")
  expect_error(simulation_design(ar1 = 0, constants = list(2)),
               "named list")
  expect_error(simulation_design(ar1 = 0, prefix = ""), "non-empty")
})

test_that("factorial scenarios run through the common benchmark engine", {
  design <- simulation_design(
    spatial_model = c("independent", "exponential"),
    ar1 = c(0, 0.4),
    constants = list(
      nrow = 3, ncol = 3, n_time = 4, spatial_rho = 0.5,
      trend_fraction = 0, constant_block = FALSE))
  methods <- list(oracle = function(series, simulation) {
    list(significant = simulation$true_signal,
         direction = simulation$true_direction)
  })

  result <- benchmark_methods(
    design, methods, n_replicates = 1, metrics = "type_i", seed = 2)
  expect_equal(nrow(result), 4L)
  expect_true(all(result$TypeI == 0))
})

test_that("simulation designs report all three timing indicators", {
  output <- capture.output(
    simulation_design(
      ar1 = c(0, 0.5), spatial_model = "independent", verbose = TRUE
    ),
    type = "message"
  )
  text <- paste(output, collapse = " ")
  expect_match(text, "progress:", fixed = TRUE)
  expect_match(text, "remaining:", fixed = TRUE)
  expect_match(text, "elapsed:", fixed = TRUE)
})
