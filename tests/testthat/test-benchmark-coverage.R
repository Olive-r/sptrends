test_that("benchmark engine rejects malformed experiment definitions", {
  method <- list(m = function(input, simulation) input)
  expect_error(benchmark_methods(list(), method), "non-empty")
  expect_error(benchmark_methods(list(list()), method), "must be named")
  expect_error(benchmark_methods(list(a = list()), list()), "methods")
  expect_error(
    benchmark_methods(list(a = list()), list(m = 1)), "methods"
  )
  expect_error(
    benchmark_methods(list(a = list()), method, n_replicates = 0),
    "positive integer"
  )
  expect_error(
    benchmark_methods(list(a = list()), method, stage = "custom"),
    "evaluator"
  )
  expect_error(
    benchmark_methods(list(a = list()), method, prepare = 1),
    "prepare"
  )
  expect_error(
    benchmark_methods(list(a = 1), method), "not an argument list"
  )
})

test_that("benchmark metadata covers absent truth and reserved names", {
  series <- terra::rast(nrows = 2, ncols = 3)
  simulation <- list(series = series)
  metadata <- sptrends:::.benchmark_scenario_data(
    list(label = "none"), simulation
  )
  expect_equal(metadata$DomainCells, 6)
  expect_true(is.na(metadata$TrueNulls))

  vector_simulation <- list(series = 1:3)
  vector_metadata <- sptrends:::.benchmark_scenario_data(
    list(label = "vector"), vector_simulation
  )
  expect_true(is.na(vector_metadata$DomainCells))
  expect_error(
    sptrends:::.benchmark_scenario_data(
      list(Method = "reserved"), simulation
    ),
    "reserved"
  )
})

test_that("benchmark score validators reject malformed outputs", {
  expect_error(
    sptrends:::.validate_benchmark_scores(list(), "m"), "data frame"
  )
  expect_error(
    sptrends:::.validate_benchmark_scores(
      data.frame(Method = c("m", "m")), "m"
    ),
    "exactly one row"
  )

  simulation <- sim_trend_stack(
    nrow = 2, ncol = 2, n_time = 4,
    constant_block = FALSE, spatial_model = "independent", seed = 2
  )
  expect_error(
    sptrends:::.benchmark_score(
      list(m = 1:3), simulation, "slope", NULL,
      character(), NULL
    ),
    "differ in length"
  )
  expect_error(
    sptrends:::.benchmark_prewhitening(
      list(m = 1:4), simulation
    ),
    "SpatRaster"
  )
  short <- simulation$series[[1:2]]
  expect_error(
    sptrends:::.benchmark_prewhitening(
      list(m = short), simulation
    ),
    "at least three"
  )

  wrong_cells <- terra::rast(nrows = 1, ncols = 5, nlyrs = 4)
  terra::values(wrong_cells) <- matrix(1, nrow = 5, ncol = 4)
  expect_error(
    sptrends:::.benchmark_prewhitening(
      list(m = wrong_cells), simulation
    ),
    "differ in cell count"
  )
})

test_that("prewhitening scores degenerate series and modification flags", {
  simulation <- sim_trend_stack(
    nrow = 2, ncol = 2, n_time = 4,
    constant_block = FALSE, spatial_model = "independent", seed = 12
  )
  transformed <- simulation$series
  terra::values(transformed) <- rbind(
    c(1, 1, 1, 1),
    c(NA, 1, 2, 3),
    c(1, 2, 1, 2),
    c(1, 3, 2, 4)
  )
  diagnostics <- transformed[[1]]
  names(diagnostics) <- "Modified"
  terra::values(diagnostics) <- c(1, 0, NA, 1)

  result <- sptrends:::.benchmark_prewhitening(
    list(method = list(
      series = transformed,
      diagnostics = diagnostics
    )),
    simulation
  )

  expect_true(is.finite(result$ResidualACF1))
  expect_equal(result$ModifiedFraction, 2 / 3)
  expect_equal(result$OutputLength, 4)
})

test_that("benchmark masks and aliases retain paired experiment semantics", {
  scenarios <- list(
    one = list(
      nrow = 2, ncol = 2, n_time = 4,
      trend_fraction = 0, constant_block = FALSE,
      spatial_model = "independent"
    )
  )
  methods <- list(oracle = function(input, simulation) {
    list(significant = simulation$true_signal,
         direction = simulation$true_direction)
  })
  mask_calls <- 0L
  mask <- function(simulation) {
    mask_calls <<- mask_calls + 1L
    rep(TRUE, terra::ncell(simulation$true_signal))
  }
  result <- benchmark_methods(
    scenarios, methods, n_replicates = 1,
    stage = "detection", evaluation_mask = mask,
    metrics = "type_i", seed = 3
  )
  expect_equal(mask_calls, 1L)
  expect_identical(attr(result, "stage"), "trend_test")
})

test_that("benchmark reports all three timing indicators", {
  simulator_verbose <- NULL
  simulator <- function(seed, verbose = TRUE) {
    simulator_verbose <<- verbose
    sim_trend_stack(
      nrow = 2, ncol = 2, n_time = 3,
      trend_fraction = 0, constant_block = FALSE,
      spatial_model = "independent", seed = seed, verbose = FALSE
    )
  }
  methods <- list(oracle = function(input, simulation) {
    list(
      significant = simulation$true_signal,
      direction = simulation$true_direction
    )
  })
  output <- capture.output(
    benchmark_methods(
      list(one = list()), methods, n_replicates = 2,
      simulator = simulator, metrics = "type_i", verbose = TRUE
    ),
    type = "message"
  )
  text <- paste(output, collapse = " ")
  expect_false(simulator_verbose)
  expect_match(text, "progress:", fixed = TRUE)
  expect_match(text, "remaining:", fixed = TRUE)
  expect_match(text, "elapsed:", fixed = TRUE)
})

test_that("benchmark summary validates structure and writes output", {
  expect_error(benchmark_summary(data.frame()), "returned by")
  malformed <- structure(
    data.frame(Scenario = "a"),
    class = c("sptrends_benchmark", "data.frame")
  )
  expect_error(benchmark_summary(malformed), "lacks")

  x <- structure(
    data.frame(
      Scenario = rep("a", 2), Method = "m", Replicate = 1:2,
      Score = c(1, 2)
    ),
    scenario_fields = character(),
    class = c("sptrends_benchmark", "sptrends", "data.frame")
  )
  path <- tempfile(fileext = ".csv")
  summary <- benchmark_summary(x, path = path)
  expect_true(file.exists(path))
  expect_equal(summary$Score_mean, 1.5)
  unlink(path)
})

test_that("benchmark summaries report all three timing indicators", {
  x <- structure(
    data.frame(
      Scenario = c("one", "one"), Method = c("a", "b"),
      Replicate = c(1L, 1L), Score = c(1, 2)
    ),
    scenario_fields = character(),
    class = c("sptrends_benchmark", "sptrends", "data.frame")
  )
  output <- capture.output(
    benchmark_summary(x, verbose = TRUE), type = "message"
  )
  text <- paste(output, collapse = " ")
  expect_match(text, "progress:", fixed = TRUE)
  expect_match(text, "remaining:", fixed = TRUE)
  expect_match(text, "elapsed:", fixed = TRUE)
})
