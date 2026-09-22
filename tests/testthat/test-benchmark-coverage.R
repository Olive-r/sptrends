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
  # The incomplete second series is outside the common evaluation domain.
  expect_equal(result$ModifiedFraction, 1)
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

test_that("benchmark_methods supports heterogeneous scenario arguments -- a real design gap found by an external audit: scenarios need not share the same argument names (a null scenario need not set spatial_rho; different spatial_model choices take different parameters), and the previous implementation only failed late, at the final rbind, after every simulation had already run", {
  scenarios <- list(
    null = list(
      nrow = 4, ncol = 4, n_time = 6, trend_fraction = 0
    ),
    spatial = list(
      nrow = 4, ncol = 4, n_time = 6, trend_fraction = 0.5,
      spatial_model = "gaussian", spatial_rho = 0.5
    )
  )

  methods <- list(
    MK = function(x, simulation) {
      result <- trend_test(x, method = "MK", report = FALSE,
                           verbose = FALSE)
      result$stats$p <= 0.05
    }
  )

  result <- benchmark_methods(
    scenarios = scenarios, methods = methods, n_replicates = 1,
    stage = "trend_test", seed = 1, verbose = FALSE
  )

  expect_s3_class(result, "sptrends_benchmark")
  expect_equal(nrow(result), 2)
  expect_true(all(c("spatial_model", "spatial_rho") %in% names(result)))

  expect_true(is.na(result$spatial_model[result$Scenario == "null"]))
  expect_true(is.na(result$spatial_rho[result$Scenario == "null"]))

  expect_equal(result$spatial_model[result$Scenario == "spatial"],
              "gaussian")
  expect_equal(result$spatial_rho[result$Scenario == "spatial"], 0.5)
})

test_that("benchmark_methods validates scenario argument names before running any simulation", {
  methods <- list(
    MK = function(x, simulation) {
      result <- trend_test(x, method = "MK", report = FALSE,
                           verbose = FALSE)
      result$stats$p <= 0.05
    }
  )

  scenarios_sin_nombre <- list(
    a = list(4, ncol = 4, n_time = 6)  # primer argumento sin nombre
  )
  expect_error(
    benchmark_methods(scenarios = scenarios_sin_nombre, methods = methods,
                      n_replicates = 1, seed = 1, verbose = FALSE),
    "uniquely named"
  )

  scenarios_duplicados <- list(
    a = list(nrow = 4, nrow = 5, ncol = 4, n_time = 6)
  )
  expect_error(
    benchmark_methods(scenarios = scenarios_duplicados, methods = methods,
                      n_replicates = 1, seed = 1, verbose = FALSE),
    "uniquely named"
  )

  scenarios_reservado <- list(
    a = list(nrow = 4, ncol = 4, n_time = 6, Method = "trampa")
  )
  expect_error(
    benchmark_methods(scenarios = scenarios_reservado, methods = methods,
                      n_replicates = 1, seed = 1, verbose = FALSE),
    "reserved"
  )
})

test_that("benchmark_methods gives a clear error, not terra's own cryptic one, when a custom simulator omits 'true_slope' for stage = 'slope' or 'prewhitening'", {
  simulador_sin_slope <- function(nrow, ncol, n_time, seed) {
    r <- terra::rast(nrow = nrow, ncol = ncol, nlyrs = n_time)
    terra::values(r) <- stats::rnorm(terra::ncell(r) * n_time)
    list(series = r, true_signal = NULL, true_direction = NULL,
         true_slope = NULL)
  }
  methods_slope <- list(m = function(x, simulation) x[[1]])
  expect_error(
    benchmark_methods(
      scenarios = list(a = list(nrow = 3, ncol = 3, n_time = 4)),
      methods = methods_slope, n_replicates = 1, stage = "slope",
      simulator = simulador_sin_slope, seed = 1, verbose = FALSE
    ),
    "true_slope"
  )

  methods_pw <- list(m = function(x, simulation) list(series = x))
  expect_error(
    benchmark_methods(
      scenarios = list(a = list(nrow = 3, ncol = 3, n_time = 4)),
      methods = methods_pw, n_replicates = 1, stage = "prewhitening",
      simulator = simulador_sin_slope, seed = 1, verbose = FALSE
    ),
    "true_slope"
  )
})

test_that(".benchmark_score(stage = 'slope') accepts a SpatRaster evaluation_mask and rejects one of mismatched length", {
  sim <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 6, seed = 1)
  outputs <- list(m = sim$true_slope)
  mascara_raster <- sim$true_slope > -Inf  # SpatRaster logico, mismo tamano

  expect_error(
    sptrends:::.benchmark_score(outputs, sim, stage = "slope",
                                evaluator = NULL, metrics = NULL,
                                evaluation_mask = mascara_raster),
    NA
  )

  mascara_corta <- terra::values(mascara_raster, mat = FALSE)[1:2]  # longitud distinta a proposito
  expect_error(
    sptrends:::.benchmark_score(outputs, sim, stage = "slope",
                                evaluator = NULL, metrics = NULL,
                                evaluation_mask = mascara_corta),
    "differ in length"
  )
})

test_that(".benchmark_score(stage = 'slope') rejects an output whose length does not match true_slope", {
  sim <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 6, seed = 2)
  outputs_mal <- list(m = terra::values(sim$true_slope, mat = FALSE)[1:3])

  expect_error(
    sptrends:::.benchmark_score(outputs_mal, sim, stage = "slope",
                                evaluator = NULL, metrics = NULL,
                                evaluation_mask = NULL),
    "differ in length"
  )
})

test_that(".benchmark_prewhitening rejects a non-SpatRaster output and one with mismatched cell count", {
  sim <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 6, seed = 3)

  outputs_no_raster <- list(m = list(series = "no es un raster"))
  expect_error(
    sptrends:::.benchmark_prewhitening(outputs_no_raster, sim,
                                       evaluation_mask = NULL),
    "must contain a SpatRaster series"
  )

  serie_corta <- sim_trend_stack(nrow = 2, ncol = 2, n_time = 6, seed = 99)$series  # menos celdas que true_slope
  outputs_pocas_celdas <- list(m = list(series = serie_corta))
  expect_error(
    sptrends:::.benchmark_prewhitening(outputs_pocas_celdas, sim,
                                       evaluation_mask = NULL),
    "differ in cell count"
  )
})

test_that(".benchmark_common_domain accepts a SpatRaster evaluation_mask, rejects mismatched length, and rejects an all-FALSE domain", {
  sim <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 6, seed = 4)
  truth <- terra::values(sim$true_slope, mat = FALSE)
  validity <- list(rep(TRUE, length(truth)))
  mascara_raster <- sim$true_slope > -Inf

  expect_error(
    sptrends:::.benchmark_common_domain(validity, truth, mascara_raster),
    NA
  )

  mascara_corta <- terra::values(mascara_raster, mat = FALSE)[1:2]
  expect_error(
    sptrends:::.benchmark_common_domain(validity, truth, mascara_corta),
    "differ in length"
  )

  mascara_vacia <- sim$true_slope * 0  # todo cero -> ningun valido
  expect_error(
    sptrends:::.benchmark_common_domain(validity, truth, mascara_vacia),
    "No valid cells"
  )
})

test_that(".benchmark_scenario_data accepts a SpatRaster evaluation_mask and rejects one of mismatched length", {
  sim <- sim_trend_stack(nrow = 3, ncol = 3, n_time = 6, seed = 5)
  mascara_raster <- sim$true_slope != 0  # SpatRaster logico, mismo tamano

  expect_error(
    sptrends:::.benchmark_scenario_data(
      arguments = list(nrow = 3, ncol = 3), simulation = sim,
      evaluation_mask = mascara_raster
    ),
    NA
  )

  mascara_corta <- terra::values(mascara_raster, mat = FALSE)[1:2]
  expect_error(
    sptrends:::.benchmark_scenario_data(
      arguments = list(nrow = 3, ncol = 3), simulation = sim,
      evaluation_mask = mascara_corta
    ),
    "differ in length"
  )
})
