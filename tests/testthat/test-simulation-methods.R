test_that("simulation objects provide unified S3 presentation", {
  simulation <- sim_trend_stack(
    nrow = 4, ncol = 5, n_time = 4,
    trend_shape = "square", signal_size = 2,
    trend_fraction = 1, constant_block = FALSE,
    spatial_model = "independent", seed = 1)

  expect_s3_class(simulation, "sptrends_simulation")
  expect_s3_class(simulation, "sptrends")
  expect_match(paste(capture.output(print(simulation)), collapse = " "),
               "sptrends simulation")
  report <- suppressMessages(summary(simulation))
  expect_true(all(c("DomainCells", "TrueSignals", "Pi0") %in%
                  names(report)))
  expect_equal(report$DomainCells, 20)
  expect_equal(report$TrueSignals, 4)
  report_path <- tempfile(fileext = ".csv")
  suppressMessages(summary(simulation, path = report_path))
  expect_true(file.exists(report_path))
  unlink(report_path)

  for (view in c("truth", "slope", "direction", "breaks", "series")) {
    expect_error(plot(simulation, which = view), NA)
  }
  plot_path <- tempfile(fileext = ".png")
  plot(simulation, path = plot_path)
  expect_true(file.exists(plot_path))
  unlink(plot_path)
})

test_that("simulation designs provide S3 inspection", {
  design <- simulation_design(
    spatial_rho = c(0.2, 0.8), ar1 = c(0, 0.5),
    constants = list(nrow = 4, ncol = 4, n_time = 5))

  expect_s3_class(design, "sptrends_simulation_design")
  expect_s3_class(design, "sptrends")
  expect_match(paste(capture.output(print(design)), collapse = " "),
               "Scenarios: 4")
  report <- suppressMessages(summary(design))
  expect_identical(report$Factor, c("spatial_rho", "ar1"))
  expect_equal(report$Levels, c(2L, 2L))
  report_path <- tempfile(fileext = ".csv")
  suppressMessages(summary(design, path = report_path))
  expect_true(file.exists(report_path))
  unlink(report_path)
  plot_path <- tempfile(fileext = ".png")
  expect_error(plot(design, path = plot_path), NA)
  expect_true(file.exists(plot_path))
  unlink(plot_path)
})

.benchmark_presentation_fixture <- function() {
  scenarios <- simulation_design(
    spatial_rho = c(0.3, 0.7), ar1 = c(0, 0.4),
    constants = list(
      nrow = 4, ncol = 4, n_time = 5,
      spatial_model = "exponential", trend_shape = "square",
      signal_size = 2, trend_fraction = 1, constant_block = FALSE))
  methods <- list(
    oracle = function(series, simulation) {
      list(significant = simulation$true_signal,
           direction = simulation$true_direction)
    },
    empty = function(series, simulation) {
      list(significant = simulation$true_signal == 2,
           direction = simulation$true_direction)
    })
  benchmark_methods(
    scenarios, methods, n_replicates = 2, seed = 4,
    metrics = c("type_i", "type_ii", "type_iii",
                "field_power", "global_power", "within_image_power"))
}

test_that("benchmarks retain explicit scenario and truth columns", {
  result <- .benchmark_presentation_fixture()

  expect_s3_class(result, "sptrends_benchmark")
  expect_s3_class(result, "sptrends")
  expect_true(all(c(
    "spatial_rho", "ar1", "DomainCells", "TrueNulls", "Pi0",
    "SignalProportion") %in% names(result)))
  expect_equal(sort(unique(result$spatial_rho)), c(0.3, 0.7))
  expect_equal(sort(unique(result$ar1)), c(0, 0.4))
  expect_true(all(result$DomainCells == 16))
  expect_true(all(result$Pi0 == 0.75))
  expect_identical(
    attr(result, "scenario_fields"),
    c("spatial_rho", "ar1", "nrow", "ncol", "n_time",
      "spatial_model", "trend_shape", "signal_size",
      "trend_fraction", "constant_block"))
  expect_identical(
    attr(result, "truth_fields"),
    c("DomainCells", "TrueNulls", "Pi0", "SignalProportion"))
})

test_that("benchmark summaries retain scenario factors", {
  result <- .benchmark_presentation_fixture()
  report <- suppressMessages(summary(result))

  expect_s3_class(report, "sptrends_benchmark_summary")
  expect_true(all(c("spatial_rho", "ar1", "Method",
                    "Pi0_mean", "TypeII_mean", "TypeII_sd") %in%
                  names(report)))
  expect_equal(nrow(report), 8L)
  expect_match(paste(capture.output(print(result)), collapse = " "),
               "scenarios: 4")
  report_path <- tempfile(fileext = ".csv")
  benchmark_summary(result, path = report_path)
  expect_true(file.exists(report_path))
  unlink(report_path)
})

test_that("benchmark S3 plots cover scenario-performance displays", {
  result <- .benchmark_presentation_fixture()
  plots <- list(
    field = list(metric = "FieldPower", scenario = "spatial_rho",
                 facet = "ar1", type = "line"),
    line = list(metric = "WithinImagePower", scenario = "spatial_rho",
                facet = "ar1", type = "line"),
    bar = list(metric = "TypeII", scenario = "spatial_rho",
               type = "bar"),
    boxplot = list(metric = "TypeII", scenario = "spatial_rho",
                   type = "boxplot"),
    heatmap = list(metric = "TypeII", scenario = "spatial_rho",
                   facet = "ar1", type = "heatmap"),
    profile = list(
      metric = c("TypeI", "TypeII", "WithinImagePower"),
      scenario = "spatial_rho", type = "profile"))

  for (arguments in plots) {
    path <- tempfile(fileext = ".png")
    arguments$x <- result
    arguments$path <- path
    expect_warning(do.call(plot, arguments), NA)
    expect_true(file.exists(path))
    expect_gt(file.info(path)$size, 0)
    unlink(path)
  }

  expect_error(plot(result), NA)
  expect_error(
    plot(result, metric = "TypeII", scenario = "spatial_rho",
         type = "profile", interval = "none"),
    NA)
  for (interval in c("se", "sd")) {
    expect_warning(
      plot(result, metric = "TypeII", scenario = "spatial_rho",
           interval = interval),
      NA)
  }
  expect_error(
    plot(
      result,
      metric = c("EmpiricalFDR", "EmpiricalFWER"),
      scenario = "spatial_rho", type = "profile"),
    NA)
})

test_that("benchmark plots draw non-degenerate uncertainty intervals", {
  result <- .benchmark_presentation_fixture()
  result$TypeII <- seq(0.05, 0.8, length.out = nrow(result))

  bar_path <- tempfile(fileext = ".png")
  expect_warning(
    plot(
      result, metric = "TypeII", scenario = "spatial_rho",
      type = "bar", interval = "sd", path = bar_path
    ),
    NA
  )
  expect_true(file.exists(bar_path))
  unlink(bar_path)

  line_path <- tempfile(fileext = ".png")
  expect_warning(
    plot(
      result, metric = "TypeII", scenario = "spatial_rho",
      type = "line", interval = "sd", path = line_path
    ),
    NA
  )
  expect_true(file.exists(line_path))
  unlink(line_path)

  boxplot_path <- tempfile(fileext = ".png")
  expect_warning(
    plot(
      result, metric = "TypeII", scenario = "spatial_rho",
      facet = "ar1", type = "boxplot", path = boxplot_path
    ),
    NA
  )
  expect_true(file.exists(boxplot_path))
  unlink(boxplot_path)
})

test_that("benchmark plots validate requested dimensions", {
  result <- .benchmark_presentation_fixture()
  expect_error(plot(result, metric = "missing"), "Unknown numerical")
  expect_error(
    plot(result, metric = c("TypeI", "TypeII"), type = "line"),
    "profile")
  expect_error(
    plot(result, metric = "TypeI", type = "heatmap"),
    "second factor")
  expect_error(plot(result, scenario = "missing"), "grouping column")
  expect_error(plot(result, level = 1), "strictly between")

  no_metric <- result[, c("Scenario", "Replicate", "Seed", "Elapsed",
                          "Method"), drop = FALSE]
  class(no_metric) <- c(
    "sptrends_benchmark", "sptrends", "data.frame")
  attr(no_metric, "scenario_fields") <- character()
  expect_error(plot(no_metric), "No numerical performance")
})

test_that("default benchmark plots adapt to non-detection stages", {
  scenarios <- simulation_design(
    ar1 = c(0, 0.4),
    constants = list(
      nrow = 3, ncol = 3, n_time = 4,
      constant_block = FALSE, spatial_model = "independent"))
  methods <- list(
    oracle = function(series, simulation) simulation$true_slope)
  result <- benchmark_methods(
    scenarios, methods, n_replicates = 1,
    stage = "slope", seed = 5)

  expect_warning(plot(result), NA)

  one_scenario <- result[result$Scenario == result$Scenario[1L], ]
  attr(one_scenario, "scenario_fields") <- character()
  class(one_scenario) <- c(
    "sptrends_benchmark", "sptrends", "data.frame")
  expect_warning(plot(one_scenario, metric = "Bias"), NA)
})

test_that("scenario metadata supports external simulator structures", {
  vector_truth <- sptrends:::.benchmark_scenario_data(
    list(setting = "vector", signal_size = c(4, 8),
         custom_mask = matrix(TRUE, 10, 10), long = 1:5),
    list(true_signal = c(1, 0, 0, NA)))
  expect_equal(vector_truth$DomainCells, 3)
  expect_equal(vector_truth$Pi0, 2 / 3)
  expect_identical(vector_truth$signal_size, "4,8")
  expect_identical(vector_truth$custom_mask, "<matrix[10x10]>")
  expect_identical(vector_truth$long, "<integer[5]>")

  raster_mask <- terra::rast(nrows = 2, ncols = 2)
  raster_metadata <- sptrends:::.benchmark_scenario_data(
    list(custom_mask = raster_mask), list(true_signal = c(0, 1, 0, 1)))
  expect_identical(
    raster_metadata$custom_mask, "<SpatRaster[4 cells]>")

  slope_truth <- sptrends:::.benchmark_scenario_data(
    list(setting = "slope"),
    list(true_slope = c(0, 0.2, -0.1, 0)))
  expect_equal(slope_truth$TrueNulls, 2)

  slope_raster <- terra::rast(nrows = 2, ncols = 2)
  terra::values(slope_raster) <- c(0, 1, 0, -1)
  slope_raster_truth <- sptrends:::.benchmark_scenario_data(
    list(setting = "slope_raster"), list(true_slope = slope_raster))
  expect_equal(slope_raster_truth$TrueNulls, 2)

  series <- terra::rast(nrows = 2, ncols = 3, nlyrs = 2)
  no_truth <- sptrends:::.benchmark_scenario_data(
    list(setting = "custom"), list(series = series))
  expect_equal(no_truth$DomainCells, 6)
  expect_true(is.na(no_truth$Pi0))

  unknown <- sptrends:::.benchmark_scenario_data(
    list(setting = "custom"), list(series = matrix(1, 2, 2)))
  expect_true(is.na(unknown$DomainCells))
  empty_truth <- sptrends:::.benchmark_scenario_data(
    list(setting = "empty"), list(true_signal = c(NA, NA)))
  expect_true(is.na(empty_truth$Pi0))
  expect_error(
    sptrends:::.benchmark_scenario_data(
      list(Scenario = "reserved"), list(true_signal = c(0, 1))),
    "reserved")
})
