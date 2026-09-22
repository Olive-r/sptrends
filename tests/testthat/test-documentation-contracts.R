.installed_vignette_lines <- function(file) {
  path <- system.file("doc", file, package = "sptrends")
  if (!nzchar(path) || !file.exists(path)) {
    testthat::skip(paste("Installed vignette source unavailable:", file))
  }
  readLines(path, warn = FALSE)
}

test_that("six introductory vignettes share the teaching architecture", {
  files <- c(
    "a-getting-started.Rmd", "b-prewhitening.Rmd", "c-trend-test.Rmd",
    "d-slope-estimation.Rmd", "e-fdr-correction.Rmd",
    "g-workflow-trends.Rmd"
  )
  headings <- c(
    "## Why this matters", "## Basic workflow",
    "## Understanding the results", "## Choosing the main options",
    "## Common mistakes", "## Next steps", "## Further details"
  )
  for (file in files) {
    lines <- .installed_vignette_lines(file)
    info <- paste("vignette:", file)
    expect_true(any(startsWith(lines, "## What ")), info = info)
    expect_true(all(headings %in% lines), info = info)
  }
})

test_that("published workflows are integrated in trend workflows", {
  workflow <- .installed_vignette_lines("g-workflow-trends.Rmd")
  expect_true(any(grepl("5. Trend workflows", workflow, fixed = TRUE)))
  expect_true(any(grepl("### Published workflows", workflow,
                        fixed = TRUE)))
  expect_true(any(grepl("methodological origin of sptrends", workflow,
                        fixed = TRUE)))
  expect_true(any(grepl("workflow_tst()", workflow, fixed = TRUE)))
  expect_true(any(grepl("workflow_rta()", workflow, fixed = TRUE)))
})

test_that("prewhitening and FDR vignettes retain teaching contracts", {
  prewhitening <- .installed_vignette_lines("b-prewhitening.Rmd")
  fdr <- .installed_vignette_lines("e-fdr-correction.Rmd")

  expect_true(all(vapply(
    c("TFPW_WS", "TFPW_Y", "TFPW_Z", "VCTFPW"),
    function(method) any(grepl(method, prewhitening, fixed = TRUE)),
    logical(1)
  )))
  expect_true(any(grepl("Classical prewhitening is deliberately",
                        prewhitening, fixed = TRUE)))
  expect_true(all(vapply(c("BH", "BKY", "BY"), function(method) {
    any(grepl(paste0("`", method, "`"), fdr, fixed = TRUE))
  }, logical(1))))
  expect_true(any(grepl(
    "The target *q* defines the false discovery rate", fdr,
    fixed = TRUE
  )))
  expect_true(any(grepl(
    "It is not an adjusted *α*", fdr,
    fixed = TRUE
  )))
})

test_that("exports resolve to installed help topics", {
  exported <- getNamespaceExports("sptrends")
  for (topic in exported) {
    help_topic <- utils::help(topic, package = "sptrends")
    expect_true(length(help_topic) > 0L, info = topic)
  }
})

test_that("local spatial autocorrelation retains its public API", {
  arguments <- names(formals(spatial_autocorrelation))
  required <- c(
    "x", "method", "scope", "connectivity", "nperm", "alternative",
    "alpha", "seed", "n_cores", "precomputed_neighbourhood",
    "report", "verbose"
  )
  expect_identical(arguments, required)
  expect_true(length(utils::help(
    "spatial_autocorrelation", package = "sptrends"
  )) > 0L)
})

.public_help_topics <- c(
  "benchmark_methods", "benchmark_summary", "compare_detections",
  "compute_anomalies",
  "example_data",
  "fdr_correction", "inspect_ts_cell", "prewhiten",
  "read_netcdf_stack", "read_ordered_stack", "sim_trend_stack",
  "simulation_design",
  "slope_estimator", "spatial_autocorrelation", "trend_test",
  "workflow_rta", "workflow_trends", "workflow_tst",
  "print.sptrends", "summary.sptrends", "plot.sptrends"
)

.rd_for_alias <- function(alias) {
  database <- tools::Rd_db("sptrends")
  matches <- vapply(database, function(topic) {
    text <- paste(as.character(topic), collapse = "")
    grepl(paste0("\\alias{", alias, "}"), text, fixed = TRUE)
  }, logical(1))
  if (sum(matches) != 1L) {
    stop("Expected exactly one help topic for alias: ", alias)
  }
  database[[which(matches)]]
}

.rd_text <- function(alias) {
  paste(as.character(.rd_for_alias(alias)), collapse = "")
}

test_that("public help follows the common documentation architecture", {
  for (topic in .public_help_topics) {
    text <- .rd_text(topic)
    info <- paste("help topic:", topic)

    expect_true(grepl("Function type:", text, fixed = TRUE), info = info)
    expect_true(grepl("\\usage", text, fixed = TRUE), info = info)
    expect_true(grepl("\\arguments", text, fixed = TRUE), info = info)
    expect_true(grepl("\\value", text, fixed = TRUE), info = info)
    expect_true(grepl("\\section{Typical use}", text, fixed = TRUE),
                info = info)
    expect_true(grepl("\\section{Methodological details}", text,
                      fixed = TRUE), info = info)
    expect_true(grepl("\\seealso", text, fixed = TRUE), info = info)
    expect_true(grepl("\\examples", text, fixed = TRUE), info = info)

    function_type <- regexpr("Function type:", text, fixed = TRUE)[[1L]]
    typical <- regexpr("\\section{Typical use}", text, fixed = TRUE)[[1L]]
    details <- regexpr(
      "\\section{Methodological details}", text, fixed = TRUE
    )[[1L]]
    expect_true(function_type < typical, info = info)
    expect_true(typical < details, info = info)
    expect_false(grepl("\\section{Function type}", text, fixed = TRUE),
                 info = info)
  }
})

test_that("Function type is never rendered as a late custom section", {
  database <- tools::Rd_db("sptrends")
  for (topic in database) {
    text <- paste(as.character(topic), collapse = "")
    expect_false(grepl("\\section{Function type}", text, fixed = TRUE))
  }
})

test_that("public help contains one copy of each common custom section", {
  count_fixed <- function(text, pattern) {
    lengths(strsplit(text, pattern, fixed = TRUE)) - 1L
  }

  for (topic in .public_help_topics) {
    text <- .rd_text(topic)
    info <- paste("help topic:", topic)
    expect_identical(
      count_fixed(text, "\\section{Typical use}"), 1L, info = info
    )
    expect_identical(
      count_fixed(text, "\\section{Methodological details}"),
      1L,
      info = info
    )
  }
})

test_that("simulation APIs retain the common methodological subsections", {
  topics <- c(
    "sim_trend_stack", "simulation_design", "benchmark_methods",
    "benchmark_summary", "compare_detections"
  )
  required <- c("Computational considerations", "Limitations",
                "Quality assurance")
  for (topic in topics) {
    text <- .rd_text(topic)
    for (heading in required) {
      expect_true(grepl(tolower(heading), tolower(text), fixed = TRUE),
                  info = paste(topic, heading))
    }
  }
})

test_that("non-simulation help examples use the bundled real dataset", {
  simulation_topics <- c(
    "sim_trend_stack", "simulation_design", "benchmark_methods",
    "benchmark_summary", "compare_detections"
  )
  topics <- setdiff(.public_help_topics, simulation_topics)
  for (topic in topics) {
    text <- .rd_text(topic)
    examples_start <- regexpr("\\examples", text, fixed = TRUE)[[1L]]
    examples <- substring(text, examples_start)
    expect_false(grepl("sim_trend_stack(", examples, fixed = TRUE),
                 info = topic)
  }
})

test_that("Usage advertises every supported method family", {
  expect_identical(
    eval(formals(prewhiten)$method),
    c("TFPW_WS", "TFPW_Y", "TFPW_Z", "VCTFPW")
  )
  expect_identical(
    eval(formals(trend_test)$method),
    c("CMK", "MK", "OLS", "MMK")
  )
  expect_identical(eval(formals(trend_test)$window_size), 3L)
  expect_identical(
    eval(formals(slope_estimator)$method),
    c("TS", "OLS", "RM")
  )
  expect_identical(
    eval(formals(fdr_correction)$method),
    c("BH", "BKY", "BY")
  )
  expect_identical(
    eval(formals(spatial_autocorrelation)$method),
    c("moran", "getis_ord")
  )
  expect_identical(
    eval(formals(spatial_autocorrelation)$scope),
    c("global", "local")
  )
  expect_identical(
    eval(formals(workflow_trends)$prewhiten_method),
    c("TFPW_WS", "TFPW_Y", "TFPW_Z", "VCTFPW", "none")
  )
  expect_identical(
    eval(formals(workflow_trends)$trend_method),
    c("CMK", "MK", "OLS", "MMK")
  )
  expect_identical(
    eval(formals(workflow_trends)$slope_method),
    c("TS", "OLS", "RM")
  )
  expect_identical(
    eval(formals(workflow_trends)$fdr_method),
    c("BKY", "BH", "BY")
  )
})
