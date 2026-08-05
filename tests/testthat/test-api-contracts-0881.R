test_that("public Usage exposes every supported method without changing defaults", {
  expect_identical(
    eval(formals(fdr_correction)$method),
    c("BH", "BKY", "BY")
  )
  expect_identical(
    eval(formals(workflow_trends)$slope_method),
    c("TS", "OLS", "RM")
  )
  expect_identical(
    eval(formals(workflow_trends)$fdr_method),
    c("BKY", "BH", "BY")
  )
  expect_identical(
    eval(formals(workflow_tst)$fdr_method),
    c("BKY", "BH", "BY")
  )

  p <- c(0.001, 0.01, 0.2, 0.8)
  result <- fdr_correction(p, report = FALSE, verbose = FALSE)
  expect_false(is.null(result$q_BH))
  expect_false(is.null(result$q_BKY))
  expect_true(is.null(result$q_BY))
})

test_that("workflow FDR choices are validated before computation", {
  expect_error(
    workflow_trends(matrix(1:4, 2), fdr_method = "bonferroni",
                    report = FALSE, verbose = FALSE),
    "'fdr_method'"
  )
  expect_error(
    workflow_trends(matrix(1:4, 2), fdr_method = character(),
                    report = FALSE, verbose = FALSE),
    "'fdr_method'"
  )
  expect_error(
    workflow_tst(matrix(1:4, 2), fdr_method = "bonferroni",
                 report = FALSE, verbose = FALSE),
    "'fdr_method'"
  )
  expect_error(
    workflow_tst(matrix(1:4, 2), fdr_method = character(),
                 report = FALSE, verbose = FALSE),
    "'fdr_method'"
  )
})

test_that("package metadata contains both author ORCIDs", {
  authors <- utils::packageDescription("sptrends", fields = "Authors@R")
  expect_true(grepl("0000-0003-2580-5465", authors, fixed = TRUE))
  expect_true(grepl("0000-0002-5514-2941", authors, fixed = TRUE))
})

test_that("ordered-stack Usage advertises common GDAL raster extensions", {
  pattern <- eval(formals(read_ordered_stack)$pattern)
  for (extension in c("tif", "tiff", "nc", "grd", "img", "vrt", "asc")) {
    expect_true(grepl(pattern, paste0("year1.", extension)),
                info = extension)
  }
})
