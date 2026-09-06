test_that("example_data() with no arguments lists the bundled files", {
  files <- example_data()
  expect_type(files, "character")
  expect_true(any(grepl("vhp_ndvi", files, fixed = TRUE)))
  expect_true(any(grepl("Readme\\.txt$", files)))
})

test_that("example_data() returns a valid path to the dataset folder", {
  dir <- example_data("vhp_ndvi")
  expect_true(dir.exists(dir))
  tifs <- list.files(dir, pattern = "\\.tif$")
  expect_length(tifs, 42)
})

test_that("example_data() errors on a non-existent path", {
  expect_error(example_data("does_not_exist.tif"))
})

test_that("the bundled dataset works end-to-end with read_ordered_stack()", {
  r <- read_ordered_stack(example_data("vhp_ndvi"), report = FALSE, verbose = FALSE)
  expect_s4_class(r, "SpatRaster")
  expect_equal(terra::nlyr(r), 42)
})
