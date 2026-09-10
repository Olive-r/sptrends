test_that("real spatial example uses one native-resolution domain", {
  series <- suppressMessages(
    read_ordered_stack(example_data("vhp_ndvi"))
  )
  expect_s4_class(series, "SpatRaster")
  expect_gt(terra::nlyr(series), 1L)

  complete <- which(stats::complete.cases(
    terra::values(series, mat = TRUE)
  ))
  centre <- terra::xyFromCell(
    series, complete[ceiling(length(complete) / 2)]
  )
  resolution <- terra::res(series)
  region <- terra::ext(c(
    centre[1] + c(-6, 6) * resolution[1],
    centre[2] + c(-6, 6) * resolution[2]
  ))
  subset <- terra::crop(series, region, snap = "near")
  expect_equal(terra::res(subset), terra::res(series))

  values <- terra::values(subset, mat = TRUE)
  ok <- stats::complete.cases(values)
  expect_gt(sum(ok), 2L)
  field <- subset[[1]]
  field_values <- terra::values(field, mat = FALSE)
  field_values[!ok] <- NA_real_
  terra::values(field) <- field_values
  neighbourhood <- sptrends:::prepare_cmk_neighbourhood(subset, ok)

  global <- spatial_autocorrelation(
    field, nperm = 3, seed = 1,
    precomputed_neighbourhood = neighbourhood,
    report = FALSE, verbose = FALSE
  )
  local <- spatial_autocorrelation(
    field, scope = "local", nperm = 3, seed = 1,
    precomputed_neighbourhood = neighbourhood,
    report = FALSE, verbose = FALSE
  )
  trend <- trend_test(
    subset, precomputed_neighbourhood = neighbourhood,
    report = FALSE, verbose = FALSE
  )

  expect_s3_class(global, "spatial_autocorrelation_global")
  expect_s3_class(local, "spatial_autocorrelation_local")
  expect_s3_class(trend, "trend_test")
})
