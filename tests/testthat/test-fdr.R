test_that("fdr_bh matches stats::p.adjust(method = 'BH')", {
  p <- c(0.001, 0.01, 0.02, 0.5, 0.8)
  expect_identical(fdr_bh(p)$q_value, stats::p.adjust(p, method = "BH"))
})

test_that("fdr_bky reproduces Example 1 of Benjamini, Krieger & Yekutieli (2006)", {
  p_paper <- c(0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298, 0.0344,
               0.0459, 0.3240, 0.4262, 0.5719, 0.6528, 0.7590, 1.0000)
  res <- fdr_bky(p_paper, q = 0.05)
  expect_identical(res$r1, 4L)
  expect_identical(res$m0_hat, 11L)
  expect_identical(sum(res$reject), 8L)
})

test_that("fdr_bh and fdr_bky preserve NA positions", {
  p <- c(0.01, NA, 0.2, 0.5, NA)
  res_bh <- fdr_bh(p)
  res_bky <- fdr_bky(p)
  expect_identical(which(is.na(res_bh$q_value)), c(2L, 5L))
  expect_identical(which(is.na(res_bky$q_value)), c(2L, 5L))
})

test_that("fdr_correction works on a plain numeric vector", {
  p <- c(0.001, 0.01, 0.02, 0.5, 0.8)
  res <- fdr_correction(p, method = c("BH", "BKY"), q = 0.05, report = FALSE, verbose = FALSE)
  expect_true(all(c("q_BH", "q_BKY", "reject_raw") %in% names(res)))
  expect_null(res$rasters)
})

test_that("fdr_correction rejects out-of-range p-value rasters", {
  r <- sim_trend_stack(nrow = 5, ncol = 5, n_time = 1, seed = 1)$series[[1]]
  expect_error(fdr_correction(r, report = FALSE, verbose = FALSE), "outside \\[0,1\\]")
})

test_that("fdr_bky's default implementation is 'multtest', unchanged from before", {
  p_paper <- c(0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298, 0.0344,
               0.0459, 0.3240, 0.4262, 0.5719, 0.6528, 0.7590, 1.0000)
  res_default <- fdr_bky(p_paper, q = 0.05)
  res_explicit <- fdr_bky(p_paper, q = 0.05, implementation = "multtest")
  expect_identical(res_default$reject, res_explicit$reject)
  expect_identical(res_default$implementation, "multtest")
})

test_that("fdr_bky's 'original' implementation shares Stage 1 with 'multtest' but is never more permissive", {
  p_paper <- c(0.0001, 0.0004, 0.0019, 0.0095, 0.0201, 0.0278, 0.0298, 0.0344,
               0.0459, 0.3240, 0.4262, 0.5719, 0.6528, 0.7590, 1.0000)
  res_multtest <- fdr_bky(p_paper, q = 0.05, implementation = "multtest")
  res_def6 <- fdr_bky(p_paper, q = 0.05, implementation = "original")

  # Stage 1 is identical regardless of implementation
  expect_identical(res_multtest$r1, res_def6$r1)
  expect_identical(res_multtest$m0_hat, res_def6$m0_hat)
  expect_identical(res_multtest$pi0_hat, res_def6$pi0_hat)

  # multtest's threshold is never stricter -> never fewer rejections
  expect_gte(sum(res_multtest$reject), sum(res_def6$reject))
})

test_that("fdr_bky's 'original' implementation handles r1 = 0 and r1 = m correctly", {
  # all p-values large -> stage 1 rejects nothing (r1 = 0)
  p_none <- c(0.6, 0.7, 0.8, 0.9, 0.95)
  res_none <- fdr_bky(p_none, q = 0.05, implementation = "original")
  expect_identical(res_none$r1, 0L)
  expect_identical(sum(res_none$reject), 0L)

  # all p-values tiny -> stage 1 rejects everything (r1 = m)
  p_all <- c(0.0001, 0.0002, 0.0003, 0.0004, 0.0005)
  res_all <- fdr_bky(p_all, q = 0.05, implementation = "original")
  expect_identical(res_all$r1, res_all$m)
  expect_true(all(res_all$reject))
})

test_that("fdr_bky's 'original' implementation preserves NA positions", {
  p <- c(0.01, NA, 0.2, 0.5, NA)
  res <- fdr_bky(p, implementation = "original")
  expect_identical(which(is.na(res$reject)), c(2L, 5L))
})

test_that("fdr_correction and workflow_tst() forward bky_implementation to fdr_bky() correctly", {
  r <- sim_trend_stack(nrow = 10, ncol = 10, n_time = 12, trend_strength = 0.3,
                        trend_fraction = 1, noise_sd = 0.3, seed = 1)$series
  trend <- trend_test(r, report = FALSE, verbose = FALSE)
  p_vals <- terra::values(trend$stats$p, mat = FALSE)

  fdr_multtest <- fdr_correction(trend$stats$p, bky_implementation = "multtest", report = FALSE, verbose = FALSE)
  fdr_def6 <- fdr_correction(trend$stats$p, bky_implementation = "original", report = FALSE, verbose = FALSE)
  expect_gte(sum(fdr_multtest$reject_BKY, na.rm = TRUE), sum(fdr_def6$reject_BKY, na.rm = TRUE))

  result <- workflow_tst(r, prewhiten = FALSE, bky_implementation = "original", report = FALSE, verbose = FALSE)
  direct <- fdr_bky(p_vals, implementation = "original")
  expect_identical(result$fdr$reject_BKY, direct$reject)
})
