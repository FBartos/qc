testthat::test_that("compute metrics", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(x, LSL = 2, USL = 18, target = 10, chains = 1, warmup = 100, iter = 200, silent = TRUE, seed = 1)

  metrics <- qc:::.compute_capability_metrics(fit, LSL = 2, USL = 18, target = 10)

  # Check that metrics were computed
  expect_true(length(metrics) > 0)

  # Test with t-distribution
  fit2 <- bpc(x, LSL = 2, USL = 18, target = 10, chains = 1, warmup = 1000, iter = 2000, silent = TRUE, seed = 1,
              prior_nu = prior("exp", list(1)))

  # Check that fit2 has metrics
  expect_true("metrics" %in% names(fit2))

})

testthat::test_that("plot posterior density of metrics", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(x, LSL = 2, USL = 18, target = 10, chains = 1, warmup = 100, iter = 200, silent = TRUE, seed = 1)

  plot_density(fit)
  plot_density(fit, what = "Cp")

  fit_t <- bpc(x, LSL = 2, USL = 18, target = 10, chains = 1, warmup = 1000, iter = 2000, silent = TRUE, seed = 1,
               distribution = "t",
             prior_nu = prior("exp",  list(1)))

  plot_density(fit_t) # TODO: fails with error? because of NaNs in one of the metrics?
  plot_density(fit_t, what = "Cp")

})

