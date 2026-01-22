testthat::test_that("integration method returns valid bpc object", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(x, LSL = 2, target = 10, USL = 18, method = "integration")

  # Check object structure
  expect_s3_class(fit, "bpc")
  expect_equal(fit$method, "integration")
  expect_true(!is.null(fit$integration_result))
  expect_true(!is.null(fit$metrics))
  expect_true(!is.null(fit$coefficients))

  # Check all metrics are present and positive
  expected_metrics <- c("Cp", "Cpk", "Cpm", "CpU", "CpL", "Cpc")
  expect_equal(names(fit$coefficients), expected_metrics)
  expect_true(all(fit$coefficients > 0))
  expect_true(all(is.finite(fit$coefficients)))
})


testthat::test_that("integration method print and summary work", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(x, LSL = 2, target = 10, USL = 18, method = "integration")

  # print should work
  expect_output(print(fit), "Bayesian Process Capability")

  # summary should work
  ss <- summary(fit)
  expect_s3_class(ss, "bpc_summary")
  expect_true(!is.null(ss$summary))
  expect_true(!is.null(ss$interval_summary))

  # print summary should work
  expect_output(print(ss), "Bayesian Process Capability")
})


testthat::test_that("integration method errors for t-distribution", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  expect_error(
    bpc(x, LSL = 2, target = 10, USL = 18, method = "integration", distribution = "t"),
    "integration method currently only supports"
  )
})


testthat::test_that("extract_samples errors for integration method", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(x, LSL = 2, target = 10, USL = 18, method = "integration")

  expect_error(
    qc:::extract_samples(fit),
    "extract_samples.*is not supported for integration method"
  )
})


testthat::test_that("integration vs MCMC agreement with Jeffreys prior", {

  skip_on_cran()  # Skip on CRAN due to long runtime

  set.seed(42)
  x <- rnorm(30, 50, 0.5)
  LSL <- 44
  USL <- 56
  target <- 50

  # Fit with both methods
  fit_int <- bpc(x, LSL = LSL, target = target, USL = USL, method = "integration")
  fit_mcmc <- bpc(x, LSL = LSL, target = target, USL = USL, method = "mcmc",
                  iter = 50000, chains = 4, silent = TRUE, seed = 42)

  # Compare posterior means - should be close
  coef_int <- fit_int$coefficients
  coef_mcmc <- fit_mcmc$coefficients

  # Check agreement within tolerance (allow 10% relative error or 0.05 absolute error)
  for (metric in names(coef_int)) {
    rel_diff <- abs(coef_int[metric] - coef_mcmc[metric]) / max(abs(coef_mcmc[metric]), 0.01)
    abs_diff <- abs(coef_int[metric] - coef_mcmc[metric])
    expect_true(
      rel_diff < 0.10 | abs_diff < 0.05,
      info = sprintf("%s: integration=%f, mcmc=%f, rel_diff=%f, abs_diff=%f",
                     metric, coef_int[metric], coef_mcmc[metric], rel_diff, abs_diff)
    )
  }
})


testthat::test_that("integration interval probabilities agree with MCMC", {

  skip_on_cran()  # Skip on CRAN due to long runtime

  set.seed(42)
  x <- rnorm(30, 50, 0.5)
  LSL <- 44
  USL <- 56
  target <- 50
  bounds <- c(1.0, 2.0)

  # Fit with both methods
  fit_int <- bpc(x, LSL = LSL, target = target, USL = USL, method = "integration")
  fit_mcmc <- bpc(x, LSL = LSL, target = target, USL = USL, method = "mcmc",
                  iter = 50000, chains = 4, silent = TRUE, seed = 42)

  # Compute interval probabilities for MCMC
  mcmc_probs <- sapply(fit_mcmc$metrics, function(values) {
    mean(values > bounds[1] & values < bounds[2])
  })

  # Get integration prior and compute probabilities
  prior_info <- qc:::.bayestools_to_integration_prior("Jeffreys_mu", "Jeffreys_sigma")
  prior <- prior_info$prior
  cached_state <- fit_int$integration_result$cached_state

  int_probs <- sapply(names(fit_int$metrics), function(m) {
    qc:::compute_cpk_prob_integration(x, LSL, USL, bounds, prior, metric = m,
                                       target = target, cached_state = cached_state)
  })

  # Check agreement within tolerance (allow 0.03 absolute difference)
  for (metric in names(int_probs)) {
    abs_diff <- abs(int_probs[metric] - mcmc_probs[metric])
    expect_true(
      abs_diff < 0.03,
      info = sprintf("%s: int_prob=%f, mcmc_prob=%f, abs_diff=%f",
                     metric, int_probs[metric], mcmc_probs[metric], abs_diff)
    )
  }
})


testthat::test_that("integration method with custom priors", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  # Test with custom normal prior on mu and gamma on sigma
  fit <- bpc(x, LSL = 2, target = 10, USL = 18,
             method = "integration",
             prior_mu = prior("normal", list(10, 5)),
             prior_sigma = prior("gamma", list(2, 1)))

  # Should still return valid results
  expect_s3_class(fit, "bpc")
  expect_equal(fit$method, "integration")
  expect_true(all(fit$coefficients > 0))
  expect_true(all(is.finite(fit$coefficients)))

  # Non-conjugate so should be slower but still work
  expect_false(fit$integration_result$is_conjugate)
})


testthat::test_that("integration vs MCMC agreement with non-conjugate priors", {

  skip_on_cran()  # Skip on CRAN due to long runtime

  set.seed(123)
  x <- rnorm(30, 10, 2)
  LSL <- 2
  USL <- 18
  target <- 10

  # Use informative normal prior on mu and gamma prior on sigma
  prior_mu <- prior("normal", list(10, 5))
  prior_sigma <- prior("gamma", list(2, 1))

  # Fit with both methods
  fit_int <- bpc(x, LSL = LSL, target = target, USL = USL, method = "integration",
                 prior_mu = prior_mu, prior_sigma = prior_sigma)
  fit_mcmc <- bpc(x, LSL = LSL, target = target, USL = USL, method = "mcmc",
                  prior_mu = prior_mu, prior_sigma = prior_sigma,
                  iter = 50000, chains = 4, silent = TRUE, seed = 123)

  # Compare posterior means - allow slightly larger tolerance for non-conjugate
  coef_int <- fit_int$coefficients
  coef_mcmc <- fit_mcmc$coefficients

  for (metric in names(coef_int)) {
    rel_diff <- abs(coef_int[metric] - coef_mcmc[metric]) / max(abs(coef_mcmc[metric]), 0.01)
    abs_diff <- abs(coef_int[metric] - coef_mcmc[metric])
    expect_true(
      rel_diff < 0.15 | abs_diff < 0.10,
      info = sprintf("%s: integration=%f, mcmc=%f, rel_diff=%f, abs_diff=%f",
                     metric, coef_int[metric], coef_mcmc[metric], rel_diff, abs_diff)
    )
  }
})


testthat::test_that("integration interval probabilities with non-conjugate priors agree with MCMC", {

  skip_on_cran()  # Skip on CRAN due to long runtime

  set.seed(123)
  x <- rnorm(30, 10, 2)
  LSL <- 2
  USL <- 18
  target <- 10
  bounds <- c(0.5, 1.5)

  prior_mu <- prior("normal", list(10, 5))
  prior_sigma <- prior("gamma", list(2, 1))

  # Fit with both methods
  fit_int <- bpc(x, LSL = LSL, target = target, USL = USL, method = "integration",
                 prior_mu = prior_mu, prior_sigma = prior_sigma)
  fit_mcmc <- bpc(x, LSL = LSL, target = target, USL = USL, method = "mcmc",
                  prior_mu = prior_mu, prior_sigma = prior_sigma,
                  iter = 50000, chains = 4, silent = TRUE, seed = 123)

  # Compute interval probabilities for MCMC
  mcmc_probs <- sapply(fit_mcmc$metrics, function(values) {
    mean(values > bounds[1] & values < bounds[2])
  })

  # Get integration prior and compute probabilities
  prior_info <- qc:::.bayestools_to_integration_prior(prior_mu, prior_sigma)
  prior <- prior_info$prior
  cached_state <- fit_int$integration_result$cached_state

  int_probs <- sapply(names(fit_int$metrics), function(m) {
    qc:::compute_cpk_prob_integration(x, LSL, USL, bounds, prior, metric = m,
                                       target = target, cached_state = cached_state)
  })

  # Check agreement within tolerance (allow 0.05 for non-conjugate)
  for (metric in names(int_probs)) {
    abs_diff <- abs(int_probs[metric] - mcmc_probs[metric])
    expect_true(
      abs_diff < 0.05,
      info = sprintf("%s: int_prob=%f, mcmc_prob=%f, abs_diff=%f",
                     metric, int_probs[metric], mcmc_probs[metric], abs_diff)
    )
  }
})


testthat::test_that("integration method speed advantage", {

  skip_on_cran()

  set.seed(1)
  x <- rnorm(50, 10, 2)

  # Time integration method
  time_int <- system.time({
    fit_int <- bpc(x, LSL = 2, target = 10, USL = 18, method = "integration")
  })["elapsed"]

  # Time MCMC method (with minimal iterations)
  time_mcmc <- system.time({
    fit_mcmc <- bpc(x, LSL = 2, target = 10, USL = 18, method = "mcmc",
                    chains = 1, iter = 1000, warmup = 500, silent = TRUE)
  })["elapsed"]

  # Integration should be faster than even minimal MCMC
  expect_true(time_int < time_mcmc,
              info = sprintf("Integration: %f sec, MCMC: %f sec", time_int, time_mcmc))
})
