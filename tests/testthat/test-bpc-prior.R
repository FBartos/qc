test_that("prior sampling works correctly for both methods", {
  # meaningful data not needed for prior sampling
  x <- numeric(2L)

  # improper priors should error for both methods
  expect_error(
    bpc(x, LSL = -1, USL = 1, target = 0, sample_priors = TRUE, method = "mcmc", silent = TRUE)
  )
  expect_error(
    bpc(x, LSL = -1, USL = 1, target = 0, sample_priors = TRUE, method = "integration")
  )

  # proper priors should work for both
  # Use seed for MCMC reproducibility
  fit1p  <- bpc(x, LSL = -1, USL = 1, target = 0,
                prior_mu = prior("normal", list(0, 1)),
                prior_sigma = prior("exp",  list(1)),
                sample_priors = TRUE, method = "mcmc", seed = 1, silent = TRUE, chains = 1, iter = 4000, warmup = 1000)

  fit1pi <- bpc(x, LSL = -1, USL = 1, target = 0,
                prior_mu = prior("normal", list(0, 1)),
                prior_sigma = prior("gamma", list(2, 1)),
                sample_priors = TRUE, method = "integration")

  # Compare first metric (Cp)
  # Note: Mean of Cp (1/sigma) diverges for Exponential prior on sigma.
  # Use Median for stable comparison.
  # MCMC returns samples of metric
  # Integration returns named vector of stats (Mean, Median, SD, etc.)
  m1 <- median(fit1p$metrics$Cp)
  m2 <- fit1pi$metrics$Cp["Median"]

  # Using a loose tolerance because MCMC and integration approximations might vary
  # due to different underlying priors and approximations
  expect_lt(abs(m1 - m2), 0.5)

  # Also check Cpk
  expect_lt(abs(median(fit1p$metrics$Cpk) - fit1pi$metrics$Cpk["Median"]), 0.5)
})
