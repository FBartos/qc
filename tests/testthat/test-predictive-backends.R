test_that("extract_predictive_samples works for integration method with conjugate prior", {
  set.seed(1)
  x   <- rnorm(50, mean = 10, sd = 2)
  uip <- create_prior_unit_information(x)

  fit     <- bpc(x, LSL = 2, target = 10, USL = 18, prior_mu = uip, method = "integration")
  samples <- extract_predictive_samples(fit, n_samples = 5000L)

  expect_true(is.numeric(samples))
  expect_length(samples, 5000L)
  expect_true(all(is.finite(samples)))
  # Predictive should be centred roughly near the data mean
  expect_lt(abs(mean(samples) - mean(x)), 1.0)
})

test_that("extract_predictive_samples works for mcmc method", {
  set.seed(2)
  x   <- rnorm(30, mean = 5, sd = 1)
  fit <- bpc(x, LSL = 1, target = 5, USL = 9,
             method = "mcmc", chains = 1, iter = 2000, warmup = 500, silent = TRUE, seed = 2)
  samples <- extract_predictive_samples(fit)

  expect_true(is.numeric(samples))
  expect_true(length(samples) > 0L)
  expect_true(all(is.finite(samples)))
  expect_lt(abs(mean(samples) - mean(x)), 1.0)
})

test_that("extract_predictive_samples: integration and MCMC posterior predictives agree", {
  # Both methods should produce predictive distributions with similar location
  # and spread. We allow generous tolerance due to Monte Carlo error.
  # Default (Jeffreys) priors are used because they are treated as conjugate NIG
  # by the integration backend, making extract_predictive_samples work for both.
  set.seed(3)
  n   <- 100
  x   <- rnorm(n, mean = 8, sd = 1.5)

  fit_int  <- bpc(x, LSL = 2, target = 8, USL = 14, method = "integration")
  fit_mcmc <- bpc(x, LSL = 2, target = 8, USL = 14,
                  method = "mcmc", chains = 2, iter = 4000, warmup = 1000,
                  silent = TRUE, seed = 3)

  s_int  <- extract_predictive_samples(fit_int,  n_samples = 20000L)
  s_mcmc <- extract_predictive_samples(fit_mcmc)

  # Means agree within 0.15 (generous MC tolerance, justified by n_samples = 20k)
  expect_lt(abs(mean(s_int) - mean(s_mcmc)), 0.15)
  # Standard deviations agree within 0.15
  expect_lt(abs(sd(s_int) - sd(s_mcmc)), 0.15)
  # 2.5% quantile agrees within 0.3
  expect_lt(abs(quantile(s_int, 0.025) - quantile(s_mcmc, 0.025)), 0.3)
  # 97.5% quantile agrees within 0.3
  expect_lt(abs(quantile(s_int, 0.975) - quantile(s_mcmc, 0.975)), 0.3)
})

test_that("extract_predictive_samples prior predictive via integration is finite", {
  set.seed(4)
  x   <- rnorm(20, mean = 0, sd = 1)
  uip <- create_prior_unit_information(x)

  fit_prior <- bpc(x, LSL = -3, target = 0, USL = 3, prior_mu = uip,
                   method = "integration", sample_priors = TRUE)
  samples <- extract_predictive_samples(fit_prior, n_samples = 2000L)

  expect_true(is.numeric(samples))
  expect_true(all(is.finite(samples)))
})

test_that("extract_predictive_samples errors informatively for non-conjugate integration", {
  set.seed(5)
  x   <- rnorm(30, mean = 5, sd = 1)
  fit <- bpc(x, LSL = 1, target = 5, USL = 9,
             prior_mu    = prior("normal", list(5, 2)),
             prior_sigma = prior("normal", list(1, 1), list(0, Inf)),
             method = "integration")

  expect_error(
    extract_predictive_samples(fit),
    regexp = "conjugate"
  )
})

test_that("extract_predictive_samples rejects improper conjugate prior-only states", {
  fit <- list(
    method = "integration",
    integration_result = list(
      prior = create_prior_conjugate(),
      cached_state = list(n = 0L, x_bar = 0, sse = 0)
    )
  )
  class(fit) <- "bpc"

  expect_error(
    extract_predictive_samples(fit, n_samples = 32L),
    regexp = "improper conjugate prior/posterior"
  )
})
