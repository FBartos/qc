test_that("extract_predictive_samples works for integration method with conjugate prior", {
  set.seed(1)
  x   <- rnorm(50, mean = 10, sd = 2)
  uip <- qc:::create_prior_unit_information(x)

  fit     <- bpc(x, LSL = 2, target = 10, USL = 18, prior = uip, method = "integration")
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
             method = "mcmc", prior = "Jeffreys",
             chains = 1, iter = 1000, warmup = 250, silent = TRUE, seed = 2)
  samples <- extract_predictive_samples(fit)

  expect_true(is.numeric(samples))
  expect_true(length(samples) > 0L)
  expect_true(all(is.finite(samples)))
  expect_lt(abs(mean(samples) - mean(x)), 1.0)
})

test_that("extract_predictive_samples: integration and MCMC posterior predictives agree", {
  set.seed(3)
  n   <- 60
  x   <- rnorm(n, mean = 8, sd = 1.5)

  fit_int  <- bpc(x, LSL = 2, target = 8, USL = 14, method = "integration")
  fit_mcmc <- bpc(x, LSL = 2, target = 8, USL = 14,
                  method = "mcmc", prior = "Jeffreys",
                  chains = 1, iter = 1500, warmup = 500,
                  silent = TRUE, seed = 3)

  s_int  <- extract_predictive_samples(fit_int,  n_samples = 5000L)
  s_mcmc <- extract_predictive_samples(fit_mcmc)

  expect_predictive_samples_close(
    s_int,
    s_mcmc,
    mean_sd_factor = 0.05,
    quantile_sd_factor = 0.12
  )
})

test_that("extract_predictive_samples prior predictive via integration is finite", {
  set.seed(4)
  x   <- rnorm(20, mean = 0, sd = 1)
  uip <- qc:::create_prior_unit_information(x)

  fit_prior <- bpc(x, LSL = -3, target = 0, USL = 3, prior = uip,
                   method = "integration", sample_priors = TRUE)
  samples <- extract_predictive_samples(fit_prior, n_samples = 2000L)

  expect_true(is.numeric(samples))
  expect_true(all(is.finite(samples)))
})

test_that("extract_predictive_samples errors informatively for non-conjugate integration", {
  fit <- list(
    method = "integration",
    distribution = "normal",
    integration_result = list(
      distribution = "normal",
      prior = qc:::create_prior_generic(function(mu, sigma) 0),
      cached_state = list(n = 1L, x_bar = 0, sse = 1)
    )
  )
  class(fit) <- "bpc"

  expect_error(
    extract_predictive_samples(fit),
    regexp = "conjugate"
  )
})

test_that("extract_predictive_samples rejects non-normal integration fits before drawing normal predictives", {
  fit <- list(
    method = "integration",
    distribution = "qc_future_distribution",
    integration_result = list(
      distribution = "qc_future_distribution",
      prior = qc:::create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 2, beta0 = 1),
      cached_state = list(n = 2L, x_bar = 0, sse = 1)
    )
  )
  class(fit) <- "bpc"

  expect_error(
    extract_predictive_samples(fit, n_samples = 32L),
    regexp = "only supported for `distribution = \"normal\"`"
  )
})

test_that("extract_predictive_samples rejects improper conjugate prior-only states", {
  fit <- list(
    method = "integration",
    distribution = "normal",
    integration_result = list(
      distribution = "normal",
      prior = qc:::create_prior_conjugate(),
      cached_state = list(n = 0L, x_bar = 0, sse = 0)
    )
  )
  class(fit) <- "bpc"

  expect_error(
    extract_predictive_samples(fit, n_samples = 32L),
    regexp = "improper"
  )
})
