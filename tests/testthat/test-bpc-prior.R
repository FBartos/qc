test_that("DCSI resolves from specification limits and ignores target", {
  state <- qc:::.as_qc_suff_stats_state(data = rnorm(5))
  resolved <- qc:::.resolve_bpc_prior(
    prior = "DCSI",
    distribution = "normal",
    method = "integration",
    LSL = 0,
    USL = 10,
    cached_state = state
  )

  expect_s3_class(resolved$specification, "PriorDCSI")
  expect_s3_class(resolved$prior, "PriorConjugate")
  expect_equal(resolved$prior$mu0, 5)
  expect_equal(resolved$prior$k0, 2)
  expect_equal(resolved$prior$alpha0, 3)
  expect_equal(resolved$prior$beta0, 6 * (10 / (6 * 1.25))^2 / 2)
})

test_that("normal bpc defaults to integration with DCSI prior", {
  metric_names <- qc:::.qc_metric_names()
  mock_metrics <- stats::setNames(as.list(rep(1, length(metric_names))), metric_names)
  mock_coefficients <- stats::setNames(rep(1, length(metric_names)), metric_names)

  testthat::local_mocked_bindings(
    .bpc_fit_integration = function(distribution, data, LSL, USL, target,
                                    prior, sigma, sample_priors = FALSE,
                                    cached_state = NULL) {
      expect_identical(distribution, "normal")
      expect_s3_class(prior, "PriorConjugate")
      expect_equal(prior$mu0, 5)
      expect_equal(prior$k0, 2)
      expect_equal(prior$alpha0, 3)
      expect_equal(prior$beta0, 6 * (10 / (6 * 1.25))^2 / 2)
      expect_identical(sample_priors, FALSE)

      list(
        metrics = mock_metrics,
        coefficients = mock_coefficients,
        cached_state = cached_state
      )
    },
    .package = "qc"
  )

  fit <- bpc(c(1, 2, 3), LSL = 0, target = 3, USL = 10)

  expect_equal(fit$method, "integration")
  expect_s3_class(fit$prior, "PriorDCSI")
  expect_s3_class(fit$prior_resolved, "PriorConjugate")
  expect_equal(fit$prior_resolved$mu0, 5)
})

test_that("DCSI summary queries re-resolve prior from requested limits", {
  state <- qc:::.as_qc_suff_stats_state(data = c(1, 2, 3))
  object <- list(
    method = "integration",
    distribution = "normal",
    prior = prior_DCSI(),
    prior_resolved = qc:::.prior_DCSI(LSL = 0, USL = 10),
    integration_result = list(cached_state = state)
  )

  resolved <- qc:::.bpc_query_prior(
    object,
    limits = list(LSL = 2, target = 10, USL = 18)
  )

  expect_s3_class(resolved, "PriorConjugate")
  expect_equal(resolved$mu0, 10)
  expect_equal(resolved$beta0, 6 * (16 / (6 * 1.25))^2 / 2)
})

test_that("unsupported named arguments are rejected at the public boundary", {
  expect_error(
    bpc(
      c(1, 2, 3),
      LSL = 0, target = 2, USL = 4,
      unsupported = TRUE
    ),
    regexp = "Unsupported argument"
  )
})

test_that("omitted priors resolve to DCSI for every method and distribution", {
  state <- qc:::.as_qc_suff_stats_state(data = c(1, 2, 3))
  for (distribution in c("normal", "t")) {
    resolved <- qc:::.resolve_bpc_prior(
      prior = eval(formals(bpc)$prior),
      distribution = distribution,
      method = "mcmc",
      LSL = 0,
      USL = 4,
      cached_state = state
    )

    expect_s3_class(resolved$specification, "PriorDCSI")
    expect_s3_class(resolved$prior, "PriorConjugate")
    expect_s3_class(resolved$stan_prior_map$mu, "qc_stan_prior")
    expect_s3_class(resolved$stan_prior_map$sigma, "qc_stan_prior")
  }
  expect_identical(resolved$stan_prior_map$nu, "uniform_nu")
})

test_that("joint priors map to the integration density for Stan", {
  state <- qc:::.as_qc_suff_stats_state(data = c(1, 2, 3))
  resolve_mcmc <- function(prior) {
    qc:::.resolve_bpc_prior(
      prior = prior, distribution = "normal", method = "mcmc",
      LSL = 0, USL = 4, cached_state = state
    )$stan_prior_map
  }

  # NIG(mu0, k0, alpha0, beta0): N(mu0, sigma^2 / k0) * sigma^-(2 alpha0 + 1) * exp(-beta0 / sigma^2)
  dcsi <- resolve_mcmc("DCSI")
  expect_equal(dcsi$mu$type, 10)
  expect_equal(dcsi$mu$parameters[1:2], c(2, 2))
  expect_equal(dcsi$sigma$type, 11)
  expect_equal(dcsi$sigma$parameters[1:2], c(7, 3 * (4 / 7.5)^2))
  expect_true(dcsi$mu$proper && dcsi$sigma$proper)

  # the flat conjugate limit is the Jeffreys prior 1 / sigma
  flat <- resolve_mcmc(prior_conjugate())
  expect_equal(flat$mu$parameters[2], 0)
  expect_equal(flat$sigma$parameters[1:2], c(1, 0))
  expect_false(flat$mu$proper || flat$sigma$proper)

  # parameter-wise priors containing conjugate components follow the integration conversion
  semi_mu <- resolve_mcmc(prior_independent(
    mu = prior_conjugate(mu0 = 1, k0 = 2, alpha0 = 3, beta0 = 4),
    sigma = BayesTools::prior("gamma", list(2, 1))
  ))
  expect_equal(semi_mu$mu$parameters[1:2], c(1, 2))
  expect_identical(semi_mu$sigma, BayesTools::prior("gamma", list(2, 1)))

  semi_sigma <- resolve_mcmc(prior_independent(
    mu = BayesTools::prior("normal", list(1, 2)),
    sigma = prior_conjugate(mu0 = 0, k0 = 0, alpha0 = 3, beta0 = 4)
  ))
  expect_identical(semi_sigma$mu, BayesTools::prior("normal", list(1, 2)))
  expect_equal(semi_sigma$sigma$parameters[1:2], c(7, 4))

  # independent BayesTools priors are passed through unchanged
  independent <- resolve_mcmc("Jeffreys")
  expect_identical(independent$mu, "Jeffreys_mu")
  expect_identical(independent$sigma, "Jeffreys_sigma")

  expect_error(
    resolve_mcmc(prior_joint(function(mu, sigma) stats::dnorm(mu, log = TRUE) - log(sigma))),
    regexp = "bayestools_priors"
  )
})

test_that("MCMC with DCSI matches the analytic Normal-Inverse-Gamma posterior", {
  skip_on_cran()
  set.seed(1)
  x <- rnorm(20, mean = 1.5, sd = 0.6)
  LSL <- 0; USL <- 4
  n <- length(x)
  prior <- qc:::.prior_DCSI(LSL, USL)
  post <- qc:::.nig_posterior(prior, n, mean(x), sum((x - mean(x))^2))

  check_moments <- function(fit, mu_n, alpha_n, beta_n) {
    samples <- rstan::extract(fit$stanfit, pars = c("mu", "sigma"))
    sigma2 <- samples$sigma^2
    # Monte Carlo standard errors from the effective sample size
    ess <- min(rstan::summary(fit$stanfit, pars = c("mu", "sigma"))$summary[, "n_eff"])
    se_mu <- stats::sd(samples$mu) / sqrt(ess)
    se_sigma2 <- stats::sd(sigma2) / sqrt(ess)
    expect_lt(abs(mean(samples$mu) - mu_n), 5 * se_mu)
    expect_lt(abs(mean(sigma2) - beta_n / (alpha_n - 1)), 5 * se_sigma2)
  }

  fit <- bpc(
    x, LSL = LSL, target = 2, USL = USL, method = "mcmc",
    chains = 2, iter = 6000, warmup = 1000, seed = 1, silent = TRUE
  )
  expect_s3_class(fit$prior, "PriorDCSI")
  check_moments(fit, post$mu_n, post$alpha_n, post$beta_n)

  prior_fit <- bpc(
    x, LSL = LSL, target = 2, USL = USL, method = "mcmc", sample_priors = TRUE,
    chains = 2, iter = 6000, warmup = 1000, seed = 1, silent = TRUE
  )
  check_moments(prior_fit, prior$mu0, prior$alpha0, prior$beta0)
})

test_that("Student-t MCMC fits accept DCSI and joint priors", {
  skip_on_cran()
  set.seed(1)
  x <- rnorm(20, mean = 1.5, sd = 0.6)
  for (prior in list("DCSI", "unit_information", prior_conjugate(mu0 = 2, k0 = 1, alpha0 = 2, beta0 = 1))) {
    fit <- bpc(
      x, LSL = 0, target = 2, USL = 4, distribution = "t",
      prior = prior, chains = 1, iter = 1000, warmup = 500, seed = 1, silent = TRUE
    )
    expect_equal(fit$method, "mcmc")
    expect_true(all(is.finite(coef(fit))))
  }
})

test_that("prior sampling rejects improper parameter-wise priors for both methods", {
  x <- numeric(2L)

  expect_error(
    bpc(
      x,
      LSL = -1, USL = 1, target = 0,
      prior = "Jeffreys",
      sample_priors = TRUE,
      method = "mcmc",
      silent = TRUE
    ),
    regexp = "Improper prior distributions"
  )
  expect_error(
    bpc(
      x,
      LSL = -1, USL = 1, target = 0,
      prior = "Jeffreys",
      sample_priors = TRUE,
      method = "integration"
    ),
    regexp = "Improper prior distributions"
  )
})

test_that("prior sampling works for MCMC with proper independent priors", {
  x <- numeric(2L)

  fit <- bpc(
    x,
    LSL = -1, USL = 1, target = 0,
    prior = prior_independent(
      mu = prior("normal", list(0, 1)),
      sigma = prior("exp", list(1))
    ),
    sample_priors = TRUE,
    method = "mcmc",
    seed = 1,
    silent = TRUE,
    chains = 1,
    iter = 1000,
    warmup = 250
  )

  expect_s3_class(fit, "bpc")
  expect_equal(fit$method, "mcmc")
  expect_s3_class(fit$prior, "PriorIndependent")
  expect_true(length(fit$metrics$Cp) > 0L)
  expect_true(length(fit$metrics$Cpk) > 0L)
  expect_true(is.finite(stats::median(fit$metrics$Cp)))
  expect_true(is.finite(stats::median(fit$metrics$Cpk)))
})

test_that("prior-only generic integration with independent BayesTools priors remains available", {
  skip_if_not_slow_tests()

  x <- numeric(2L)
  fit <- bpc(
    x,
    LSL = -1, USL = 1, target = 0,
    prior = prior_independent(
      mu = prior("normal", list(0, 1)),
      sigma = prior("gamma", list(2, 1))
    ),
    sample_priors = TRUE,
    method = "integration"
  )

  expect_s3_class(fit, "bpc")
  expect_equal(fit$method, "integration")
  expect_true(is.finite(unname(fit$metrics$Cp["Median"])))
  expect_true(is.finite(unname(fit$metrics$Cpk["Median"])))
})

test_that("prior-only integration rejects improper direct generic priors", {
  skip_if_not_installed("cubature")

  improper_prior <- qc::prior_joint(function(mu, sigma) {
    ifelse(sigma <= 0, -Inf, 0)
  })

  expect_error(
    bpc(
      NULL,
      LSL = -1, USL = 1, target = 0,
      prior = improper_prior,
      sample_priors = TRUE,
      method = "integration"
    ),
    regexp = "improper.*arbitrary finite box"
  )
})

test_that("prior-only sampling rejects improper conjugate priors", {
  expect_error(
    bpc(
      NULL,
      LSL = -1, USL = 1, target = 0,
      prior = prior_conjugate(),
      sample_priors = TRUE,
      method = "integration"
    ),
    regexp = "Improper prior distributions cannot be sampled from with `sample_priors = TRUE`"
  )
})

test_that("proper conjugate priors remain valid for prior-only integration", {
  fit <- expect_no_error(
    bpc(
      NULL,
      LSL = -1, USL = 1, target = 0,
      prior = prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 1, beta0 = 1),
      sample_priors = TRUE,
      method = "integration"
    )
  )

  expect_s3_class(fit, "bpc")
  pred <- expect_no_error(extract_predictive_samples(fit, n_samples = 128L))
  expect_true(all(is.finite(pred)))
})

test_that("BPC prior maps drive Stan prior encoding", {
  mu_prior <- prior("normal", list(0, 1))
  sigma_prior <- prior("gamma", list(2, 1))
  extra_prior <- prior("normal", list(0, 1))
  prior_map <- list(mu = mu_prior, sigma = sigma_prior, nu = extra_prior)

  mapped <- qc:::.bpc_priors(
    "normal",
    prior_map = prior_map
  )

  expect_true("prior_type_mu" %in% names(mapped))
  expect_true("prior_type_sigma" %in% names(mapped))
  expect_false(any(grepl("_nu$", names(mapped))))
  expect_identical(mapped$sample_priors, 0)

  theta_prior <- prior("normal", list(0, 1))
  testthat::local_mocked_bindings(
    .qc_distribution_parameter_names = function(distribution, type = c("sample", "prior")) {
      type <- match.arg(type)
      if (identical(type, "prior")) {
        return(c("mu", "sigma", "theta"))
      }
      c("mu", "sigma")
    },
    .package = "qc"
  )

  extended <- qc:::.bpc_priors(
    "normal",
    prior_map = list(mu = mu_prior, sigma = sigma_prior, theta = theta_prior)
  )
  expect_true("prior_type_theta" %in% names(extended))
  expect_true("estimate_theta" %in% names(extended))
})

test_that("prior-map validation is generic across registered prior parameters", {
  proper_map <- list(
    mu = prior("normal", list(0, 1)),
    sigma = prior("gamma", list(2, 1)),
    nu = prior("exp", list(rate = 1 / 30), list(lower = 2, upper = Inf))
  )

  expect_silent(
    qc:::.validate_bpc_prior_configuration(
      method = "mcmc",
      distribution = "t",
      prior_info = list(prior_map = proper_map, prior = prior_independent(), kind = "independent"),
      sample_priors = TRUE
    )
  )

  improper_map <- proper_map
  improper_map$nu <- "uniform_nu"
  expect_error(
    qc:::.validate_bpc_prior_configuration(
      method = "mcmc",
      distribution = "t",
      prior_info = list(prior_map = improper_map, prior = prior_independent(), kind = "independent"),
      sample_priors = TRUE
    ),
    regexp = "prior\\$nu"
  )

  expect_error(
    bpc(
      c(1, 2, 3),
      LSL = 0, target = 2, USL = 4,
      method = "mcmc",
      prior = prior_independent(shape = prior("normal", list(0, 1)))
    ),
    regexp = "not a prior parameter"
  )
})

test_that("BPC fit objects preserve resolved unified priors", {
  expected_mu <- prior("normal", list(0, 1))
  expected_sigma <- prior("gamma", list(2, 1))
  metric_names <- qc:::.qc_metric_names()
  mock_metrics <- stats::setNames(as.list(rep(1, length(metric_names))), metric_names)
  mock_coefficients <- stats::setNames(rep(1, length(metric_names)), metric_names)

  testthat::local_mocked_bindings(
    .bpc_fit_integration = function(distribution, data, LSL, USL, target,
                                    prior, sigma,
                                    sample_priors = FALSE,
                                    cached_state = NULL) {
      expect_s3_class(prior, "PriorGeneric")
      expect_identical(prior$bayestools_priors$mu, expected_mu)
      expect_identical(prior$bayestools_priors$sigma, expected_sigma)
      list(
        metrics = mock_metrics,
        coefficients = mock_coefficients,
        cached_state = cached_state
      )
    },
    .package = "qc"
  )

  fit <- bpc(
    c(1, 2, 3),
    LSL = 0, target = 2, USL = 4,
    prior = prior_independent(mu = expected_mu, sigma = expected_sigma),
    method = "integration"
  )

  expect_s3_class(fit$prior, "PriorIndependent")
  expect_s3_class(fit$prior_resolved, "PriorGeneric")
  expect_identical(fit$prior_map, list(mu = expected_mu, sigma = expected_sigma))
})
