test_that("DCSI resolves from specification limits and ignores target", {
  state <- qc:::.as_qc_suff_stats_state(data = rnorm(5))
  resolved <- qc:::.resolve_bpc_prior(
    prior = "DCSI",
    prior_missing = FALSE,
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

test_that("DCSI and other joint normal priors are integration-only until Stan supports joint priors", {
  expect_error(
    bpc(
      c(1, 2, 3),
      LSL = 0, target = 2, USL = 4,
      method = "mcmc",
      prior = "DCSI",
      chains = 1,
      iter = 100,
      warmup = 50,
      silent = TRUE
    ),
    regexp = "DCSI.*integration"
  )

  expect_error(
    bpc(
      c(1, 2, 3),
      LSL = 0, target = 2, USL = 4,
      method = "mcmc",
      prior = prior_conjugate(mu0 = 2, k0 = 1, alpha0 = 2, beta0 = 1),
      chains = 1,
      iter = 100,
      warmup = 50,
      silent = TRUE
    ),
    regexp = "Joint normal likelihood priors"
  )
})

test_that("omitted priors resolve to Jeffreys for explicit MCMC", {
  state <- qc:::.as_qc_suff_stats_state(data = c(1, 2, 3))
  resolved <- qc:::.resolve_bpc_prior(
    prior = "DCSI",
    prior_missing = TRUE,
    distribution = "normal",
    method = "mcmc",
    LSL = 0,
    USL = 4,
    cached_state = state
  )

  expect_s3_class(resolved$specification, "PriorIndependent")
  expect_identical(resolved$prior_map$mu, "Jeffreys_mu")
  expect_identical(resolved$prior_map$sigma, "Jeffreys_sigma")
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
