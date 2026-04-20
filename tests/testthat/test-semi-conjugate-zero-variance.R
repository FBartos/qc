semi_sigma_generic_reference <- function(mu_prior, alpha0, beta0) {
  log_prior_mu <- qc:::.make_prior_log_dens_fn(mu_prior)
  log_prior_sigma <- function(sigma) {
    ifelse(
      sigma <= 0,
      -Inf,
      log(2) + alpha0 * log(beta0) - lgamma(alpha0) -
        (2 * alpha0 + 1) * log(sigma) - beta0 / sigma^2
    )
  }

  qc:::create_prior_generic(function(mu, sigma) {
    log_prior_mu(mu) + log_prior_sigma(sigma)
  })
}

test_that("constant-data Jeffreys semi-conjugate sigma posterior is rejected explicitly", {
  skip_if_not_installed("BayesTools")

  expect_error(
    qc::bpc(
      rep(5, 10),
      LSL = 0, USL = 10, target = 5,
      prior_mu = BayesTools::prior("normal", list(5, 2)),
      prior_sigma = "Jeffreys_sigma",
      method = "integration"
    ),
    regexp = "improper because beta_n reaches zero inside the mu support"
  )
})

test_that("constant-data semi-conjugate sigma moments match an equivalent generic prior", {
  skip_if_not_installed("BayesTools")

  x <- rep(5, 10)
  mu_prior <- BayesTools::prior("normal", list(5, 2))
  sigma_prior <- qc:::create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 2, beta0 = 1)

  semi_prior <- qc:::.bayestools_to_integration_prior(mu_prior, sigma_prior)$prior
  generic_prior <- semi_sigma_generic_reference(mu_prior, alpha0 = 2, beta0 = 1)

  semi_cp <- qc:::compute_metric_moments(
    data = x, LSL = 0, USL = 10, prior = semi_prior,
    metric = "Cp", target = 5
  )
  generic_cp <- qc:::compute_metric_moments(
    data = x, LSL = 0, USL = 10, prior = generic_prior,
    metric = "Cp", target = 5
  )
  semi_cpk <- qc:::compute_metric_moments(
    data = x, LSL = 0, USL = 10, prior = semi_prior,
    metric = "Cpk", target = 5
  )
  generic_cpk <- qc:::compute_metric_moments(
    data = x, LSL = 0, USL = 10, prior = generic_prior,
    metric = "Cpk", target = 5
  )

  expect_equal(semi_cp$mean, generic_cp$mean, tolerance = 0.05)
  expect_equal(semi_cp$sd, generic_cp$sd, tolerance = 0.05)
  expect_equal(semi_cpk$mean, generic_cpk$mean, tolerance = 0.05)
  expect_equal(semi_cpk$sd, generic_cpk$sd, tolerance = 0.05)

  fit <- expect_no_error(
    qc::bpc(
      x,
      LSL = 0, USL = 10, target = 5,
      prior_mu = mu_prior,
      prior_sigma = sigma_prior,
      method = "integration"
    )
  )
  expect_true(all(is.finite(fit$coefficients[c("Cp", "Cpk")])))
  expect_gt(unname(fit$coefficients["Cp"]), 0)
  expect_gt(unname(fit$coefficients["Cpk"]), 0)
})

test_that("prior-only semi-conjugate sigma uses prior-centered domains", {
  skip_if_not_installed("BayesTools")

  mu_prior <- BayesTools::prior("normal", list(5, 2))
  sigma_prior <- qc:::create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 2, beta0 = 1)

  semi_prior <- qc:::.bayestools_to_integration_prior(mu_prior, sigma_prior)$prior
  generic_prior <- semi_sigma_generic_reference(mu_prior, alpha0 = 2, beta0 = 1)

  semi_cp <- qc:::compute_metric_moments(
    data = numeric(0), LSL = 0, USL = 10, prior = semi_prior,
    metric = "Cp", target = 5
  )
  generic_cp <- qc:::compute_metric_moments(
    data = numeric(0), LSL = 0, USL = 10, prior = generic_prior,
    metric = "Cp", target = 5
  )

  expect_equal(semi_cp$mean, generic_cp$mean, tolerance = 0.05)
  expect_equal(semi_cp$sd, generic_cp$sd, tolerance = 0.05)

  fit <- expect_no_error(
    qc::bpc(
      NULL,
      LSL = 0, USL = 10, target = 5,
      prior_mu = mu_prior,
      prior_sigma = sigma_prior,
      sample_priors = TRUE,
      method = "integration"
    )
  )
  expect_true(is.finite(unname(fit$coefficients["Cp"])))
  expect_gt(unname(fit$coefficients["Cp"]), 0)
})
