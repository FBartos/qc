test_that("semi-conjugate mu moments honor sigma support far from the data scale", {
  skip_if_not_installed("BayesTools")

  set.seed(1)
  x <- rnorm(10)
  sigma_prior <- BayesTools::prior(
    "uniform",
    parameters = list(a = 1000, b = 1001)
  )
  prior <- qc:::.bayestools_to_integration_prior("Jeffreys_mu", sigma_prior)$prior

  cp_moments <- qc:::compute_metric_moments(
    x, -3, 3, prior,
    metric = "Cp",
    target = 0
  )

  expect_true(all(is.finite(unlist(cp_moments))))
  expect_gt(cp_moments$mean, 0)

  fit <- qc::bpc(
    x,
    LSL = -3,
    target = 0,
    USL = 3,
    method = "integration",
    prior = qc::prior_independent(mu = "Jeffreys_mu", sigma = sigma_prior)
  )

  expect_true(is.finite(fit$coefficients[["Cp"]]))
  expect_gt(fit$coefficients[["Cp"]], 0)
})

test_that("semi-conjugate mu rejects sigma tails that cannot be normalized", {
  prior <- qc:::create_prior_semi_mu(
    mu0 = 0,
    k0 = 1,
    log_dens_sigma = function(sigma) rep(0, length(sigma))
  )

  expect_error(
    qc:::make_solver(numeric(0), -1, 1, prior, metric = "Cp"),
    "semi-conjugate mu posterior.*improper"
  )
  expect_error(
    qc:::make_density_solver(numeric(0), -1, 1, prior, metric = "Cp"),
    "semi-conjugate mu posterior.*improper"
  )
})

test_that("semi-conjugate sigma with alpha0 = 0 keeps moments and summaries finite", {
  skip_if_not_installed("BayesTools")

  set.seed(2)
  x <- rnorm(12, mean = 0.3, sd = 0.8)
  prior <- qc:::.bayestools_to_integration_prior(
    BayesTools::prior("normal", parameters = list(mean = 0, sd = 1)),
    qc:::create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 0, beta0 = 1)
  )$prior

  moment_result <- qc:::compute_metric_moments(
    x, -3, 3, prior,
    metric = "Cp"
  )
  density_result <- qc:::analyze_capability_integration(
    x, -3, 3, prior,
    metric = "Cp",
    use_density_solver = TRUE
  )
  survival_result <- qc:::analyze_capability_integration(
    x, -3, 3, prior,
    metric = "Cp",
    use_density_solver = FALSE
  )

  expect_true(all(is.finite(unlist(moment_result))))
  expect_gt(moment_result$mean, 0)
  expect_true(all(is.finite(unname(density_result$stats[c("Mean", "SD")]))))
  expect_true(all(is.finite(unname(survival_result$stats[c("Mean", "SD")]))))
  expect_equal(
    moment_result$mean,
    unname(density_result$stats["Mean"]),
    tolerance = 0.03
  )
  expect_equal(
    moment_result$mean,
    unname(survival_result$stats["Mean"]),
    tolerance = 0.03
  )
})
