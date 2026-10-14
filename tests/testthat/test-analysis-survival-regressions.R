test_that("semi-conjugate sigma survival ignores optional mu metadata", {
  skip_if_not_installed("BayesTools")

  log_dens_mu <- function(mu) stats::dnorm(mu, mean = 0, sd = 20, log = TRUE)
  mu_prior <- BayesTools::prior("normal", list(0, 20))
  prior_plain <- qc:::create_prior_semi_sigma(
    alpha0 = 2,
    beta0 = 1,
    log_dens_mu = log_dens_mu
  )
  prior_meta <- qc:::create_prior_semi_sigma(
    alpha0 = 2,
    beta0 = 1,
    log_dens_mu = log_dens_mu,
    bayestools_priors = list(mu = mu_prior)
  )

  solver_plain <- qc:::make_solver(
    numeric(0), 0, 10, prior_plain,
    metric = "Cpu", target = 5
  )
  solver_meta <- qc:::make_solver(
    numeric(0), 0, 10, prior_meta,
    metric = "Cpu", target = 5
  )

  expect_equal(solver_plain(0.5), solver_meta(0.5), tolerance = 1e-8)
  expect_equal(solver_plain(10), solver_meta(10), tolerance = 1e-8)
})

test_that("semi-conjugate sigma moments ignore optional mu metadata", {
  skip_if_not_installed("BayesTools")

  log_dens_mu <- function(mu) stats::dnorm(mu, mean = 0, sd = 20, log = TRUE)
  mu_prior <- BayesTools::prior("normal", list(0, 20))
  prior_plain <- qc:::create_prior_semi_sigma(
    alpha0 = 2,
    beta0 = 1,
    log_dens_mu = log_dens_mu
  )
  prior_meta <- qc:::create_prior_semi_sigma(
    alpha0 = 2,
    beta0 = 1,
    log_dens_mu = log_dens_mu,
    bayestools_priors = list(mu = mu_prior)
  )

  moments_plain <- qc:::compute_metric_moments(
    numeric(0), 0, 10, prior_plain,
    metric = "Cpu", target = 5
  )
  moments_meta <- qc:::compute_metric_moments(
    numeric(0), 0, 10, prior_meta,
    metric = "Cpu", target = 5
  )

  expect_equal(moments_plain$mean, moments_meta$mean, tolerance = 1e-8)
  expect_equal(moments_plain$sd, moments_meta$sd, tolerance = 1e-8)
})

test_that("prior-only Monte Carlo analysis accepts PriorConjugate sigma metadata", {
  skip_if_not_installed("BayesTools")
  skip_if_not_installed("cubature")

  mu_prior <- BayesTools::prior("normal", list(0, 1))
  sigma_prior <- qc:::create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 2, beta0 = 1)
  log_dens_mu <- qc:::.make_prior_log_dens_fn(mu_prior)
  log_dens_sigma <- function(sigma) {
    stats::dgamma(1 / sigma^2, shape = 2, rate = 1, log = TRUE) +
      log(2) - 3 * log(sigma)
  }

  prior <- qc:::create_prior_generic(
    function(mu, sigma) log_dens_mu(mu) + log_dens_sigma(sigma),
    bayestools_priors = list(mu = mu_prior, sigma = sigma_prior)
  )

  mu_samples <- c(-0.2, 0.1, 0.4)
  sig_samples <- c(0.3, 0.5, 0.9)
  result <- qc:::analyze_capability_integration(
    numeric(0), -1, 1, prior,
    metric = "Cp",
    target = 0,
    mc_samples = list(mu = mu_samples, sig = sig_samples)
  )

  expected <- 1 / (3 * sig_samples)

  expect_equal(result$samples, expected)
  expect_equal(unname(result$stats["Mean"]), mean(expected))
  expect_equal(unname(result$stats["Median"]), stats::median(expected))
})

test_that("analysis falls back to survival when the density solver errors", {
  fake_survival <- function(c) pmax(0, 1 - c / 10)

  testthat::local_mocked_bindings(
    make_density_solver = function(...) stop("density boom", call. = FALSE),
    make_solver = function(...) fake_survival,
    .package = "qc"
  )

  result <- expect_no_error(
    qc:::analyze_capability_integration(
      numeric(0), 0, 1,
      qc:::create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 2, beta0 = 1),
      metric = "Cp",
      target = 0,
      use_density_solver = TRUE
    )
  )

  expect_true(all(is.finite(unname(result$stats[c("Mean", "Median", "SD")]))))
  expect_gt(result$area, 0)
})

test_that("interval probabilities fall back to survival when the density solver errors", {
  fake_survival <- function(c) pmax(0, 1 - c / 10)

  testthat::local_mocked_bindings(
    make_density_solver = function(...) stop("density boom", call. = FALSE),
    make_solver = function(...) fake_survival,
    .package = "qc"
  )

  prob <- expect_no_error(
    qc:::compute_cpk_prob_integration(
      numeric(0), 0, 1, c(1, 2),
      qc:::create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 2, beta0 = 1),
      metric = "Cp",
      target = 0
    )
  )

  expect_equal(prob, fake_survival(1) - fake_survival(2))
})

test_that("divergent negative prior-only means preserve the negative sign", {
  skip_if_not_installed("BayesTools")

  prior <- qc:::.bayestools_to_integration_prior(
    BayesTools::prior("normal", list(2, 0.1)),
    BayesTools::prior("exp", list(1))
  )$prior

  result <- qc:::analyze_capability_integration(
    numeric(0), -1, 1, prior,
    metric = "Cpu",
    target = 0,
    mc_samples = list(
      mu = c(1.4, 1.5, 1.6),
      sig = c(0.1, 0.2, 0.3)
    ),
    divergence_info = list(
      mean_divergent = TRUE,
      sd_divergent = TRUE,
      alpha = 1,
      reason = "test"
    )
  )

  expect_true(all(result$samples < 0))
  expect_identical(unname(result$stats["Mean"]), -Inf)
  expect_identical(unname(result$stats["SD"]), Inf)
})

test_that("heavy-tailed Cpm and Cpc survival agrees with independent importance sampling", {
  skip_if_not_installed("BayesTools")

  x <- c(0.05, 0.12)
  prior <- qc:::.bayestools_to_integration_prior(
    BayesTools::prior("cauchy", list(0, 5)),
    qc:::create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 1, beta0 = 1)
  )$prior

  # Under a flat mu prior, this two-observation posterior has precision
  # Gamma(1.5, 1.001225) and mu|sigma ~ Normal(0.085, sigma/sqrt(2)).
  # Draw from that proposal and weight by the specified Cauchy mu prior.
  # No package posterior, constraint, or solver helpers enter this reference.
  set.seed(20261014)
  sigma <- 1 / sqrt(stats::rgamma(500000, shape = 1.5, rate = 1.001225))
  mu <- stats::rnorm(length(sigma), 0.085, sigma / sqrt(2))
  weight <- stats::dcauchy(mu, 0, 5)
  weight <- weight / sum(weight)
  mean_abs <- 2 * sigma * stats::dnorm(mu / sigma) +
    mu * (2 * stats::pnorm(mu / sigma) - 1)
  indices <- list(
    Cpm = 1 / (3 * sqrt(sigma^2 + mu^2)),
    Cpc = 1 / (3 * sqrt(pi / 2) * mean_abs)
  )

  for (metric in c("Cpm", "Cpc")) {
    event <- indices[[metric]] > 0.5
    ref_prob <- sum(weight * event)
    mc_se <- sqrt(sum(weight^2 * (event - ref_prob)^2))
    solver <- qc:::make_solver(
      x, -1, 1, prior,
      metric = metric,
      target = 0
    )

    expect_lt(mc_se, 0.001)
    expect_true(abs(solver(0.5) - ref_prob) < 5 * mc_se, info = metric)
  }
})

test_that("semi-conjugate sigma Cpm and Cpc survival stay monotone at tiny positive thresholds", {
  prior <- qc:::create_prior_semi_sigma(
    alpha0 = 2,
    beta0 = 1,
    log_dens_mu = function(mu) stats::dnorm(mu, 0, 1, log = TRUE)
  )

  for (metric in c("Cpm", "Cpc")) {
    solver <- qc:::make_solver(
      numeric(0), -1, 1, prior,
      metric = metric,
      target = 0
    )
    probs <- vapply(c(0, 1e-5, 1e-4, 1e-3), solver, numeric(1))

    expect_equal(probs[1], 1, tolerance = 1e-12, info = metric)
    expect_true(all(diff(probs) <= 1e-10), info = metric)
    expect_true(all(probs >= -1e-12 & probs <= 1 + 1e-12), info = metric)
  }
})

test_that("semi-conjugate sigma solvers default target to the midpoint for Cpm and Cpc", {
  prior <- qc:::create_prior_semi_sigma(
    alpha0 = 2,
    beta0 = 1,
    log_dens_mu = function(mu) stats::dnorm(mu, 0, 1, log = TRUE)
  )

  midpoint <- 0

  for (metric in c("Cpm", "Cpc")) {
    solver_null <- qc:::make_solver(
      numeric(0), -1, 1, prior,
      metric = metric,
      target = NULL
    )
    solver_mid <- qc:::make_solver(
      numeric(0), -1, 1, prior,
      metric = metric,
      target = midpoint
    )

    expect_equal(solver_null(1), solver_mid(1), tolerance = 1e-10, info = metric)
  }
})

test_that("semi-conjugate sigma Cpm and Cpc avoid the unstable density backend", {
  skip_if_not_installed("BayesTools")

  x <- c(0.05, 0.12)
  prior <- qc:::.bayestools_to_integration_prior(
    BayesTools::prior("cauchy", list(0, 5)),
    qc:::create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 1, beta0 = 1)
  )$prior

  for (metric in c("Cpm", "Cpc")) {
    expect_false(qc:::.integration_can_use_density(prior, metric), info = metric)

    solver <- qc:::make_solver(
      x, -1, 1, prior,
      metric = metric,
      target = 0
    )

    expect_equal(
      qc:::compute_cpk_prob_integration(
        x, -1, 1, c(0.5, Inf),
        prior,
        metric = metric,
        target = 0
      ),
      solver(0.5),
      tolerance = 1e-6,
      info = metric
    )

    fit_default <- qc:::analyze_capability_integration(
      x, -1, 1, prior,
      metric = metric,
      target = 0,
      n_grid = 128,
      use_density_solver = TRUE
    )
    fit_survival <- qc:::analyze_capability_integration(
      x, -1, 1, prior,
      metric = metric,
      target = 0,
      n_grid = 128,
      use_density_solver = FALSE
    )

    expect_equal(
      unname(fit_default$stats[c("Mean", "Median", "SD")]),
      unname(fit_survival$stats[c("Mean", "Median", "SD")]),
      tolerance = 1e-4,
      info = metric
    )

    pdf_fn <- qc:::make_density_solver(
      x, -1, 1, prior,
      metric = metric,
      target = 0
    )
    density_mass <- stats::integrate(
      Vectorize(pdf_fn),
      0.5,
      0.8,
      rel.tol = 1e-5,
      subdivisions = 400
    )$value

    expect_equal(
      density_mass,
      solver(0.5) - solver(0.8),
      tolerance = 0.01,
      info = metric
    )
  }
})
