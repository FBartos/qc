sigma_log_dens_reference <- function(sigma_prior) {
  alpha0 <- sigma_prior$alpha0
  beta0 <- sigma_prior$beta0
  log_const <- log(2) + alpha0 * log(beta0) - lgamma(alpha0)

  function(sigma) {
    ifelse(
      sigma <= 0,
      -Inf,
      log_const - (2 * alpha0 + 1) * log(sigma) - beta0 / sigma^2
    )
  }
}

log_sigma_interval_prob_reference <- function(lower, upper, shape, rate) {
  log_prob <- rep(-Inf, length(lower))

  full_region <- lower <= 0 & is.infinite(upper) & upper > 0
  log_prob[full_region] <- 0

  upper_only <- lower <= 0 & is.finite(upper) & upper > 0
  if (any(upper_only)) {
    log_prob[upper_only] <- stats::pgamma(
      1 / upper[upper_only]^2,
      shape = shape,
      rate = rate[upper_only],
      lower.tail = FALSE,
      log.p = TRUE
    )
  }

  lower_only <- is.finite(lower) & lower > 0 & is.infinite(upper) & upper > 0
  if (any(lower_only)) {
    log_prob[lower_only] <- stats::pgamma(
      1 / lower[lower_only]^2,
      shape = shape,
      rate = rate[lower_only],
      lower.tail = TRUE,
      log.p = TRUE
    )
  }

  finite_band <- is.finite(lower) & lower > 0 &
    is.finite(upper) & upper > lower
  if (any(finite_band)) {
    lower_band <- lower[finite_band]
    upper_band <- upper[finite_band]
    rate_band <- rate[finite_band]

    log_cdf_lower <- stats::pgamma(
      1 / lower_band^2,
      shape = shape,
      rate = rate_band,
      lower.tail = TRUE,
      log.p = TRUE
    )
    log_cdf_upper <- stats::pgamma(
      1 / upper_band^2,
      shape = shape,
      rate = rate_band,
      lower.tail = TRUE,
      log.p = TRUE
    )
    log_surv_upper <- stats::pgamma(
      1 / upper_band^2,
      shape = shape,
      rate = rate_band,
      lower.tail = FALSE,
      log.p = TRUE
    )
    log_surv_lower <- stats::pgamma(
      1 / lower_band^2,
      shape = shape,
      rate = rate_band,
      lower.tail = FALSE,
      log.p = TRUE
    )

    log_mass_cdf <- qc:::log_diff_exp(log_cdf_lower, log_cdf_upper)
    log_mass_surv <- qc:::log_diff_exp(log_surv_upper, log_surv_lower)

    log_band <- log_mass_cdf
    use_surv <- !is.finite(log_band) | (is.finite(log_mass_surv) & log_mass_surv > log_band)
    log_band[use_surv] <- log_mass_surv[use_surv]
    log_prob[finite_band] <- log_band
  }

  log_prob
}

semi_sigma_survival_reference <- function(data, prior, metric, threshold,
                                          LSL, USL, target,
                                          sigma_level = 3) {
  state <- qc:::.prepare_semi_sigma_state(
    data, prior, cached_state = NULL,
    context = "The semi-conjugate sigma posterior for test reference"
  )

  log_Z <- qc:::.compute_semi_sigma_log_Z_survival(
    state$n, state$x_bar, state$sse,
    state$alpha0, state$beta0,
    prior$log_dens_mu,
    jeffreys_adj = state$jeffreys_adj,
    mu_domain = state$mu_domain,
    context = "The semi-conjugate sigma posterior for test reference"
  )

  log_integrand <- function(mu) {
    sse_mu <- state$sse + state$n * (mu - state$x_bar)^2
    beta_n <- state$beta0 + sse_mu / 2
    log_weight_mu <- prior$log_dens_mu(mu) - state$alpha_n_eff * log(beta_n)

    sigma_region <- qc:::.metric_sigma_region_from_mu(
      metric, mu, threshold, LSL, USL, target,
      sigma_level = sigma_level
    )
    log_prob_sigma <- log_sigma_interval_prob_reference(
      sigma_region$lower, sigma_region$upper,
      shape = state$alpha_n_eff,
      rate = beta_n
    )

    vals <- log_weight_mu + log_prob_sigma - log_Z
    vals[!is.finite(vals)] <- -Inf
    vals
  }

  mu_low <- -Inf
  mu_high <- Inf
  if (metric == "Cpm") {
    mu_radius <- qc:::.cpm_spec_distance(LSL, USL, target) / (sigma_level * threshold)
    mu_low <- target - mu_radius
    mu_high <- target + mu_radius
  } else if (metric == "Cpc") {
    mu_radius <- (USL - LSL) / ((2 * sigma_level) * sqrt(pi / 2) * threshold)
    mu_low <- target - mu_radius
    mu_high <- target + mu_radius
  }
  if (is.na(mu_low) || is.na(mu_high) ||
      (is.finite(mu_low) && is.finite(mu_high) && mu_high <= mu_low)) {
    return(0)
  }

  qc:::.semi_sigma_integrate_log_kernel(
    log_integrand,
    state$mu_domain,
    lower = mu_low,
    upper = mu_high,
    rel.tol = 1e-6,
    subdivisions = 400
  )$value
}

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
  log_dens_sigma <- sigma_log_dens_reference(sigma_prior)

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

  expected <- qc:::compute_metric_value(
    mu = mu_samples,
    sigma = sig_samples,
    LSL = -1,
    USL = 1,
    target = 0,
    metric = "Cp"
  )

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

test_that("semi-conjugate sigma survival stays finite for heavy-tailed Cpm and Cpc", {
  skip_if_not_installed("BayesTools")

  x <- c(0.05, 0.12)
  prior <- qc:::.bayestools_to_integration_prior(
    BayesTools::prior("cauchy", list(0, 5)),
    qc:::create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 1, beta0 = 1)
  )$prior

  for (metric in c("Cpm", "Cpc")) {
    ref_prob <- semi_sigma_survival_reference(
      data = x,
      prior = prior,
      metric = metric,
      threshold = 0.5,
      LSL = -1,
      USL = 1,
      target = 0
    )
    solver <- qc:::make_solver(
      x, -1, 1, prior,
      metric = metric,
      target = 0
    )

    expect_true(ref_prob > 0, info = metric)
    expect_equal(solver(0.5), ref_prob, tolerance = 1e-6, info = metric)
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
