test_that("Cpm uses the nearest target-to-spec distance off the midpoint", {
  LSL <- 0
  USL <- 10
  target <- 1
  sigma_level <- 3
  threshold <- 0.25
  cpm_dist <- min(USL - target, target - LSL)

  # The feasible sigma radius and direct metric evaluation must both use the
  # nearest target-to-spec distance, not the full tolerance width.
  constraints <- qc:::get_metric_constraints(
    "Cpm", threshold, LSL, USL, target,
    sigma_level = sigma_level
  )
  expect_equal(constraints$s_max_fn(), cpm_dist / (sigma_level * threshold))

  expected_value <- cpm_dist / (sigma_level * sqrt(1 + (0.5 - target)^2))
  expect_equal(
    qc:::compute_metric_value(
      mu = 0.5, sigma = 1,
      LSL = LSL, USL = USL, target = target,
      metric = "Cpm", sigma_level = sigma_level
    ),
    expected_value
  )

  degenerate <- qc:::.degenerate_conjugate_metric_distribution(
    mu_n = 0.5, k_n = 2,
    LSL = LSL, USL = USL, target = target,
    metric = "Cpm", sigma_level = sigma_level
  )
  expect_equal(degenerate$type, "point")
  expect_equal(degenerate$value, cpm_dist / (sigma_level * abs(0.5 - target)))
})

test_that("off-center Cpm density and survival solvers stay aligned", {
  set.seed(123)

  # This target is intentionally off the midpoint so the test fails if any
  # duplicated contour code still uses tolerance / 2.
  x <- rnorm(25, mean = 0.8, sd = 0.3)
  LSL <- 0
  USL <- 10
  target <- 1
  prior <- qc:::create_prior_conjugate()

  density_result <- qc:::analyze_capability_integration(
    x, LSL, USL, prior,
    metric = "Cpm",
    target = target,
    use_density_solver = TRUE
  )
  survival_result <- qc:::analyze_capability_integration(
    x, LSL, USL, prior,
    metric = "Cpm",
    target = target,
    use_density_solver = FALSE
  )

  expect_equal(
    unname(density_result$stats["Mean"]),
    unname(survival_result$stats["Mean"]),
    tolerance = 0.02
  )
  expect_equal(
    unname(density_result$stats["Median"]),
    unname(survival_result$stats["Median"]),
    tolerance = 0.03
  )
})

test_that("prior-only Monte Carlo analysis reuses shared metric geometry", {
  skip_if_not_installed("BayesTools")

  prior <- qc:::.bayestools_to_integration_prior(
    BayesTools::prior("normal", list(0, 1)),
    BayesTools::prior("gamma", list(2, 1))
  )$prior

  mu_samples <- c(-1.3, -0.1, 0.4, 1.2)
  sig_samples <- c(0.35, 0.8, 1.4, 2.1)
  LSL <- -1
  USL <- 1
  target <- 0.6

  for (metric in c("Cp", "Cpu", "Cpl", "Cpk", "Cpm", "Cpc")) {
    result <- qc:::analyze_capability_integration(
      numeric(0), LSL, USL, prior,
      metric = metric,
      target = target,
      mc_samples = list(mu = mu_samples, sig = sig_samples)
    )

    expected <- qc:::compute_metric_value(
      mu = mu_samples,
      sigma = sig_samples,
      LSL = LSL,
      USL = USL,
      target = target,
      metric = metric
    )

    expect_equal(result$samples, expected, info = metric)
    expect_equal(unname(result$stats["Mean"]), mean(expected), info = metric)
    expect_equal(unname(result$stats["Median"]), stats::median(expected), info = metric)
    expect_equal(result$area, 1, info = metric)

    if (metric %in% c("Cpu", "Cpl", "Cpk")) {
      expect_true(any(result$samples < 0), info = metric)
    }
  }
})

test_that("semi-conjugate mu moments reuse cached sufficient statistics", {
  skip_if_not_installed("BayesTools")

  x <- c(4.2, 5.1, 5.6, 6.4, 7.0)
  cached_state <- list(
    n = length(x),
    x_bar = mean(x),
    sse = sum((x - mean(x))^2)
  )
  prior <- qc:::.bayestools_to_integration_prior(
    "Jeffreys_mu",
    BayesTools::prior("gamma", list(2, 1))
  )$prior

  # The stats-only path should produce the same moments as the raw-data path.
  raw_result <- qc:::compute_metric_moments(
    x, 4, 8, prior,
    metric = "Cpk",
    target = 6
  )
  cached_result <- qc:::compute_metric_moments(
    numeric(0), 4, 8, prior,
    metric = "Cpk",
    target = 6,
    cached_state = cached_state
  )

  expect_equal(cached_result, raw_result, tolerance = 1e-6)
})

test_that("partial cached_state fails fast with a clear suff-stat error", {
  prior <- qc:::create_prior_conjugate()

  expect_error(
    qc:::compute_metric_moments(
      c(1, 2, 3), 0, 4, prior,
      metric = "Cp",
      cached_state = list(n = 3L, x_bar = 2)
    ),
    "cached_state must contain `n`, `x_bar`, and `sse`"
  )
})

test_that("semi-conjugate sigma moments match density and survival summaries", {
  skip_if_not_installed("BayesTools")

  prior <- qc:::.bayestools_to_integration_prior(
    BayesTools::prior("normal", list(0, 2)),
    "Jeffreys_sigma"
  )$prior
  x <- c(0, 1)

  # The Jeffreys-adjusted shape parameter must be used consistently across the
  # moments, density, and survival implementations.
  moment_result <- qc:::compute_metric_moments(
    x, -1, 1, prior,
    metric = "Cp"
  )
  density_result <- qc:::analyze_capability_integration(
    x, -1, 1, prior,
    metric = "Cp",
    use_density_solver = TRUE
  )
  survival_result <- qc:::analyze_capability_integration(
    x, -1, 1, prior,
    metric = "Cp",
    use_density_solver = FALSE
  )

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

test_that("diffuse conjugate numerical moments stay aligned with analytic and Monte Carlo references", {
  set.seed(42)

  prior <- qc:::create_prior_conjugate(
    mu0 = -2,
    k0 = 0.1,
    alpha0 = 0.51,
    beta0 = 1e-4
  )
  cached_state <- list(n = 0L, x_bar = 0, sse = 0)

  cpu_analytic <- qc:::compute_metric_moments(
    numeric(0), -1, 1, prior,
    metric = "Cpu",
    target = 0,
    use_analytic = TRUE,
    cached_state = cached_state
  )
  cpu_numerical <- qc:::compute_metric_moments(
    numeric(0), -1, 1, prior,
    metric = "Cpu",
    target = 0,
    use_analytic = FALSE,
    cached_state = cached_state
  )

  expect_equal(cpu_numerical, cpu_analytic, tolerance = 1e-3)

  n_mc <- 200000L
  sigma2 <- 1 / stats::rgamma(n_mc, shape = prior$alpha0, rate = prior$beta0)
  sigma <- sqrt(sigma2)
  mu <- stats::rnorm(n_mc, mean = prior$mu0, sd = sigma / sqrt(prior$k0))

  for (metric in c("Cpm", "Cpc")) {
    mc_vals <- qc:::compute_metric_value(mu, sigma, -1, 1, 0, metric)
    result <- qc:::compute_metric_moments(
      numeric(0), -1, 1, prior,
      metric = metric,
      target = 0,
      use_analytic = TRUE,
      cached_state = cached_state
    )

    expect_true(abs(result$mean - mean(mc_vals)) < 0.002, info = metric)
    expect_true(abs(result$sd - stats::sd(mc_vals)) < 0.0015, info = metric)
  }
})

test_that("conjugate smooth-metric moments use the numerical reference path", {
  prior <- qc:::create_prior_conjugate(
    mu0 = 0.7,
    k0 = 0.01,
    alpha0 = 1.5,
    beta0 = 5
  )
  cached_state <- list(n = 0L, x_bar = 0, sse = 0)

  cpm_analytic <- qc:::compute_metric_moments(
    numeric(0), -1, 1, prior,
    metric = "Cpm",
    target = 0.8,
    use_analytic = TRUE,
    cached_state = cached_state
  )
  cpm_numerical <- qc:::compute_metric_moments(
    numeric(0), -1, 1, prior,
    metric = "Cpm",
    target = 0.8,
    use_analytic = FALSE,
    cached_state = cached_state
  )
  cpc_analytic <- qc:::compute_metric_moments(
    numeric(0), -1, 1, prior,
    metric = "Cpc",
    target = 0.8,
    use_analytic = TRUE,
    cached_state = cached_state
  )
  cpc_numerical <- qc:::compute_metric_moments(
    numeric(0), -1, 1, prior,
    metric = "Cpc",
    target = 0.8,
    use_analytic = FALSE,
    cached_state = cached_state
  )

  expect_equal(cpm_analytic, cpm_numerical, tolerance = 1e-6)
  expect_equal(cpc_analytic, cpc_numerical, tolerance = 1e-6)
})

test_that("semi-conjugate sigma diffuse moments agree with generic references", {
  skip_if_not_installed("BayesTools")
  skip_if_not_installed("cubature")

  prior <- qc:::.bayestools_to_integration_prior(
    BayesTools::prior("normal", list(0, 10)),
    qc:::create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 0.5001, beta0 = 1e-8)
  )$prior
  generic_prior <- qc:::.semi_sigma_as_generic_prior(prior)

  analytic_cp <- qc:::compute_metric_moments(
    numeric(0), -1, 1, prior,
    metric = "Cp",
    target = 0,
    use_analytic = TRUE
  )
  exact_cp <- qc:::compute_metric_moments(
    numeric(0), -1, 1, prior,
    metric = "Cp",
    target = 0,
    use_analytic = FALSE
  )
  generic_cp <- qc:::compute_metric_moments(
    numeric(0), -1, 1,
    generic_prior,
    metric = "Cp",
    target = 0
  )
  analytic_cpk <- qc:::compute_metric_moments(
    numeric(0), -1, 1, prior,
    metric = "Cpk",
    target = 0,
    use_analytic = TRUE
  )
  generic_cpk <- qc:::compute_metric_moments(
    numeric(0), -1, 1, generic_prior,
    metric = "Cpk",
    target = 0
  )
  analytic_cpm <- qc:::compute_metric_moments(
    numeric(0), -1, 1, prior,
    metric = "Cpm",
    target = 0.1,
    use_analytic = TRUE
  )
  exact_cpm <- qc:::compute_metric_moments(
    numeric(0), -1, 1, prior,
    metric = "Cpm",
    target = 0.1,
    use_analytic = FALSE
  )
  generic_cpm <- qc:::compute_metric_moments(
    numeric(0), -1, 1, generic_prior,
    metric = "Cpm",
    target = 0.1
  )

  expect_equal(analytic_cp$mean, generic_cp$mean, tolerance = 0.02)
  expect_equal(analytic_cp$sd, generic_cp$sd, tolerance = 0.02)
  expect_equal(analytic_cpk$mean, generic_cpk$mean, tolerance = 0.02)
  expect_equal(analytic_cpk$sd, generic_cpk$sd, tolerance = 0.02)
  expect_equal(analytic_cpm$mean, generic_cpm$mean, tolerance = 0.02)
  expect_equal(analytic_cpm$sd, generic_cpm$sd, tolerance = 0.02)
  expect_equal(exact_cp, generic_cp, tolerance = 1e-6)
  expect_equal(exact_cpm, generic_cpm, tolerance = 1e-6)
})

test_that("semi-conjugate mu smooth-metric moments match an equivalent generic prior", {
  skip_if_not_installed("BayesTools")
  skip_if_not_installed("cubature")

  sigma_prior <- BayesTools::prior("lognormal", list(0, 1))
  prior <- qc:::.bayestools_to_integration_prior("Jeffreys_mu", sigma_prior)$prior
  generic_prior <- qc:::.semi_mu_as_generic_prior(prior)
  cached_state <- list(n = 5L, x_bar = 0.7, sse = 1.2)

  semi_cpm <- qc:::compute_metric_moments(
    numeric(0), -1, 1, prior,
    metric = "Cpm",
    target = 0.8,
    use_analytic = TRUE,
    cached_state = cached_state
  )
  generic_cpm <- qc:::compute_metric_moments(
    numeric(0), -1, 1, generic_prior,
    metric = "Cpm",
    target = 0.8,
    cached_state = cached_state
  )
  semi_cpc <- qc:::compute_metric_moments(
    numeric(0), -1, 1, prior,
    metric = "Cpc",
    target = 0.8,
    use_analytic = FALSE,
    cached_state = cached_state
  )
  generic_cpc <- qc:::compute_metric_moments(
    numeric(0), -1, 1, generic_prior,
    metric = "Cpc",
    target = 0.8,
    cached_state = cached_state
  )

  expect_equal(semi_cpm, generic_cpm, tolerance = 1e-6)
  expect_equal(semi_cpc, generic_cpc, tolerance = 1e-6)
})

test_that("semi-conjugate mu rejects prior-only improper flat-mu posteriors explicitly", {
  skip_if_not_installed("BayesTools")

  prior <- qc:::.bayestools_to_integration_prior(
    "Jeffreys_mu",
    BayesTools::prior("lognormal", list(0, 0.5))
  )$prior

  expect_error(
    qc:::compute_metric_moments(
      numeric(0), 0, 10, prior,
      metric = "Cp",
      cached_state = list(n = 0L, x_bar = 0, sse = 0)
    ),
    "semi-conjugate mu posterior.*improper"
  )
})

test_that("semi-conjugate sigma smooth-metric fallback rejects prior-only improper flat-mu posteriors", {
  skip_if_not_installed("cubature")

  prior <- qc:::create_prior_semi_sigma(
    alpha0 = 2,
    beta0 = 1,
    log_dens_mu = function(mu) rep(0, length(mu))
  )

  expect_error(
    qc:::compute_metric_moments(
      numeric(0), 0, 10, prior,
      metric = "Cpm",
      target = 5
    ),
    "improper"
  )
})

test_that("generic Cp moments share the same state bounds as solver normalization", {
  skip_if_not_installed("BayesTools")
  skip_if_not_installed("cubature")

  make_prior <- function(sd_mu) {
    qc:::.bayestools_to_integration_prior(
      BayesTools::prior("normal", list(0, sd_mu)),
      BayesTools::prior("gamma", list(5, 5))
    )$prior
  }

  prior_narrow <- make_prior(1)
  prior_wide <- make_prior(50)

  # Cp depends only on sigma, so changing the width of the mu prior must not
  # change the posterior mean once the same integration box is used for Z.
  narrow_moments <- qc:::compute_metric_moments(
    numeric(0), -1, 1, prior_narrow,
    metric = "Cp"
  )
  wide_moments <- qc:::compute_metric_moments(
    numeric(0), -1, 1, prior_wide,
    metric = "Cp"
  )
  analysis_result <- qc:::analyze_capability_integration(
    numeric(0), -1, 1, prior_narrow,
    metric = "Cp"
  )
  solver <- qc:::make_solver(
    numeric(0), -1, 1, prior_narrow,
    metric = "Cp"
  )

  expect_equal(narrow_moments$mean, wide_moments$mean, tolerance = 0.01)
  expect_equal(
    narrow_moments$mean,
    unname(analysis_result$stats["Mean"]),
    tolerance = 0.01
  )
  expect_true(is.finite(solver(0.2)))
  expect_gte(solver(0.2), 0)
  expect_lte(solver(0.2), 1)
})

test_that("generic Cp density matches survival mass for truncated informative mu priors", {
  skip_if_not_installed("BayesTools")
  skip_if_not_installed("cubature")

  x <- c(0.02)
  LSL <- -1
  USL <- 1
  sigma_level <- 3
  prior <- qc:::.bayestools_to_integration_prior(
    BayesTools::prior("normal", list(1.2, 0.03), list(1.15, 1.25)),
    BayesTools::prior("uniform", list(0.03, 1.0))
  )$prior
  cached_state <- qc:::precompute_generic_state(x, prior)

  pdf_fn <- qc:::make_density_solver(
    x, LSL, USL, prior,
    metric = "Cp",
    cached_state = cached_state,
    sigma_level = sigma_level
  )
  S <- qc:::make_solver(
    x, LSL, USL, prior,
    metric = "Cp",
    cached_state = cached_state,
    sigma_level = sigma_level
  )

  cp_support <- c(
    (USL - LSL) / ((2 * sigma_level) * cached_state$sigma_upper),
    (USL - LSL) / ((2 * sigma_level) * cached_state$sigma_lower)
  )
  bounds <- c(
    cp_support[1] * 1.05,
    min(cp_support[2] * 0.95, cp_support[1] + 1)
  )

  density_mass <- stats::integrate(
    Vectorize(pdf_fn),
    bounds[1], bounds[2],
    rel.tol = 1e-5,
    subdivisions = 600
  )$value
  solver_mass <- S(bounds[1]) - S(bounds[2])

  expect_gt(abs(cached_state$map_mu - cached_state$x_bar), 1)
  expect_gt(density_mass, 0.8)
  expect_equal(density_mass, solver_mass, tolerance = 0.01)
})

test_that("prior conversion rejects point priors for integration and handles truncation tails in log space", {
  skip_if_not_installed("BayesTools")

  point_fn <- qc:::.make_prior_log_dens_fn(
    BayesTools::prior("point", list(location = 1.5))
  )
  expect_equal(unname(point_fn(1.5)), 0)
  expect_equal(unname(point_fn(c(1.5, 1.6))), c(0, -Inf))

  expect_error(
    qc:::.bayestools_to_integration_prior(
      BayesTools::prior("point", list(location = 1.5)),
      BayesTools::prior("lognormal", list(0, 0.5))
    ),
    "Point priors on `prior\\$mu`"
  )
  expect_error(
    qc:::.bayestools_to_integration_prior(
      BayesTools::prior("normal", list(0, 1)),
      BayesTools::prior("point", list(location = 1.5))
    ),
    "Point priors on `prior\\$sigma`"
  )

  extreme_log_norm <- qc:::.compute_truncation_norm(
    "normal",
    list(mean = 0, sd = 1),
    lower = 10,
    upper = 10.1
  )
  expect_true(is.finite(extreme_log_norm))

  truncated_fn <- qc:::.make_prior_log_dens_fn(
    BayesTools::prior("normal", list(0, 1), list(10, 10.1))
  )
  expect_true(is.finite(unname(truncated_fn(10.05))))
})
