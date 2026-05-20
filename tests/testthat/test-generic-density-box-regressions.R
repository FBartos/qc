generic_box_prior <- function(mu_mean, mu_sd, sigma_shape, sigma_rate) {
  qc:::create_prior_generic(function(mu, sigma) {
    out <- rep(-Inf, length(mu))
    valid <- is.finite(mu) & is.finite(sigma) & sigma > 0
    if (any(valid)) {
      out[valid] <- stats::dnorm(mu[valid], mu_mean, mu_sd, log = TRUE) +
        stats::dgamma(
          sigma[valid],
          shape = sigma_shape,
          rate = sigma_rate,
          log = TRUE
        )
    }
    out
  })
}

expect_generic_density_mass_matches_survival <- function(
  x, LSL, USL, prior, metric, bounds,
  target = NULL, sigma_level = 3,
  cached_state = NULL, tolerance = 0.01
) {
  if (is.null(cached_state)) {
    cached_state <- qc:::precompute_generic_state(x, prior)
  }

  pdf_fn <- qc:::make_density_solver(
    x, LSL, USL, prior,
    metric = metric,
    target = target,
    cached_state = cached_state,
    sigma_level = sigma_level
  )
  S_fn <- qc:::make_solver(
    x, LSL, USL, prior,
    metric = metric,
    target = target,
    cached_state = cached_state,
    sigma_level = sigma_level
  )

  density_mass <- stats::integrate(
    Vectorize(pdf_fn),
    bounds[1], bounds[2],
    rel.tol = 1e-5,
    subdivisions = 1200
  )$value
  solver_mass <- as.numeric(S_fn(bounds[1]) - S_fn(bounds[2]))

  testthat::expect_equal(
    density_mass,
    solver_mass,
    tolerance = tolerance,
    info = sprintf(
      "%s bounds=[%.6f, %.6f]: density=%.8f survival=%.8f",
      metric, bounds[1], bounds[2], density_mass, solver_mass
    )
  )
}

testthat::test_that("generic Cpc bypasses the direct density backend", {
  prior <- generic_box_prior(mu_mean = 0, mu_sd = 1, sigma_shape = 2, sigma_rate = 1)

  testthat::expect_false(qc:::.integration_can_use_density(prior, "Cpc"))
  testthat::expect_true(qc:::.integration_can_use_density(prior, "Cpm"))
})

testthat::test_that("prior-only raw generic state rejects box-dependent normalization", {
  testthat::skip_if_not_installed("cubature")

  improper_prior <- qc:::create_prior_generic(function(mu, sigma) {
    ifelse(
      sigma <= 0,
      -Inf,
      stats::dgamma(sigma, shape = 2, rate = 1, log = TRUE)
    )
  })

  testthat::expect_error(
    qc:::precompute_generic_state(numeric(0), improper_prior),
    regexp = "improper.*arbitrary finite box"
  )
})

testthat::test_that("prior-only generic state still rejects box-dependent normalization with partial metadata", {
  testthat::skip_if_not_installed("BayesTools")
  testthat::skip_if_not_installed("cubature")

  improper_prior <- qc:::create_prior_generic(
    function(mu, sigma) {
      ifelse(sigma <= 0, -Inf, stats::dnorm(mu, 0, 1, log = TRUE))
    },
    bayestools_priors = list(
      mu = BayesTools::prior("normal", list(0, 1))
    )
  )

  testthat::expect_error(
    qc:::precompute_generic_state(numeric(0), improper_prior),
    regexp = "improper.*arbitrary finite box"
  )
})

testthat::test_that("prior-only generic state still rejects box-dependent normalization with sigma-only metadata", {
  testthat::skip_if_not_installed("cubature")

  sigma_prior <- qc:::create_prior_conjugate(
    mu0 = 0, k0 = 1,
    alpha0 = 2, beta0 = 1
  )
  improper_prior <- qc:::create_prior_generic(
    function(mu, sigma) {
      ifelse(
        sigma <= 0,
        -Inf,
        qc:::.integration_make_sigma_log_dens_fn(sigma_prior)(sigma)
      )
    },
    bayestools_priors = list(sigma = sigma_prior)
  )

  testthat::expect_error(
    qc:::precompute_generic_state(numeric(0), improper_prior),
    regexp = "improper.*arbitrary finite box"
  )
})

testthat::test_that("generic linear contour densities use the cached normalization box", {
  testthat::skip_if_not_installed("cubature")

  LSL <- -1
  USL <- 1
  sigma_level <- 3

  cases <- list(
    Cpu = list(
      x = c(rep(-0.1, 28), rep(-0.08, 2)),
      prior = generic_box_prior(mu_mean = 0.95, mu_sd = 0.03, sigma_shape = 120, sigma_rate = 2400),
      bounds = c(15, 20)
    ),
    Cpl = list(
      x = c(rep(0.1, 28), rep(0.12, 2)),
      prior = generic_box_prior(mu_mean = -0.95, mu_sd = 0.03, sigma_shape = 120, sigma_rate = 2400),
      bounds = c(15, 20)
    )
  )

  for (metric in names(cases)) {
    case <- cases[[metric]]
    cached_state <- qc:::precompute_generic_state(case$x, case$prior)
    limit <- if (metric == "Cpu") {
      (USL - cached_state[["mu_lower"]]) /
        (sigma_level * cached_state[["sigma_lower"]])
    } else {
      (cached_state[["mu_upper"]] - LSL) /
        (sigma_level * cached_state[["sigma_lower"]])
    }

    # These intervals sit entirely above the largest positive value reachable
    # inside the cached normalization box, so the correct density mass is zero.
    testthat::expect_lt(limit, case$bounds[1])
    testthat::expect_gt(abs(cached_state[["map_mu"]] - cached_state[["x_bar"]]), 0.5)

    expect_generic_density_mass_matches_survival(
      case$x, LSL, USL, case$prior,
      metric = metric,
      bounds = case$bounds,
      cached_state = cached_state,
      sigma_level = sigma_level
    )
  }
})

testthat::test_that("generic smooth contour densities drop contour segments outside the cached box", {
  testthat::skip_if_not_installed("cubature")

  x <- rep(0.02, 30)
  LSL <- -1
  USL <- 1
  target <- 0.1
  prior <- generic_box_prior(mu_mean = 0.9, mu_sd = 0.05, sigma_shape = 20, sigma_rate = 80)
  cached_state <- qc:::precompute_generic_state(x, prior)

  # The target-centered contours live outside the cached mu box for this
  # posterior, so any positive density mass here indicates inconsistent bounds.
  testthat::expect_lt(cached_state[["mu_upper"]], target)

  cases <- list(
    Cpm = c(0.2870954, 0.44714877),
    Cpc = c(0.3091027, 0.46915610)
  )

  for (metric in names(cases)) {
    expect_generic_density_mass_matches_survival(
      x, LSL, USL, prior,
      metric = metric,
      bounds = cases[[metric]],
      target = target,
      cached_state = cached_state
    )
  }
})
