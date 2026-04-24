testthat::test_that("conjugate Cpk density keeps the low-c sigma tail", {
  prior <- qc:::create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 0.6, beta0 = 0.001)
  bounds <- c(0.01, 0.1)

  pdf_fn <- qc:::make_density_solver(
    numeric(0), -1, 1, prior,
    metric = "Cpk"
  )
  solver <- qc:::make_solver(
    numeric(0), -1, 1, prior,
    metric = "Cpk"
  )

  density_mass <- stats::integrate(
    Vectorize(pdf_fn),
    bounds[1],
    bounds[2],
    rel.tol = 1e-5,
    subdivisions = 1200
  )$value
  solver_mass <- solver(bounds[1]) - solver(bounds[2])

  testthat::expect_gt(density_mass, 1e-3)
  testthat::expect_equal(density_mass, solver_mass, tolerance = 1e-3)
})

testthat::test_that("raw generic priors keep sigma boxes wide enough to normalize", {
  testthat::skip_if_not_installed("cubature")

  prior <- qc:::create_prior_generic(function(mu, sigma) {
    ifelse(
      sigma <= 0,
      -Inf,
      stats::dnorm(mu, 0, 1, log = TRUE) +
        stats::dnorm(sigma, 500, 5, log = TRUE)
    )
  })

  state <- qc:::precompute_generic_state(numeric(0), prior)
  fit <- qc:::analyze_capability_integration(
    numeric(0), -1, 1, prior,
    metric = "Cp"
  )

  testthat::expect_equal(state$map_sig, 500, tolerance = 1)
  testthat::expect_gt(state$sigma_upper, state$map_sig)
  testthat::expect_gt(state$Z, 0)
  testthat::expect_true(all(is.finite(unname(fit$stats[c("Mean", "Median", "SD")]))))
})

testthat::test_that("survival-derived densities use a one-sided derivative at the lower support boundary", {
  S <- function(c) ifelse(c <= 0, 1, exp(-c))
  pdf_fn <- qc:::.density_from_survival_fn(S, support_lower = 0)

  testthat::expect_equal(pdf_fn(0), 1, tolerance = 1e-4)
  testthat::expect_equal(pdf_fn(1e-6), exp(-1e-6), tolerance = 1e-4)
  testthat::expect_equal(pdf_fn(1), exp(-1), tolerance = 1e-4)
})

testthat::test_that("semi-conjugate sigma Cp density stays aligned with survival on broad mu priors", {
  log_dens_mu <- function(mu) stats::dnorm(mu, mean = 0, sd = 20, log = TRUE)
  prior <- qc:::create_prior_semi_sigma(
    alpha0 = 2,
    beta0 = 1,
    log_dens_mu = log_dens_mu
  )

  pdf_fn <- qc:::make_density_solver(
    numeric(0), 0, 10, prior,
    metric = "Cp",
    target = 5
  )
  solver <- qc:::make_solver(
    numeric(0), 0, 10, prior,
    metric = "Cp",
    target = 5
  )

  bounds <- c(0.05, 1)
  density_mass <- stats::integrate(
    Vectorize(pdf_fn),
    bounds[1],
    bounds[2],
    rel.tol = 1e-5,
    subdivisions = 400
  )$value
  solver_mass <- solver(bounds[1]) - solver(bounds[2])

  testthat::expect_equal(density_mass, solver_mass, tolerance = 1e-6)
})

testthat::test_that("semi-conjugate sigma Cpm density ignores optional mu metadata", {
  testthat::skip_if_not_installed("BayesTools")

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

  pdf_plain <- qc:::make_density_solver(
    numeric(0), 0, 10, prior_plain,
    metric = "Cpm",
    target = 5
  )
  pdf_meta <- qc:::make_density_solver(
    numeric(0), 0, 10, prior_meta,
    metric = "Cpm",
    target = 5
  )

  testthat::expect_equal(pdf_plain(0.5), pdf_meta(0.5), tolerance = 1e-10)
})

testthat::test_that("generic negative-support densities cover Cpu, Cpl, and Cpk", {
  testthat::skip_if_not_installed("cubature")

  cases <- list(
    Cpu = 2,
    Cpl = -2,
    Cpk = 2
  )

  for (metric in names(cases)) {
    prior <- qc:::create_prior_generic(function(mu, sigma) {
      ifelse(
        sigma <= 0,
        -Inf,
        stats::dnorm(mu, cases[[metric]], 1, log = TRUE) +
          stats::dgamma(sigma, shape = 2, rate = 1, log = TRUE)
      )
    })

    pdf_fn <- qc:::make_density_solver(
      numeric(0), -1, 1, prior,
      metric = metric
    )
    solver <- qc:::make_solver(
      numeric(0), -1, 1, prior,
      metric = metric
    )

    c0 <- -0.5
    h <- max(abs(c0) * 1e-4, 1e-5)
    approx <- (solver(c0 - h) - solver(c0 + h)) / (2 * h)

    testthat::expect_gt(approx, 0)
    testthat::expect_equal(pdf_fn(c0), approx, tolerance = 1e-4, info = metric)
  }
})

testthat::test_that("semi-conjugate mu negative-support densities cover Cpu, Cpl, and Cpk", {
  cases <- list(
    Cpu = 2,
    Cpl = -2,
    Cpk = 2
  )

  for (metric in names(cases)) {
    prior <- qc:::create_prior_semi_mu(
      mu0 = cases[[metric]],
      k0 = 1,
      log_dens_sigma = function(sigma) {
        stats::dgamma(sigma, shape = 2, rate = 1, log = TRUE)
      }
    )

    pdf_fn <- qc:::make_density_solver(
      numeric(0), -1, 1, prior,
      metric = metric
    )
    solver <- qc:::make_solver(
      numeric(0), -1, 1, prior,
      metric = metric
    )

    c0 <- -0.5
    h <- max(abs(c0) * 1e-4, 1e-5)
    approx <- (solver(c0 - h) - solver(c0 + h)) / (2 * h)

    testthat::expect_gt(approx, 0)
    testthat::expect_equal(pdf_fn(c0), approx, tolerance = 1e-4, info = metric)
  }
})

testthat::test_that("generic Cpm density defaults target to the midpoint", {
  testthat::skip_if_not_installed("cubature")

  prior <- qc:::create_prior_generic(function(mu, sigma) {
    ifelse(
      sigma <= 0,
      -Inf,
      stats::dnorm(mu, 0, 1, log = TRUE) +
        stats::dgamma(sigma, shape = 2, rate = 1, log = TRUE)
    )
  })

  pdf_null <- qc:::make_density_solver(
    numeric(0), -1, 1, prior,
    metric = "Cpm",
    target = NULL
  )
  pdf_mid <- qc:::make_density_solver(
    numeric(0), -1, 1, prior,
    metric = "Cpm",
    target = 0
  )

  testthat::expect_gt(pdf_mid(0.5), 0)
  testthat::expect_equal(pdf_null(0.5), pdf_mid(0.5), tolerance = 1e-8)
})
