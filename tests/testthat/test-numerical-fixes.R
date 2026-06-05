# Tests for numerical fixes: #5 (HDI on log grids), #8 (Jensen bias), #9 (GH order)

testthat::test_that("#5: HDI is correct on log-spaced grids (prior-only, heavy tail)", {
  # With a broad sigma prior (gamma(2,1)), the prior-only density of Cp is

  # heavy-tailed and uses a log-spaced grid. The HDI should cover ~95% of mass
  # and its width should be reasonable (not inflated by constant dx[1]).
  skip_if_not_installed("BayesTools")

  LSL <- 0; USL <- 6; target <- 3
  mu_prior    <- BayesTools::prior("normal", list(3, 1))
  sigma_prior <- BayesTools::prior("gamma",  list(2, 1))

  conv <- qc:::.bayestools_to_integration_prior(mu_prior, sigma_prior)
  result <- qc:::analyze_capability_integration(
    numeric(0), LSL, USL, conv$prior,
    metric = "Cp", target = target, n_grid = 1024L
  )

  hdi_lo <- result$stats["HDI_Lo"]
  hdi_hi <- result$stats["HDI_Hi"]
  hdi_width <- hdi_hi - hdi_lo

  # The HDI must be a finite positive interval
  expect_true(is.finite(hdi_lo) && is.finite(hdi_hi))
  expect_true(hdi_width > 0)

  # Verify HDI covers ~95% mass using the density grid
  g <- result$grid
  in_hdi <- g$x >= hdi_lo & g$x <= hdi_hi
  dx <- diff(c(g$x[1] - (g$x[2] - g$x[1]) / 2,
               (g$x[-length(g$x)] + g$x[-1]) / 2,
               g$x[length(g$x)] + (g$x[length(g$x)] - g$x[length(g$x) - 1]) / 2))
  mass_in_hdi <- sum(g$density[in_hdi] * dx[in_hdi])
  total_mass  <- sum(g$density * dx)
  coverage <- mass_in_hdi / total_mass

  expect_true(coverage > 0.90,
              info = sprintf("HDI coverage = %.3f, expected >= 0.90", coverage))
})

testthat::test_that("#8: semi-conjugate sigma moments match 2D numerical integration", {
  # compute_metric_moments.PriorSemiConjugateSigma previously used a
  # plug-in E[sigma|mu] which introduced Jensen-inequality bias for 1/sigma
  # metrics. Compare against a brute-force 2D numerical integral.
  skip_if_not_installed("BayesTools")
  skip_if_not_installed("cubature")

  set.seed(99)
  x <- rnorm(30, 5, 1)
  LSL <- 2; USL <- 8; target <- 5

  # Build semi-conjugate sigma prior (Jeffreys sigma, Normal mu)
  mu_prior    <- BayesTools::prior("normal", list(5, 2))
  sigma_prior <- "Jeffreys_sigma"

  conv <- qc:::.bayestools_to_integration_prior(mu_prior, sigma_prior)
  prior <- conv$prior

  analytic <- qc:::compute_metric_moments(
    x, LSL, USL, prior,
    metric = "Cp", target = target,
    use_analytic = TRUE
  )
  exact <- qc:::compute_metric_moments(
    x, LSL, USL, prior,
    metric = "Cp", target = target,
    use_analytic = FALSE
  )

  # Reference: full 2D numerical integration (no Jensen approximation)
  n <- length(x); x_bar <- mean(x); sse <- sum((x - x_bar)^2)
  alpha0 <- prior$alpha0; beta0 <- prior$beta0
  alpha_n <- alpha0 + n / 2
  log_dens_mu <- qc:::.make_prior_log_dens_fn(mu_prior)
  alpha_0_times_logbeta0 <- if (beta0 == 0) 0 else alpha0 * log(beta0)

  log_joint <- function(mu, sigma) {
    if (sigma <= 0) return(-Inf)
    sse_mu <- sse + n * (mu - x_bar)^2
    beta_n <- beta0 + sse_mu / 2
    log_marg <- lgamma(alpha_n) - lgamma(alpha0) + alpha_0_times_logbeta0 - alpha_n * log(beta_n)
    log_dens_mu(mu) - log(sigma)   # sigma kernel from InvGamma(alpha0=-0.5, beta0=0)
  }

  sd_data <- sqrt(sse / (n - 1))
  mu_lo <- x_bar - 10 * sd_data; mu_hi <- x_bar + 10 * sd_data
  sig_hi <- 5 * sd_data
  tol <- USL - LSL

  log_post_2d <- function(mu, sigma) {
    if (sigma <= 0) return(-Inf)
    sse_mu <- sse + n * (mu - x_bar)^2
    -n * log(sigma) - sse_mu / (2 * sigma^2) + log_dens_mu(mu) - log(sigma)
  }

  h_max <- -optim(c(x_bar, sd_data), function(p) -log_post_2d(p[1], p[2]),
                   method = "L-BFGS-B", lower = c(-100, 1e-6))$value

  integrand_Z <- function(x) {
    mu <- x[1, ]; sigma <- x[2, ]
    v <- vapply(seq_along(mu), function(i) log_post_2d(mu[i], sigma[i]), numeric(1))
    matrix(exp(v - h_max), nrow = 1)
  }
  Z_ref <- cubature::pcubature(integrand_Z,
                                c(mu_lo, 1e-6), c(mu_hi, sig_hi),
                                tol = 1e-5, vectorInterface = TRUE)$integral

  integrand_E1 <- function(x) {
    mu <- x[1, ]; sigma <- x[2, ]
    v <- vapply(seq_along(mu), function(i) log_post_2d(mu[i], sigma[i]), numeric(1))
    m <- tol / (6 * sigma)
    matrix(exp(v - h_max) * m, nrow = 1)
  }
  E1_ref <- cubature::pcubature(integrand_E1,
                                 c(mu_lo, 1e-6), c(mu_hi, sig_hi),
                                 tol = 1e-5, vectorInterface = TRUE)$integral / Z_ref

  # Both the semi-analytic and generic-reference paths should stay close to the
  # same 2D posterior integral.
  expect_equal(analytic$mean, E1_ref, tolerance = 0.02,
               info = sprintf("semi-conj mean=%.4f, 2D ref=%.4f", analytic$mean, E1_ref))
  expect_equal(exact$mean, E1_ref, tolerance = 0.02,
               info = sprintf("generic-fallback mean=%.4f, 2D ref=%.4f", exact$mean, E1_ref))
  expect_equal(analytic$mean, exact$mean, tolerance = 0.01)
})

testthat::test_that("#9: Gauss-Hermite quadrature is accurate for Cpk", {
  # Cpk = min(Cpu, Cpl) has a kink at mu = midpoint. Low-order GH can miss this.
  # Compare compute_metric_moments.PriorSemiConjugateMu against a brute-force 2D reference.
  skip_if_not_installed("BayesTools")
  skip_if_not_installed("cubature")

  set.seed(42)
  x <- rnorm(20, 5.5, 1)  # Slightly off-center to stress the kink

  LSL <- 2; USL <- 8; target <- 5

  mu_prior_bt    <- "Jeffreys_mu"
  sigma_prior_bt <- BayesTools::prior("gamma", list(2, 1))

  conv <- qc:::.bayestools_to_integration_prior(mu_prior_bt, sigma_prior_bt)
  prior <- conv$prior
  stopifnot(inherits(prior, "PriorSemiConjugateMu"))

  result <- qc:::compute_metric_moments(x, LSL, USL, prior, metric = "Cpk", target = target)

  # Reference: direct 2D integration
  n <- length(x); x_bar <- mean(x); sse <- sum((x - x_bar)^2)
  k_n <- prior$k0 + n
  mu_n <- if (k_n > 0) (prior$k0 * prior$mu0 + n * x_bar) / k_n else x_bar
  sse_n <- sse + prior$k0 * n * (x_bar - prior$mu0)^2 / max(k_n, 1)
  log_dens_sigma <- prior$log_dens_sigma

  log_post_2d <- function(mu, sigma) {
    if (sigma <= 0) return(-Inf)
    -n * log(sigma) - (sse + n * (mu - x_bar)^2) / (2 * sigma^2) + log_dens_sigma(sigma)
  }
  sd_hat <- sqrt(sse / (n - 1))
  h_max <- -optim(c(x_bar, sd_hat), function(p) -log_post_2d(p[1], p[2]),
                   method = "L-BFGS-B", lower = c(-100, 1e-6))$value

  mu_lo <- x_bar - 10 * sd_hat; mu_hi <- x_bar + 10 * sd_hat; sig_hi <- 5 * sd_hat

  mk_integrand <- function(power) {
    function(x) {
      mu <- x[1, ]; sigma <- x[2, ]
      v <- vapply(seq_along(mu), function(i) log_post_2d(mu[i], sigma[i]), numeric(1))
      m <- qc:::compute_metric_value(mu, sigma, LSL, USL, target, "Cpk")
      matrix(exp(v - h_max) * m^power, nrow = 1)
    }
  }

  Z_ref <- cubature::pcubature(mk_integrand(0),
                                c(mu_lo, 1e-6), c(mu_hi, sig_hi),
                                tol = 1e-5, vectorInterface = TRUE)$integral
  E1_ref <- cubature::pcubature(mk_integrand(1),
                                 c(mu_lo, 1e-6), c(mu_hi, sig_hi),
                                 tol = 1e-5, vectorInterface = TRUE)$integral / Z_ref

  expect_equal(result$mean, E1_ref, tolerance = 0.02,
               info = sprintf("GH Cpk mean=%.4f, 2D ref=%.4f", result$mean, E1_ref))
})

testthat::test_that("conjugate Cpk integration is stable for high-information off-center data", {
  LSL <- 6
  USL <- 9
  target <- 7.5
  prior <- qc:::.prior_DCSI(LSL, USL)
  state <- qc:::.as_qc_suff_stats_state(
    cached_state = list(
      n = 1566,
      x_bar = 7.452067,
      sse = (1566 - 1) * 0.5162514^2
    )
  )

  cpk_moments <- qc:::compute_metric_moments(
    numeric(0), LSL, USL, prior,
    metric = "Cpk", target = target,
    cached_state = state
  )
  cp_moments <- qc:::compute_metric_moments(
    numeric(0), LSL, USL, prior,
    metric = "Cp", target = target,
    cached_state = state
  )

  expect_lt(abs(cpk_moments$mean - 0.93847), 1e-4)
  expect_lt(abs(cpk_moments$sd - 0.01874), 1e-4)
  expect_lt(cpk_moments$mean, cp_moments$mean)

  request <- qc:::.new_qc_integration_request(
    data = numeric(0),
    LSL = LSL,
    USL = USL,
    prior = prior,
    metric = "Cpk",
    target = target,
    cached_state = state,
    sigma_level = 3
  )
  S <- qc:::.integration_make_solver(request = request)
  thresholds <- c(0.70, 0.75, 0.80, 0.88, 0.90, 0.94, 0.97, 1.00)
  survival <- S(thresholds)

  expect_true(all(diff(survival) <= 1e-8))
  expect_gt(survival[thresholds == 0.75], 0.999)
  expect_equal(survival[thresholds == 0.94], 0.467, tolerance = 0.005)

  result <- qc:::analyze_capability_integration(
    request = request,
    n_grid = 512L
  )
  mode_x <- result$grid$x[which.max(result$grid$density)]

  expect_gt(min(result$grid$x), 0.85)
  expect_lt(max(result$grid$x), 1.01)
  expect_equal(mode_x, 0.938, tolerance = 0.01)
})

testthat::test_that("conjugate Cpk numerical moments use finite high-information bounds", {
  LSL <- 48
  USL <- 52
  target <- 50
  prior <- qc:::.bayestools_to_integration_prior(
    "Jeffreys_mu",
    "Jeffreys_sigma"
  )$prior
  state <- qc:::.as_qc_suff_stats_state(
    cached_state = list(
      n = 100,
      x_bar = 50,
      sse = (100 - 1) * 0.001^2
    )
  )

  cp_moments <- qc:::compute_metric_moments(
    numeric(0), LSL, USL, prior,
    metric = "Cp", target = target,
    cached_state = state
  )
  cpk_moments <- qc:::compute_metric_moments(
    numeric(0), LSL, USL, prior,
    metric = "Cpk", target = target,
    cached_state = state
  )

  expect_gt(cpk_moments$mean, 600)
  expect_gt(cpk_moments$sd, 40)
  expect_lt(abs(cp_moments$mean - cpk_moments$mean), 0.1)
})

testthat::test_that("survival-derived densities support scalar survival solvers", {
  scalar_survival <- function(c) {
    if (length(c) != 1L) {
      stop("scalar only")
    }
    exp(-c)
  }
  density <- qc:::.density_from_survival_fn(
    scalar_survival,
    rel_step = 1e-5,
    abs_step = 1e-6,
    support_lower = 0
  )
  x <- c(0.1, 1, 2)

  expect_equal(density(x), exp(-x), tolerance = 1e-4)
})

testthat::test_that("#1: gamma overflow is fixed for large n", {
  set.seed(1)
  x <- rnorm(500, 10, 2)
  prior <- qc:::create_prior_conjugate()
  result <- qc:::compute_metric_moments(x, 4, 16, prior, metric = "Cp", use_analytic = TRUE)
  expect_true(is.finite(result$mean), info = "E[Cp] should be finite for n=500")
  expect_true(is.finite(result$sd),   info = "SD[Cp] should be finite for n=500")
  expect_true(result$mean > 0)
})
