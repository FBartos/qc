# Tests for numerical fixes: #5 (HDI on log grids), #8 (Jensen bias), #9 (GH order)

testthat::test_that("#5: HDI is correct on log-spaced grids (prior-only, heavy tail)", {
  # With a broad sigma prior (gamma(2,1)), the prior-only density of Cp is

  # heavy-tailed and uses a log-spaced grid. The HDI should cover ~95% of mass
  # and its width should be reasonable (not inflated by constant dx[1]).
  skip_if_not_installed("BayesTools")

  LSL <- 0; USL <- 6; target <- 3
  prior_mu    <- BayesTools::prior("normal", list(3, 1))
  prior_sigma <- BayesTools::prior("gamma",  list(2, 1))

  conv <- qc:::.bayestools_to_integration_prior(prior_mu, prior_sigma)
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
  prior_mu    <- BayesTools::prior("normal", list(5, 2))
  prior_sigma <- "Jeffreys_sigma"

  conv <- qc:::.bayestools_to_integration_prior(prior_mu, prior_sigma)
  prior <- conv$prior

  # Method under test
  result <- qc:::compute_metric_moments(x, LSL, USL, prior,
                                         metric = "Cp", target = target)

  # Reference: full 2D numerical integration (no Jensen approximation)
  n <- length(x); x_bar <- mean(x); sse <- sum((x - x_bar)^2)
  alpha0 <- prior$alpha0; beta0 <- prior$beta0
  alpha_n <- alpha0 + n / 2
  log_dens_mu <- qc:::.make_prior_log_dens_fn(prior_mu)
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

  # The semi-conjugate method should be close to the 2D reference
  expect_equal(result$mean, E1_ref, tolerance = 0.02,
               info = sprintf("semi-conj mean=%.4f, 2D ref=%.4f", result$mean, E1_ref))
})

testthat::test_that("#9: Gauss-Hermite quadrature is accurate for Cpk", {
  # Cpk = min(Cpu, Cpl) has a kink at mu = midpoint. Low-order GH can miss this.
  # Compare compute_metric_moments.PriorSemiConjugateMu against a brute-force 2D reference.
  skip_if_not_installed("BayesTools")
  skip_if_not_installed("cubature")

  set.seed(42)
  x <- rnorm(20, 5.5, 1)  # Slightly off-center to stress the kink

  LSL <- 2; USL <- 8; target <- 5

  prior_mu_bt    <- "Jeffreys_mu"
  prior_sigma_bt <- BayesTools::prior("gamma", list(2, 1))

  conv <- qc:::.bayestools_to_integration_prior(prior_mu_bt, prior_sigma_bt)
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

testthat::test_that("#1: gamma overflow is fixed for large n", {
  set.seed(1)
  x <- rnorm(500, 10, 2)
  prior <- qc:::create_prior_conjugate()
  result <- qc:::compute_metric_moments(x, 4, 16, prior, metric = "Cp", use_analytic = TRUE)
  expect_true(is.finite(result$mean), info = "E[Cp] should be finite for n=500")
  expect_true(is.finite(result$sd),   info = "SD[Cp] should be finite for n=500")
  expect_true(result$mean > 0)
})
