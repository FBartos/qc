testthat::test_that("request-only analyses retain the data, limits, and sigma level", {
  # x = (2, 4, 6), prior alpha = 3, beta = 1, mu0 = 4:
  # posterior precision is Gamma(shape = 4.5, rate = 5).
  # Cp = 12/(2*4*sigma), so its moments are scaled Gamma half-moments.
  request <- qc:::.new_qc_integration_request(
    data = c(2, 4, 6), LSL = -2, USL = 10,
    prior = qc::prior_conjugate(mu0 = 4, k0 = 2, alpha0 = 3, beta0 = 1),
    metric = "Cp", sigma_level = 4
  )
  result <- qc:::analyze_capability_integration(request = request)
  scale <- 12 / (8 * sqrt(5))
  gamma_mean <- exp(lgamma(5) - lgamma(4.5))

  testthat::expect_equal(unname(result$stats["Mean"]), scale * gamma_mean,
                        tolerance = 1e-8)
  testthat::expect_equal(unname(result$stats["SD"]),
                        scale * sqrt(4.5 - gamma_mean^2), tolerance = 1e-8)
  testthat::expect_equal(unname(result$stats["Median"]),
                        scale * sqrt(stats::qgamma(0.5, shape = 4.5)),
                        tolerance = 1e-3)
})

testthat::test_that("integration analyses default the Cpm target to the specification midpoint", {
  # Jeffreys posterior: sigma^2 ~ InvGamma((n - 1) / 2, sse / 2) and
  # mu | sigma ~ Normal(x_bar, sigma^2 / n). With LSL = 2 and USL = 10 the
  # midpoint target is 6, so Cpm = 4 / (3 * sqrt(sigma^2 + (mu - 6)^2)).
  x <- c(4.4, 4.9, 5.2, 5.6, 5.1, 4.7)
  n <- length(x)
  set.seed(11)
  sigma <- sqrt((sum((x - mean(x))^2) / 2) / stats::rgamma(4e5, shape = (n - 1) / 2))
  mu <- stats::rnorm(4e5, mean(x), sigma / sqrt(n))
  cpm <- 4 / (3 * sqrt(sigma^2 + (mu - 6)^2))

  result <- qc:::analyze_capability_integration(
    x, LSL = 2, USL = 10, prior = qc::prior_conjugate(),
    metric = "Cpm", target = NULL
  )

  testthat::expect_equal(unname(result$stats["Median"]), stats::median(cpm),
                         tolerance = 5e-3)
  testthat::expect_equal(unname(result$stats["Q2.5"]),
                         unname(stats::quantile(cpm, 0.025)), tolerance = 5e-3)
})
