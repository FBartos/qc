testthat::test_that("conjugate survival probabilities match independent Gamma integrals", {
  # Hand-calculated posterior for x = (4, 6, 8) and this prior:
  # mu_n = 5.75, k_n = 4, alpha_n = 5.5, beta_n = 7.375.
  # These constants deliberately do not use the package's posterior helpers.
  x <- c(4, 6, 8)
  prior <- qc::prior_conjugate(mu0 = 5, k0 = 1, alpha0 = 4, beta0 = 3)
  points <- c(-1, 0, 0.25, 0.5, 1, 2.5)

  for (sigma_level in c(3, 4)) {
    cp <- qc:::make_solver(x, 2, 8, prior, metric = "Cp", sigma_level = sigma_level)
    expected_cp <- ifelse(points <= 0, 1, stats::pgamma(
      (2 * sigma_level * points / 6)^2,
      shape = 5.5, rate = 7.375, lower.tail = FALSE
    ))
    testthat::expect_equal(vapply(points, cp, numeric(1)), expected_cp, tolerance = 1e-8)

    for (metric in c("Cpu", "Cpl")) {
      solver <- qc:::make_solver(x, 2, 8, prior, metric = metric,
                                sigma_level = sigma_level)
      distance <- if (metric == "Cpu") 8 - 5.75 else 5.75 - 2
      expected <- vapply(points, function(c) {
        stats::integrate(function(precision) {
          stats::pnorm(2 * distance * sqrt(precision) - 2 * sigma_level * c) *
            stats::dgamma(precision, shape = 5.5, rate = 7.375)
        }, 0, Inf, rel.tol = 1e-8)$value
      }, numeric(1))
      testthat::expect_equal(vapply(points, solver, numeric(1)), expected,
                            tolerance = 1e-5, info = metric)
    }
  }
})

testthat::test_that("direct Cpu/Cpl densities agree with survival probabilities", {
  set.seed(42)
  x <- rnorm(30, mean = 5, sd = 1)
  prior <- qc::prior_conjugate()

  for (metric in c("Cpu", "Cpl")) {
    S <- qc:::make_solver(x, 2, 8, prior, metric = metric, target = 5)
    pdf <- qc:::make_density_solver(x, 2, 8, prior, metric = metric, target = 5)
    for (c in c(0.5, 1, 1.5, 2, 2.5)) {
      mass <- stats::integrate(pdf, c, Inf, subdivisions = 200)$value
      testthat::expect_equal(mass, S(c), tolerance = 0.02, info = metric)
    }
  }
})
