testthat::test_that(".integration_compute_metric_moments(request=...) matches the old signature for a conjugate case", {
  x <- c(4.8, 5.1, 5.4, 4.9, 5.2, 5.0)
  prior <- qc:::create_prior_conjugate(mu0 = 5, k0 = 1, alpha0 = 2, beta0 = 1)
  request <- qc:::.new_qc_integration_request(
    data = x,
    LSL = 2,
    USL = 8,
    prior = prior,
    metric = "Cp",
    target = NULL,
    sigma_level = 3
  )

  direct <- qc:::compute_metric_moments(
    x, 2, 8, prior,
    metric = "Cp",
    target = NULL,
    sigma_level = 3
  )
  via_request <- qc:::.integration_compute_metric_moments(request = request)

  testthat::expect_equal(via_request, direct, tolerance = 1e-12)
})

testthat::test_that(".integration_make_solver(request=...) matches the old signature for a conjugate case", {
  x <- c(9.6, 10.1, 10.4, 9.9, 10.2, 10.0)
  prior <- qc:::create_prior_conjugate(mu0 = 10, k0 = 1, alpha0 = 2, beta0 = 1)
  request <- qc:::.new_qc_integration_request(
    data = x,
    LSL = 4,
    USL = 16,
    prior = prior,
    metric = "Cpk",
    target = NULL,
    sigma_level = 3
  )

  direct <- qc:::make_solver(
    x, 4, 16, prior,
    metric = "Cpk",
    target = NULL,
    sigma_level = 3
  )
  via_request <- qc:::.integration_make_solver(request = request)

  points <- c(0.5, 1, 1.5, 2)
  testthat::expect_equal(vapply(points, direct, numeric(1)), vapply(points, via_request, numeric(1)), tolerance = 1e-12)
})

testthat::test_that(".integration_make_density_solver(request=...) matches the old signature for a semi-conjugate fallback", {
  x <- c(0.2, -0.1, 0.05, 0.15)
  prior <- qc:::create_prior_semi_mu(
    mu0 = 0,
    k0 = 1,
    log_dens_sigma = function(sigma) {
      stats::dgamma(sigma, shape = 2, rate = 1, log = TRUE)
    }
  )
  request <- qc:::.new_qc_integration_request(
    data = x,
    LSL = -1,
    USL = 1,
    prior = prior,
    metric = "Cpk",
    target = NULL,
    sigma_level = 3
  )

  direct <- qc:::make_density_solver(
    x, -1, 1, prior,
    metric = "Cpk",
    target = NULL,
    sigma_level = 3
  )
  via_request <- qc:::.integration_make_density_solver(request = request)

  points <- c(0.25, 0.5, 1, 2)
  testthat::expect_equal(vapply(points, direct, numeric(1)), vapply(points, via_request, numeric(1)), tolerance = 1e-10)
})
