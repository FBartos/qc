testthat::test_that("semi-conjugate mu Cp density normalizes on narrow sigma support", {
  testthat::skip_if_not_installed("BayesTools")

  x <- rep(5, 6)
  LSL <- 4.5
  USL <- 5.5
  target <- 5
  sigma_bounds <- c(0.09, 0.11)
  prior <- qc:::.bayestools_to_integration_prior(
    "Jeffreys_mu",
    BayesTools::prior("uniform", as.list(sigma_bounds))
  )$prior

  pdf_fn <- qc:::make_density_solver(x, LSL, USL, prior, metric = "Cp", target = target)
  S_fn <- qc:::make_solver(x, LSL, USL, prior, metric = "Cp", target = target)

  cp_support <- sort((USL - LSL) / (2 * 3 * rev(sigma_bounds)))
  total_mass <- integrate(Vectorize(pdf_fn), cp_support[1], cp_support[2],
                          rel.tol = 1e-5, subdivisions = 400)$value
  mid_c <- mean(cp_support)
  tail_mass <- integrate(Vectorize(pdf_fn), mid_c, cp_support[2],
                         rel.tol = 1e-5, subdivisions = 400)$value

  testthat::expect_gt(pdf_fn(mid_c), 0)
  testthat::expect_equal(total_mass, 1, tolerance = 1e-3)
  testthat::expect_equal(tail_mass, S_fn(mid_c), tolerance = 0.01)
})

testthat::test_that("semi-conjugate mu density and survival agree for narrow truncated sigma priors", {
  testthat::skip_if_not_installed("BayesTools")

  x <- c(5, 5.000001, 4.999999, 5.000002, 4.999998)
  LSL <- 4.5
  USL <- 5.5
  target <- 5
  bounds <- c(1.6, 1.8)
  prior <- qc:::.bayestools_to_integration_prior(
    "Jeffreys_mu",
    BayesTools::prior("gamma", list(20, 200), list(0.08, 0.12))
  )$prior

  S_fn <- qc:::make_solver(x, LSL, USL, prior, metric = "Cpu", target = target)
  pdf_fn <- qc:::make_density_solver(x, LSL, USL, prior, metric = "Cpu", target = target)

  density_mass <- integrate(Vectorize(pdf_fn), bounds[1], bounds[2],
                            rel.tol = 1e-5, subdivisions = 400)$value
  solver_mass <- S_fn(bounds[1]) - S_fn(bounds[2])

  testthat::expect_gt(pdf_fn(bounds[1]), 0)
  testthat::expect_gt(density_mass, 0.1)
  testthat::expect_equal(density_mass, solver_mass, tolerance = 0.01)
})
