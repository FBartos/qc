# Independent reference for the Jeffreys conjugate posterior:
# sigma^2 ~ InvGamma((n - 1) / 2, sse / 2) and mu | sigma ~ Normal(x_bar, sigma^2 / n).
# Cpu <= c iff mu >= USL - 3 c sigma; Cpl <= c iff mu <= LSL + 3 c sigma.
jeffreys_one_sided_prob <- function(x, bound, c, side = c("upper", "lower")) {
  side <- match.arg(side)
  n <- length(x)
  x_bar <- mean(x)
  sse <- sum((x - x_bar)^2)
  s <- sqrt(sse / (n - 1))
  shape <- (n - 1) / 2
  rate <- sse / 2

  integrand <- function(z) {
    sigma <- s * z
    # Density of sigma when sigma^2 ~ InvGamma(shape, rate), times ds / dz.
    log_dens <- log(2) + shape * log(rate) - lgamma(shape) -
      (2 * shape + 1) * log(sigma) - rate / sigma^2 + log(s)
    zscore <- if (side == "upper") {
      (bound - 3 * c * sigma - x_bar) * sqrt(n) / sigma
    } else {
      (bound + 3 * c * sigma - x_bar) * sqrt(n) / sigma
    }
    cond <- if (side == "upper") {
      stats::pnorm(zscore, lower.tail = FALSE)
    } else {
      stats::pnorm(zscore)
    }
    exp(log_dens) * cond
  }
  stats::integrate(integrand, 0.5, 2, rel.tol = 1e-10)$value
}

testthat::test_that("narrow posteriors on wide grids keep the far interval empty", {
  set.seed(1)
  x <- rnorm(150, 25.003006666666668, 0.015062953651489967)
  LSL <- 24.95
  USL <- 25.05
  fit <- qc::bpc(
    x, LSL = LSL, USL = USL, target = 25,
    prior = qc::prior_independent(mu = "Jeffreys_mu", sigma = "Jeffreys_sigma"),
    method = "integration"
  )
  probs <- as.data.frame(summary(fit)$interval_summary)

  # Every metric's posterior sits near 1.2 with sd ~0.07, so (2, Inf] is
  # ~10 sd out. Grid quadrature error used to land there (~2e-4 for Cpu).
  far <- stats::setNames(probs[["(2, Inf]"]], probs$metric)
  testthat::expect_true(all(far < 1e-8), info = paste(names(far), far, collapse = ", "))

  # The coarse shoulder cells near 1 (about -2 sd for Cpu) also biased the
  # finite bins: relative errors were ~7e-3 for Cpu and ~1.5e-2 for Cpl.
  rel_error <- function(metric, reference) {
    abs(probs[probs$metric == metric, "[-Inf,1]"] / reference - 1)
  }
  testthat::expect_lt(rel_error("Cpu", jeffreys_one_sided_prob(x, USL, 1, "upper")), 1e-3)
  testthat::expect_lt(rel_error("Cpl", jeffreys_one_sided_prob(x, LSL, 1, "lower")), 5e-3)
})

testthat::test_that("density tail masses measure the right tail directly", {
  pdf_fn <- function(x) stats::dnorm(x, 1.2, 0.07)
  grid_x <- seq(0.9, 1.5, length.out = 101)
  true_grid <- stats::pnorm(1.5, 1.2, 0.07) - stats::pnorm(0.9, 1.2, 0.07)

  # A small midpoint-rule deficit belongs to the grid, not the right tail.
  tails <- qc:::.integration_density_tail_masses(
    pdf_fn, grid_x, grid_mass = true_grid - 4e-4
  )
  testthat::expect_equal(tails$right_tail_mass,
                         stats::pnorm(1.5, 1.2, 0.07, lower.tail = FALSE),
                         tolerance = 1e-4)
  testthat::expect_equal(tails$left_tail_mass, stats::pnorm(0.9, 1.2, 0.07),
                         tolerance = 1e-4)
  testthat::expect_equal(
    tails$mass_nonpositive + tails$left_tail_mass + tails$grid_mass +
      tails$right_tail_mass,
    1
  )

  # A large gap is not quadrature error, so it is kept as the remainder.
  half_pdf <- function(x) 0.5 * pdf_fn(x)
  tails <- qc:::.integration_density_tail_masses(
    half_pdf, grid_x, grid_mass = 0.5 * true_grid
  )
  testthat::expect_equal(tails$grid_mass, 0.5 * true_grid)
  testthat::expect_equal(tails$right_tail_mass,
                         1 - tails$left_tail_mass - 0.5 * true_grid)
})
