# Independent references for the Jeffreys conjugate prior: the posterior is
# sigma^2 ~ InvGamma((n - 1) / 2, sse / 2), mu | sigma ~ Normal(x_bar, sigma^2 / n),
# so the marginal of mu is x_bar + s / sqrt(n) * t_{n - 1}.
narrow_jeffreys_draws <- function(x, n_draws, seed) {
  n <- length(x)
  sse <- sum((x - mean(x))^2)
  set.seed(seed)
  sigma <- sqrt((sse / 2) / stats::rgamma(n_draws, shape = (n - 1) / 2))
  mu <- stats::rnorm(n_draws, mean(x), sigma / sqrt(n))
  list(mu = mu, sigma = sigma)
}

# Cpk <= 0 exactly when mu lies outside [LSL, USL].
nonpositive_cpk_prob <- function(x, LSL, USL) {
  n <- length(x)
  se <- stats::sd(x) / sqrt(n)
  stats::pt((LSL - mean(x)) / se, n - 1) +
    stats::pt((mean(x) - USL) / se, n - 1)
}

testthat::test_that("a narrow Cpk posterior keeps its median and interval", {
  # n = 50000 gives a posterior sd of ~0.0036, far narrower than the
  # coarse grid spacing used to locate the peak.
  x <- 25.003 + 0.015 * stats::qnorm(stats::ppoints(50000))
  LSL <- 24.95
  USL <- 25.05
  result <- qc:::analyze_capability_integration(
    x, LSL, USL, qc::prior_conjugate(), metric = "Cpk"
  )

  draws <- narrow_jeffreys_draws(x, 4e5, seed = 1)
  cpk <- pmin(USL - draws$mu, draws$mu - LSL) / (3 * draws$sigma)
  expected <- stats::quantile(cpk, c(0.025, 0.5, 0.975), names = FALSE)
  tolerance <- 0.1 * stats::sd(cpk)

  testthat::expect_lt(abs(result$stats[["Q2.5"]] - expected[1]), tolerance)
  testthat::expect_lt(abs(result$stats[["Median"]] - expected[2]), tolerance)
  testthat::expect_lt(abs(result$stats[["Q97.5"]] - expected[3]), tolerance)

  # The density must resolve the peak: a Normal with this sd peaks at
  # 1 / (sd * sqrt(2 * pi)); a collapsed grid flattened it by ~40%.
  peak <- 1 / (stats::sd(cpk) * sqrt(2 * pi))
  testthat::expect_equal(max(result$grid$density), peak, tolerance = 0.05)
})

testthat::test_that("density tail masses fall back to the survival function", {
  x <- c(8.4, 9.0, 9.3, 9.6, 9.8, 10.1, 10.5, 10.9)
  LSL <- 0
  USL <- 10
  request <- qc:::.new_qc_integration_request(
    data = x, LSL = LSL, USL = USL,
    prior = qc::prior_conjugate(), metric = "Cpk"
  )
  pdf <- qc:::.integration_make_density_solver(request = request)
  make_S <- function() qc:::.integration_make_solver(request = request)

  # A grid that starts above zero, so P(Cpk <= 0) must come from the tails.
  edges <- seq(0.1, 3, length.out = 201)
  mids <- (edges[-1] + edges[-201]) / 2
  dens <- pdf(mids)
  area <- sum(dens * diff(edges))
  grid_df <- data.frame(
    x = mids, density = dens / area,
    x_left = edges[-201], x_right = edges[-1], width = diff(edges)
  )
  expected <- nonpositive_cpk_prob(x, LSL, USL)

  entry_for <- function(pdf_fn, make_S) {
    qc:::.integration_grid_entry(
      metric = "Cpk", grid_df = grid_df, area = area,
      backend = list(mode = "density", pdf_fn = pdf_fn, S = NULL),
      metric_can_be_negative = TRUE, make_S = make_S
    )
  }

  direct <- entry_for(pdf, make_S)
  testthat::expect_equal(direct$mass_nonpositive, expected, tolerance = 1e-3)

  # A density that cannot be evaluated below zero must not move that mass
  # into the right tail.
  broken <- function(c) {
    vals <- pdf(c)
    vals[c < 0] <- NA_real_
    vals
  }
  testthat::expect_no_warning(fallback <- entry_for(broken, make_S))
  testthat::expect_equal(fallback$mass_nonpositive, expected, tolerance = 1e-3)
  testthat::expect_equal(fallback$right_tail_mass, direct$right_tail_mass,
                         tolerance = 1e-3)

  testthat::expect_warning(entry_for(broken, NULL), "could not be integrated")
})

testthat::test_that("analysis falls back to the survival grid when the density fails", {
  x <- c(8.4, 9.0, 9.3, 9.6, 9.8, 10.1, 10.5, 10.9)
  LSL <- 0
  USL <- 10
  real_solver <- qc:::.integration_make_density_solver
  testthat::local_mocked_bindings(
    .integration_make_density_solver = function(...) {
      pdf <- real_solver(...)
      function(c) {
        vals <- pdf(c)
        vals[c < 0] <- NA_real_
        vals
      }
    },
    .package = "qc"
  )

  result <- qc:::analyze_capability_integration(
    x, LSL, USL, qc::prior_conjugate(), metric = "Cpk"
  )
  probs <- qc:::.integration_direct_interval_probs(result, c(-Inf, 0, Inf))
  survival <- qc:::analyze_capability_integration(
    x, LSL, USL, qc::prior_conjugate(), metric = "Cpk",
    use_density_solver = FALSE
  )
  survival_probs <- qc:::.integration_direct_interval_probs(survival, c(-Inf, 0, Inf))

  testthat::expect_true(all(is.finite(result$grid$density)))
  # The survival grid omits mass left of its first point, so it is ~1% low.
  testthat::expect_equal(probs[[1]], nonpositive_cpk_prob(x, LSL, USL),
                         tolerance = 0.02)
  testthat::expect_equal(probs, survival_probs, tolerance = 1e-6)
})
