# Semi-conjugate sigma prior: mu ~ Normal(25, 2), sigma^2 ~ InvGamma(2, 0.001).
# Measurements with sd ~ 0.015 give a mu posterior far narrower than the
# prior (sd ~ 0.0012 at n = 150), which the integrator must resolve.
semi_sigma_narrow_prior <- function(k = 1) {
  qc::prior_semi_conjugate(
    "sigma", alpha0 = 2, beta0 = 0.001 * k^2,
    log_dens_mu = function(mu) stats::dnorm(mu, 25 * k, 2 * k, log = TRUE)
  )
}

semi_sigma_narrow_data <- function() {
  set.seed(1)
  list(
    small = stats::rnorm(8, 25.01, 0.02),
    issue5 = stats::rnorm(150, 25.003006666666668, 0.015062953651489967)
  )
}

# Independent reference: the joint (mu, sigma) posterior evaluated on a
# brute-force midpoint grid over (mu, log sigma), from the Normal likelihood
# and the two priors written out directly. Each cell carries its posterior
# mass; metric summaries are weighted means and interpolated weighted
# quantiles over the cells.
semi_sigma_grid_reference <- function(x, LSL = 24.95, USL = 25.05,
                                      target = 25, alpha0 = 2, beta0 = 0.001,
                                      n_mu = 1000L, n_sigma = 1000L) {
  n <- length(x)
  x_bar <- mean(x)
  sse <- sum((x - x_bar)^2)
  s <- stats::sd(x)

  mu_half <- 12 * s / sqrt(n)
  log_sigma_half <- 10 / sqrt(2 * n)
  mu_edges <- seq(x_bar - mu_half, x_bar + mu_half, length.out = n_mu + 1L)
  ls_edges <- seq(log(s) - log_sigma_half, log(s) + log_sigma_half,
                  length.out = n_sigma + 1L)
  mu_mid <- (mu_edges[-1L] + mu_edges[-(n_mu + 1L)]) / 2
  ls_mid <- (ls_edges[-1L] + ls_edges[-(n_sigma + 1L)]) / 2

  mu <- rep(mu_mid, times = n_sigma)
  sigma <- exp(rep(ls_mid, each = n_mu))

  # Normal log likelihood + log Normal(25, 2) for mu + log InvGamma(alpha0,
  # beta0) for sigma^2 expressed on log(sigma) (Jacobian 2 * sigma^2).
  log_post <- -n * log(sigma) - (sse + n * (mu - x_bar)^2) / (2 * sigma^2) +
    stats::dnorm(mu, 25, 2, log = TRUE) +
    (-(alpha0 + 1) * log(sigma^2) - beta0 / sigma^2) + 2 * log(sigma)
  w <- exp(log_post - max(log_post))
  w <- w / sum(w)

  on_edge <- rep(seq_len(n_mu) %in% c(1L, n_mu), times = n_sigma) |
    rep(seq_len(n_sigma) %in% c(1L, n_sigma), each = n_mu)

  half_tol <- (USL - LSL) / 2
  d_target <- min(USL - target, target - LSL)
  delta <- abs(mu - target)
  mean_abs_dev <- sigma * sqrt(2 / pi) * exp(-delta^2 / (2 * sigma^2)) +
    delta * (1 - 2 * stats::pnorm(-delta / sigma))
  metrics <- list(
    Cp = half_tol / (3 * sigma),
    Cpu = (USL - mu) / (3 * sigma),
    Cpl = (mu - LSL) / (3 * sigma),
    Cpk = pmin(USL - mu, mu - LSL) / (3 * sigma),
    Cpm = d_target / (3 * sqrt(sigma^2 + (mu - target)^2)),
    Cpc = half_tol / (3 * sqrt(pi / 2) * mean_abs_dev)
  )

  summaries <- lapply(metrics, function(value) {
    ord <- order(value)
    cdf <- cumsum(w[ord]) - w[ord] / 2
    quant <- stats::approx(cdf, value[ord], c(0.025, 0.5, 0.975),
                           ties = "ordered")$y
    m <- sum(w * value)
    c(Mean = m, SD = sqrt(sum(w * (value - m)^2)),
      Q2.5 = quant[1], Median = quant[2], Q97.5 = quant[3])
  })

  list(summaries = summaries, edge_mass = sum(w[on_edge]))
}

testthat::test_that("semi-conjugate sigma analysis matches a brute-force posterior grid", {
  data <- semi_sigma_narrow_data()
  prior <- semi_sigma_narrow_prior()
  metrics <- c("Cp", "Cpu", "Cpl", "Cpk", "Cpm", "Cpc")

  for (name in names(data)) {
    reference <- semi_sigma_grid_reference(data[[name]])
    # The reference grid itself must contain the posterior.
    testthat::expect_lt(reference$edge_mass, 1e-6)

    for (metric in metrics) {
      label <- paste(name, metric)
      expected <- reference$summaries[[metric]]
      result <- qc:::analyze_capability_integration(
        data[[name]], 24.95, 25.05, prior,
        metric = metric, target = 25
      )
      stats <- result$stats
      scale <- expected[["SD"]]

      testthat::expect_gt(result$area, 0.99, label = label)
      testthat::expect_lt(abs(stats[["Mean"]] - expected[["Mean"]]), 0.02 * scale,
                          label = paste(label, "mean"))
      testthat::expect_lt(abs(stats[["SD"]] - expected[["SD"]]), 0.02 * scale,
                          label = paste(label, "sd"))
      for (q in c("Q2.5", "Median", "Q97.5")) {
        testthat::expect_lt(abs(stats[[q]] - expected[[q]]), 0.05 * scale,
                            label = paste(label, q))
      }
    }
  }
})

testthat::test_that("semi-conjugate sigma results do not depend on measurement units", {
  # Capability indices are dimensionless: rescaling data, limits, target and
  # the priors (mu prior by k, the sigma^2 rate by k^2) must not move them.
  x <- semi_sigma_narrow_data()$issue5
  analyze <- function(k, metric) {
    qc:::analyze_capability_integration(
      x * k, 24.95 * k, 25.05 * k, semi_sigma_narrow_prior(k),
      metric = metric, target = 25 * k
    )$stats[c("Mean", "SD", "Q2.5", "Median", "Q97.5")]
  }

  for (metric in c("Cpk", "Cpm")) {
    base <- analyze(1, metric)
    testthat::expect_equal(analyze(1e-3, metric), base, tolerance = 1e-4,
                           label = paste(metric, "k = 1e-3"))
    testthat::expect_equal(analyze(1e3, metric), base, tolerance = 1e-4,
                           label = paste(metric, "k = 1e3"))
  }
})

testthat::test_that("semi-conjugate sigma density solver resolves a narrow posterior", {
  x <- semi_sigma_narrow_data()$issue5
  prior <- semi_sigma_narrow_prior()
  reference <- semi_sigma_grid_reference(x)

  for (metric in c("Cp", "Cpm")) {
    pdf <- qc:::make_density_solver(x, 24.95, 25.05, prior,
                                    metric = metric, target = 25)
    expected <- reference$summaries[[metric]]

    # Midpoint rule over the reference central 95% interval.
    edges <- seq(expected[["Q2.5"]], expected[["Q97.5"]], length.out = 201L)
    mids <- (edges[-1L] + edges[-201L]) / 2
    mass <- sum(vapply(mids, pdf, numeric(1))) * diff(edges[1:2])

    testthat::expect_equal(mass, 0.95, tolerance = 0.005, label = metric)
  }
})
