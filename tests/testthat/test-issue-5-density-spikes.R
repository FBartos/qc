# Independent reference: integrate the conditional Normal density against the
# one-dimensional sigma posterior, using sigma / sd(x) rather than log(sigma).
issue5_density_reference <- function(x, prior_name) {
  n <- length(x)
  x_bar <- mean(x)
  sse <- sum((x - x_bar)^2)
  sigma_ref <- stats::sd(x)

  if (prior_name == "DCSI") {
    k_n <- n + 2
    mu_n <- (n * x_bar + 2 * 25) / k_n
    alpha_n <- n / 2 + 3
    beta_n <- 3 * (0.1 / 7.5)^2 + sse / 2 +
      n * (x_bar - 25)^2 / k_n
    log_sigma_weight <- function(sigma) {
      -(2 * alpha_n + 1) * log(sigma) - beta_n / sigma^2
    }
    conditional_mu <- function(sigma) {
      list(mean = mu_n, sd = sigma / sqrt(k_n))
    }
    sigma_upper <- Inf
  } else {
    # mu ~ Normal(25, 1), sigma ~ Uniform(0, 1); integrate mu out exactly.
    log_sigma_weight <- function(sigma) {
      -(n - 1) * log(sigma) - sse / (2 * sigma^2) -
        log(n + sigma^2) / 2 - n * (x_bar - 25)^2 / (2 * (n + sigma^2))
    }
    conditional_mu <- function(sigma) {
      variance <- sigma^2 / (n + sigma^2)
      list(mean = (n * x_bar + 25 * sigma^2) / (n + sigma^2),
           sd = sqrt(variance))
    }
    sigma_upper <- 1
  }

  offset <- log_sigma_weight(sigma_ref)
  weight <- function(sigma) exp(log_sigma_weight(sigma) - offset)
  integrate_scaled <- function(fn, upper) {
    # At n = 150, sigma > 4*sd(x) has negligible posterior mass. Keeping
    # this reference interval compact also resolves tail contour peaks.
    upper <- min(upper, 4)
    split <- min(1, upper)
    first <- stats::integrate(fn, 0, split, rel.tol = 1e-8)$value
    second <- if (upper > split) {
      stats::integrate(fn, split, upper, rel.tol = 1e-8)$value
    } else 0
    first + second
  }
  Z <- integrate_scaled(function(z) weight(sigma_ref * z), sigma_upper / sigma_ref)

  function(c, metric) {
    vapply(c, function(ci) {
      upper <- sigma_upper
      if (metric == "Cpk" && ci > 0) upper <- min(upper, 0.1 / (6 * ci))
      integrate_scaled(function(z) {
        sigma <- sigma_ref * z
        cond <- conditional_mu(sigma)
        upper_pdf <- stats::dnorm(25.05 - 3 * ci * sigma, cond$mean, cond$sd)
        lower_pdf <- stats::dnorm(24.95 + 3 * ci * sigma, cond$mean, cond$sd)
        mu_pdf <- switch(metric, Cpu = upper_pdf, Cpl = lower_pdf,
                         Cpk = upper_pdf + lower_pdf)
        3 * sigma * mu_pdf * weight(sigma) / Z
      }, upper / sigma_ref)
    }, numeric(1))
  }
}

for (prior_name in c("DCSI", "independent")) {
  testthat::test_that(paste("issue #5 has smooth, accurate densities with", prior_name), {
    set.seed(1)
    x <- rnorm(150, mean = 25.003006666666668, sd = 0.015062953651489967)
    prior <- if (prior_name == "DCSI") "DCSI" else qc::prior_independent(
      mu = qc::prior("normal", list(mean = 25, sd = 1)),
      sigma = qc::prior("uniform", list(a = 0, b = 1))
    )
    # Keep the issue's arguments, including its implicit integration default.
    fit <- qc::bpc(x, chains = 1, warmup = 1000, iter = 5000, silent = TRUE,
                   seed = 1, target = 25, LSL = 24.95, USL = 25.05,
                   prior = prior, sample_priors = FALSE)
    s <- summary(fit, ci.level = 0.95)
    reference <- issue5_density_reference(x, prior_name)
    distributions <- attr(s$metrics, "distributions")
    for (metric in c("Cpu", "Cpl", "Cpk")) {
      entry <- distributions[[metric]]$entry
      g <- entry$grid
      testthat::expect_gt(entry$area, 0.99)

      # Compare every point, normalizing the reference on the same finite
      # plotting grid (as the package does).
      expected <- reference(g$x, metric)
      expected <- expected / sum(expected * g$width)
      testthat::expect_equal(g$density, expected, tolerance = 1e-5, info = metric)

      built <- ggplot2::ggplot_build(qc::plot_density(s, what = metric))
      plotted <- built$data[[1]]$y[match(g$x, built$data[[1]]$x)]
      testthat::expect_equal(plotted, expected, tolerance = 1e-5, info = metric)
    }
  })
}
