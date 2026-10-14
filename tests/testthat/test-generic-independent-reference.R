# Regression tests for the generic (independent-prior) integration backend on
# narrow posteriors (GitHub issue #5): mu ~ Normal(25, 1), sigma ~ Uniform(0, 1),
# at n = 150 and n = 5000.
#
# The reference integrates mu out analytically (Normal prior x Normal
# likelihood), leaving a one-dimensional sigma posterior handled by
# stats::integrate. Metric formulas are written out here rather than taken
# from the package.

narrow_reference <- function(x, LSL = 24.95, USL = 25.05, target = 25) {
  n <- length(x)
  x_bar <- mean(x)
  sse <- sum((x - x_bar)^2)
  s_ref <- stats::sd(x)
  d <- (USL - LSL) / 2
  M <- (USL + LSL) / 2

  log_w <- function(sigma) {
    -(n - 1) * log(sigma) - sse / (2 * sigma^2) -
      log(n + sigma^2) / 2 - n * (x_bar - 25)^2 / (2 * (n + sigma^2))
  }
  # mu | sigma, data ~ Normal(cond_mean, cond_sd)
  cond_mean <- function(sigma) (n * x_bar + 25 * sigma^2) / (n + sigma^2)
  cond_sd <- function(sigma) sqrt(sigma^2 / (n + sigma^2))

  # The sigma posterior has sd ~ s_ref / sqrt(2n); +-12 of those leaves
  # negligible mass outside (capped at s_ref * [0.5, 2] for small n).
  half_width <- min(0.5, 12 / sqrt(2 * n))
  lo <- s_ref * (1 - half_width)
  hi <- s_ref * (1 + 2 * half_width)
  offset <- log_w(s_ref)
  w <- function(sigma) exp(log_w(sigma) - offset)
  # `breaks` split the range at sharp features of fn so that QUADPACK
  # cannot step over them.
  integrate_sigma <- function(fn, upper = hi, breaks = numeric(0)) {
    upper <- min(upper, hi)
    if (upper <= lo) return(0)
    edges <- sort(unique(c(lo, breaks[breaks > lo & breaks < upper], upper)))
    sum(vapply(seq_len(length(edges) - 1L), function(i) {
      stats::integrate(function(s) fn(s) * w(s), edges[i], edges[i + 1L],
                       rel.tol = 1e-10, subdivisions = 1000L)$value
    }, numeric(1)))
  }
  Z <- integrate_sigma(function(s) 1)
  post_mean <- function(fn, upper = hi, breaks = numeric(0)) {
    integrate_sigma(fn, upper, breaks) / Z
  }

  # Conditional raw moments E[metric^k | sigma] for the linear metrics.
  folded_mean <- function(s) {
    delta <- cond_mean(s) - M
    v <- cond_sd(s)
    v * sqrt(2 / pi) * exp(-delta^2 / (2 * v^2)) +
      delta * (1 - 2 * stats::pnorm(-delta / v))
  }
  conditional <- list(
    Cp = list(m1 = function(s) (USL - LSL) / (6 * s),
              m2 = function(s) ((USL - LSL) / (6 * s))^2),
    Cpu = list(m1 = function(s) (USL - cond_mean(s)) / (3 * s),
               m2 = function(s) ((USL - cond_mean(s))^2 + cond_sd(s)^2) / (9 * s^2)),
    Cpl = list(m1 = function(s) (cond_mean(s) - LSL) / (3 * s),
               m2 = function(s) ((cond_mean(s) - LSL)^2 + cond_sd(s)^2) / (9 * s^2)),
    Cpk = list(m1 = function(s) (d - folded_mean(s)) / (3 * s),
               m2 = function(s) {
                 (d^2 - 2 * d * folded_mean(s) +
                    (cond_mean(s) - M)^2 + cond_sd(s)^2) / (9 * s^2)
               })
  )

  moments <- function(metric) {
    m1 <- post_mean(conditional[[metric]]$m1)
    m2 <- post_mean(conditional[[metric]]$m2)
    c(mean = m1, sd = sqrt(m2 - m1^2))
  }

  # P(metric > c), with the mu constraint handled by the Normal CDF.
  survival <- function(c, metric) {
    breaks <- c((USL - x_bar) / (3 * c), (x_bar - LSL) / (3 * c))
    switch(metric,
      Cpu = post_mean(function(s) {
        stats::pnorm(USL - 3 * c * s, cond_mean(s), cond_sd(s))
      }, breaks = breaks),
      Cpl = post_mean(function(s) {
        stats::pnorm(LSL + 3 * c * s, cond_mean(s), cond_sd(s), lower.tail = FALSE)
      }, breaks = breaks),
      Cpk = post_mean(function(s) {
        stats::pnorm(USL - 3 * c * s, cond_mean(s), cond_sd(s)) -
          stats::pnorm(LSL + 3 * c * s, cond_mean(s), cond_sd(s))
      }, upper = (USL - LSL) / (6 * c), breaks = breaks)
    )
  }

  # Density of the metric at c > 0. Cp: change of variables
  # sigma = (USL - LSL) / (6c). Cpu/Cpl/Cpk: mu = USL - 3c sigma (or
  # LSL + 3c sigma) has Jacobian 3 sigma; Cpk takes either branch while the
  # other index exceeds c, i.e. sigma < (USL - LSL) / (6c).
  density <- function(c, metric) {
    if (metric == "Cp") {
      s <- (USL - LSL) / (6 * c)
      return(w(s) * s / c / Z)
    }
    vapply(c, function(ci) {
      upper <- if (metric == "Cpk") (USL - LSL) / (6 * ci) else hi
      # The integrand peaks where USL - 3c sigma (or LSL + 3c sigma) meets
      # the conditional mean of mu; for large n that peak is very narrow.
      breaks <- c((USL - x_bar) / (3 * ci), (x_bar - LSL) / (3 * ci))
      post_mean(function(s) {
        upper_pdf <- stats::dnorm(USL - 3 * ci * s, cond_mean(s), cond_sd(s))
        lower_pdf <- stats::dnorm(LSL + 3 * ci * s, cond_mean(s), cond_sd(s))
        3 * s * switch(metric, Cpu = upper_pdf, Cpl = lower_pdf,
                       Cpk = upper_pdf + lower_pdf)
      }, upper = upper, breaks = breaks)
    }, numeric(1))
  }

  # Exact posterior draws (inverse CDF of the sigma marginal on a fine grid,
  # then mu | sigma) for the non-linear metrics Cpm and Cpc.
  draws <- function(n_draws = 2e5, seed = 123) {
    if (exists(".Random.seed", envir = globalenv())) {
      old_seed <- get(".Random.seed", envir = globalenv())
      on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
    }
    set.seed(seed)
    grid <- seq(lo, hi, length.out = 1e5 + 1)
    cdf <- cumsum(w(grid))
    cdf <- cdf / cdf[length(cdf)]
    sigma <- stats::approx(cdf, grid, stats::runif(n_draws), ties = "ordered", rule = 2)$y
    mu <- stats::rnorm(n_draws, cond_mean(sigma), cond_sd(sigma))
    delta <- mu - target
    abs_dev <- sigma * sqrt(2 / pi) * exp(-(delta / sigma)^2 / 2) +
      delta * (1 - 2 * stats::pnorm(-delta / sigma))
    list(
      Cpm = (USL - LSL) / (6 * sqrt(sigma^2 + delta^2)),
      Cpc = (USL - LSL) / (6 * sqrt(pi / 2) * abs_dev)
    )
  }

  list(moments = moments, survival = survival, density = density, draws = draws)
}

fit_narrow_independent <- function(n) {
  set.seed(1)
  x <- rnorm(n, 25.003, 0.015)
  prior <- qc::prior_independent(
    mu = qc::prior("normal", list(mean = 25, sd = 1)),
    sigma = qc::prior("uniform", list(a = 0, b = 1))
  )
  fit <- qc::bpc(x, chains = 1, warmup = 200, iter = 400, silent = TRUE,
                 seed = 1, target = 25, LSL = 24.95, USL = 25.05,
                 prior = prior, sample_priors = FALSE)
  list(x = x, summary = summary(fit, ci.level = 0.95))
}

summary_value <- function(table, metric, columns) {
  row <- table[unlist(table[[1]]) == metric, , drop = FALSE]
  unname(vapply(columns, function(col) as.numeric(row[[col]]), numeric(1)))
}

check_narrow_posterior <- function(n) {
  fit <- fit_narrow_independent(n)
  s <- fit$summary
  ref <- narrow_reference(fit$x)
  dists <- attr(s$metrics, "distributions")

  # Moments: mean and SD agree with the analytic sigma-marginal reference,
  # and the SD is not collapsed to zero.
  for (metric in c("Cp", "Cpu", "Cpl", "Cpk")) {
    expected <- ref$moments(metric)
    expect_equal(summary_value(s$summary, metric, "mean"), expected[["mean"]],
                 tolerance = 1e-4, info = paste(n, metric, "mean"))
    expect_equal(summary_value(s$summary, metric, "sd"), expected[["sd"]],
                 tolerance = 1e-3, info = paste(n, metric, "sd"))
  }

  # Interval probabilities (read off the finite density grid, hence the
  # absolute tolerance).
  for (metric in c("Cpu", "Cpl", "Cpk")) {
    expected <- ref$survival(1, metric) - ref$survival(1.33, metric)
    actual <- summary_value(s$interval_summary, metric, "(1,1.33]")
    expect_lt(abs(actual - expected), 2e-3, label = paste(n, metric, "interval error"))
  }

  # Every metric has a normalizable density grid.
  for (metric in c("Cp", "Cpu", "Cpl", "Cpk", "Cpm", "Cpc")) {
    expect_gt(dists[[metric]]$entry$area, 0.99)
  }

  # Direct density solvers: compare pointwise, normalizing the reference on
  # the same finite grid (as the package does).
  for (metric in c("Cp", "Cpu", "Cpl", "Cpk")) {
    g <- dists[[metric]]$entry$grid
    expected <- ref$density(g$x, metric)
    expected <- expected / sum(expected * g$width)
    expect_equal(g$density, expected, tolerance = 1e-3,
                 info = paste(n, metric, "density"))
  }

  # Non-linear metrics against exact posterior draws.
  draws <- ref$draws()
  for (metric in c("Cpm", "Cpc")) {
    sample_sd <- stats::sd(draws[[metric]])
    expected <- c(
      mean = mean(draws[[metric]]),
      stats::quantile(draws[[metric]], c(0.025, 0.5, 0.975))
    )
    actual <- summary_value(s$summary, metric, c("mean", "lower", "median", "upper"))
    expect_lt(max(abs(actual - expected)) / sample_sd, 0.05,
              label = paste(n, metric, "location error / sd"))
    expect_equal(summary_value(s$summary, metric, "sd"), sample_sd,
                 tolerance = 0.01, info = paste(n, metric, "sd"))
  }
}

test_that("independent prior is accurate for a narrow posterior (n = 150)", {
  check_narrow_posterior(150)
})

test_that("independent prior normalizes a very narrow posterior (n = 5000)", {
  check_narrow_posterior(5000)
})
