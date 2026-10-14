# Generic (PriorGeneric) posteriors for data with a small sd. The posterior is
# much narrower than the prior, so integration must be anchored on the
# posterior itself rather than on a fixed (mu, sigma) box.

# mu ~ Normal(25 k, 2 k); p(sigma) proportional to sigma^-5 exp(-b / sigma^2)
# with b = 0.001 k^2 (an inverse-gamma on sigma^2, written in sigma).
generic_narrow_log_prior <- function(k = 1) {
  b <- 0.001 * k^2
  function(mu, sigma) {
    stats::dnorm(mu, 25 * k, 2 * k, log = TRUE) +
      ifelse(sigma > 0,
             log(2) + 2 * log(b) - lgamma(2) - 5 * log(sigma) - b / sigma^2,
             -Inf)
  }
}

generic_narrow_data <- function() {
  set.seed(1)
  small <- stats::rnorm(8, 25.01, 0.02)
  issue5 <- stats::rnorm(150, 25.003006666666668, 0.015062953651489967)
  list(small = small, issue5 = issue5)
}

# Independent reference: brute-force midpoint grid over (mu, log sigma) built
# from the sufficient statistics, wide enough (+-15 and +-12 local sds) and fine
# enough (800 x 800) to resolve the posterior. Returns posterior summaries of
# each metric and log(integral of p(data | mu, sigma) p(mu, sigma) dmu dsigma)
# up to the constant the likelihood drops (-n/2 log(2 pi)).
generic_narrow_reference <- function(x, log_prior, LSL, USL, target,
                                     n_grid = 800L) {
  n <- length(x)
  x_bar <- mean(x)
  sse <- sum((x - x_bar)^2)
  s <- sqrt(sse / (n - 1))
  mu <- x_bar + seq(-15, 15, length.out = n_grid) * s / sqrt(n)
  t <- log(s) + seq(-12, 12, length.out = n_grid) / sqrt(2 * n)
  grid <- expand.grid(mu = mu, t = t)
  sigma <- exp(grid$t)
  log_w <- -n * grid$t - (sse + n * (grid$mu - x_bar)^2) / (2 * sigma^2) +
    log_prior(grid$mu, sigma) + grid$t
  log_w_max <- max(log_w)
  w <- exp(log_w - log_w_max)
  log_norm <- log(sum(w) * diff(mu[1:2]) * diff(t[1:2])) + log_w_max
  w <- w / sum(w)

  d <- grid$mu - target
  e_abs_dev <- sigma * sqrt(2 / pi) * exp(-d^2 / (2 * sigma^2)) +
    d * (1 - 2 * stats::pnorm(-d / sigma))
  values <- list(
    Cp = (USL - LSL) / (6 * sigma),
    Cpu = (USL - grid$mu) / (3 * sigma),
    Cpl = (grid$mu - LSL) / (3 * sigma),
    Cpk = pmin(USL - grid$mu, grid$mu - LSL) / (3 * sigma),
    Cpm = (USL - LSL) / (6 * sqrt(sigma^2 + d^2)),
    Cpc = (USL - LSL) / (6 * sqrt(pi / 2) * e_abs_dev)
  )
  stats <- t(vapply(values, function(v) {
    ord <- order(v)
    cum_w <- cumsum(w[ord])
    m <- sum(w * v)
    q <- vapply(c(0.025, 0.5, 0.975),
                function(p) v[ord][which(cum_w >= p)[1]], numeric(1))
    c(Mean = m, SD = sqrt(sum(w * (v - m)^2)),
      Q2.5 = q[1], Median = q[2], Q97.5 = q[3])
  }, numeric(5)))
  list(stats = stats, log_norm = log_norm)
}

generic_narrow_stats <- function(x, k, metric) {
  prior <- prior_joint(generic_narrow_log_prior(k))
  res <- analyze_capability_integration(
    x * k, 24.95 * k, 25.05 * k, prior,
    metric = metric, target = 25 * k
  )
  res$stats[c("Mean", "SD", "Q2.5", "Median", "Q97.5")]
}

test_that("generic posterior normalization matches a brute-force grid for narrow data", {
  data <- generic_narrow_data()
  for (name in names(data)) {
    x <- data[[name]]
    ref <- generic_narrow_reference(
      x, generic_narrow_log_prior(), 24.95, 25.05, 25
    )
    state <- precompute_generic_state(x, prior_joint(generic_narrow_log_prior()))
    # Absolute error on the log scale, i.e. relative error of the normalizer.
    expect_lt(abs(log(state$Z) + state$h_max - ref$log_norm), 1e-4)
  }
})

test_that("generic posterior moments and quantiles match a brute-force grid for narrow data", {
  data <- generic_narrow_data()
  for (name in names(data)) {
    x <- data[[name]]
    ref <- generic_narrow_reference(
      x, generic_narrow_log_prior(), 24.95, 25.05, 25
    )$stats
    for (metric in c("Cp", "Cpu", "Cpl", "Cpk", "Cpm", "Cpc")) {
      moments <- compute_metric_moments(
        x, 24.95, 25.05, prior_joint(generic_narrow_log_prior()),
        metric = metric, target = 25
      )
      label <- paste(name, metric)
      expect_equal(moments$mean, ref[metric, "Mean"], tolerance = 1e-3,
                   label = paste(label, "moment mean"))
      expect_equal(moments$sd, ref[metric, "SD"], tolerance = 5e-3,
                   label = paste(label, "moment sd"))

      stats <- generic_narrow_stats(x, 1, metric)
      expect_equal(unname(stats[["Mean"]]), ref[metric, "Mean"],
                   tolerance = 1e-3, label = paste(label, "mean"))
      expect_equal(unname(stats[["SD"]]), ref[metric, "SD"],
                   tolerance = 5e-3, label = paste(label, "sd"))
      # Grid quantiles are resolved to about 1e-3 in the metric.
      for (q in c("Q2.5", "Median", "Q97.5")) {
        expect_lt(abs(stats[[q]] - ref[metric, q]), 5e-3)
      }
    }
  }
})

test_that("generic narrow-posterior Cpc survival quantiles match the reported reference", {
  x <- generic_narrow_data()$issue5
  stats <- generic_narrow_stats(x, 1, "Cpc")
  expect_lt(abs(stats[["Median"]] - 1.162), 2e-3)
  expect_lt(abs(stats[["Q2.5"]] - 1.034), 2e-3)
  expect_lt(abs(stats[["Q97.5"]] - 1.296), 2e-3)
})

test_that("generic posterior capability indices are invariant to the measurement scale", {
  x <- generic_narrow_data()$issue5
  for (metric in c("Cp", "Cpu", "Cpl", "Cpk", "Cpm", "Cpc")) {
    base <- generic_narrow_stats(x, 1, metric)
    for (k in c(1e-3, 1e3)) {
      expect_equal(generic_narrow_stats(x, k, metric), base,
                   tolerance = 1e-5,
                   label = sprintf("%s at scale k = %g", metric, k))
    }
  }
})
