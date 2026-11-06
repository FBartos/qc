testthat::test_that("adaptive density grids evaluate vectorized density batches", {
  call_lengths <- integer()
  pdf_fn <- function(x) {
    call_lengths <<- c(call_lengths, length(x))
    stats::dnorm(x)
  }

  result <- qc:::.adaptive_density_grid(pdf_fn, 0, 3, n_grid = 32L)

  testthat::expect_true(any(call_lengths > 1L))
  testthat::expect_lt(length(call_lengths), length(result$mid_x))
  testthat::expect_true(all(is.finite(result$pdf_vals)))
})

testthat::test_that("survival-derived densities evaluate vectorized survival batches", {
  call_lengths <- integer()
  S <- function(c) {
    call_lengths <<- c(call_lengths, length(c))
    ifelse(c <= 0, 1, exp(-c))
  }

  pdf_fn <- qc:::.density_from_survival_fn(S, support_lower = 0)
  vals <- pdf_fn(c(0, 0.5, 1))

  testthat::expect_true(any(call_lengths > 1L))
  testthat::expect_equal(length(vals), 3L)
  testthat::expect_true(all(is.finite(vals) & vals >= 0))
})

testthat::test_that("density backends used by adaptive grids accept vector inputs", {
  testthat::skip_if_not_installed("BayesTools")
  testthat::skip_if_not_installed("cubature")

  set.seed(11)
  x <- stats::rnorm(8)
  cases <- list(
    conjugate = list(
      prior = qc:::create_prior_conjugate(),
      metrics = c("Cp", "Cpm", "Cpc")
    ),
    semi_mu = list(
      prior = qc:::.bayestools_to_integration_prior(
        "Jeffreys_mu",
        BayesTools::prior("gamma", list(2, 1))
      )$prior,
      metrics = c("Cp", "Cpm", "Cpc")
    ),
    generic = list(
      prior = qc:::.bayestools_to_integration_prior(
        BayesTools::prior("normal", list(0, 1)),
        BayesTools::prior("gamma", list(2, 1))
      )$prior,
      metrics = c("Cp", "Cpm")
    )
  )

  points <- c(0.25, 0.5, 1)
  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    for (metric in case$metrics) {
      pdf_fn <- qc:::.integration_make_density_solver(
        data = x,
        LSL = -3,
        USL = 3,
        target = 0,
        prior = case$prior,
        metric = metric
      )
      vals <- pdf_fn(points)
      testthat::expect_equal(length(vals), length(points), info = paste(case_name, metric))
      testthat::expect_true(all(is.finite(vals) & vals >= 0), info = paste(case_name, metric))
    }
  }
})

testthat::test_that("conjugate smooth-metric survival solvers accept vector thresholds", {
  set.seed(12)
  x <- stats::rnorm(8)
  prior <- qc:::create_prior_conjugate()

  for (metric in c("Cpm", "Cpc")) {
    solver <- qc:::.integration_make_solver(
      data = x,
      LSL = -3,
      USL = 3,
      target = 0,
      prior = prior,
      metric = metric
    )
    vals <- solver(c(0.25, 0.5, 1))
    testthat::expect_equal(length(vals), 3L, info = metric)
    testthat::expect_true(all(is.finite(vals)), info = metric)
    testthat::expect_true(all(diff(vals) <= 1e-8), info = metric)
  }
})

conjugate_log_joint <- function(post) {
  log_Z <- lgamma(post$alpha_n) - log(2) - post$alpha_n * log(post$beta_n)
  function(mu, sigma) {
    -(2 * post$alpha_n + 1) * log(sigma) - post$beta_n / sigma^2 - log_Z +
      stats::dnorm(mu, post$mu_n, sigma / sqrt(post$k_n), log = TRUE)
  }
}

testthat::test_that("linear contour densities evaluate the kernel in batches", {
  post <- list(mu_n = 10, k_n = 30, alpha_n = 15, beta_n = 15)
  log_joint <- conjugate_log_joint(post)
  n_calls <- 0L
  pdf_fn <- qc:::.linear_metric_density(
    function(mu, sigma) {
      n_calls <<- n_calls + 1L
      log_joint(mu, sigma)
    },
    LSL = 6, USL = 14, metric = "Cpk", sigma_level = 3, sigma_ref = 1
  )

  pdf_fn(1.3)
  calls_single <- n_calls
  n_calls <- 0L
  vals <- pdf_fn(seq(0.8, 1.9, length.out = 400L))

  testthat::expect_lte(n_calls, calls_single)
  testthat::expect_equal(vals[200], pdf_fn(seq(0.8, 1.9, length.out = 400L)[200]))
})

testthat::test_that("linear contour densities match survival differences", {
  set.seed(21)
  x <- stats::rnorm(30, 11, 1.2)
  prior <- qc:::create_prior_conjugate(10, 2, 3, 3)
  state <- qc:::.as_qc_suff_stats_state(data = x)
  S <- qc:::make_solver(x, 6, 14, prior, metric = "Cpk", cached_state = state)
  pdf_fn <- qc:::make_density_solver(
    x, 6, 14, prior, metric = "Cpk", cached_state = state
  )

  edges <- seq(0.3, 1.6, length.out = 501L)
  mids <- (edges[-1L] + edges[-501L]) / 2
  expected <- (S(edges[-501L]) - S(edges[-1L])) / diff(edges)

  testthat::expect_equal(pdf_fn(mids), expected, tolerance = 1e-3)
})

testthat::test_that("linear metric summaries survive very large samples", {
  n <- 1e6L
  state <- qc:::.new_qc_suff_stats_state(n = n, x_bar = 10, sse = n)
  prior <- qc:::create_prior_conjugate(10, 1, 3, 3)

  for (metric in c("Cpu", "Cpk")) {
    request <- qc:::.new_qc_integration_request(
      numeric(0), 6, 14, prior, metric, 10, state, 3
    )
    entry <- qc:::analyze_capability_integration(request = request)

    posterior_sd <- sqrt(1 / (9 * n) + (4 / 3)^2 / (2 * n))
    testthat::expect_equal(entry$area, 1, tolerance = 1e-3)
    testthat::expect_lt(abs(unname(entry$stats["Median"]) - 4 / 3), 3 * posterior_sd)
    testthat::expect_lt(
      abs(unname(entry$stats["Q97.5"]) - (4 / 3 + 1.96 * posterior_sd)),
      posterior_sd
    )
  }
})

testthat::test_that("contour densities respect a hard truncation of the mean prior", {
  testthat::skip_if_not_installed("BayesTools")

  set.seed(22)
  x <- stats::rnorm(100, 10, 1)
  prior <- qc:::.bayestools_to_integration_prior(
    BayesTools::prior("normal", list(10, 1), list(9.95, 10.05)),
    "Jeffreys_sigma"
  )$prior
  state <- qc:::.as_qc_suff_stats_state(data = x)
  request <- qc:::.new_qc_integration_request(x, 6, 14, prior, "Cpu", 10, state, 3)
  entry <- qc:::analyze_capability_integration(request = request)
  S <- qc:::.integration_make_solver(request = request)

  testthat::expect_equal(entry$area, 1, tolerance = 1e-3)
  testthat::expect_equal(S(unname(entry$stats["Median"])), 0.5, tolerance = 1e-3)
})

testthat::test_that("conjugate Cpm and Cpc moments match Monte Carlo draws", {
  draws <- 2e6L
  cases <- list(
    list(mu_n = 10.4, k_n = 30, alpha_n = 15, beta_n = 18, target = 10),
    list(mu_n = 12, k_n = 0.05, alpha_n = 0.6, beta_n = 0.5, target = 9),
    list(mu_n = 10, k_n = 0.5, alpha_n = 1.2, beta_n = 0.05, target = 10)
  )

  set.seed(23)
  for (case in cases) {
    sigma <- sqrt(1 / stats::rgamma(draws, case$alpha_n, rate = case$beta_n))
    mu <- stats::rnorm(draws, case$mu_n, sigma / sqrt(case$k_n))

    for (metric in c("Cpm", "Cpc")) {
      moments <- qc:::.conjugate_smooth_metric_moments(
        case$mu_n, case$k_n, case$alpha_n, case$beta_n,
        6, 14, case$target, metric
      )
      samples <- qc:::compute_metric_value(mu, sigma, 6, 14, case$target, metric)
      mc_error <- stats::sd(samples) / sqrt(draws)

      testthat::expect_lt(abs(moments$mean - mean(samples)), 5 * mc_error)
      testthat::expect_lt(abs(moments$sd / stats::sd(samples) - 1), 2e-3)
    }
  }
})

testthat::test_that("contour densities resolve heavy tails and second modes", {
  mixture <- function(mu) {
    log(0.5 * stats::dnorm(mu, 9.5, 0.02) + 0.5 * stats::dnorm(mu, 10.5, 0.02))
  }
  set.seed(9)
  x <- stats::rnorm(4, 10, 0.8)
  cases <- list(
    heavy_tail = list(
      data = numeric(0),
      prior = qc:::create_prior_conjugate(10, 0.1, 0.05, 0.05)
    ),
    second_mode = list(
      data = x,
      prior = qc:::create_prior_semi_sigma(2, 1, mixture)
    )
  )

  for (case in cases) {
    state <- qc:::.as_qc_suff_stats_state(data = case$data)
    request <- qc:::.new_qc_integration_request(
      case$data, 6, 14, case$prior, "Cpl", 10, state, 3
    )
    entry <- qc:::analyze_capability_integration(request = request)
    S <- qc:::.integration_make_solver(request = request)

    testthat::expect_equal(S(unname(entry$stats["Median"])), 0.5, tolerance = 2e-3)
    testthat::expect_equal(S(unname(entry$stats["Q97.5"])), 0.025, tolerance = 0.05)
  }
})

testthat::test_that("contour quadrature failures are reported as missing density", {
  post <- list(mu_n = 10, k_n = 30, alpha_n = 15, beta_n = 15)
  log_joint <- conjugate_log_joint(post)
  pdf_fn <- qc:::.linear_metric_density(
    function(mu, sigma) {
      spike <- abs(log(sigma) - 0.4) < 1e-3
      log_joint(mu, sigma) + ifelse(spike, 1e4, 0) + 40 * sin(sigma * 1e4)
    },
    LSL = 6, USL = 14, metric = "Cpu", sigma_level = 3, sigma_ref = 1
  )

  testthat::expect_true(any(is.na(pdf_fn(seq(1, 1.5, length.out = 5L)))))
})

expect_entry_matches_draws <- function(entry, draws, tolerance = 0.02) {
  expected <- stats::quantile(draws, c(0.025, 0.5, 0.975), names = FALSE)
  width <- expected[3] - expected[1]
  actual <- unname(entry$stats[c("Q2.5", "Median", "Q97.5")])
  testthat::expect_lt(max(abs(actual - expected)) / width, tolerance)
  testthat::expect_lt(abs(unname(entry$stats["Mean"]) - mean(draws)) / width, tolerance)
  testthat::expect_lt(abs(unname(entry$stats["SD"]) - stats::sd(draws)) / width, tolerance)
}

analyze_from_state <- function(state, prior, metric, target = 10) {
  qc:::analyze_capability_integration(
    request = qc:::.new_qc_integration_request(
      numeric(0), 6, 14, prior, metric, target, state, 3
    )
  )
}

conjugate_draws <- function(state, prior, metric, target = 10, n_draws = 5e5L) {
  post <- qc:::.compute_validated_conjugate_posterior(
    prior, numeric(0), state
  )$post
  sigma <- sqrt(1 / stats::rgamma(n_draws, post$alpha_n, rate = post$beta_n))
  mu <- stats::rnorm(n_draws, post$mu_n, sigma / sqrt(post$k_n))
  qc:::compute_metric_value(mu, sigma, 6, 14, target, metric)
}

testthat::test_that("conjugate Cpm and Cpc are right far from the target", {
  prior <- qc:::create_prior_conjugate(10, 2, 3, 3)
  cases <- list(
    out_of_spec = c(x_bar = 15, sd = 1.5, n = 40),
    off_target = c(x_bar = 12, sd = 0.5, n = 1e4),
    far_narrow = c(x_bar = 15, sd = 0.1, n = 1e5),
    centred = c(x_bar = 10, sd = 1, n = 1e6),
    tiny_sd = c(x_bar = 10, sd = 0.001, n = 30)
  )

  set.seed(31)
  for (case in cases) {
    state <- qc:::.new_qc_suff_stats_state(
      n = case[["n"]], x_bar = case[["x_bar"]],
      sse = (case[["n"]] - 1) * case[["sd"]]^2
    )
    for (metric in c("Cpm", "Cpc")) {
      entry <- analyze_from_state(state, prior, metric)
      testthat::expect_equal(entry$area, 1, tolerance = 1e-3)
      expect_entry_matches_draws(entry, conjugate_draws(state, prior, metric))
    }
  }
})

testthat::test_that("radial densities match the survival function", {
  set.seed(32)
  x <- stats::rnorm(40, 12.5, 1)
  state <- qc:::.as_qc_suff_stats_state(data = x)
  priors <- list(
    conjugate = qc:::create_prior_conjugate(10, 2, 3, 3),
    semi_mu = qc:::.bayestools_to_integration_prior(
      "Jeffreys_mu", BayesTools::prior("gamma", list(2, 1))
    )$prior
  )

  for (prior in priors) {
    for (metric in c("Cpm", "Cpc")) {
      request <- qc:::.new_qc_integration_request(x, 6, 14, prior, metric, 10, state, 3)
      S <- qc:::.integration_make_solver(request = request)
      pdf_fn <- qc:::.integration_make_density_solver(request = request)
      moments <- qc:::.integration_compute_metric_moments(request = request)

      edges <- moments$mean + seq(-3, 3, length.out = 61L) * moments$sd
      survival <- vapply(edges, S, numeric(1))
      expected <- -diff(survival) / diff(edges)
      mids <- (edges[-1L] + edges[-61L]) / 2

      testthat::expect_equal(pdf_fn(mids), expected, tolerance = 2e-2)
      testthat::expect_equal(
        pdf_fn(mids[30:31]), pdf_fn(mids)[30:31]
      )
    }
  }
})

testthat::test_that("generic Cpm is right for narrow and off-target posteriors", {
  testthat::skip_if_not_installed("BayesTools")
  prior <- qc:::.bayestools_to_integration_prior(
    BayesTools::prior("normal", list(10, 5)),
    BayesTools::prior("lognormal", list(0, 1))
  )$prior

  for (case in list(c(x_bar = 5, sd = 2, n = 1e3), c(x_bar = 12, sd = 0.7, n = 1e5))) {
    base_state <- qc:::.new_qc_suff_stats_state(
      n = case[["n"]], x_bar = case[["x_bar"]],
      sse = (case[["n"]] - 1) * case[["sd"]]^2
    )
    state <- qc:::precompute_generic_state(numeric(0), prior, cached_state = base_state)
    entry <- analyze_from_state(state, prior, "Cpm")

    # The likelihood dominates, so the sampling distribution is the reference.
    plug_in <- qc:::compute_metric_value(
      case[["x_bar"]], case[["sd"]], 6, 14, 10, "Cpm"
    )
    testthat::expect_equal(entry$area, 1, tolerance = 1e-3)
    testthat::expect_equal(unname(entry$stats["Median"]), plug_in, tolerance = 0.02)
  }
})

testthat::test_that("semi-conjugate mu posteriors are normalized at large n", {
  testthat::skip_if_not_installed("BayesTools")
  prior <- qc:::.bayestools_to_integration_prior(
    "Jeffreys_mu", BayesTools::prior("gamma", list(2, 1))
  )$prior

  for (case in list(c(x_bar = 12, sd = 0.5, n = 1e4), c(x_bar = 15, sd = 0.1, n = 1e5))) {
    state <- qc:::.new_qc_suff_stats_state(
      n = case[["n"]], x_bar = case[["x_bar"]],
      sse = (case[["n"]] - 1) * case[["sd"]]^2
    )
    # The sigma prior is negligible here, so Jeffreys is an exact-enough reference.
    reference <- qc:::create_prior_conjugate(0, 0, -0.5, 0)

    for (metric in c("Cp", "Cpu", "Cpk", "Cpm", "Cpc")) {
      entry <- analyze_from_state(state, prior, metric)
      expected <- analyze_from_state(state, reference, metric)
      width <- unname(expected$stats["Q97.5"] - expected$stats["Q2.5"])

      testthat::expect_lt(
        max(abs(entry$stats[c("Mean", "Q2.5", "Median", "Q97.5")] -
                  expected$stats[c("Mean", "Q2.5", "Median", "Q97.5")])) / width,
        0.05
      )
    }
  }
})

testthat::test_that("scalar and batched conjugate survival functions agree", {
  prior <- qc:::create_prior_conjugate(0, 0, -0.5, 0)
  state <- qc:::.new_qc_suff_stats_state(n = 1000L, x_bar = 10, sse = 999e-6)

  for (metric in c("Cpm", "Cpc")) {
    S <- qc:::make_solver(numeric(0), 6, 14, prior, metric = metric,
                          target = 10, cached_state = state)
    testthat::expect_equal(S(1000), 1, tolerance = 1e-6)
    testthat::expect_equal(S(1300), S(c(1300, 1400))[1])

    entry <- analyze_from_state(state, prior, metric)
    testthat::expect_equal(entry$area, 1, tolerance = 1e-3)
    expect_entry_matches_draws(entry, conjugate_draws(state, prior, metric))
  }
})

testthat::test_that("very diffuse priors keep their mass near zero on the grid", {
  state <- qc:::.new_qc_suff_stats_state(n = 0L, x_bar = 0, sse = 0)
  priors <- list(
    centred = qc:::create_prior_conjugate(10, 0.1, 0.05, 0.05),
    tiny_scale = qc:::create_prior_conjugate(10, 0.1, 0.05, 5e-5),
    # Cpm and Cpc pile up in a spike against their bound.
    off_target = qc:::create_prior_conjugate(25, 0.1, 0.05, 1e-5)
  )

  set.seed(33)
  for (prior in priors) {
    for (metric in c("Cp", "Cpm", "Cpc", "Cpk")) {
      entry <- analyze_from_state(state, prior, metric)
      draws <- conjugate_draws(state, prior, metric, n_draws = 1e6L)
      expected <- stats::quantile(draws, c(0.25, 0.5, 0.75), names = FALSE)

      testthat::expect_lt(
        abs(unname(entry$stats["Median"]) - expected[2]) / (expected[3] - expected[1]),
        0.02
      )
      testthat::expect_lt(
        abs(unname(entry$stats["Q97.5"]) / stats::quantile(draws, 0.975) - 1),
        0.05
      )
    }
  }
})

testthat::test_that("adaptive grids resolve the shoulders of a zoomed narrow peak", {
  pdf_fn <- function(x) stats::dnorm(x, 4 / 3, 1e-3)

  result <- qc:::.adaptive_density_grid(pdf_fn, 0, 2, n_grid = 512L)

  testthat::expect_equal(result$area, 1, tolerance = 1e-5)
})

testthat::test_that("densities that lose mass fall back to the survival function", {
  near_modes <- function(mu) {
    log(0.5 * stats::dnorm(mu, 9.9, 0.01) + 0.5 * stats::dnorm(mu, 10.1, 0.01))
  }
  set.seed(9)
  x <- stats::rnorm(4, 10, 0.8)
  prior <- qc:::create_prior_semi_sigma(2, 1, near_modes)
  state <- qc:::.as_qc_suff_stats_state(data = x)

  for (metric in c("Cpu", "Cpl")) {
    request <- qc:::.new_qc_integration_request(x, 6, 14, prior, metric, 10, state, 3)
    entry <- qc:::analyze_capability_integration(request = request)
    S <- qc:::.integration_make_solver(request = request)

    testthat::expect_equal(S(unname(entry$stats["Median"])), 0.5, tolerance = 5e-3)
    testthat::expect_equal(S(unname(entry$stats["Q97.5"])), 0.025, tolerance = 0.1)
  }
})

semi_sigma_draws <- function(prior, state, n_draws = 5e5L) {
  frame <- qc:::.prepare_semi_sigma_state(numeric(0), prior, state)
  grid <- seq(frame$mu_frame$lower, frame$mu_frame$upper, length.out = 20001L)
  log_weight <- frame$log_weight(grid)
  mu <- sample(grid, n_draws, TRUE, exp(log_weight - max(log_weight))) +
    stats::runif(n_draws, -0.5, 0.5) * (grid[2] - grid[1])
  beta_n <- frame$beta0 + (frame$sse + frame$n * (mu - frame$x_bar)^2) / 2
  list(
    mu = mu,
    sigma = sqrt(1 / stats::rgamma(n_draws, frame$alpha_n_eff, rate = beta_n))
  )
}

testthat::test_that("semi-conjugate sigma Cp, Cpm and Cpc match draws", {
  testthat::skip_if_not_installed("BayesTools")

  priors <- list(
    student = qc:::.bayestools_to_integration_prior(
      BayesTools::prior("t", list(10, 5, 3)), "Jeffreys_sigma"
    )$prior,
    truncated = qc:::.bayestools_to_integration_prior(
      BayesTools::prior("normal", list(10, 2), list(9.9, Inf)), "Jeffreys_sigma"
    )$prior
  )
  cases <- list(
    small = c(x_bar = 10.2, sd = 1, n = 5),
    off_target = c(x_bar = 12, sd = 0.7, n = 100),
    out_of_spec = c(x_bar = 15, sd = 1.5, n = 40),
    narrow = c(x_bar = 10.3, sd = 0.02, n = 1e4),
    large = c(x_bar = 10.2, sd = 1, n = 1e5)
  )

  set.seed(41)
  for (prior in priors) {
    for (case in cases) {
      state <- qc:::.new_qc_suff_stats_state(
        n = case[["n"]], x_bar = case[["x_bar"]],
        sse = (case[["n"]] - 1) * case[["sd"]]^2
      )
      draws <- semi_sigma_draws(prior, state)
      for (metric in c("Cp", "Cpm", "Cpc")) {
        expect_entry_matches_draws(
          analyze_from_state(state, prior, metric),
          qc:::compute_metric_value(draws$mu, draws$sigma, 6, 14, 10, metric)
        )
      }
    }
  }
})

testthat::test_that("semi-conjugate sigma densities are evaluated in batches", {
  calls <- 0L
  prior <- qc:::create_prior_semi_sigma(2, 1, function(mu) {
    calls <<- calls + 1L
    stats::dt((mu - 10) / 5, 3, log = TRUE)
  })
  set.seed(42)
  x <- stats::rnorm(100, 10, 1)
  points <- seq(1, 2, length.out = 200L)

  for (metric in c("Cp", "Cpm", "Cpc")) {
    pdf_fn <- qc:::make_density_solver(x, 6, 14, prior, metric = metric, target = 10)
    calls <- 0L
    vals <- pdf_fn(points)
    testthat::expect_true(all(is.finite(vals)), info = metric)
    testthat::expect_lt(calls, length(points))
    testthat::expect_equal(pdf_fn(points[100]), vals[100], tolerance = 1e-6, info = metric)
  }
})

testthat::test_that("generic Cpc matches draws and the survival function", {
  testthat::skip_if_not_installed("BayesTools")

  prior <- qc:::.bayestools_to_integration_prior(
    BayesTools::prior("normal", list(10, 5)),
    BayesTools::prior("lognormal", list(0, 1))
  )$prior
  cases <- list(
    centred = c(x_bar = 10.2, sd = 1, n = 30),
    off_target = c(x_bar = 12, sd = 0.7, n = 100),
    narrow = c(x_bar = 10.3, sd = 0.02, n = 1e4)
  )

  log_prior <- function(mu, sigma) {
    stats::dnorm(mu, 10, 5, log = TRUE) + stats::dlnorm(sigma, 0, 1, log = TRUE)
  }

  set.seed(43)
  for (case in cases) {
    n <- case[["n"]]
    x_bar <- case[["x_bar"]]
    s <- case[["sd"]]
    base_state <- qc:::.new_qc_suff_stats_state(
      n = n, x_bar = x_bar, sse = (n - 1) * s^2
    )
    state <- qc:::precompute_generic_state(numeric(0), prior, cached_state = base_state)
    entry <- analyze_from_state(state, prior, "Cpc")

    # Sample the posterior from a fine (mu, log(sigma)) grid.
    mu_grid <- x_bar + seq(-8, 8, length.out = 601L) * s / sqrt(n)
    t_grid <- log(s) + seq(-8, 8, length.out = 401L) / sqrt(2 * n)
    grid <- expand.grid(mu = mu_grid, t = t_grid)
    log_post <- log_prior(grid$mu, exp(grid$t)) + (1 - n) * grid$t -
      ((n - 1) * s^2 + n * (grid$mu - x_bar)^2) / (2 * exp(2 * grid$t))
    pick <- sample(nrow(grid), 5e5L, TRUE, exp(log_post - max(log_post)))
    mu <- grid$mu[pick] + stats::runif(5e5L, -0.5, 0.5) * (mu_grid[2] - mu_grid[1])
    sigma <- exp(grid$t[pick] + stats::runif(5e5L, -0.5, 0.5) * (t_grid[2] - t_grid[1]))
    expect_entry_matches_draws(
      entry, qc:::compute_metric_value(mu, sigma, 6, 14, 10, "Cpc")
    )

    solver <- qc:::make_solver(
      numeric(0), 6, 14, prior, metric = "Cpc", target = 10, cached_state = state
    )
    median_c <- unname(entry$stats["Median"])
    testthat::expect_equal(solver(median_c), 0.5, tolerance = 0.01)
  }
})

testthat::test_that("long-tailed posteriors get a grid that resolves their peak", {
  testthat::skip_if_not_installed("BayesTools")

  # The prior holds sigma near 1 while two nearly equal points pull it to 0.
  prior <- qc:::.bayestools_to_integration_prior(
    "Jeffreys_mu",
    BayesTools::prior("gamma", list(2, 1))
  )$prior
  state <- qc:::.new_qc_suff_stats_state(n = 2, x_bar = 10, sse = 1e-6)

  set.seed(44)
  n_draws <- 6e6L
  t <- stats::runif(n_draws, log(1e-3) - 6, log(1e-3) + 9)
  sigma <- exp(t)
  mu <- stats::rnorm(n_draws, 10, sigma / sqrt(2))
  log_weight <- stats::dgamma(sigma, 2, 1, log = TRUE) - 1e-6 / (2 * sigma^2)
  pick <- sample(n_draws, 2e6L, TRUE, exp(log_weight - max(log_weight)))

  for (metric in c("Cpu", "Cpc")) {
    entry <- analyze_from_state(state, prior, metric)
    draws <- qc:::compute_metric_value(mu[pick], sigma[pick], 6, 14, 10, metric)
    expected <- stats::quantile(draws, c(0.025, 0.5, 0.975), names = FALSE)
    actual <- unname(entry$stats[c("Q2.5", "Median", "Q97.5")])
    testthat::expect_equal(actual, expected, tolerance = 0.02, info = metric)
  }
})

testthat::test_that("a step in the density does not trigger the survival fallback", {
  prior <- qc:::create_prior_conjugate(0, 0, -0.5, 0)
  state <- qc:::.new_qc_suff_stats_state(n = 2, x_bar = 10, sse = 1e-6)
  entry <- analyze_from_state(state, prior, "Cpk")

  # Finite differences of the survival function ripple on this flat density.
  rising <- diff(entry$grid$density) > 0
  testthat::expect_lte(sum(diff(rising) != 0), 2L)
  draws <- conjugate_draws(state, prior, "Cpk")
  expect_entry_matches_draws(entry, draws)
})

testthat::test_that("conjugate Cpm and Cpc moments are finite for very diffuse priors", {
  # With mu0 = target and k0 = 1, Cpm = (4 / 3) sqrt(G / 10) / sqrt(1 + Z^2)
  # for independent G ~ Gamma(alpha, 1) and Z ~ N(0, 1); Cpc replaces
  # sqrt(1 + Z^2) by sqrt(pi / 2) g(|Z|).
  spread <- list(
    Cpm = function(z) sqrt(1 + z^2),
    Cpc = function(z) sqrt(pi / 2) * qc:::.cpc_g(abs(z))
  )
  for (alpha in c(0.01, 0.001)) {
    prior <- qc:::create_prior_conjugate(10, 1, alpha, 10)
    for (metric in names(spread)) {
      z_moment <- function(power) {
        stats::integrate(function(z) 2 * stats::dnorm(z) / spread[[metric]](z)^power,
                         0, Inf, rel.tol = 1e-10)$value
      }
      mean_ref <- 4 / 3 / sqrt(10) *
        exp(lgamma(alpha + 0.5) - lgamma(alpha)) * z_moment(1)
      sd_ref <- sqrt((4 / 3)^2 / 10 * alpha * z_moment(2) - mean_ref^2)

      moments <- qc:::compute_metric_moments(
        numeric(0), 6, 14, prior, metric, target = 10
      )
      info <- paste(alpha, metric)
      testthat::expect_equal(moments$mean, mean_ref, tolerance = 1e-6, info = info)
      testthat::expect_equal(moments$sd, sd_ref, tolerance = 1e-6, info = info)
    }
  }
})

testthat::test_that("radial densities keep a narrow mean support away from the target", {
  # mu ~ normal(12.2, 1) truncated to [12, 12.1] and sigma ~ inverse gamma
  # (3, scale 2), so |mu - target| >= 2 and Cpm <= 4 / (3 * 2).
  log_joint_on <- function(lower, upper) {
    mu_mass <- stats::pnorm(upper, 12.2, 1) - stats::pnorm(lower, 12.2, 1)
    function(mu, sigma) {
      vals <- stats::dnorm(mu, 12.2, 1, log = TRUE) - log(mu_mass) +
        stats::dgamma(1 / sigma, 3, rate = 2, log = TRUE) - 2 * log(sigma)
      vals[mu < lower | mu > upper] <- -Inf
      vals
    }
  }
  radial_pdf <- function(metric, lower, upper) {
    qc:::.radial_metric_density(
      log_joint_on(lower, upper), 6, 14, target = 10, metric = metric,
      sigma_level = 3, sigma_ref = 1, mu_ref = (lower + upper) / 2,
      mu_lower = lower, mu_upper = upper
    )
  }

  # On a Cpm contour sigma = sqrt(R^2 - delta^2) with R = 4 / (3 c), and
  # |d sigma / d c| = R^2 / (c sigma).
  log_joint <- log_joint_on(12, 12.1)
  reference <- function(c) {
    stats::integrate(function(mu) {
      R <- 4 / (3 * c)
      sigma <- sqrt(pmax(R^2 - (mu - 10)^2, 0))
      dens <- exp(log_joint(mu, sigma)) * R^2 / (c * sigma)
      dens[!is.finite(dens)] <- 0
      dens
    }, 12, 12.1, rel.tol = 1e-10)$value
  }
  pdf_fn <- radial_pdf("Cpm", 12, 12.1)
  c_values <- c(0.35, 0.45, 0.5, 0.55, 0.6, 0.64)
  testthat::expect_equal(pdf_fn(c_values), vapply(c_values, reference, numeric(1)),
                         tolerance = 1e-6)
  testthat::expect_equal(pdf_fn(0.7), 0)

  # A support far thinner than the Cpc lookup table's resolution on log(z).
  for (metric in c("Cpm", "Cpc")) {
    pdf_fn <- radial_pdf(metric, 12, 12 + 1e-8)
    mass <- stats::integrate(pdf_fn, 0, 0.7, subdivisions = 400L)$value
    testthat::expect_equal(mass, 1, tolerance = 1e-4, info = metric)
  }
})

testthat::test_that("truncated mean priors keep the radial density on its support", {
  testthat::skip_if_not_installed("BayesTools")

  set.seed(5)
  x <- stats::rnorm(20, 12.05, 1.4)
  sigma_priors <- list(
    semi_sigma = "Jeffreys_sigma",
    # The generic posterior box is far wider than the truncation here.
    generic = BayesTools::prior("lognormal", list(0, 1))
  )
  # |mu - target| >= 2 bounds Cpm by 0.667 and Cpc by about 0.532.
  c_values <- seq(0.0025, 0.7, by = 0.005)

  for (case in names(sigma_priors)) {
    prior <- qc:::.bayestools_to_integration_prior(
      BayesTools::prior("normal", list(12.2, 1), list(12, 12.1)),
      sigma_priors[[case]]
    )$prior
    for (metric in c("Cpm", "Cpc")) {
      pdf_fn <- qc:::make_density_solver(x, 6, 14, prior, metric = metric, target = 10)
      vals <- pdf_fn(c_values)
      on_support <- range(which(vals > 0))
      info <- paste(case, metric)
      # The contours sweep the support continuously, so it has no holes.
      testthat::expect_true(all(vals[on_support[1]:on_support[2]] > 0), info = info)
      # The generic box normalization is good to about 0.3% under a hard
      # truncation of mu.
      testthat::expect_equal(sum(vals) * 0.005, 1, tolerance = 5e-3, info = info)
    }
  }
})
