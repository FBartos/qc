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
