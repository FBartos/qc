test_that("Conjugate backends agree for all metrics under matched Jeffreys priors", {
  skip_if_not_slow_tests()

  pair <- get_backend_test_fixture("backend-slow-conjugate", function() {
    set.seed(42)
    x <- rnorm(30, mean = 50, sd = 0.5)

    fit_backend_pair(
      x = x,
      LSL = 44, target = 50, USL = 56,
      prior = "Jeffreys",
      mcmc_chains = 2,
      mcmc_iter = 4000,
      mcmc_warmup = 1000,
      seed = 42
    )
  })

  summaries <- expect_fit_summaries_close(
    pair,
    metrics = qc:::.qc_metric_names(),
    columns = c("mean", "median", "lower", "upper"),
    abs_tol = 0.03,
    rel_tol = 0.02
  )

  expect_interval_summary_close(
    summary_int = summaries$integration,
    summary_mcmc = summaries$mcmc,
    metrics = qc:::.qc_metric_names(),
    abs_tol = 0.02
  )
})

test_that("Semi-conjugate mu backends agree on summaries and interval probabilities", {
  skip_if_not_slow_tests()
  skip_if_not_installed("BayesTools")

  pair <- get_backend_test_fixture("backend-slow-semi-mu", function() {
    set.seed(123)
    x <- rnorm(30, mean = 10, sd = 2)

    fit_backend_pair(
      x = x,
      LSL = 2, target = 10, USL = 18,
      prior = qc::prior_independent(
        mu = "Jeffreys_mu",
        sigma = BayesTools::prior("gamma", list(2, 1))
      ),
      mcmc_chains = 2,
      mcmc_iter = 4000,
      mcmc_warmup = 1000,
      seed = 123
    )
  })

  summaries <- expect_fit_summaries_close(
    pair,
    metrics = qc:::.qc_metric_names(),
    columns = c("mean", "median", "lower", "upper"),
    abs_tol = 0.05,
    rel_tol = 0.05,
    interval_probability = c(0, 1, 4 / 3, 1.5, 2)
  )

  expect_interval_summary_close(
    summary_int = summaries$integration,
    summary_mcmc = summaries$mcmc,
    metrics = c("Cpu", "Cpl", "Cpk"),
    abs_tol = 0.03
  )
})

test_that("Semi-conjugate sigma backends agree on summaries and interval probabilities", {
  skip_if_not_slow_tests()
  skip_if_not_installed("BayesTools")

  pair <- get_backend_test_fixture("backend-slow-semi-sigma", function() {
    set.seed(124)
    x <- rnorm(30, mean = 10, sd = 2)

    fit_backend_pair(
      x = x,
      LSL = 2, target = 10, USL = 18,
      prior = qc::prior_independent(
        mu = BayesTools::prior("normal", list(10, 5)),
        sigma = "Jeffreys_sigma"
      ),
      mcmc_chains = 2,
      mcmc_iter = 4000,
      mcmc_warmup = 1000,
      seed = 124
    )
  })

  summaries <- expect_fit_summaries_close(
    pair,
    metrics = qc:::.qc_metric_names(),
    columns = c("mean", "median", "lower", "upper"),
    abs_tol = 0.05,
    rel_tol = 0.05,
    interval_probability = c(0, 1, 4 / 3, 1.5, 2)
  )

  expect_interval_summary_close(
    summary_int = summaries$integration,
    summary_mcmc = summaries$mcmc,
    metrics = c("Cpu", "Cpl", "Cpk"),
    abs_tol = 0.03
  )
})

test_that("Non-conjugate matched priors agree on summaries and interval probabilities", {
  skip_if_not_slow_tests()
  skip_if_not_installed("BayesTools")

  pair <- get_backend_test_fixture("backend-slow-generic", function() {
    set.seed(125)
    x <- rnorm(30, mean = 10, sd = 2)
    prior_spec <- qc::prior_independent(
      mu = BayesTools::prior("normal", list(10, 5)),
      sigma = BayesTools::prior("gamma", list(2, 1))
    )

    fit_backend_pair(
      x = x,
      LSL = 2, target = 10, USL = 18,
      prior = prior_spec,
      mcmc_chains = 2,
      mcmc_iter = 4000,
      mcmc_warmup = 1000,
      seed = 125
    )
  })

  summaries <- expect_fit_summaries_close(
    pair,
    metrics = qc:::.qc_metric_names(),
    columns = c("mean", "median", "lower", "upper"),
    abs_tol = 0.08,
    rel_tol = 0.05
  )

  expect_interval_summary_close(
    summary_int = summaries$integration,
    summary_mcmc = summaries$mcmc,
    metrics = qc:::.qc_metric_names(),
    abs_tol = 0.03
  )
})

test_that("Prior-only proper matched priors agree on finite summaries and interval probabilities", {
  skip_if_not_slow_tests()
  skip_if_not_installed("BayesTools")

  pair <- get_backend_test_fixture("backend-slow-prior-only-proper", function() {
    fit_backend_pair(
      x = numeric(2L),
      LSL = -1, target = 0, USL = 1,
      prior = qc::prior_independent(
        mu = BayesTools::prior("normal", list(0, 1)),
        sigma = BayesTools::prior("gamma", list(3, 1))
      ),
      sample_priors = TRUE,
      mcmc_chains = 2,
      mcmc_iter = 2500,
      mcmc_warmup = 500,
      seed = 126
    )
  })

  summaries <- expect_fit_summaries_close(
    pair,
    metrics = c("Cp", "Cpk", "Cpm"),
    columns = c("median", "lower", "upper"),
    abs_tol = 0.08,
    rel_tol = 0.05,
    interval_probability = c(1, 4 / 3, 1.5, 2)
  )

  expect_interval_summary_close(
    summary_int = summaries$integration,
    summary_mcmc = summaries$mcmc,
    metrics = c("Cp", "Cpk", "Cpm"),
    abs_tol = 0.03
  )
})
