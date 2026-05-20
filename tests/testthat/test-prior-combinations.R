# Selected breadth coverage across prior families. This is a slow, opt-in
# parity sweep rather than the primary backend regression guard.
testthat::test_that("selected prior-family combinations keep matched-prior backend parity", {
  skip_if_not_slow_tests()
  skip_if_not_installed("BayesTools")

  # Generate test data
  set.seed(123)
  x <- rnorm(30, mean = 10, sd = 2)

  scenarios <- list(
    list(
      name = "with-data jeffreys/jeffreys",
      x = x,
      LSL = 2, target = 10, USL = 18,
      prior = qc::prior_independent(mu = "Jeffreys_mu", sigma = "Jeffreys_sigma"),
      sample_priors = FALSE,
      metrics = qc:::.qc_metric_names(),
      columns = c("mean", "median", "lower", "upper"),
      summary_abs_tol = 0.03,
      summary_rel_tol = 0.02,
      interval_breaks = c(1.00, 1.33, 1.50, 2.00),
      interval_abs_tol = 0.02,
      expect_failure = FALSE
    ),
    list(
      name = "with-data normal/invgamma",
      x = x,
      LSL = 2, target = 10, USL = 18,
      prior = qc::prior_independent(
        mu = BayesTools::prior("normal", list(10, 5)),
        sigma = BayesTools::prior("invgamma", list(2, 1))
      ),
      sample_priors = FALSE,
      metrics = c("Cp", "Cpk", "Cpm"),
      columns = c("mean", "median", "lower", "upper"),
      summary_abs_tol = 0.05,
      summary_rel_tol = 0.05,
      interval_breaks = c(1.00, 1.33, 1.50, 2.00),
      interval_abs_tol = 0.03,
      expect_failure = FALSE
    ),
    list(
      name = "with-data truncated-normal/gamma",
      x = x,
      LSL = 2, target = 10, USL = 18,
      prior = qc::prior_independent(
        mu = BayesTools::prior("normal", list(10, 5), list(5, 15)),
        sigma = BayesTools::prior("gamma", list(2, 1))
      ),
      sample_priors = FALSE,
      metrics = c("Cp", "Cpk", "Cpm"),
      columns = c("mean", "median", "lower", "upper"),
      summary_abs_tol = 0.08,
      summary_rel_tol = 0.05,
      interval_breaks = c(1.00, 1.33, 1.50, 2.00),
      interval_abs_tol = 0.03,
      expect_failure = FALSE
    ),
    list(
      name = "with-data student-t/lognormal",
      x = x,
      LSL = 2, target = 10, USL = 18,
      prior = qc::prior_independent(
        mu = BayesTools::prior("t", list(10, 3, 5)),
        sigma = BayesTools::prior("lognormal", list(0, 1))
      ),
      sample_priors = FALSE,
      metrics = c("Cp", "Cpk", "Cpm"),
      columns = c("mean", "median", "lower", "upper"),
      summary_abs_tol = 0.08,
      summary_rel_tol = 0.05,
      interval_breaks = c(1.00, 1.33, 1.50, 2.00),
      interval_abs_tol = 0.03,
      expect_failure = FALSE
    ),
    list(
      name = "with-data uniform/exponential",
      x = x,
      LSL = 2, target = 10, USL = 18,
      prior = qc::prior_independent(
        mu = BayesTools::prior("uniform", list(0, 20)),
        sigma = BayesTools::prior("exp", list(1))
      ),
      sample_priors = FALSE,
      metrics = c("Cp", "Cpk", "Cpm"),
      columns = c("mean", "median", "lower", "upper"),
      summary_abs_tol = 0.08,
      summary_rel_tol = 0.05,
      interval_breaks = c(1.00, 1.33, 1.50, 2.00),
      interval_abs_tol = 0.03,
      expect_failure = FALSE
    ),
    list(
      name = "prior-only proper normal/gamma",
      x = numeric(2L),
      LSL = -1, target = 0, USL = 1,
      prior = qc::prior_independent(
        mu = BayesTools::prior("normal", list(0, 1)),
        sigma = BayesTools::prior("gamma", list(3, 1))
      ),
      sample_priors = TRUE,
      metrics = c("Cp", "Cpk", "Cpm"),
      columns = c("median", "lower", "upper"),
      summary_abs_tol = 0.08,
      summary_rel_tol = 0.05,
      interval_breaks = c(1.00, 1.33, 1.50, 2.00),
      interval_abs_tol = 0.03,
      expect_failure = FALSE
    ),
    list(
      name = "prior-only improper jeffreys/gamma",
      x = numeric(2L),
      LSL = -1, target = 0, USL = 1,
      prior = qc::prior_independent(
        mu = "Jeffreys_mu",
        sigma = BayesTools::prior("gamma", list(3, 1))
      ),
      sample_priors = TRUE,
      expect_failure = TRUE
    ),
    list(
      name = "prior-only improper normal/jeffreys",
      x = numeric(2L),
      LSL = -1, target = 0, USL = 1,
      prior = qc::prior_independent(
        mu = BayesTools::prior("normal", list(0, 1)),
        sigma = "Jeffreys_sigma"
      ),
      sample_priors = TRUE,
      expect_failure = TRUE
    )
  )

  for (scenario in scenarios) {
    fit_args <- list(
      x = scenario$x,
      LSL = scenario$LSL,
      target = scenario$target,
      USL = scenario$USL,
      prior = scenario$prior,
      sample_priors = scenario$sample_priors
    )

    fit_int <- tryCatch(
      do.call(qc::bpc, c(fit_args, list(method = "integration"))),
      error = function(e) e
    )
    fit_mcmc <- tryCatch(
      do.call(
        qc::bpc,
        c(
          fit_args,
          list(
            method = "mcmc",
            chains = 2,
            iter = 2500,
            warmup = 500,
            cores = 1,
            silent = TRUE,
            seed = 123
          )
        )
      ),
      error = function(e) e
    )

    if (scenario$expect_failure) {
      expect_s3_class(fit_int, "error", info = scenario$name)
      expect_s3_class(fit_mcmc, "error", info = scenario$name)
      next
    }

    if (inherits(fit_int, "error") || inherits(fit_mcmc, "error")) {
      testthat::fail(sprintf(
        "%s unexpected failure: integration=%s; mcmc=%s",
        scenario$name,
        if (inherits(fit_int, "error")) fit_int$message else "<ok>",
        if (inherits(fit_mcmc, "error")) fit_mcmc$message else "<ok>"
      ))
    }

    pair <- list(integration = fit_int, mcmc = fit_mcmc)
    summaries <- expect_fit_summaries_close(
      pair,
      metrics = scenario$metrics,
      columns = scenario$columns,
      abs_tol = scenario$summary_abs_tol,
      rel_tol = scenario$summary_rel_tol,
      interval_probability = scenario$interval_breaks
    )

    expect_interval_summary_close(
      summary_int = summaries$integration,
      summary_mcmc = summaries$mcmc,
      metrics = scenario$metrics,
      abs_tol = scenario$interval_abs_tol
    )
  }
})
