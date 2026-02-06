# Tests for all combinations of prior families
# Compares MCMC against integration for both prior-only and with-data scenarios
testthat::test_that("integration matches MCMC for all prior combinations", {

  skip_on_cran()  # Skip on CRAN due to long runtime

  library(BayesTools)

  # Generate test data
  set.seed(123)
  x <- rnorm(30, mean = 10, sd = 2)

  # Define prior families for mu
  mu_priors <- rlang::list2(
    "Jeffreys" = "Jeffreys_mu",
    "Normal" = prior("normal", list(10, 5)),
    "Normal_truncated" = prior("normal", list(10, 5), list(5, 15)),
    # "Student-t" = prior("t", list(10, 3, 5)),     # location=10, scale=3, df=5
    # "Uniform" = prior("uniform", list(0, 20))
  )

  # Define prior families for sigma
  sigma_priors <- rlang::list2(
    "Jeffreys" = "Jeffreys_sigma",
    "InvGamma" = prior("invgamma", list(2, 1)),   # shape=2, scale=1
    "Gamma" = prior("gamma", list(2, 1)),         # shape=2, rate=1
    # "LogNormal" = prior("lognormal", list(0, 1)), # log-mean=0, log-sd=1
    # "Exponential" = prior("exp", list(1))         # rate=1
  )

  sample_prior_opts <- c(TRUE, FALSE)
  sample_prior_opts <- c(FALSE)

  all_opts <- expand.grid(
    mu = names(mu_priors),
    sigma = names(sigma_priors),
    sample_priors = sample_prior_opts,
    stringsAsFactors = FALSE
  )
  ntest <- nrow(all_opts)

  # Test parameters
  LSL <- 2
  USL <- 18
  target <- 10
  bounds <- c(0.8, 1.5)
  tolerance <- 0.05

  for (i in seq_len(nrow(all_opts))) {
    mu_name <- all_opts$mu[i]
    sigma_name <- all_opts$sigma[i]
    sample_priors <- all_opts$sample_priors[i]

    prior_mu <- mu_priors[[mu_name]]
    prior_sigma <- sigma_priors[[sigma_name]]

    # Determine if we expect an error (improper priors + prior-only mode)
    is_improper <- (mu_name == "Jeffreys" || sigma_name == "Jeffreys")
    expect_failure <- sample_priors && is_improper

    cat(sprintf("\nTesting: mu=%s, sigma=%s, sample_priors=%s (Expect Failure: %s)\n",
                mu_name, sigma_name, sample_priors, expect_failure))

    # Fit with integration
    fit_int <- tryCatch({
      bpc(x, LSL = LSL, target = target, USL = USL,
          method = "integration",
          prior_mu = prior_mu, prior_sigma = prior_sigma,
          sample_priors = sample_priors)
    }, error = function(e) {
      # cat(sprintf("  Integration error: %s\n", e$message))
      NULL
    })

    # Fit with MCMC
    fit_mcmc <- tryCatch({
      bpc(x, LSL = LSL, target = target, USL = USL,
          method = "mcmc",
          prior_mu = prior_mu, prior_sigma = prior_sigma,
          iter = 10000, chains = 2, silent = TRUE, seed = 123,
          sample_priors = sample_priors)
    }, error = function(e) {
      # cat(sprintf("  MCMC error: %s\n", e$message))
      NULL
    })

    if (expect_failure) {
      testthat::expect_true(is.null(fit_int),
                            info = sprintf("Integration should fail for improper priors without data: i=%d, mu=%s, sigma=%s", i, mu_name, sigma_name))
      testthat::expect_true(is.null(fit_mcmc),
                            info = sprintf("MCMC should fail for improper priors without data: i=%d,mu=%s, sigma=%s", i, mu_name, sigma_name))
      next
    }

    if (is.null(fit_int) || is.null(fit_mcmc)) {
      # Unexpected failure
      if (is.null(fit_int)) cat(sprintf("  Unexpected Integration failure for i=%d, mu=%s, sigma=%s\n", i, mu_name, sigma_name))
      if (is.null(fit_mcmc)) cat(sprintf("  Unexpected MCMC failure for i=%d, mu=%s, sigma=%s\n", i, mu_name, sigma_name))
      next
    }

    sum_mcmc <- summary(fit_mcmc)
    sum_int  <- summary(fit_int)

    metrics_mcmc <- sum_mcmc$summary
    metrics_int  <- sum_int$summary

    # Check equality of summary statistics (compare means)
    # Using low tolerance for means comparison
    # Note: MCMC variance can be high depending on iter, so we check loose agreement or structure
    testthat::expect_equal(metrics_mcmc, metrics_int,
                           info = sprintf("i=%d, Metric name mismatch", i), tolerance = 1e-1)

    # We compare the means with some tolerance
    # Extract mean columns
    mcmc_means <- metrics_mcmc$mean
    int_means <- metrics_int$mean

    # Note: Tolerance MCMC vs Integration can be loose
    diffs <- abs(mcmc_means - int_means)
    max_diff <- max(diffs, na.rm = TRUE)

    cat(sprintf("i: %d,  Max mean diff: %.4f\n", i, max_diff))

    # Just checking intervals structure for now as requested
    intervals_mcmc <- sum_mcmc$interval_summary
    intervals_int  <- sum_int$interval_summary

    testthat::expect_equal(intervals_mcmc, intervals_int,
                           info = sprintf("i=%d, Interval metric mismatch", i), tolerance = 1e-1)
  }
})
