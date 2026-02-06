
testthat::test_that("integration method works when data is missing from scope", {

  # Create a local environment to fit the model and then discard the data
  fit_model <- function() {
    set.seed(123)
    # Create dataset locally
    local_data <- rnorm(30, mean = 10, sd = 1)

    # Fit model with integration method
    # This should now cache the sufficient statistics
    fit <- bpc(
      local_data,
      LSL = 7,
      target = 10,
      USL = 13,
      method = "integration"
    )
    return(fit)
  }

  fit <- fit_model()

  # Ensure local_data is definitely not in the global environment or current scope
  expect_false(exists("local_data", envir = environment()))

  # This summary call would previously fail with "object 'local_data' not found"
  # because it tried to eval(call$x)
  # It should now succeed using the cached sufficient statistics
  summ <- summary(fit)

  # Check results are valid
  expect_s3_class(summ, "bpc_summary")
  expect_true(!is.null(summ$summary))

  # Check that stats are reasonable (mean approx 10)
  cpk_idx <- which(summ$summary$metric == "Cpk")
  expect_equal(summ$summary$mean[cpk_idx], 1.0, tolerance = 0.5)
})

testthat::test_that("integration method works with BayesTools priors when data is missing from scope", {

  # Create a local environment to fit the model and then discard the data
  fit_model_priors <- function() {
    set.seed(123)
    local_data <- rnorm(30, mean = 10, sd = 1)

    # Use BayesTools priors (PriorGeneric path)
    fit <- bpc(
      local_data,
      LSL = 7,
      target = 10,
      USL = 13,
      method = "integration",
      prior_mu = BayesTools::prior("normal", list(10, 5)),
      prior_sigma = BayesTools::prior("gamma", list(1, 1))
    )
    return(fit)
  }

  fit <- fit_model_priors()

  # Ensure local_data is definitely not in the global environment or current scope
  expect_false(exists("local_data", envir = environment()))

  # This covers the PriorGeneric path which also needs cached_state
  summ <- summary(fit)

  expect_s3_class(summ, "bpc_summary")
  expect_true(!is.null(summ$summary))
})
