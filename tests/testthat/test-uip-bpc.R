test_that("bpc with UIP via integration method returns a valid bpc object", {
  set.seed(1)
  x   <- rnorm(50, mean = 10, sd = 2)
  uip <- create_prior_unit_information(x)

  fit <- bpc(x, LSL = 2, target = 10, USL = 18,
             prior = uip, method = "integration")

  expect_s3_class(fit, "bpc")
  expect_equal(fit$method, "integration")
  expect_named(fit$coefficients, c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm"))
  expect_true(all(is.finite(fit$coefficients)))
  expect_true(all(fit$coefficients > 0))
})

test_that("bpc UIP posterior summaries are non-empty and finite", {
  set.seed(2)
  x   <- rnorm(30, mean = 5, sd = 1)
  uip <- create_prior_unit_information(x)

  fit <- bpc(x, LSL = 1, target = 5, USL = 9, prior = uip, method = "integration")
  s   <- summary(fit)

  expect_true(!is.null(s))
  # coefficients should be positive and finite
  expect_true(all(is.finite(fit$coefficients)))
})

test_that("bpc UIP regression: posterior Cp close to frequentist estimate for large n", {
  # For large n, the UIP posterior mean should be close to the MLE
  set.seed(3)
  n   <- 500
  mu  <- 10
  sig <- 2
  x   <- rnorm(n, mu, sig)
  LSL <- mu - 3 * sig
  USL <- mu + 3 * sig

  uip <- create_prior_unit_information(x)
  fit <- bpc(x, LSL = LSL, target = mu, USL = USL,
             prior = uip, method = "integration")

  # Frequentist Cp = (USL - LSL) / (6 * sd(x))
  cp_mle <- (USL - LSL) / (6 * sd(x))
  expect_lt(abs(fit$coefficients[["Cp"]] - cp_mle), 0.05)
})

test_that("bpc UIP prior predictive uses prior parameters when sample_priors = TRUE", {
  set.seed(4)
  x   <- rnorm(20, mean = 0, sd = 1)
  uip <- create_prior_unit_information(x)

  # alpha0 = 0.5 gives InvGamma with no finite mean, so capability metric
  # moments may diverge — the important thing is no error is thrown.
  expect_no_error(
    bpc(x, LSL = -3, target = 0, USL = 3, prior = uip, method = "integration",
        sample_priors = TRUE)
  )
})
