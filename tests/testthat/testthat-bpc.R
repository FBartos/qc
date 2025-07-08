
testthat::test_that("default settings", {

  set.seed(1)
  x <- rnorm(100, 0, 1)

  fit <- bpc(x, LSL = -3, target = 0, USL = 3)

  fit
  summary(fit)
})


testthat::test_that("custom priors", {

  set.seed(1)
  x <- rnorm(100, 0, 1)

  fit <- bpc(
    x, LSL = -3, target = 0, USL = 3,
    prior_mu    = prior("normal",  list(15, 10), list(10, Inf)),
    prior_sigma = prior("uniform", list(5, 10)))

  fit
  summary(fit)
  summary(fit, LSL = -2.5, target = 0, USL = 2.5)
})


testthat::test_that("fit control", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(x, LSL = -3, target = 0, USL = 3, chains = 1, warmup = 100, iter = 200, silent = TRUE, seed = 1)

  fit
  summary(fit)
  summary(fit, LSL = -2, target = 0, USL = 2)

})


testthat::test_that("fit with t-distribution", {

  set.seed(1)
  x <- rnorm(100, 0, 1)

  fit <- bpc(x, LSL = -3, target = 0, USL = 3, chains = 1, warmup = 100, iter = 200, silent = TRUE, seed = 1)

  # debugonce(qc:::.bpc_compute_capability_metrics)
  fit <- bpc(
    x, LSL = -3, target = 0, USL = 3#,
    # prior_nu    = prior("exp",  list(1))
  )
  ss <- summary(fit)

  plot_density(fit, ci = "HPD")

  # TODO: this confuses me, discuss this?
  lapply(fit$metrics, var)
  lapply(fit$metrics, var)
  other_cpc <- qc:::.bpc_compute_E_abs_dev.default(fit, 0)
  var(other_cpc)
  freq_cpc <- 6 / (6 * sqrt(pi /  2) * mean(abs(x)))
  mean(fit$metrics$Cpc) # looks better!
  mean(other_cpc)

  fit2 <- bpc(
    x, LSL = -3, target = 0, USL = 3,
    prior_nu    = prior("uniform", list(0, 100))
  )

  fit
  summary(fit)
  summary(fit, LSL = -2, target = 0, USL = 2)

  fit2
  summary(fit2)
  summary(fit2, LSL = -2, target = 0, USL = 2)
})
