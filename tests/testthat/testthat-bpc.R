
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


testthat::test_that("summary-statistics inputs are marked correctly for Stan", {

  prepared <- qc:::.bpc_data(x = NULL, mean = 5, sd = 2, N = 30)

  expect_equal(prepared$is_ss, 1L)
  expect_equal(prepared$N, 30L)
  expect_equal(prepared$ss_mean, as.array(c(5)))
  expect_equal(prepared$ss_sd, as.array(c(2)))
  expect_length(prepared$x, 0)
})


testthat::test_that("MCMC accepts summary-statistics inputs", {

  set.seed(1)
  x <- rnorm(20, mean = 5, sd = 1)

  fit <- bpc(
    NULL,
    LSL = 2, target = 5, USL = 8,
    mean = mean(x), sd = stats::sd(x), N = length(x),
    method = "mcmc",
    chains = 1, iter = 200, warmup = 100, cores = 1,
    silent = TRUE, seed = 1
  )

  expect_s3_class(fit, "bpc")
  expect_true(all(is.finite(fit$coefficients[c("Cp", "Cpk")])))
})


testthat::test_that("t-distribution rejects summary-statistics inputs up front", {

  expect_error(
    bpc(
      NULL,
      LSL = -3, target = 0, USL = 3,
      distribution = "t",
      mean = 0, sd = 1, N = 20,
      prior_nu = prior("exponential", list(1 / 30), list(2, Inf)),
      method = "mcmc",
      chains = 1, iter = 100, warmup = 50, cores = 1,
      silent = TRUE, seed = 1
    ),
    regexp = "Summary-statistics inputs .* `distribution = \"normal\"`"
  )

  expect_error(
    bpc(
      NULL,
      LSL = -3, target = 0, USL = 3,
      distribution = "t",
      prior_nu = prior("exponential", list(1 / 30), list(2, Inf)),
      method = "mcmc",
      chains = 1, iter = 100, warmup = 50, cores = 1,
      silent = TRUE, seed = 1
    ),
    regexp = "supply raw observations in `x`"
  )
})


testthat::test_that("fit with t-distribution", {

  set.seed(1)
  x <- rnorm(100, 0, 1)


  fit <- bpc(
    x, LSL = -3, target = 0, USL = 3, prior_nu = prior("exponential", list(1/30), list(2, Inf)),
    distribution = "t"
  )
  ss <- summary(fit)

  plot_density(fit, ci = "HPD")

  skip()
  # TODO: this confuses me, discuss this?
  lapply(fit$metrics, var)
  lapply(fit$metrics, var)
  other_cpc <- qc:::compute_E_abs_dev.default(fit, 0) # this function was renamed
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
