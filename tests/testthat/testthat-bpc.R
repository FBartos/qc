
testthat::test_that("default settings", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(x)

  fit$stanfit
})


testthat::test_that("custom priors", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(
    x,
    prior_mu    = prior("normal",  list(15, 10), list(10, Inf)),
    prior_sigma = prior("uniform", list(5, 10)))

  fit$stanfit
})


testthat::test_that("fit control", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(x, chains = 1, warmup = 100, iter = 200, silent = TRUE, seed = 1)

  fit$stanfit
})
