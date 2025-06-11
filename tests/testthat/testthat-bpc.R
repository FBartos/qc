
testthat::test_that("default settings", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(x)

})


testthat::test_that("custom priors", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(
    x,
    prior_mu    = prior("normal",  list(15, 10), list(10, Inf)),
    prior_sigma = prior("uniform", list(5, 10)))

  summary(fit$stanfit)
})
