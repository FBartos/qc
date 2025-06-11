
testthat::test_that("default bpc", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(x)

})
