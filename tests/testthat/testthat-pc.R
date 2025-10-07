
testthat::test_that("default settings", {

  set.seed(1)
  x <- rnorm(100, 0, 1)

  fit <- pc(x, LSL = -3, target = 0, USL = 3)

  fit
  summary(fit)
})


testthat::test_that("fit with t-distribution", {

  set.seed(1)
  x <- rnorm(100, 0, 1)


  fit <- pc(
    x, LSL = -3, target = 0, USL = 3,
    distribution = "t"
  )
  # TODO: how do we deal with NAs from incomplete model fits? They happen within the bootstrap
  ss <- summary(fit)


  fit <- pc(
    x, LSL = -3, target = 0, USL = 3,
    distribution = "t", parallel = TRUE
  )
})
