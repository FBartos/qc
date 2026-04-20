testthat::test_that("integration handles generic sigma priors whose support excludes the sample sd", {
  set.seed(1)
  x <- rnorm(40, 0, 1)

  fit <- bpc(
    x, LSL = -3, target = 0, USL = 3,
    prior_mu = prior("normal", list(0, 1)),
    prior_sigma = prior("uniform", list(0, 0.15)),
    method = "integration"
  )

  testthat::expect_s3_class(fit, "bpc")
  testthat::expect_equal(fit$method, "integration")
  testthat::expect_true(all(is.finite(unname(fit$coefficients))))
})
