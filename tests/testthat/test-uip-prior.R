test_that("create_prior_unit_information returns a PriorConjugate with correct parameters", {
  set.seed(1)
  x    <- rnorm(20, mean = 5, sd = 2)
  xbar <- mean(x)
  s2   <- var(x)

  uip <- qc:::create_prior_unit_information(x)

  expect_s3_class(uip, "PriorConjugate")
  expect_equal(uip$mu0,   xbar)
  expect_equal(uip$k0,    1)
  expect_equal(uip$alpha0, 0.5)
  expect_equal(uip$beta0,  s2 / 2)
})

test_that("create_prior_unit_information is deterministic given the same data", {
  x   <- c(1.2, 3.4, 2.1, 4.5, 2.8)
  u1  <- qc:::create_prior_unit_information(x)
  u2  <- qc:::create_prior_unit_information(x)
  expect_identical(u1, u2)
})

test_that("create_prior_unit_information strips NAs before computing parameters", {
  x_clean <- c(1, 2, 4, 5)
  x_na    <- c(NA, 1, 2, NA, 4, 5, NA)
  expect_equal(qc:::create_prior_unit_information(x_clean),
               qc:::create_prior_unit_information(x_na))
})

test_that("create_prior_unit_information errors on fewer than 2 finite observations", {
  expect_error(qc:::create_prior_unit_information(c(1)))
  expect_error(qc:::create_prior_unit_information(c(NA, 1)))
  expect_error(qc:::create_prior_unit_information(numeric(0)))
})

test_that("create_prior_unit_information handles edge case of near-zero variance", {
  x <- c(1.0, 1.0 + 1e-10, 1.0 - 1e-10)
  expect_no_error(qc:::create_prior_unit_information(x))
  uip <- qc:::create_prior_unit_information(x)
  expect_true(is.finite(uip$beta0))
  expect_true(uip$beta0 >= 0)
})
