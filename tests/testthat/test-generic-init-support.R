testthat::test_that("generic integration finds feasible initials when sigma support excludes the sample sd", {
  x <- c(-0.12, 0.12)
  x_bar <- mean(x)
  sse <- sum((x - x_bar)^2)
  init_sd <- max(sqrt(sse / max(length(x) - 1, 1)), 1e-6)

  prior_obj <- qc:::.bayestools_to_integration_prior(
    prior("normal", list(0, 1)),
    prior("uniform", list(0, 0.10))
  )$prior

  log_post <- function(mu, sigma) {
    if (sigma <= 0) {
      return(-Inf)
    }

    log_lik <- -length(x) * log(sigma) -
      (sse + length(x) * (mu - x_bar)^2) / (2 * sigma^2)

    log_lik + prior_obj$log_dens(mu, sigma)
  }

  init <- qc:::.find_feasible_generic_init(
    init_mu = x_bar,
    init_sd = init_sd,
    prior = prior_obj,
    log_post = log_post
  )

  testthat::expect_gt(stats::sd(x), 0.10)
  testthat::expect_true(is.finite(init$mu))
  testthat::expect_true(is.finite(init$sigma))
  testthat::expect_gt(init$sigma, 0)
  testthat::expect_lte(init$sigma, 0.10)
  testthat::expect_true(is.finite(log_post(init$mu, init$sigma)))
})
