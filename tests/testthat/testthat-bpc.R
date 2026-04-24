testthat::test_that("mcmc bpc normal fits expose stable posterior summaries", {
  set.seed(1)
  x <- rnorm(80, mean = 10, sd = 2)

  fit <- bpc(
    x,
    LSL = 4,
    target = 10,
    USL = 16,
    method = "mcmc",
    prior = "Jeffreys",
    chains = 1,
    iter = 600,
    warmup = 200,
    cores = 1,
    silent = TRUE,
    seed = 1
  )

  expect_s3_class(fit, "bpc")
  expect_equal(fit$method, "mcmc")
  expect_named(fit$coefficients, c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm"))
  expect_true(all(is.finite(fit$coefficients)))
  expect_true(all(vapply(fit$metrics, function(draws) length(draws) > 0L, logical(1))))

  summary_fit <- summary(fit)
  expect_s3_class(summary_fit, "bpc_summary")
  expect_true(all(is.finite(summary_fit$summary$mean)))
  expect_true(all(is.finite(summary_fit$summary$median)))
  expect_true(all(is.finite(summary_fit$summary$sd)))
})


testthat::test_that("mcmc summary-statistics inputs stay aligned with raw-data fits", {
  set.seed(2)
  x <- rnorm(40, mean = 5, sd = 1)

  fit_raw <- bpc(
    x,
    LSL = 2,
    target = 5,
    USL = 8,
    method = "mcmc",
    prior = "Jeffreys",
    chains = 1,
    iter = 800,
    warmup = 200,
    cores = 1,
    silent = TRUE,
    seed = 2
  )
  fit_ss <- bpc(
    NULL,
    LSL = 2,
    target = 5,
    USL = 8,
    mean = mean(x),
    sd = stats::sd(x),
    N = length(x),
    method = "mcmc",
    prior = "Jeffreys",
    chains = 1,
    iter = 800,
    warmup = 200,
    cores = 1,
    silent = TRUE,
    seed = 2
  )

  expect_equal(fit_raw$stan_data$is_ss, 0L)
  expect_equal(fit_ss$stan_data$is_ss, 1L)

  sum_raw <- summary(fit_raw)
  sum_ss <- summary(fit_ss)
  rows_raw <- sum_raw$summary[match(c("Cp", "Cpk", "Cpm"), sum_raw$summary$metric), ]
  rows_ss <- sum_ss$summary[match(c("Cp", "Cpk", "Cpm"), sum_ss$summary$metric), ]

  expect_equal(rows_raw$mean, rows_ss$mean, tolerance = 0.12)
  expect_equal(rows_raw$median, rows_ss$median, tolerance = 0.12)
})


testthat::test_that("custom mcmc priors are preserved on the fit and give finite summaries", {
  set.seed(3)
  x <- rnorm(60, mean = 10, sd = 2)
  mu_prior <- prior("normal", list(mean = 10, sd = 4), list(lower = 0, upper = 20))
  sigma_prior <- prior("lognormal", list(meanlog = 0, sdlog = 0.5))

  fit <- bpc(
    x,
    LSL = 4,
    target = 10,
    USL = 16,
    method = "mcmc",
    prior = prior_independent(mu = mu_prior, sigma = sigma_prior),
    chains = 1,
    iter = 600,
    warmup = 200,
    cores = 1,
    silent = TRUE,
    seed = 3
  )

  expect_identical(fit$prior_map$mu, mu_prior)
  expect_identical(fit$prior_map$sigma, sigma_prior)
  expect_true(all(is.finite(fit$coefficients[c("Cp", "Cpk", "Cpm", "Cpc")])))

  summary_fit <- summary(fit)
  expect_true(all(is.finite(summary_fit$summary$lower)))
  expect_true(all(is.finite(summary_fit$summary$upper)))
})


testthat::test_that("t-distribution rejects summary-statistics inputs up front", {
  expect_error(
    bpc(
      NULL,
      LSL = -3,
      target = 0,
      USL = 3,
      distribution = "t",
      mean = 0,
      sd = 1,
      N = 20,
      prior = prior_independent(nu = prior("exponential", list(1 / 30), list(2, Inf))),
      method = "mcmc",
      chains = 1,
      iter = 100,
      warmup = 50,
      cores = 1,
      silent = TRUE,
      seed = 1
    ),
    regexp = "Summary-statistics inputs .* `distribution = \"normal\"`"
  )

  expect_error(
    bpc(
      NULL,
      LSL = -3,
      target = 0,
      USL = 3,
      distribution = "t",
      prior = prior_independent(nu = prior("exponential", list(1 / 30), list(2, Inf))),
      method = "mcmc",
      chains = 1,
      iter = 100,
      warmup = 50,
      cores = 1,
      silent = TRUE,
      seed = 1
    ),
    regexp = "supply raw observations in `x`"
  )
})


testthat::test_that("student-t mcmc fits support densities, plotting, and predictive samples", {
  set.seed(4)
  x <- 0.25 + 1.1 * stats::rt(120, df = 7)

  fit <- bpc(
    x,
    LSL = -5,
    target = 0.25,
    USL = 6,
    distribution = "t",
    prior = prior_independent(nu = prior("exp", list(rate = 1 / 30), list(lower = 2, upper = Inf))),
    method = "mcmc",
    chains = 1,
    iter = 800,
    warmup = 200,
    cores = 1,
    silent = TRUE,
    seed = 4
  )

  expect_true(all(vapply(fit$metrics, function(draws) all(is.finite(draws)), logical(1))))

  summary_fit <- summary(fit)
  expect_true(all(is.finite(summary_fit$summary$mean)))
  expect_true(all(is.finite(summary_fit$summary$sd)))

  density_data <- extract_density_data(fit, what = c("Cp", "Cpc"))
  expect_s3_class(density_data, "tbl_df")
  expect_equal(sort(unique(as.character(density_data$metric))), c("Cp", "Cpc"))

  plot_obj <- plot_density(fit, what = "Cpc", ci = "HPD")
  expect_s3_class(plot_obj, "ggplot")
  expect_false(inherits(try(ggplot2::ggplot_build(plot_obj), silent = TRUE), "try-error"))

  predictive <- extract_predictive_samples(fit)
  expect_true(is.numeric(predictive))
  expect_true(length(predictive) > 0L)
  expect_true(all(is.finite(predictive)))
})
