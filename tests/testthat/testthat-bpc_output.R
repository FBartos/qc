testthat::test_that("bpc print methods expose the expected section headers", {
  set.seed(1)
  x <- rnorm(60, 10, 2)
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

  expect_output(print(fit), "Bayesian Process Capability")
  expect_output(print(fit), "Cp")

  summary_fit <- summary(fit)
  expect_output(print(summary_fit), "Interval Probability")
  expect_output(print(summary_fit), "Cpk")
})


testthat::test_that("custom and support intervals pass through the plotting helpers exactly", {
  set.seed(2)
  x <- 0.5 + stats::rt(120, df = 8)
  fit <- bpc(
    x,
    LSL = -5,
    target = 0.5,
    USL = 6,
    distribution = "t",
    prior = prior_independent(nu = prior("exp", list(rate = 1 / 30), list(lower = 2, upper = Inf))),
    method = "mcmc",
    chains = 1,
    iter = 800,
    warmup = 200,
    cores = 1,
    silent = TRUE,
    seed = 2
  )

  df_cp <- extract_density_data(fit, what = "Cp")
  custom_ci <- extract_ci_data(
    fit,
    what = "Cp",
    ci = "custom",
    ci_level = 0.95,
    dfDensity = df_cp,
    ci_custom_left = 0.8,
    ci_custom_right = 1.2
  )
  expect_equal(custom_ci$dfCi$xmin, 0.8)
  expect_equal(custom_ci$dfCi$xmax, 1.2)

  df_cpc <- extract_density_data(fit, what = "Cpc")
  support_ci <- extract_ci_data(
    fit,
    what = "Cpc",
    ci = "support",
    ci_level = 0.95,
    dfDensity = df_cpc,
    bf_support = list(lower = 0.9, upper = 1.4)
  )
  expect_equal(support_ci$dfCi$xmin, 0.9)
  expect_equal(support_ci$dfCi$xmax, 1.4)

  plot_custom <- plot_density(
    fit,
    what = "Cp",
    ci = "custom",
    ci_custom_left = 0.8,
    ci_custom_right = 1.2
  )
  plot_support <- plot_density(
    fit,
    what = "Cpc",
    ci = "support",
    bf_support = list(lower = 0.9, upper = 1.4)
  )

  expect_s3_class(plot_custom, "ggplot")
  expect_s3_class(plot_support, "ggplot")

  built_custom <- ggplot2::ggplot_build(plot_custom)
  built_support <- ggplot2::ggplot_build(plot_support)

  custom_layers <- attr(built_custom, "data")
  support_layers <- attr(built_support, "data")
  expect_length(custom_layers, 4L)
  expect_length(support_layers, 4L)
  expect_equal(custom_layers[[3]]$xmax[[1]], 1.2)
  expect_equal(support_layers[[3]]$xmax[[1]], 1.4)
  expect_equal(custom_layers[[4]]$label[[1]], "Custom CI [0.800, 1.200]")
  expect_equal(support_layers[[4]]$label[[1]], "Support [0.900, 1.400]")
})
