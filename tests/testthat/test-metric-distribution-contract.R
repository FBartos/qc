invisible(lapply(distribution_contract_cases(supported = TRUE), function(case) {
  test_that(case$label, {
    expect_distribution_contract_support(case)
  })
}))

invisible(lapply(distribution_contract_cases(supported = FALSE), function(case) {
  test_that(case$label, {
    expect_distribution_contract_rejection(case)
  })
}))

test_that("pc supports point-estimate-only fits when bootstrap is disabled", {
  fit <- distribution_contract_fit("pc-normal-raw")

  expect_s3_class(fit, "pc")
  expect_null(fit$metrics_boot)

  ss <- summary(fit)
  expect_s3_class(ss, "pc_summary")
  expect_false(attr(ss, "has_bootstrap"))
  expect_true(all(is.na(ss$summary$median)))
  expect_true(all(is.na(ss$summary$sd)))
  expect_true(all(is.na(ss$summary$lower)))
  expect_true(all(is.na(ss$summary$upper)))
  expect_error(extract_density_data(fit, what = "Cp"), "require bootstrap draws")
  expect_error(plot_density(fit, what = "Cp"), "require bootstrap draws")
})

test_that("pc bootstrap fits support shared density and plotting APIs", {
  fit <- distribution_contract_pc_bootstrap_fit()

  df_fit <- extract_density_data(fit, what = c("Cp", "Cpk"))
  ss <- summary(fit)
  df_summary <- extract_density_data(ss, what = c("Cp", "Cpk"))

  expect_s3_class(df_fit, "tbl_df")
  expect_equal(df_fit, df_summary)

  plt_fit <- plot_density(fit, what = c("Cp", "Cpk"))
  plt_summary <- plot_density(ss, what = c("Cp", "Cpk"))

  expect_s3_class(plt_fit, "ggplot")
  expect_s3_class(plt_summary, "ggplot")
})

test_that("integration results are first-class readers for density, intervals, and plotting", {
  fit <- distribution_contract_fit("bpc-integration-normal-raw")

  df_from_result <- extract_density_data(fit$integration_result, what = c("Cp", "Cpk"))
  df_from_metrics <- extract_density_data(fit$metrics, what = c("Cp", "Cpk"))
  expect_equal(df_from_result, df_from_metrics)

  df_points <- qc:::extract_point_estimates(
    obj = fit$integration_result,
    what = c("Cp", "Cpk"),
    point_estimate = "mean",
    dfDensity = df_from_result
  )
  expect_s3_class(df_points, "tbl_df")
  expect_equal(nrow(df_points), 2)

  ci_data <- qc:::extract_ci_data(
    obj = fit$integration_result,
    what = c("Cp", "Cpk"),
    ci = "HPD",
    ci_level = 0.95,
    dfDensity = df_from_result
  )
  expect_type(ci_data, "list")
  expect_true(all(c("dfCi", "dfArea") %in% names(ci_data)))

  plt <- plot_density(fit$integration_result, what = c("Cp", "Cpk"))
  expect_s3_class(plt, "ggplot")
})

test_that("default E_abs_dev fallback computes one conditional expectation per posterior draw", {
  samples <- structure(list(theta = c(1, 2, 3)), class = "mock_distribution")
  calls <- 0L

  testthat::local_mocked_bindings(
    samples_to_posterior_predictives = function(x) {
      calls <<- calls + 1L
      expect_equal(length(unique(x$theta)), 1L)
      x$theta
    },
    .package = "qc"
  )

  out <- qc:::samples_to_E_abs_dev.default(samples, target = 0, n_inner = 4L)

  expect_equal(out, c(1, 2, 3))
  expect_equal(calls, 3L)
})

test_that("default E_abs_dev fallback validates aligned posterior draw lengths", {
  bad_samples <- structure(
    list(mu = c(0, 1, 2), sigma = c(1, 2)),
    class = "bad_distribution"
  )

  expect_error(
    qc:::samples_to_E_abs_dev.default(bad_samples, target = 0),
    regexp = "same length"
  )
})
