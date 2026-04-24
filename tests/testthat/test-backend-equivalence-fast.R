test_that("MCMC backend supports the shared density, interval, and plotting contracts", {
  fixture <- get_backend_test_fixture("backend-fast-mcmc-contract", function() {
    set.seed(101)
    x <- rnorm(30, mean = 10, sd = 2)

    qc::bpc(
      x,
      LSL = 4, target = 10, USL = 16,
      method = "mcmc",
      prior = "Jeffreys",
      chains = 1,
      iter = 1000,
      warmup = 250,
      cores = 1,
      silent = TRUE,
      seed = 101
    )
  })

  ss <- summary(fixture)
  df_density <- extract_density_data(fixture, what = c("Cp", "Cpk"))
  df_points <- extract_point_estimates(
    obj = fixture,
    what = c("Cp", "Cpk"),
    point_estimate = "median",
    dfDensity = df_density
  )
  ci_95 <- extract_ci_data(
    obj = fixture,
    what = c("Cp", "Cpk"),
    ci = "central",
    ci_level = 0.95,
    dfDensity = df_density
  )
  ci_50 <- extract_ci_data(
    obj = fixture,
    what = c("Cp", "Cpk"),
    ci = "central",
    ci_level = 0.50,
    dfDensity = df_density
  )

  expect_s3_class(ss, "bpc_summary")
  expect_s3_class(df_density, "tbl_df")
  expect_equal(names(df_density), c("x", "density", "metric"))
  expect_equal(sort(unique(as.character(df_density$metric))), c("Cp", "Cpk"))
  expect_true(all(df_density$density >= 0))

  expect_s3_class(df_points, "tbl_df")
  expect_equal(names(df_points), c("x", "y", "metric"))
  expect_equal(nrow(df_points), 2L)
  expect_true(all(is.finite(df_points$x)))

  expect_type(ci_95, "list")
  expect_equal(names(ci_95), c("dfCi", "dfArea"))
  expect_equal(nrow(ci_95$dfCi), 2L)
  expect_equal(nrow(ci_50$dfCi), 2L)
  expect_true(all(ci_50$dfCi$xmin >= ci_95$dfCi$xmin))
  expect_true(all(ci_50$dfCi$xmax <= ci_95$dfCi$xmax))

  expect_s3_class(plot_density(fixture, what = c("Cp", "Cpk")), "ggplot")
  expect_s3_class(plot_density(ss, what = "Cp"), "ggplot")
})

test_that("Matched conjugate backends agree for raw and summary-statistics inputs", {
  raw_pair <- get_backend_test_fixture("backend-fast-conjugate-raw", function() {
    set.seed(11)
    x <- rnorm(30, mean = 10, sd = 2)

    fit_backend_pair(
      x = x,
      LSL = 4, target = 10, USL = 16,
      mcmc_chains = 1,
      mcmc_iter = 1500,
      mcmc_warmup = 500,
      seed = 11
    )
  })
  ss_pair <- get_backend_test_fixture("backend-fast-conjugate-ss", function() {
    set.seed(11)
    x <- rnorm(30, mean = 10, sd = 2)

    fit_backend_pair(
      x = NULL,
      mean = mean(x),
      sd = stats::sd(x),
      N = length(x),
      LSL = 4, target = 10, USL = 16,
      mcmc_chains = 1,
      mcmc_iter = 1500,
      mcmc_warmup = 500,
      seed = 11
    )
  })

  sum_raw <- expect_fit_summaries_close(
    raw_pair,
    metrics = c("Cp", "Cpk", "Cpm"),
    columns = c("mean", "median", "lower", "upper"),
    abs_tol = 0.15,
    rel_tol = 0.03
  )
  sum_ss <- expect_fit_summaries_close(
    ss_pair,
    metrics = c("Cp", "Cpk", "Cpm"),
    columns = c("mean", "median", "lower", "upper"),
    abs_tol = 0.15,
    rel_tol = 0.03
  )

  expect_summary_columns_close(
    summary_int = sum_raw$integration,
    summary_mcmc = sum_ss$integration,
    metrics = c("Cp", "Cpk", "Cpm"),
    columns = c("mean", "median", "lower", "upper"),
    abs_tol = 1e-8
  )
  expect_summary_columns_close(
    summary_int = sum_raw$mcmc,
    summary_mcmc = sum_ss$mcmc,
    metrics = c("Cp", "Cpk", "Cpm"),
    columns = c("mean", "median", "lower", "upper"),
    abs_tol = 0.12,
    rel_tol = 0.03
  )
})

test_that("Matched conjugate backends agree under sigma scaling", {
  pair <- get_backend_test_fixture("backend-fast-conjugate-raw", function() {
    set.seed(11)
    x <- rnorm(30, mean = 10, sd = 2)

    fit_backend_pair(
      x = x,
      LSL = 4, target = 10, USL = 16,
      mcmc_chains = 1,
      mcmc_iter = 1500,
      mcmc_warmup = 500,
      seed = 11
    )
  })

  sum_scaled <- expect_fit_summaries_close(
    pair,
    metrics = c("Cp", "Cpk", "Cpm"),
    columns = c("mean", "median", "lower", "upper"),
    abs_tol = 0.20,
    rel_tol = 0.03,
    sigma = 4
  )

  expect_interval_summary_close(
    summary_int = sum_scaled$integration,
    summary_mcmc = sum_scaled$mcmc,
    metrics = c("Cp", "Cpk", "Cpm"),
    abs_tol = 0.05
  )
})

test_that("Matched conjugate backends agree on negative-support upper-tail events", {
  pair <- get_backend_test_fixture("backend-fast-negative-upper", function() {
    set.seed(21)
    x <- rnorm(40, mean = 20, sd = 1)

    fit_backend_pair(
      x = x,
      LSL = 0, target = 5, USL = 10,
      mcmc_chains = 1,
      mcmc_iter = 1500,
      mcmc_warmup = 500,
      seed = 21
    )
  })

  summaries <- expect_fit_summaries_close(
    pair,
    metrics = c("Cpu", "Cpk"),
    columns = c("median"),
    abs_tol = 0.15,
    rel_tol = 0.03,
    interval_probability = c(0)
  )

  expect_true(all(summary_metric_rows(summaries$integration$summary, c("Cpu", "Cpk"))$median < 0))
  expect_true(all(summary_metric_rows(summaries$mcmc$summary, c("Cpu", "Cpk"))$median < 0))
  expect_interval_summary_close(
    summary_int = summaries$integration,
    summary_mcmc = summaries$mcmc,
    metrics = c("Cpu", "Cpk"),
    abs_tol = 0.05
  )
})

test_that("Matched conjugate backends agree on negative-support lower-tail events", {
  pair <- get_backend_test_fixture("backend-fast-negative-lower", function() {
    set.seed(22)
    x <- rnorm(40, mean = -10, sd = 1)

    fit_backend_pair(
      x = x,
      LSL = 0, target = 5, USL = 10,
      mcmc_chains = 1,
      mcmc_iter = 1500,
      mcmc_warmup = 500,
      seed = 22
    )
  })

  summaries <- expect_fit_summaries_close(
    pair,
    metrics = c("Cpl", "Cpk"),
    columns = c("median"),
    abs_tol = 0.15,
    rel_tol = 0.03,
    interval_probability = c(0)
  )

  expect_true(all(summary_metric_rows(summaries$integration$summary, c("Cpl", "Cpk"))$median < 0))
  expect_true(all(summary_metric_rows(summaries$mcmc$summary, c("Cpl", "Cpk"))$median < 0))
  expect_interval_summary_close(
    summary_int = summaries$integration,
    summary_mcmc = summaries$mcmc,
    metrics = c("Cpl", "Cpk"),
    abs_tol = 0.05
  )
})

test_that("Matched conjugate backends agree for target-sensitive metrics", {
  pair <- get_backend_test_fixture("backend-fast-target-sensitive", function() {
    set.seed(23)
    x <- rnorm(40, mean = 10.8, sd = 1.2)

    fit_backend_pair(
      x = x,
      LSL = 4, target = 10, USL = 16,
      mcmc_chains = 1,
      mcmc_iter = 1500,
      mcmc_warmup = 500,
      seed = 23
    )
  })

  summaries <- expect_fit_summaries_close(
    pair,
    metrics = c("Cpc", "Cpm"),
    columns = c("mean", "median", "lower", "upper"),
    abs_tol = 0.15,
    rel_tol = 0.03
  )

  expect_interval_summary_close(
    summary_int = summaries$integration,
    summary_mcmc = summaries$mcmc,
    metrics = c("Cpc", "Cpm"),
    abs_tol = 0.05
  )
})
