test_that("constant-data conjugate posterior yields infinite interior capability metrics", {
  x <- rep(5, 10)
  LSL <- 0
  USL <- 10
  target <- 5
  prior <- qc:::create_prior_conjugate()

  metrics <- c("Cp", "Cpu", "Cpl", "Cpk", "Cpm", "Cpc")
  moments <- lapply(metrics, function(metric) {
    qc:::compute_metric_moments(x, LSL, USL, prior, metric = metric, target = target)
  })
  names(moments) <- metrics

  for (metric in metrics) {
    expect_true(is.infinite(moments[[metric]]$mean), info = metric)
    expect_true(is.infinite(moments[[metric]]$sd), info = metric)
  }

  fit <- expect_no_error(
    qc::bpc(x, LSL = LSL, target = target, USL = USL,
            method = "integration", prior = "Jeffreys")
  )
  expect_true(all(is.infinite(fit$coefficients)))

  summary_fit <- expect_no_error(summary(fit))
  expect_true(all(is.infinite(summary_fit$summary$mean)))

  pred <- qc::extract_predictive_samples(fit, n_samples = 128L)
  expect_equal(pred, rep(5, 128))
})

test_that("single-observation Jeffreys conjugate posteriors are rejected before degeneracy handling", {
  x <- 1
  LSL <- 0
  USL <- 2
  target <- 1
  prior <- qc:::create_prior_conjugate()

  expect_error(
    qc:::compute_metric_moments(x, LSL, USL, prior, metric = "Cp", target = target),
    regexp = "improper"
  )

  expect_error(
    qc:::make_solver(x, LSL, USL, prior, metric = "Cp", target = target),
    regexp = "improper"
  )

  expect_error(
    qc:::make_density_solver(x, LSL, USL, prior, metric = "Cp", target = target),
    regexp = "improper"
  )

  expect_error(
    qc:::compute_cpk_prob_integration(
      x, LSL, USL, c(0, 1), prior,
      metric = "Cp", target = target
    ),
    regexp = "improper"
  )

  expect_error(
    qc:::analyze_capability_integration(x, LSL, USL, prior, metric = "Cp", target = target),
    regexp = "improper"
  )

  expect_error(
    qc::bpc(x, LSL = LSL, target = target, USL = USL,
            method = "integration", prior = "Jeffreys"),
    regexp = "improper"
  )
})

test_that("degenerate +Inf integration summaries keep mass in the open upper tail", {
  fit <- qc::bpc(rep(5, 10), LSL = 0, target = 5, USL = 10,
                 method = "integration", prior = "Jeffreys")
  ss <- summary(fit, interval_probability = c(1.00, 1.33, 1.50, 2.00))

  expect_true(all(vapply(fit$integration_result$results, function(x) {
    identical(x$degenerate$type, "pos_inf")
  }, logical(1))))

  expect_true(all(ss$interval_summary$`(2, Inf]` == 1))
  expect_true(all(ss$interval_summary$`[-Inf,1]` == 0))
  expect_true(all(ss$interval_summary$`(1,1.33]` == 0))
  expect_true(all(ss$interval_summary$`(1.33,1.5]` == 0))
  expect_true(all(ss$interval_summary$`(1.5,2]` == 0))
})

test_that("degenerate conjugate Cpm and Cpc collapse to finite point masses off target", {
  x <- rep(5, 10)
  LSL <- 0
  USL <- 10
  target <- 3
  prior <- qc:::create_prior_conjugate()

  expected_cpm <- qc:::.cpm_spec_distance(LSL, USL, target) / (3 * abs(5 - target))
  expected_cpc <- (USL - LSL) / ((2 * 3) * sqrt(pi / 2) * abs(5 - target))

  cpm_mom <- qc:::compute_metric_moments(x, LSL, USL, prior, metric = "Cpm", target = target)
  cpc_mom <- qc:::compute_metric_moments(x, LSL, USL, prior, metric = "Cpc", target = target)

  expect_equal(cpm_mom$mean, expected_cpm)
  expect_equal(cpm_mom$sd, 0)
  expect_equal(cpc_mom$mean, expected_cpc)
  expect_equal(cpc_mom$sd, 0)

  cpm_result <- qc:::analyze_capability_integration(
    x, LSL, USL, prior, metric = "Cpm", target = target
  )
  expect_equal(unname(cpm_result$stats["Mean"]), expected_cpm)
  expect_equal(unname(cpm_result$stats["Median"]), expected_cpm)
  expect_equal(unname(cpm_result$stats["Q2.5"]), expected_cpm)
  expect_equal(unname(cpm_result$stats["Q97.5"]), expected_cpm)

  expect_equal(
    qc:::compute_cpk_prob_integration(
      x, LSL, USL,
      c(expected_cpm - 0.01, expected_cpm + 0.01),
      prior, metric = "Cpm", target = target
    ),
    1
  )
  expect_equal(
    qc:::compute_cpk_prob_integration(
      x, LSL, USL,
      c(expected_cpm + 0.01, expected_cpm + 0.1),
      prior, metric = "Cpm", target = target
    ),
    0
  )
})

test_that("extract_point_estimates keeps compatibility with stats-only degenerate entries", {
  cpm_fit <- qc::bpc(rep(5, 10), LSL = 0, target = 3, USL = 10,
                     method = "integration", prior = "Jeffreys")
  cpm_result <- cpm_fit$integration_result$results$Cpm
  cpm_density <- qc::extract_density_data(cpm_fit, what = "Cpm")
  expected_cpm <- qc:::.cpm_spec_distance(0, 10, 3) / (3 * abs(5 - 3))

  cpm_mode_full <- qc:::extract_point_estimates(
    obj = NULL,
    what = "Cpm",
    point_estimate = "mode",
    dfDensity = cpm_density,
    stats_list = list(Cpm = cpm_result)
  )
  cpm_mode_stats <- qc:::extract_point_estimates(
    obj = NULL,
    what = "Cpm",
    point_estimate = "mode",
    dfDensity = cpm_density,
    stats_list = list(Cpm = cpm_result$stats)
  )

  expect_equal(cpm_mode_full$x, expected_cpm)
  expect_equal(cpm_mode_stats$x, expected_cpm)

  cpu_fit <- qc::bpc(rep(5, 10), LSL = 0, target = 5, USL = 10,
                     method = "integration", prior = "Jeffreys")
  cpu_result <- cpu_fit$integration_result$results$Cpu
  cpu_density <- qc::extract_density_data(cpu_fit, what = "Cpu")

  cpu_mode_full <- qc:::extract_point_estimates(
    obj = NULL,
    what = "Cpu",
    point_estimate = "mode",
    dfDensity = cpu_density,
    stats_list = list(Cpu = cpu_result)
  )
  cpu_mode_stats <- qc:::extract_point_estimates(
    obj = NULL,
    what = "Cpu",
    point_estimate = "mode",
    dfDensity = cpu_density,
    stats_list = list(Cpu = cpu_result$stats)
  )

  expect_equal(nrow(cpu_mode_full), 0)
  expect_equal(nrow(cpu_mode_stats), 0)
})

test_that("degenerate conjugate boundary case yields exact normal limit for Cpu and Cpk", {
  x <- rep(10, 10)
  LSL <- 0
  USL <- 10
  target <- 5
  prior <- qc:::create_prior_conjugate()
  expected_sd <- 1 / (3 * sqrt(length(x)))

  cpu_mom <- qc:::compute_metric_moments(x, LSL, USL, prior, metric = "Cpu", target = target)
  cpk_mom <- qc:::compute_metric_moments(x, LSL, USL, prior, metric = "Cpk", target = target)

  expect_equal(cpu_mom$mean, 0, tolerance = 1e-12)
  expect_equal(cpu_mom$sd, expected_sd, tolerance = 1e-12)
  expect_equal(cpk_mom$mean, 0, tolerance = 1e-12)
  expect_equal(cpk_mom$sd, expected_sd, tolerance = 1e-12)

  cpu_solver <- qc:::make_solver(x, LSL, USL, prior, metric = "Cpu", target = target)
  cpk_solver <- qc:::make_solver(x, LSL, USL, prior, metric = "Cpk", target = target)

  expect_equal(cpu_solver(0), 0.5, tolerance = 1e-12)
  expect_equal(cpk_solver(0), 0.5, tolerance = 1e-12)
  expect_equal(
    cpu_solver(0.1),
    stats::pnorm(0.1, mean = 0, sd = expected_sd, lower.tail = FALSE),
    tolerance = 1e-12
  )
  expect_equal(
    qc:::compute_cpk_prob_integration(
      x, LSL, USL, c(-Inf, 0), prior, metric = "Cpu", target = target
    ),
    0.5,
    tolerance = 1e-12
  )
})
