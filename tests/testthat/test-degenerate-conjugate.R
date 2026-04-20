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
    qc::bpc(x, LSL = LSL, target = target, USL = USL, method = "integration")
  )
  expect_true(all(is.infinite(fit$coefficients)))

  summary_fit <- expect_no_error(summary(fit))
  expect_true(all(is.infinite(summary_fit$summary$mean)))

  pred <- qc::extract_predictive_samples(fit, n_samples = 128L)
  expect_equal(pred, rep(5, 128))
})

test_that("degenerate conjugate Cpm and Cpc collapse to finite point masses off target", {
  x <- rep(5, 10)
  LSL <- 0
  USL <- 10
  target <- 3
  prior <- qc:::create_prior_conjugate()

  expected_cpm <- (USL - LSL) / ((2 * 3) * abs(5 - target))
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
