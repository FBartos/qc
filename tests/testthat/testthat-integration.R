testthat::test_that("integration method returns valid bpc object", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(x, LSL = 2, target = 10, USL = 18, method = "integration")

  # Check object structure
  expect_s3_class(fit, "bpc")
  expect_equal(fit$method, "integration")
  expect_true(!is.null(fit$integration_result))
  expect_true(!is.null(fit$metrics))
  expect_true(!is.null(fit$coefficients))

  # Check all metrics are present and positive
  expected_metrics <- c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm")
  expect_equal(names(fit$coefficients), expected_metrics)
  expect_true(all(fit$coefficients > 0))
  expect_true(all(is.finite(fit$coefficients)))
})


testthat::test_that("integration method uses the summary-statistics path", {

  set.seed(123)
  x <- rnorm(30, mean = 5, sd = 1)

  fit_raw <- bpc(x, LSL = 2, target = 5, USL = 8, method = "integration")
  fit_ss <- bpc(
    NULL,
    LSL = 2, target = 5, USL = 8,
    mean = mean(x), sd = stats::sd(x), N = length(x),
    method = "integration"
  )

  expect_equal(fit_ss$coefficients, fit_raw$coefficients, tolerance = 1e-10)
  expect_equal(
    fit_ss$integration_result$cached_state[c("n", "x_bar", "sse")],
    fit_raw$integration_result$cached_state[c("n", "x_bar", "sse")],
    tolerance = 1e-10
  )
})


testthat::test_that("integration method respects sigma scaling", {

  set.seed(123)
  x <- rnorm(40, mean = 10, sd = 1.5)

  fit_sigma3 <- bpc(x, LSL = 6, target = 10, USL = 14,
                    method = "integration", sigma = 3)
  fit_sigma4 <- bpc(x, LSL = 6, target = 10, USL = 14,
                    method = "integration", sigma = 4)

  expect_equal(fit_sigma4$integration_result$sigma, 4)

  for (metric in c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm")) {
    mean3 <- unname(fit_sigma3$metrics[[metric]]["Mean"])
    mean4 <- unname(fit_sigma4$metrics[[metric]]["Mean"])

    expect_false(isTRUE(all.equal(mean3, mean4)))
    expect_equal(mean4 / mean3, 3 / 4, tolerance = 0.03,
                 info = metric)
  }

  summary_default <- summary(fit_sigma4, LSL = 5.5, target = 10, USL = 14.5)
  summary_explicit <- summary(fit_sigma4, LSL = 5.5, target = 10, USL = 14.5,
                              sigma = 4)

  expect_equal(summary_default$summary, summary_explicit$summary)
  expect_equal(summary_default$interval_summary, summary_explicit$interval_summary)
})


testthat::test_that("density solver vs survival function methods agree", {

  set.seed(42)
  x <- rnorm(30, mean = 5, sd = 1)
  LSL <- 2; USL <- 8
  target <- 5
  prior <- qc:::create_prior_conjugate()

  # Test each supported metric
  for (metric in c("Cpk", "Cp", "Cpu", "Cpl", "Cpm", "Cpc")) {

    # Use density solver (new default)
    result_density <- qc:::analyze_capability_integration(x, LSL, USL, prior,
                                                           metric = metric,
                                                          target = target,
                                                           use_density_solver = TRUE)

    # Use survival function + finite diff (old method)
    result_survival <- qc:::analyze_capability_integration(x, LSL, USL, prior,
                                                            metric = metric,
                                                            use_density_solver = FALSE)

    # Stats should agree within tolerance
    for (stat in c("Mean", "Median", "SD", "Q2.5", "Q97.5")) {
      diff <- abs(result_density$stats[stat] - result_survival$stats[stat])
      rel_diff <- diff / max(abs(result_survival$stats[stat]), 0.01)
      expect_true(
        rel_diff < 0.05 || diff < 0.02,
        info = sprintf("%s/%s: density=%.4f, survival=%.4f, diff=%.4f",
                       metric, stat, result_density$stats[stat],
                       result_survival$stats[stat], diff)
      )
    }
  }
})


testthat::test_that("Cpc depends on target off the midpoint", {

  set.seed(42)
  x <- rnorm(30, mean = 10, sd = 1)
  LSL <- 4
  USL <- 16
  prior <- qc:::create_prior_conjugate()

  density_t10 <- qc:::analyze_capability_integration(
    x, LSL, USL, prior,
    metric = "Cpc",
    target = 10,
    use_density_solver = TRUE
  )
  density_t11 <- qc:::analyze_capability_integration(
    x, LSL, USL, prior,
    metric = "Cpc",
    target = 11,
    use_density_solver = TRUE
  )
  survival_t10 <- qc:::analyze_capability_integration(
    x, LSL, USL, prior,
    metric = "Cpc",
    target = 10,
    use_density_solver = FALSE
  )
  survival_t11 <- qc:::analyze_capability_integration(
    x, LSL, USL, prior,
    metric = "Cpc",
    target = 11,
    use_density_solver = FALSE
  )

  expect_gt(abs(density_t10$stats["Mean"] - density_t11$stats["Mean"]), 0.10)
  expect_gt(abs(survival_t10$stats["Mean"] - survival_t11$stats["Mean"]), 0.10)

  density_vs_survival <- abs(density_t11$stats["Mean"] - survival_t11$stats["Mean"])
  rel_diff <- density_vs_survival / max(abs(survival_t11$stats["Mean"]), 0.01)
  expect_true(
    rel_diff < 0.05 || density_vs_survival < 0.02,
    info = sprintf(
      "Cpc target=11: density=%.4f, survival=%.4f, diff=%.4f",
      density_t11$stats["Mean"],
      survival_t11$stats["Mean"],
      density_vs_survival
    )
  )
})


testthat::test_that("density solver produces normalized PDFs", {

  set.seed(42)
  x <- rnorm(30, mean = 5, sd = 1)
  LSL <- 2; USL <- 8
  target <- 5
  prior <- qc:::create_prior_conjugate()

  for (metric in c("Cpk", "Cp", "Cpu", "Cpl", "Cpm", "Cpc")) {
    pdf_fn <- qc:::make_density_solver(x, LSL, USL, prior, metric = metric, target = target)
    pdf_fn_vec <- Vectorize(pdf_fn)

    # PDF should integrate to ~1
    total_mass <- integrate(pdf_fn_vec, 0, 3, subdivisions = 100)$value
    expect_true(
      abs(total_mass - 1) < 0.01,
      info = sprintf("%s: total_mass = %.4f", metric, total_mass)
    )
  }
})


testthat::test_that("integration handles semi-conjugate mu priors with Cp density", {

  set.seed(1)
  x <- rnorm(30, mean = 10, sd = 2)

  fit <- bpc(
    x,
    LSL = 2, target = 10, USL = 18,
    method = "integration",
    prior = prior_independent(
      mu = "Jeffreys_mu",
      sigma = BayesTools::prior("gamma", list(2, 1))
    )
  )

  expect_s3_class(fit, "bpc")

  ss <- summary(fit)
  expect_s3_class(ss, "bpc_summary")
  expect_true("Cp" %in% ss$summary$metric)
  expect_true(all(is.finite(ss$summary$mean)))
})


testthat::test_that("semi-conjugate survival solvers agree with density mass", {
  skip_if_not_installed("BayesTools")

  set.seed(123)
  x <- rnorm(30, mean = 10, sd = 2)
  LSL <- 2
  USL <- 18
  target <- 10
  bounds <- c(0.5, 1.5)

  prior_cases <- list(
    semi_mu = qc:::.bayestools_to_integration_prior(
      "Jeffreys_mu",
      BayesTools::prior("gamma", list(2, 1))
    )$prior,
    semi_sigma = qc:::.bayestools_to_integration_prior(
      BayesTools::prior("normal", list(10, 5)),
      "Jeffreys_sigma"
    )$prior
  )

  for (case_name in names(prior_cases)) {
    prior <- prior_cases[[case_name]]
    for (metric in c("Cpu", "Cpl", "Cpk")) {
      S <- qc:::make_solver(x, LSL, USL, prior, metric = metric, target = target)
      pdf_fn <- qc:::make_density_solver(x, LSL, USL, prior,
                                         metric = metric, target = target)
      ref_mass <- integrate(
        Vectorize(pdf_fn),
        bounds[1], bounds[2],
        rel.tol = 1e-4,
        subdivisions = if (metric == "Cpk") 1000 else 200
      )$value
      solver_mass <- S(bounds[1]) - S(bounds[2])

      expect_equal(
        solver_mass,
        ref_mass,
        tolerance = 0.02,
        info = sprintf(
          "%s/%s: solver=%.6f ref=%.6f",
          case_name, metric, solver_mass, ref_mass
        )
      )
    }
  }
})


testthat::test_that("semi-conjugate positive interval probabilities match density mass", {
  skip_if_not_installed("BayesTools")

  set.seed(123)
  x <- rnorm(30, mean = 10, sd = 2)
  LSL <- 2
  USL <- 18
  target <- 10
  bounds <- c(0.5, 1.5)

  prior_cases <- list(
    semi_mu = qc:::.bayestools_to_integration_prior(
      "Jeffreys_mu",
      BayesTools::prior("gamma", list(2, 1))
    )$prior,
    semi_sigma = qc:::.bayestools_to_integration_prior(
      BayesTools::prior("normal", list(10, 5)),
      "Jeffreys_sigma"
    )$prior
  )

  for (case_name in names(prior_cases)) {
    prior <- prior_cases[[case_name]]
    for (metric in c("Cpu", "Cpl", "Cpk")) {
      pdf_fn <- qc:::make_density_solver(x, LSL, USL, prior,
                                         metric = metric, target = target)
      ref_mass <- integrate(
        Vectorize(pdf_fn),
        bounds[1], bounds[2],
        rel.tol = 1e-4,
        subdivisions = if (metric == "Cpk") 1000 else 200
      )$value
      interval_prob <- qc:::compute_cpk_prob_integration(
        x, LSL, USL, bounds, prior,
        metric = metric, target = target
      )

      expect_equal(
        interval_prob,
        ref_mass,
        tolerance = 0.02,
        info = sprintf(
          "%s/%s: interval=%.6f ref=%.6f",
          case_name, metric, interval_prob, ref_mass
        )
      )
    }
  }
})


testthat::test_that("integration method print and summary work", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(x, LSL = 2, target = 10, USL = 18, method = "integration")

  # print should work
  expect_output(print(fit), "Bayesian Process Capability")

  # summary should work
  ss <- summary(fit)
  expect_s3_class(ss, "bpc_summary")
  expect_true(!is.null(ss$summary))
  expect_true(!is.null(ss$interval_summary))

  # print summary should work
  expect_output(print(ss), "Bayesian Process Capability")
})


testthat::test_that("integration summary respects ci.level", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(x, LSL = 2, target = 10, USL = 18, method = "integration")

  summary_95 <- summary(fit, ci.level = 0.95)$summary
  summary_50 <- summary(fit, ci.level = 0.50)$summary

  expect_false(isTRUE(all.equal(summary_95$lower, summary_50$lower)))
  expect_false(isTRUE(all.equal(summary_95$upper, summary_50$upper)))
  expect_true(all(summary_50$lower >= summary_95$lower))
  expect_true(all(summary_50$upper <= summary_95$upper))

  expected_95 <- vapply(fit$integration_result$results, function(result) {
    result$stats[c("Q2.5", "Q97.5")]
  }, numeric(2L))

  # Summary() routes through interval extraction on the stored integration
  # densities, so the public summary can differ slightly from the raw cached
  # quantiles due to grid interpolation.
  expect_equal(summary_95$lower, unname(expected_95[1, ]), tolerance = 0.002)
  expect_equal(summary_95$upper, unname(expected_95[2, ]), tolerance = 0.002)
})


testthat::test_that("integration method errors for t-distribution", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  expect_error(
    bpc(x, LSL = 2, target = 10, USL = 18, method = "integration", distribution = "t"),
    "integration method currently only supports"
  )
})


testthat::test_that("extract_samples errors for integration method", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(x, LSL = 2, target = 10, USL = 18, method = "integration")

  expect_error(
    qc:::extract_samples(fit),
    "extract_samples.*is not supported for integration method"
  )
})


testthat::test_that("integration preserves negative Cpu and Cpk support with data", {

  set.seed(1)
  x <- rnorm(100, mean = 20, sd = 1)
  LSL <- 0
  USL <- 10
  target <- 5

  fit <- bpc(x, LSL = LSL, target = target, USL = USL,
             method = "integration", prior = "Jeffreys")
  prior <- qc:::.bayestools_to_integration_prior("Jeffreys_mu", "Jeffreys_sigma")$prior

  cpu_mom <- qc:::compute_metric_moments(x, LSL, USL, prior, metric = "Cpu", target = target)
  cpk_mom <- qc:::compute_metric_moments(x, LSL, USL, prior, metric = "Cpk", target = target)

  expect_lt(unname(fit$metrics$Cpu["Mean"]), 0)
  expect_lt(unname(fit$metrics$Cpk["Mean"]), 0)
  expect_equal(unname(fit$metrics$Cpu["Mean"]), cpu_mom$mean, tolerance = 0.05)
  expect_equal(unname(fit$metrics$Cpk["Mean"]), cpk_mom$mean, tolerance = 0.05)

  prob_cpu_nonpos <- qc:::compute_cpk_prob_integration(
    x, LSL, USL, c(-Inf, 0), prior,
    metric = "Cpu", target = target,
    cached_state = fit$integration_result$cached_state
  )
  prob_cpk_nonpos <- qc:::compute_cpk_prob_integration(
    x, LSL, USL, c(-Inf, 0), prior,
    metric = "Cpk", target = target,
    cached_state = fit$integration_result$cached_state
  )

  expect_gt(prob_cpu_nonpos, 0.99)
  expect_gt(prob_cpk_nonpos, 0.99)
})


testthat::test_that("integration preserves negative Cpl support with data", {

  set.seed(2)
  x <- rnorm(100, mean = -10, sd = 1)
  LSL <- 0
  USL <- 10
  target <- 5

  fit <- bpc(x, LSL = LSL, target = target, USL = USL,
             method = "integration", prior = "Jeffreys")
  prior <- qc:::.bayestools_to_integration_prior("Jeffreys_mu", "Jeffreys_sigma")$prior

  cpl_mom <- qc:::compute_metric_moments(x, LSL, USL, prior, metric = "Cpl", target = target)
  cpk_mom <- qc:::compute_metric_moments(x, LSL, USL, prior, metric = "Cpk", target = target)

  expect_lt(unname(fit$metrics$Cpl["Mean"]), 0)
  expect_lt(unname(fit$metrics$Cpk["Mean"]), 0)
  expect_equal(unname(fit$metrics$Cpl["Mean"]), cpl_mom$mean, tolerance = 0.05)
  expect_equal(unname(fit$metrics$Cpk["Mean"]), cpk_mom$mean, tolerance = 0.05)

  prob_cpl_nonpos <- qc:::compute_cpk_prob_integration(
    x, LSL, USL, c(-Inf, 0), prior,
    metric = "Cpl", target = target,
    cached_state = fit$integration_result$cached_state
  )

  expect_gt(prob_cpl_nonpos, 0.99)
})


testthat::test_that("integration tail probabilities keep open upper bounds", {

  x <- rep(50, 30) + seq(-0.05, 0.05, length.out = 30)
  LSL <- 48
  USL <- 52
  target <- 50

  fit <- bpc(x, LSL = LSL, target = target, USL = USL,
             method = "integration", prior = "Jeffreys")
  prior <- qc:::.bayestools_to_integration_prior("Jeffreys_mu", "Jeffreys_sigma")$prior
  cached_state <- fit$integration_result$cached_state

  prob_gt_10 <- qc:::compute_cpk_prob_integration(
    x, LSL, USL, c(10, Inf), prior,
    metric = "Cp", target = target,
    cached_state = cached_state
  )
  S_cp <- qc:::make_solver(
    x, LSL, USL, prior,
    metric = "Cp", target = target,
    cached_state = cached_state
  )

  expect_gt(unname(fit$coefficients["Cp"]), 20)
  expect_gt(prob_gt_10, 0.99)
  expect_equal(prob_gt_10, S_cp(10), tolerance = 1e-6)
})


testthat::test_that("integration method with custom priors", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  # Test with custom normal prior on mu and gamma on sigma
  fit <- bpc(x, LSL = 2, target = 10, USL = 18,
             method = "integration",
             prior = prior_independent(
               mu = prior("normal", list(10, 5)),
               sigma = prior("gamma", list(2, 1))
             ))

  # Should still return valid results
  expect_s3_class(fit, "bpc")
  expect_equal(fit$method, "integration")
  expect_true(all(fit$coefficients > 0))
  expect_true(all(is.finite(fit$coefficients)))

  # Non-conjugate so should be slower but still work
  expect_false(fit$integration_result$is_conjugate)
})


testthat::test_that("integration summary recomputes with new specification limits", {

  set.seed(42)
  x <- rnorm(50, mean = 10, sd = 2)

  # Fit with original limits
  fit <- bpc(x, LSL = 4, target = 10, USL = 16,
             method = "integration", prior = "Jeffreys")
  original_coef <- fit$coefficients

  # Get summary with new limits (wider tolerance)
  ss_new <- summary(fit, LSL = 2, target = 10, USL = 18)
  fit_new <- bpc(x, LSL = 2, target = 10, USL = 18,
                 method = "integration", prior = "Jeffreys")
  ss_refit <- summary(fit_new)

  # New limits are wider, so capability indices should be higher
  expect_true(ss_new$summary$mean[ss_new$summary$metric == "Cp"] >
              original_coef["Cp"])
  expect_true(ss_new$summary$mean[ss_new$summary$metric == "Cpk"] >
              original_coef["Cpk"])

  expect_equal(attr(ss_new$metrics, "LSL"), 2)
  expect_equal(attr(ss_new$metrics, "target"), 10)
  expect_equal(attr(ss_new$metrics, "USL"), 18)

  for (metric in names(ss_new$metrics)) {
    summary_mean <- ss_new$summary$mean[ss_new$summary$metric == metric]
    expect_equal(summary_mean, unname(ss_new$metrics[[metric]]["Mean"]), tolerance = 1e-8,
                 info = sprintf("%s: summary=%f, metrics=%f",
                               metric, summary_mean, ss_new$metrics[[metric]]["Mean"]))
  }

  expect_equal(ss_new$interval_summary, ss_refit$interval_summary, tolerance = 1e-6)

  # Get summary with original limits (should match original coefficients)
  ss_orig <- summary(fit, LSL = 4, target = 10, USL = 16)

  for (metric in names(original_coef)) {
    orig_mean <- ss_orig$summary$mean[ss_orig$summary$metric == metric]
    expect_equal(orig_mean, unname(original_coef[metric]), tolerance = 0.001,
                 info = sprintf("%s: summary=%f, original=%f",
                               metric, orig_mean, original_coef[metric]))
  }

  # Summary with narrower limits should give lower capability
  ss_narrow <- summary(fit, LSL = 6, target = 10, USL = 14)
  expect_true(ss_narrow$summary$mean[ss_narrow$summary$metric == "Cp"] <
              original_coef["Cp"])
})


testthat::test_that("plot_density works for integration method", {

  set.seed(42)
  x <- rnorm(50, 10, 2)

  fit <- bpc(x, LSL = 4, target = 10, USL = 16, method = "integration")

  # Basic plot should work
  p <- plot_density(fit)
  expect_s3_class(p, "ggplot")

  # Plot with options should work
  p2 <- plot_density(fit, point_estimate = "mean", ci = "central")
  expect_s3_class(p2, "ggplot")

  # Summary plot should work
  ss <- summary(fit)
  p3 <- plot_density(ss)
  expect_s3_class(p3, "ggplot")

  # Summary with new limits should also work for plotting
  ss_new <- summary(fit, LSL = 2, target = 10, USL = 18)
  p4 <- plot_density(ss_new)
  expect_s3_class(p4, "ggplot")
})


testthat::test_that("plot_density does not recycle finite point labels onto infinite integration facets", {

  metrics <- c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm")
  fit <- bpc(rep(5, 10), LSL = 0, target = 3, USL = 10,
             method = "integration", prior = "Jeffreys")

  p <- plot_density(fit, what = metrics, point_estimate = "mean", ci = "none")
  layers <- p[["layers"]]
  text_idx <- which(vapply(layers, function(layer) inherits(layer[["geom"]], "GeomRichText"), logical(1)))

  expect_length(text_idx, 1)

  text_data <- layers[[text_idx]][["data"]]
  expected_cpc <- (10 - 0) / ((2 * 3) * sqrt(pi / 2) * abs(5 - 3))
  expected_cpm <- min(10 - 3, 3 - 0) / (3 * abs(5 - 3))

  expect_equal(levels(text_data$metric), metrics)
  expect_equal(as.character(text_data$metric), c("Cpc", "Cpm"))
  expect_equal(
    text_data$labels,
    c(sprintf("Mean = %.3f", expected_cpc), sprintf("Mean = %.3f", expected_cpm))
  )
})


testthat::test_that("plot_density handles sample-backed prior-only integration results", {

  fit <- bpc(
    NULL,
    LSL = 2, target = 10, USL = 18,
    method = "integration",
    prior = prior_independent(
      mu = prior("normal", list(10, 5)),
      sigma = prior("gamma", list(2, 1))
    ),
    sample_priors = TRUE
  )

  df <- extract_density_data(fit, what = c("Cp", "Cpk"))
  expect_s3_class(df, "tbl_df")
  expect_equal(names(df), c("x", "density", "metric"))
  expect_equal(unique(as.character(df$metric)), c("Cp", "Cpk"))
  expect_true(all(is.finite(df$x)))
  expect_true(all(is.finite(df$density)))

  p <- plot_density(fit, what = c("Cp", "Cpk"))
  expect_s3_class(p, "ggplot")
  build_fit <- try(ggplot2::ggplot_build(p), silent = TRUE)
  expect_false(inherits(build_fit, "try-error"))

  ss <- summary(fit)
  p_summary <- plot_density(ss, what = "Cp")
  expect_s3_class(p_summary, "ggplot")
  build_summary <- suppressMessages(try(ggplot2::ggplot_build(p_summary), silent = TRUE))
  expect_false(inherits(build_summary, "try-error"))
})


testthat::test_that("extract_density_data works for integration method", {

  set.seed(42)
  x <- rnorm(50, 10, 2)

  fit <- bpc(x, LSL = 4, target = 10, USL = 16, method = "integration")

  # Extract from bpc object
  df <- extract_density_data(fit)
  expect_s3_class(df, "tbl_df")
  expect_equal(names(df), c("x", "density", "metric"))
  expect_equal(levels(df$metric), c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm"))
  expect_true(all(df$density >= 0))

  # Extract single metric
  df_cpk <- extract_density_data(fit, what = "Cpk")
  expect_equal(unique(as.character(df_cpk$metric)), "Cpk")

  # Extract from summary
  ss <- summary(fit)
  df_ss <- extract_density_data(ss)
  expect_s3_class(df_ss, "tbl_df")
  expect_equal(names(df_ss), c("x", "density", "metric"))
})


testthat::test_that("extract_point_estimates works for both methods", {

  set.seed(42)
  x <- rnorm(50, 10, 2)

  fit_int <- bpc(x, LSL = 4, target = 10, USL = 16, method = "integration")

  # Extract density data first
  dfDensity <- extract_density_data(fit_int, what = c("Cp", "Cpk"))

  # Extract point estimates
  dfPoints <- qc:::extract_point_estimates(
    obj = fit_int$integration_result,
    what = c("Cp", "Cpk"),
    point_estimate = "mean",
    dfDensity = dfDensity
  )

  expect_s3_class(dfPoints, "tbl_df")
  expect_equal(names(dfPoints), c("x", "y", "metric"))
  expect_equal(nrow(dfPoints), 2)

  # Mode should also work
  dfMode <- qc:::extract_point_estimates(
    obj = fit_int$integration_result,
    what = c("Cp"),
    point_estimate = "mode",
    dfDensity = dfDensity
  )
  expect_equal(nrow(dfMode), 1)

  dfModeFromEntries <- qc:::extract_point_estimates(
    obj = NULL,
    what = c("Cp"),
    point_estimate = "mode",
    dfDensity = dfDensity,
    stats_list = list(Cp = fit_int$integration_result$results$Cp)
  )
  expect_equal(nrow(dfModeFromEntries), 1)
  expect_equal(dfMode$x, dfModeFromEntries$x)
})


testthat::test_that("extract_point_estimates preserves metric levels when infinite integration estimates are skipped", {

  metrics <- c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm")
  fit <- bpc(rep(5, 10), LSL = 0, target = 3, USL = 10,
             method = "integration", prior = "Jeffreys")
  dfDensity <- extract_density_data(fit, what = metrics)

  dfPoints <- qc:::extract_point_estimates(
    obj = fit$integration_result,
    what = metrics,
    point_estimate = "mean",
    dfDensity = dfDensity
  )

  expect_equal(levels(dfPoints$metric), metrics)
  expect_equal(as.character(dfPoints$metric), c("Cpc", "Cpm"))
})


testthat::test_that("extract_ci_data works for both methods", {

  set.seed(42)
  x <- rnorm(50, 10, 2)

  fit_int <- bpc(x, LSL = 4, target = 10, USL = 16, method = "integration")

  # Extract density data first
  dfDensity <- extract_density_data(fit_int, what = c("Cp", "Cpk"))

  # Extract CI data
  ci_data <- qc:::extract_ci_data(
    obj = fit_int$integration_result,
    what = c("Cp", "Cpk"),
    ci = "HPD",
    ci_level = 0.95,
    dfDensity = dfDensity
  )

  expect_type(ci_data, "list")
  expect_true("dfCi" %in% names(ci_data))
  expect_true("dfArea" %in% names(ci_data))
  expect_equal(nrow(ci_data$dfCi), 2)
})


testthat::test_that("extract_ci_data respects ci_level for integration results", {

  set.seed(42)
  x <- rnorm(50, 10, 2)

  fit_int <- bpc(x, LSL = 4, target = 10, USL = 16, method = "integration")
  dfDensity <- extract_density_data(fit_int, what = c("Cp", "Cpk"))

  ci_central_95 <- qc:::extract_ci_data(
    obj = fit_int$integration_result,
    what = c("Cp", "Cpk"),
    ci = "central",
    ci_level = 0.95,
    dfDensity = dfDensity
  )$dfCi

  ci_central_50 <- qc:::extract_ci_data(
    obj = fit_int$integration_result,
    what = c("Cp", "Cpk"),
    ci = "central",
    ci_level = 0.50,
    dfDensity = dfDensity
  )$dfCi

  expect_false(isTRUE(all.equal(ci_central_95$xmin, ci_central_50$xmin)))
  expect_false(isTRUE(all.equal(ci_central_95$xmax, ci_central_50$xmax)))
  expect_true(all(ci_central_50$xmin >= ci_central_95$xmin))
  expect_true(all(ci_central_50$xmax <= ci_central_95$xmax))

  ci_hpd_95 <- qc:::extract_ci_data(
    obj = fit_int$integration_result,
    what = c("Cp", "Cpk"),
    ci = "HPD",
    ci_level = 0.95,
    dfDensity = dfDensity
  )$dfCi

  ci_hpd_50 <- qc:::extract_ci_data(
    obj = fit_int$integration_result,
    what = c("Cp", "Cpk"),
    ci = "HPD",
    ci_level = 0.50,
    dfDensity = dfDensity
  )$dfCi

  width_95 <- ci_hpd_95$xmax - ci_hpd_95$xmin
  width_50 <- ci_hpd_50$xmax - ci_hpd_50$xmin

  expect_false(isTRUE(all.equal(width_95, width_50)))
  expect_true(all(width_50 <= width_95))
})
