divergence_metrics <- qc:::.qc_metric_names()

make_divergence_map <- function(alpha, label) {
  stats::setNames(
    lapply(divergence_metrics, function(metric) {
      qc:::.check_moment_divergence(metric, alpha, label)
    }),
    divergence_metrics
  )
}

make_divergence_distributions <- function(alpha, label) {
  divergence_map <- make_divergence_map(alpha, label)
  context <- qc:::.new_qc_metric_context(divergence = divergence_map)

  stats::setNames(
    lapply(divergence_metrics, function(metric) {
      info <- divergence_map[[metric]]
      qc:::.new_qc_metric_distribution(
        metric = metric,
        stats = c(
          Mean = if (info$mean_divergent) Inf else 1.25,
          Median = 1.00,
          SD = if (info$sd_divergent) Inf else 0.20,
          Q2.5 = 0.75,
          Q97.5 = 1.50
        ),
        context = context
      )
    }),
    divergence_metrics
  )
}

make_divergence_summary <- function(alpha, label) {
  divergence_map <- make_divergence_map(alpha, label)
  distributions <- make_divergence_distributions(alpha, label)
  metrics <- stats::setNames(lapply(distributions, function(dist) dist$stats), divergence_metrics)

  summary_tbl <- data.frame(
    metric = divergence_metrics,
    mean = vapply(distributions, function(dist) unname(dist$stats["Mean"]), numeric(1)),
    median = vapply(distributions, function(dist) unname(dist$stats["Median"]), numeric(1)),
    sd = vapply(distributions, function(dist) unname(dist$stats["SD"]), numeric(1)),
    lower = vapply(distributions, function(dist) unname(dist$stats["Q2.5"]), numeric(1)),
    upper = vapply(distributions, function(dist) unname(dist$stats["Q97.5"]), numeric(1)),
    stringsAsFactors = FALSE
  )

  qc:::.new_bpc_summary(
    call = quote(summary(bpc(...))),
    metrics = qc:::.new_capability_metrics(
      metrics = metrics,
      LSL = -1,
      USL = 1,
      target = 0,
      method = "integration",
      distributions = distributions,
      divergence = divergence_map
    ),
    summary = summary_tbl,
    interval_summary = data.frame(metric = divergence_metrics, stringsAsFactors = FALSE),
    divergence_diagnostics = qc:::.qc_metric_distributions_divergence_table(distributions)
  )
}

test_that(".extract_alpha_parameter classifies prior families correctly", {
  skip_if_not_installed("BayesTools")

  p_exp <- BayesTools::prior("exp", list(1))
  expect_equal(qc:::.extract_alpha_parameter(p_exp), 1)

  p_gamma2 <- BayesTools::prior("gamma", list(2, 1))
  expect_equal(qc:::.extract_alpha_parameter(p_gamma2), 2)

  p_gamma05 <- BayesTools::prior("gamma", list(0.5, 1))
  expect_equal(qc:::.extract_alpha_parameter(p_gamma05), 0.5)

  p_invgamma <- BayesTools::prior("invgamma", list(2, 1))
  expect_equal(qc:::.extract_alpha_parameter(p_invgamma), Inf)

  p_lnorm <- BayesTools::prior("lognormal", list(0, 1))
  expect_equal(qc:::.extract_alpha_parameter(p_lnorm), Inf)

  p_cauchy <- BayesTools::prior("cauchy", list(0, 2.5), list(0, Inf))
  expect_equal(qc:::.extract_alpha_parameter(p_cauchy), 1)

  p_t <- BayesTools::prior("t", list(0, 2.5, 5), list(0, Inf))
  expect_equal(qc:::.extract_alpha_parameter(p_t), 1)

  p_normal <- BayesTools::prior("normal", list(0, 1), list(0, Inf))
  expect_equal(qc:::.extract_alpha_parameter(p_normal), 1)

  p_unif <- BayesTools::prior("uniform", list(0, 10))
  expect_equal(qc:::.extract_alpha_parameter(p_unif), 1)

  p_unif_pos <- BayesTools::prior("uniform", list(0.1, 10))
  expect_equal(qc:::.extract_alpha_parameter(p_unif_pos), Inf)

  p_exp_trunc <- BayesTools::prior("exp", list(1), list(0.1, Inf))
  expect_equal(qc:::.extract_alpha_parameter(p_exp_trunc), Inf)

  p_beta <- BayesTools::prior("beta", list(alpha = 2, beta = 4))
  expect_equal(qc:::.extract_alpha_parameter(p_beta), 2)

  p_normal_trunc_away <- BayesTools::prior("normal", list(0, 1), list(0.1, Inf))
  expect_equal(qc:::.extract_alpha_parameter(p_normal_trunc_away), Inf)
})

test_that(".check_moment_divergence returns correct flags for Cp-family", {
  d1 <- qc:::.check_moment_divergence("Cp", 1, "Exp(1)")
  expect_true(d1$mean_divergent)
  expect_true(d1$sd_divergent)

  d15 <- qc:::.check_moment_divergence("Cp", 1.5, "Gamma(1.5, 1)")
  expect_false(d15$mean_divergent)
  expect_true(d15$sd_divergent)

  d3 <- qc:::.check_moment_divergence("Cp", 3, "Gamma(3, 1)")
  expect_false(d3$mean_divergent)
  expect_false(d3$sd_divergent)

  d_inf <- qc:::.check_moment_divergence("Cp", Inf, "InvGamma(2, 1)")
  expect_false(d_inf$mean_divergent)
  expect_false(d_inf$sd_divergent)

  for (metric in c("Cpu", "Cpl", "Cpk")) {
    d <- qc:::.check_moment_divergence(metric, 1, "Exp(1)")
    expect_true(d$mean_divergent, info = paste(metric, "alpha=1 mean"))
    expect_true(d$sd_divergent, info = paste(metric, "alpha=1 sd"))
  }
})

test_that(".check_moment_divergence returns correct flags for Cpm-family", {
  d05 <- qc:::.check_moment_divergence("Cpm", 0.5, "Gamma(0.5, 1)")
  expect_false(d05$mean_divergent)
  expect_true(d05$sd_divergent)

  d1 <- qc:::.check_moment_divergence("Cpm", 1, "Exp(1)")
  expect_false(d1$mean_divergent)
  expect_true(d1$sd_divergent)

  d2 <- qc:::.check_moment_divergence("Cpm", 2, "Gamma(2, 1)")
  expect_false(d2$mean_divergent)
  expect_false(d2$sd_divergent)

  d_inf <- qc:::.check_moment_divergence("Cpm", Inf, "InvGamma(2, 1)")
  expect_false(d_inf$mean_divergent)
  expect_false(d_inf$sd_divergent)
})

test_that("divergence tables encode the expected alpha=1 prior-only behavior", {
  diagnostics <- qc:::.qc_metric_distributions_divergence_table(
    make_divergence_distributions(1, "Exp(1)")
  )

  expect_s3_class(diagnostics, "data.frame")
  expect_equal(
    diagnostics$metric,
    c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm")
  )
  expect_true(all(diagnostics$sd_divergent))

  cp_rows <- diagnostics$metric %in% c("Cp", "Cpu", "Cpl", "Cpk")
  expect_true(all(diagnostics$mean_divergent[cp_rows]))
  expect_true(all(!diagnostics$mean_divergent[!cp_rows]))
})

test_that("synthetic bpc summaries expose analytic flags without a full integration fit", {
  s <- make_divergence_summary(1, "Exp(1)")

  for (metric in c("Cp", "Cpu", "Cpl", "Cpk")) {
    row <- s$summary[s$summary$metric == metric, , drop = FALSE]
    expect_true(is.infinite(row$mean), info = paste(metric, "mean should be Inf"))
    expect_true(is.infinite(row$sd), info = paste(metric, "sd should be Inf"))
    expect_true(is.finite(row$median), info = paste(metric, "median should be finite"))
    expect_true(is_analytic(s, metric, "mean"))
    expect_true(is_analytic(s, metric, "sd"))
  }

  for (metric in c("Cpm", "Cpc")) {
    row <- s$summary[s$summary$metric == metric, , drop = FALSE]
    expect_true(is.finite(row$mean), info = paste(metric, "mean should be finite"))
    expect_true(is.infinite(row$sd), info = paste(metric, "sd should be Inf"))
    expect_false(is_analytic(s, metric, "mean"))
    expect_true(is_analytic(s, metric, "sd"))
  }

  flags <- get_analytic_flags(s)
  expect_s3_class(flags, "data.frame")
  expect_true(all(c("metric", "mean_divergent", "sd_divergent", "reason") %in% names(flags)))
  expect_equal(nrow(flags), 6L)
})

test_that("non-divergent summaries return no analytic flags", {
  s <- make_divergence_summary(Inf, "InvGamma(2, 1)")

  expect_false(isTRUE(attr(s, "has_divergent_moments")))

  for (metric in divergence_metrics) {
    expect_false(is_analytic(s, metric, "mean"), info = paste(metric, "mean"))
    expect_false(is_analytic(s, metric, "sd"), info = paste(metric, "sd"))
  }

  expect_null(get_analytic_flags(s))
})

test_that("prior-only UIP integration does not flag proper conjugate sigma priors as divergent", {
  x <- c(-1.5, -0.5, 0.25, 0.75, 1.5)
  uip <- qc:::create_prior_unit_information(x)

  fit <- qc::bpc(
    x,
    LSL = -3,
    USL = 3,
    target = 0,
    prior = uip,
    sample_priors = TRUE,
    method = "integration"
  )
  s <- summary(fit)

  expect_false(isTRUE(attr(s, "has_divergent_moments")))
  expect_null(get_analytic_flags(s))
  expect_true(all(is.finite(s$summary$mean)))
  expect_true(all(is.finite(s$summary$sd)))
})
