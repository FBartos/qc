test_that(".extract_alpha_parameter classifies prior families correctly", {
  # Exponential: alpha = 1
  p_exp <- BayesTools::prior("exp", list(1))
  expect_equal(qc:::.extract_alpha_parameter(p_exp), 1)

  # Gamma(shape, rate): alpha = shape
  p_gamma2 <- BayesTools::prior("gamma", list(2, 1))
  expect_equal(qc:::.extract_alpha_parameter(p_gamma2), 2)

  p_gamma05 <- BayesTools::prior("gamma", list(0.5, 1))
  expect_equal(qc:::.extract_alpha_parameter(p_gamma05), 0.5)

  # Inverse-gamma: always alpha = Inf (all negative moments finite)
  p_invgamma <- BayesTools::prior("invgamma", list(2, 1))
  expect_equal(qc:::.extract_alpha_parameter(p_invgamma), Inf)

  # Log-normal: alpha = Inf
  p_lnorm <- BayesTools::prior("lognormal", list(0, 1))
  expect_equal(qc:::.extract_alpha_parameter(p_lnorm), Inf)

  # Cauchy (half): alpha = 1
  p_cauchy <- BayesTools::prior("cauchy", list(0, 2.5), list(0, Inf))
  expect_equal(qc:::.extract_alpha_parameter(p_cauchy), 1)

  # t (half): alpha = 1
  p_t <- BayesTools::prior("t", list(0, 2.5, 5), list(0, Inf))
  expect_equal(qc:::.extract_alpha_parameter(p_t), 1)

  # Normal (half): alpha = Inf
  p_normal <- BayesTools::prior("normal", list(0, 1), list(0, Inf))
  expect_equal(qc:::.extract_alpha_parameter(p_normal), Inf)

  # Uniform(0, b): alpha = 1 (linear near zero)
  p_unif <- BayesTools::prior("uniform", list(0, 10))
  expect_equal(qc:::.extract_alpha_parameter(p_unif), 1)

  # Uniform(a > 0, b): alpha = Inf (bounded away from zero)
  p_unif_pos <- BayesTools::prior("uniform", list(0.1, 10))
  expect_equal(qc:::.extract_alpha_parameter(p_unif_pos), Inf)

  # Truncation with lower > 0 overrides to Inf
  p_exp_trunc <- BayesTools::prior("exp", list(1), list(0.1, Inf))
  expect_equal(qc:::.extract_alpha_parameter(p_exp_trunc), Inf)
})

test_that(".check_moment_divergence returns correct flags for Cp-family", {
  # Cp requires alpha > 1 for mean, alpha > 2 for variance
  # alpha = 1: both diverge
  d1 <- qc:::.check_moment_divergence("Cp", 1, "Exp(1)")
  expect_true(d1$mean_divergent)
  expect_true(d1$sd_divergent)

  # alpha = 1.5: mean finite, variance diverges
  d15 <- qc:::.check_moment_divergence("Cp", 1.5, "Gamma(1.5, 1)")
  expect_false(d15$mean_divergent)
  expect_true(d15$sd_divergent)

  # alpha = 3: both finite
  d3 <- qc:::.check_moment_divergence("Cp", 3, "Gamma(3, 1)")
  expect_false(d3$mean_divergent)
  expect_false(d3$sd_divergent)

  # alpha = Inf: both finite
  dInf <- qc:::.check_moment_divergence("Cp", Inf, "InvGamma(2, 1)")
  expect_false(dInf$mean_divergent)
  expect_false(dInf$sd_divergent)

  # Cpu, Cpl, Cpk behave identical to Cp
  for (m in c("Cpu", "Cpl", "Cpk")) {
    d <- qc:::.check_moment_divergence(m, 1, "Exp(1)")
    expect_true(d$mean_divergent, info = paste(m, "alpha=1 mean"))
    expect_true(d$sd_divergent, info = paste(m, "alpha=1 sd"))
  }
})

test_that(".check_moment_divergence returns correct flags for Cpm-family", {
  # Cpm/Cpc require alpha > 0 for mean, alpha > 1 for variance
  # alpha = 0.5: mean finite, variance diverges
  d05 <- qc:::.check_moment_divergence("Cpm", 0.5, "Gamma(0.5, 1)")
  expect_false(d05$mean_divergent)
  expect_true(d05$sd_divergent)

  # alpha = 1: mean finite, variance diverges
  d1 <- qc:::.check_moment_divergence("Cpm", 1, "Exp(1)")
  expect_false(d1$mean_divergent)
  expect_true(d1$sd_divergent)

  # alpha = 2: both finite
  d2 <- qc:::.check_moment_divergence("Cpm", 2, "Gamma(2, 1)")
  expect_false(d2$mean_divergent)
  expect_false(d2$sd_divergent)

  # alpha = Inf: both finite
  dInf <- qc:::.check_moment_divergence("Cpm", Inf, "InvGamma(2, 1)")
  expect_false(dInf$mean_divergent)
  expect_false(dInf$sd_divergent)
})

test_that("bpc with Exp(1) prior: all Cp-family means are Inf", {
  x <- numeric(2)
  fit <- bpc(x, LSL = -1, USL = 1, target = 0,
             prior_mu    = BayesTools::prior("normal", list(0, 1)),
             prior_sigma = BayesTools::prior("exp", list(1)),
             sample_priors = TRUE, method = "integration")
  s <- summary(fit)

  cp_family <- c("Cp", "Cpu", "Cpl", "Cpk")
  for (m in cp_family) {
    row <- s$summary[s$summary$metric == m, ]
    expect_true(is.infinite(row$mean), info = paste(m, "mean should be Inf"))
    expect_true(is.infinite(row$sd),   info = paste(m, "sd should be Inf"))
    expect_true(is.finite(row$median), info = paste(m, "median should be finite"))
    expect_true(is_analytic(s, m, "mean"))
    expect_true(is_analytic(s, m, "sd"))
  }

  # Cpm/Cpc: mean finite (alpha=1 > 0), sd infinite (alpha=1 <= 1)
  for (m in c("Cpm", "Cpc")) {
    row <- s$summary[s$summary$metric == m, ]
    expect_true(is.finite(row$mean),   info = paste(m, "mean should be finite"))
    expect_true(is.infinite(row$sd),   info = paste(m, "sd should be Inf"))
    expect_false(is_analytic(s, m, "mean"))
    expect_true(is_analytic(s, m, "sd"))
  }
})

test_that("bpc with Gamma(2,1) prior: Cp mean finite, sd diverges", {
  x <- numeric(2)
  fit <- bpc(x, LSL = -1, USL = 1, target = 0,
             prior_mu    = BayesTools::prior("normal", list(0, 1)),
             prior_sigma = BayesTools::prior("gamma", list(2, 1)),
             sample_priors = TRUE, method = "integration")
  s <- summary(fit)

  row_cp <- s$summary[s$summary$metric == "Cp", ]
  expect_true(is.finite(row_cp$mean))
  expect_true(is.infinite(row_cp$sd))
  expect_false(is_analytic(s, "Cp", "mean"))
  expect_true(is_analytic(s, "Cp", "sd"))

  # Cpm: alpha=2 > 1 for variance => both finite
  row_cpm <- s$summary[s$summary$metric == "Cpm", ]
  expect_true(is.finite(row_cpm$mean))
  expect_true(is.finite(row_cpm$sd))
  expect_false(is_analytic(s, "Cpm", "mean"))
  expect_false(is_analytic(s, "Cpm", "sd"))
})

test_that("bpc with InvGamma(2,1) prior: no divergence", {
  x <- numeric(2)
  fit <- bpc(x, LSL = -1, USL = 1, target = 0,
             prior_mu    = BayesTools::prior("normal", list(0, 1)),
             prior_sigma = BayesTools::prior("invgamma", list(2, 1)),
             sample_priors = TRUE, method = "integration")
  s <- summary(fit)

  expect_false(isTRUE(attr(s, "has_divergent_moments")))
  for (m in c("Cp", "Cpu", "Cpl", "Cpk", "Cpm", "Cpc")) {
    expect_false(is_analytic(s, m, "mean"), info = paste(m, "mean"))
    expect_false(is_analytic(s, m, "sd"),   info = paste(m, "sd"))
  }
})

test_that("bpc with Uniform(0, b) prior: treated as alpha=1", {
  x <- numeric(2)
  fit <- bpc(x, LSL = -1, USL = 1, target = 0,
             prior_mu    = BayesTools::prior("normal", list(0, 1)),
             prior_sigma = BayesTools::prior("uniform", list(0, 10)),
             sample_priors = TRUE, method = "integration")
  s <- summary(fit)

  expect_true(is.infinite(s$summary$mean[s$summary$metric == "Cp"]))
  expect_true(is_analytic(s, "Cp", "mean"))
})

test_that("bpc with Half-Cauchy prior: treated as alpha=1", {
  x <- numeric(2)
  fit <- bpc(x, LSL = -1, USL = 1, target = 0,
             prior_mu    = BayesTools::prior("normal", list(0, 1)),
             prior_sigma = BayesTools::prior("cauchy", list(0, 2.5), list(0, Inf)),
             sample_priors = TRUE, method = "integration")
  s <- summary(fit)

  expect_true(is.infinite(s$summary$mean[s$summary$metric == "Cp"]))
  expect_true(is_analytic(s, "Cp", "mean"))
})

test_that("bpc with data: likelihood regularizes, no divergence", {
  set.seed(42)
  x <- rnorm(30, mean = 5, sd = 1)
  fit <- bpc(x, LSL = 2, USL = 8, target = 5,
             prior_mu    = BayesTools::prior("normal", list(5, 10)),
             prior_sigma = BayesTools::prior("exp", list(1)),
             method = "integration")
  s <- summary(fit)

  expect_false(isTRUE(attr(s, "has_divergent_moments")))
  expect_true(is.finite(s$summary$mean[s$summary$metric == "Cp"]))
})

test_that("bpc with truncated Exp(1) at lower=0.1: alpha=Inf, no divergence", {
  x <- numeric(2)
  fit <- bpc(x, LSL = -1, USL = 1, target = 0,
             prior_mu    = BayesTools::prior("normal", list(0, 1)),
             prior_sigma = BayesTools::prior("exp", list(1), list(0.1, Inf)),
             sample_priors = TRUE, method = "integration")
  s <- summary(fit)

  expect_false(isTRUE(attr(s, "has_divergent_moments")))
  expect_true(is.finite(s$summary$mean[s$summary$metric == "Cp"]))
})

test_that("get_analytic_flags returns correct structure", {
  x <- numeric(2)
  fit <- bpc(x, LSL = -1, USL = 1, target = 0,
             prior_mu    = BayesTools::prior("normal", list(0, 1)),
             prior_sigma = BayesTools::prior("exp", list(1)),
             sample_priors = TRUE, method = "integration")
  s <- summary(fit)

  flags <- get_analytic_flags(s)
  expect_s3_class(flags, "data.frame")
  expect_true("metric" %in% names(flags))
  expect_true("mean_divergent" %in% names(flags))
  expect_true("sd_divergent" %in% names(flags))
  expect_true("reason" %in% names(flags))
  expect_equal(nrow(flags), 6L)
})

test_that("is_analytic returns FALSE for non-divergent fit", {
  x <- numeric(2)
  fit <- bpc(x, LSL = -1, USL = 1, target = 0,
             prior_mu    = BayesTools::prior("normal", list(0, 1)),
             prior_sigma = BayesTools::prior("invgamma", list(2, 1)),
             sample_priors = TRUE, method = "integration")
  s <- summary(fit)

  expect_false(is_analytic(s, "Cp", "mean"))
  expect_false(is_analytic(s, "Cp", "sd"))

  flags <- get_analytic_flags(s)
  expect_null(flags)
})
