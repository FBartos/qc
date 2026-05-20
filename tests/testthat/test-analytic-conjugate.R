testthat::test_that("analytic Cp solver matches numerical integration", {
  set.seed(123)
  prior <- qc:::create_prior_conjugate(mu0 = 0, k0 = 2, alpha0 = 3, beta0 = 2)
  data <- rnorm(30, mean = 5, sd = 1)
  LSL <- 2; USL <- 8

  S_analytic <- qc:::make_solver(data, LSL, USL, prior, metric = "Cp")

  # Reference: numerical integration via the generic solver path
  ss <- qc:::.extract_suff_stats(data, NULL)
  post <- qc:::.nig_posterior(prior, ss$n, ss$x_bar, ss$SS)
  tol <- USL - LSL; beta_n <- post$beta_n; alpha_n <- post$alpha_n
  df_p <- 2 * alpha_n

  S_numerical <- function(c) {
    if (c <= 0) return(1.0)
    constr <- qc:::get_metric_constraints("Cp", c, LSL, USL, NULL)
    s_max <- constr$s_max_fn()
    if (!is.infinite(s_max) && s_max <= 0) return(0.0)
    y_min <- if (is.infinite(s_max)) 0 else (2 * beta_n) / (s_max^2)
    k_n <- post$k_n; mu_n <- post$mu_n
    y_mode <- max(df_p - 2, 1e-6)
    h_max <- stats::dchisq(y_mode, df_p, log = TRUE)
    integrand <- function(y) {
      sigma <- sqrt((2 * beta_n) / y)
      sd_mu <- sigma / sqrt(k_n)
      mb <- constr$mu_b_fn_vec(sigma)
      z_U <- ifelse(is.infinite(mb$upper), Inf, (mb$upper - mu_n) / sd_mu)
      z_L <- ifelse(is.infinite(mb$lower), -Inf, (mb$lower - mu_n) / sd_mu)
      log_prob <- qc:::log_diff_exp(pnorm(z_U, log.p = TRUE), pnorm(z_L, log.p = TRUE))
      exp(log_prob + dchisq(y, df_p, log = TRUE) - h_max)
    }
    res <- integrate(integrand, y_min, Inf)$value
    if (res <= 0) return(0.0)
    exp(h_max + log(res))
  }

  c_vals <- seq(0.3, 2.5, by = 0.1)
  for (c_val in c_vals) {
    expect_equal(S_analytic(c_val), S_numerical(c_val), tolerance = 1e-4,
                 info = sprintf("Cp at c=%.1f", c_val))
  }
})

testthat::test_that("analytic Cpu solver matches numerical integration", {
  set.seed(123)
  prior <- qc:::create_prior_conjugate(mu0 = 5, k0 = 1, alpha0 = 4, beta0 = 3)
  data <- rnorm(30, mean = 5, sd = 1)
  LSL <- 2; USL <- 8

  S_analytic <- qc:::make_solver(data, LSL, USL, prior, metric = "Cpu")

  ss <- qc:::.extract_suff_stats(data, NULL)
  post <- qc:::.nig_posterior(prior, ss$n, ss$x_bar, ss$SS)
  beta_n <- post$beta_n; alpha_n <- post$alpha_n; k_n <- post$k_n; mu_n <- post$mu_n
  df_p <- 2 * alpha_n
  y_mode <- max(df_p - 2, 1e-6)
  h_max <- stats::dchisq(y_mode, df_p, log = TRUE)

  S_numerical <- function(c) {
    constr <- qc:::get_metric_constraints("Cpu", c, LSL, USL, NULL)
    s_max <- constr$s_max_fn()
    if (!is.infinite(s_max) && s_max <= 0) return(0.0)
    y_min <- if (is.infinite(s_max)) 0 else (2 * beta_n) / (s_max^2)
    integrand <- function(y) {
      sigma <- sqrt((2 * beta_n) / y)
      sd_mu <- sigma / sqrt(k_n)
      mb <- constr$mu_b_fn_vec(sigma)
      z_U <- ifelse(is.infinite(mb$upper), Inf, (mb$upper - mu_n) / sd_mu)
      z_L <- ifelse(is.infinite(mb$lower), -Inf, (mb$lower - mu_n) / sd_mu)
      log_prob <- qc:::log_diff_exp(pnorm(z_U, log.p = TRUE), pnorm(z_L, log.p = TRUE))
      exp(log_prob + dchisq(y, df_p, log = TRUE) - h_max)
    }
    res <- integrate(integrand, y_min, Inf)$value
    if (res <= 0) return(0.0)
    exp(h_max + log(res))
  }

  c_vals <- seq(-1.5, 2.5, by = 0.1)
  for (c_val in c_vals) {
    expect_equal(S_analytic(c_val), S_numerical(c_val), tolerance = 1e-4,
                 info = sprintf("Cpu at c=%.1f", c_val))
  }
})

testthat::test_that("analytic Cpl solver matches numerical integration", {
  set.seed(123)
  prior <- qc:::create_prior_conjugate(mu0 = 5, k0 = 1, alpha0 = 4, beta0 = 3)
  data <- rnorm(30, mean = 5, sd = 1)
  LSL <- 2; USL <- 8

  S_analytic <- qc:::make_solver(data, LSL, USL, prior, metric = "Cpl")

  ss <- qc:::.extract_suff_stats(data, NULL)
  post <- qc:::.nig_posterior(prior, ss$n, ss$x_bar, ss$SS)
  beta_n <- post$beta_n; alpha_n <- post$alpha_n; k_n <- post$k_n; mu_n <- post$mu_n
  df_p <- 2 * alpha_n
  y_mode <- max(df_p - 2, 1e-6)
  h_max <- stats::dchisq(y_mode, df_p, log = TRUE)

  S_numerical <- function(c) {
    constr <- qc:::get_metric_constraints("Cpl", c, LSL, USL, NULL)
    s_max <- constr$s_max_fn()
    if (!is.infinite(s_max) && s_max <= 0) return(0.0)
    y_min <- if (is.infinite(s_max)) 0 else (2 * beta_n) / (s_max^2)
    integrand <- function(y) {
      sigma <- sqrt((2 * beta_n) / y)
      sd_mu <- sigma / sqrt(k_n)
      mb <- constr$mu_b_fn_vec(sigma)
      z_U <- ifelse(is.infinite(mb$upper), Inf, (mb$upper - mu_n) / sd_mu)
      z_L <- ifelse(is.infinite(mb$lower), -Inf, (mb$lower - mu_n) / sd_mu)
      log_prob <- qc:::log_diff_exp(pnorm(z_U, log.p = TRUE), pnorm(z_L, log.p = TRUE))
      exp(log_prob + dchisq(y, df_p, log = TRUE) - h_max)
    }
    res <- integrate(integrand, y_min, Inf)$value
    if (res <= 0) return(0.0)
    exp(h_max + log(res))
  }

  c_vals <- seq(-1.5, 2.5, by = 0.1)
  for (c_val in c_vals) {
    expect_equal(S_analytic(c_val), S_numerical(c_val), tolerance = 1e-4,
                 info = sprintf("Cpl at c=%.1f", c_val))
  }
})

testthat::test_that("analytic density solver for Cpu/Cpl produces normalized PDFs", {
  set.seed(42)
  data <- rnorm(30, mean = 5, sd = 1)
  LSL <- 2; USL <- 8; target <- 5
  prior <- qc:::create_prior_conjugate()

  for (metric in c("Cpu", "Cpl")) {
    pdf_fn <- qc:::make_density_solver(data, LSL, USL, prior, metric = metric, target = target)
    pdf_fn_vec <- Vectorize(pdf_fn)
    total_mass <- integrate(pdf_fn_vec, 0, 5, subdivisions = 200)$value
    expect_true(
      abs(total_mass - 1) < 0.02,
      info = sprintf("%s: total_mass = %.4f", metric, total_mass)
    )
  }
})

testthat::test_that("analytic density solver agrees with survival function for Cpu/Cpl", {
  set.seed(42)
  data <- rnorm(30, mean = 5, sd = 1)
  LSL <- 2; USL <- 8; target <- 5
  prior <- qc:::create_prior_conjugate()

  for (metric in c("Cpu", "Cpl")) {
    S_fn <- qc:::make_solver(data, LSL, USL, prior, metric = metric, target = target)
    pdf_fn <- qc:::make_density_solver(data, LSL, USL, prior, metric = metric, target = target)
    pdf_fn_vec <- Vectorize(pdf_fn)

    c_vals <- seq(0.5, 2.5, by = 0.5)
    for (c_val in c_vals) {
      S_val <- S_fn(c_val)
      cdf_from_pdf <- integrate(pdf_fn_vec, c_val, 5, subdivisions = 200)$value
      expect_equal(cdf_from_pdf, S_val, tolerance = 0.02,
                   info = sprintf("%s at c=%.1f: S=%.4f, cdf_from_pdf=%.4f",
                                  metric, c_val, S_val, cdf_from_pdf))
    }
  }
})

testthat::test_that("analytic solvers work in analyze_capability_integration", {
  set.seed(42)
  data <- rnorm(30, mean = 5, sd = 1)
  LSL <- 2; USL <- 8; target <- 5
  prior <- qc:::create_prior_conjugate()

  for (metric in c("Cp", "Cpu", "Cpl")) {
    result <- qc:::analyze_capability_integration(data, LSL, USL, prior,
                                                   metric = metric, target = target)
    expect_true(result$stats["Mean"] > 0,
                info = sprintf("%s: mean should be positive", metric))
    expect_true(result$stats["SD"] > 0,
                info = sprintf("%s: sd should be positive", metric))
    expect_true(result$stats["Q2.5"] < result$stats["Q97.5"],
                info = sprintf("%s: Q2.5 should be less than Q97.5", metric))
  }
})
