metric_registry_stats <- function(mean = 1, median = mean, sd = 0.1) {
  c(
    Mean = mean,
    Median = median,
    SD = sd,
    Q2.5 = mean - 0.2,
    Q97.5 = mean + 0.2,
    HDI_Lo = mean - 0.15,
    HDI_Hi = mean + 0.15
  )
}

testthat::test_that("metric distribution helpers validate requested metrics and coercion inputs", {
  testthat::expect_equal(
    qc:::.qc_metric_names(),
    c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm")
  )
  testthat::expect_equal(qc:::.qc_match_metrics(c("Cpk", "Cp")), c("Cpk", "Cp"))

  testthat::expect_error(
    qc:::.qc_match_metrics("not-a-metric"),
    regexp = "should be one of"
  )
  testthat::expect_error(
    qc:::.coerce_qc_metric_distributions(),
    regexp = "Provide either"
  )
  testthat::expect_error(
    qc:::.coerce_qc_metric_distributions(
      what = c("Cp", "Cpk"),
      stats_list = list(Cp = metric_registry_stats())
    ),
    regexp = "Missing distribution input for metric 'Cpk'"
  )
})

testthat::test_that("metric distribution coercion wraps numeric stats and preserves prepared distributions", {
  context <- qc:::.new_qc_metric_context(LSL = 0, USL = 10, target = 5, sigma_level = 3)
  ready <- qc:::.new_qc_metric_distribution(
    metric = "Cpk",
    stats = metric_registry_stats(mean = 1.25),
    context = context
  )

  distributions <- qc:::.coerce_qc_metric_distributions(
    what = c("Cp", "Cpk"),
    stats_list = list(
      Cp = metric_registry_stats(mean = 1.1),
      Cpk = ready
    ),
    context = context
  )

  testthat::expect_true(all(vapply(distributions, inherits, logical(1), "qc_metric_distribution")))
  testthat::expect_identical(distributions$Cpk, ready)
  testthat::expect_identical(distributions$Cp$metric, "Cp")
  testthat::expect_equal(unname(distributions$Cp$stats["Mean"]), 1.1)
  testthat::expect_identical(distributions$Cp$context, context)
})

testthat::test_that("integration-backed capability metrics rebuild distributions from attached results metadata", {
  metric_names <- qc:::.qc_metric_names()
  results <- stats::setNames(lapply(seq_along(metric_names), function(i) {
    list(stats = metric_registry_stats(mean = 1 + i / 10))
  }), metric_names)
  prior <- qc::create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 2, beta0 = 1)
  cached_state <- list(n = 8L, x_bar = 0.25, sse = 3.5)
  divergence <- stats::setNames(rep(list(list(
    mean_divergent = FALSE,
    sd_divergent = FALSE,
    alpha = Inf,
    reason = NULL
  )), length(metric_names)), metric_names)

  metrics <- qc:::.new_capability_metrics(
    metrics = stats::setNames(lapply(metric_names, function(metric) results[[metric]]$stats), metric_names),
    LSL = 0,
    USL = 10,
    target = 5,
    sigma = 3,
    method = "integration",
    results = results,
    prior = prior,
    cached_state = cached_state,
    divergence = divergence
  )

  distributions <- qc:::.as_qc_metric_distributions.capability_metrics(
    metrics,
    what = c("Cpk", "Cp")
  )

  testthat::expect_named(distributions, c("Cpk", "Cp"))
  testthat::expect_identical(distributions$Cpk$entry, results$Cpk)
  testthat::expect_identical(distributions$Cpk$stats, results$Cpk$stats)
  testthat::expect_identical(distributions$Cpk$context$prior, prior)
  testthat::expect_identical(distributions$Cpk$context$cached_state, cached_state)
  testthat::expect_equal(distributions$Cpk$context$LSL, 0)
  testthat::expect_equal(distributions$Cpk$context$USL, 10)
  testthat::expect_equal(distributions$Cpk$context$target, 5)
  testthat::expect_equal(distributions$Cpk$context$sigma_level, 3)
  testthat::expect_identical(distributions$Cpk$context$divergence, divergence)

  metrics_without_results <- qc:::.new_capability_metrics(
    metrics = stats::setNames(lapply(metric_names, function(metric) results[[metric]]$stats), metric_names),
    LSL = 0,
    USL = 10,
    target = 5,
    method = "integration"
  )

  testthat::expect_error(
    qc:::.as_qc_metric_distributions.capability_metrics(
      metrics_without_results,
      what = "Cp"
    ),
    regexp = "require attached integration results"
  )
})

testthat::test_that("distribution registry centralizes supported methods and parameters", {
  testthat::expect_equal(qc:::.qc_distribution_names(), c("normal", "t"))
  testthat::expect_equal(qc:::.qc_distribution_names(method = "pc"), c("normal", "t"))
  testthat::expect_equal(qc:::.qc_distribution_names(method = "mcmc"), c("normal", "t"))
  testthat::expect_equal(qc:::.qc_distribution_names(method = "integration"), "normal")

  testthat::expect_true(qc:::.qc_distribution_supports_summary_statistics("normal"))
  testthat::expect_false(qc:::.qc_distribution_supports_summary_statistics("t"))

  testthat::expect_equal(
    qc:::.qc_distribution_parameter_names("normal", type = "sample"),
    c("mu", "sigma")
  )
  testthat::expect_equal(
    qc:::.qc_distribution_parameter_names("t", type = "sample"),
    c("mu", "scale", "nu")
  )
  testthat::expect_equal(
    qc:::.qc_distribution_parameter_names("t", type = "prior"),
    c("mu", "sigma", "nu")
  )
  testthat::expect_equal(
    qc:::.qc_distribution_prior_defaults("normal"),
    list(mu = "Jeffreys_mu", sigma = "Jeffreys_sigma")
  )
  testthat::expect_equal(
    qc:::.qc_distribution_prior_defaults("t"),
    list(mu = "Jeffreys_mu", sigma = "Jeffreys_sigma", nu = "uniform_nu")
  )
  testthat::expect_true(qc:::.qc_distribution_supports_parameter("t", "nu", type = "prior"))
  testthat::expect_false(qc:::.qc_distribution_supports_parameter("normal", "nu", type = "prior"))

  testthat::expect_equal(qc:::.qc_distribution_stan_model("normal"), "normal")
  testthat::expect_equal(qc:::.qc_distribution_stan_model("t"), "t")
  testthat::expect_equal(
    qc:::.qc_distribution_integration_backend("normal"),
    ".bpc_fit_integration_normal"
  )
  testthat::expect_true(qc:::.qc_distribution_has_registered_integration_backend("normal"))
})

testthat::test_that("distribution registry normalizes independent prior maps from defaults and parameters", {
  t_defaults <- qc:::.qc_distribution_prior_map("t")
  testthat::expect_equal(
    t_defaults,
    list(mu = "Jeffreys_mu", sigma = "Jeffreys_sigma", nu = "uniform_nu")
  )

  nu_prior <- prior("exp", list(rate = 1 / 30), list(lower = 2, upper = Inf))
  t_map <- qc:::.qc_distribution_prior_map(
    "t",
    parameters = list(nu = nu_prior)
  )
  testthat::expect_identical(t_map$nu, nu_prior)

  theta_prior <- prior("normal", list(0, 1))
  testthat::local_mocked_bindings(
    .qc_distribution_parameter_names = function(distribution, type = c("sample", "prior")) {
      type <- match.arg(type)
      if (identical(type, "prior")) {
        return(c("mu", "sigma", "theta"))
      }
      c("mu", "sigma")
    },
    .qc_distribution_prior_defaults = function(distribution) {
      list(mu = "Jeffreys_mu", sigma = "Jeffreys_sigma", theta = NULL)
    },
    .package = "qc"
  )
  extended_map <- qc:::.qc_distribution_prior_map(
    "normal",
    parameters = list(theta = theta_prior)
  )
  testthat::expect_identical(extended_map$theta, theta_prior)

  testthat::expect_error(
    qc:::.qc_distribution_prior_map(
      "normal",
      parameters = list(shape = theta_prior)
    ),
    regexp = "not a prior parameter"
  )
})

testthat::test_that("distribution registry exposes normalized adapters and contracts", {
  adapter <- qc:::.qc_distribution_adapter("normal")

  testthat::expect_s3_class(adapter, "qc_distribution_adapter")
  testthat::expect_identical(adapter$name, "normal")
  testthat::expect_identical(adapter$sample_class, "normal")
  testthat::expect_equal(adapter$parameters$sample, c("mu", "sigma"))
  testthat::expect_equal(adapter$parameters$prior, c("mu", "sigma"))
  testthat::expect_equal(
    adapter$prior_defaults,
    list(mu = "Jeffreys_mu", sigma = "Jeffreys_sigma")
  )
  testthat::expect_equal(adapter$mcmc$stan_model_name, "normal")
  testthat::expect_true(adapter$mcmc$available)
  testthat::expect_equal(adapter$integration$backend_name, ".bpc_fit_integration_normal")
  testthat::expect_true(adapter$integration$available)
  testthat::expect_true(adapter$s3$percentiles$available)
  testthat::expect_true(adapter$s3$posterior_predictives$available)
  testthat::expect_true(adapter$s3$E_abs_dev$available)
  testthat::expect_true(adapter$s3$pc_fit$available)

  contract <- qc:::.qc_distribution_contract("normal")
  testthat::expect_s3_class(contract, "data.frame")
  testthat::expect_true(all(c(
    "capability", "kind", "symbol", "required", "available"
  ) %in% names(contract)))
  testthat::expect_true(all(contract$available[contract$required]))
  testthat::expect_true(".bpc_fit_integration_normal" %in% contract$symbol)

  pc_contract <- qc:::.qc_distribution_contract("t", method = "pc")
  testthat::expect_true("pc_fit_distribution.t" %in% pc_contract$symbol)
  testthat::expect_false(any(grepl("stanmodels", pc_contract$symbol, fixed = TRUE)))

  mcmc_contract <- qc:::.qc_distribution_contract("t", method = "mcmc")
  testthat::expect_true("stanmodels[[\"t\"]]" %in% mcmc_contract$symbol)
  testthat::expect_false(any(mcmc_contract$symbol == ".bpc_fit_integration_normal"))
})

testthat::test_that("distribution registry drives compatibility errors", {
  testthat::expect_error(
    qc:::.qc_distribution_spec("t", method = "integration"),
    regexp = "The integration method currently only supports distribution = 'normal'"
  )

  testthat::expect_error(
    qc:::.qc_distribution_spec("lognormal"),
    regexp = "not recognized by the 'distribution' argument"
  )
})

testthat::test_that("registered distributions satisfy the declared adapter contract", {
  for (distribution in qc:::.qc_distribution_names()) {
    testthat::expect_equal(
      qc:::.qc_distribution_missing_requirements(distribution),
      character(),
      info = distribution
    )
    for (method in c("pc", "mcmc")) {
      if (distribution %in% qc:::.qc_distribution_names(method = method)) {
        testthat::expect_equal(
          qc:::.qc_distribution_missing_requirements(distribution, method = method),
          character(),
          info = paste(distribution, method)
        )
      }
    }
  }

  for (distribution in qc:::.qc_distribution_names(method = "integration")) {
    testthat::expect_equal(
      qc:::.qc_distribution_missing_requirements(distribution, method = "integration"),
      character(),
      info = paste(distribution, "integration")
    )
  }
})
