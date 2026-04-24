testthat::test_that("capability entry points reject targets outside the specification limits", {
  regexp <- paste0(
    "`target` must lie within \\[`LSL`, `USL`\\] ",
    "so Cpm/Cpc are well-defined"
  )

  testthat::expect_no_error(qc:::.validate_LSL_USL_target(LSL = 0, USL = 1, target = 0))
  testthat::expect_no_error(qc:::.validate_LSL_USL_target(LSL = 0, USL = 1, target = 1))

  testthat::expect_error(
    qc::bpc(rnorm(8), LSL = 0, target = 2, USL = 1, method = "integration"),
    regexp = regexp
  )
  testthat::expect_error(
    qc::pc(rnorm(8), LSL = 0, target = -0.25, USL = 1, bootstrap = FALSE),
    regexp = regexp
  )
})

testthat::test_that("capability entry points reject nonpositive sigma levels", {
  testthat::expect_error(
    qc::bpc(rnorm(8), LSL = 0, target = 0.5, USL = 1, method = "integration", sigma = 0),
    regexp = "`sigma` must be greater than 0"
  )
  testthat::expect_error(
    qc::pc(rnorm(8), LSL = 0, target = 0.5, USL = 1, bootstrap = FALSE, sigma = 0),
    regexp = "`sigma` must be greater than 0"
  )
})

testthat::test_that("data preparation centralizes summary-statistics validation by distribution", {
  prepared <- qc:::.bpc_prepare_data(
    distribution = "normal",
    x = NULL,
    mean = 1.5,
    sd = 2,
    N = 4L
  )

  testthat::expect_equal(prepared$raw_data, numeric())
  testthat::expect_equal(prepared$stan_data$N, 4L)
  testthat::expect_equal(prepared$stan_data$is_ss, 1L)
  testthat::expect_equal(prepared$cached_state$n, 4L)
  testthat::expect_equal(prepared$cached_state$x_bar, 1.5)
  testthat::expect_equal(prepared$cached_state$sse, 12)

  testthat::expect_error(
    qc:::.bpc_prepare_data(
      distribution = "t",
      x = NULL,
      mean = 1.5,
      sd = 2,
      N = 4L
    ),
    regexp = "Summary-statistics inputs"
  )

  testthat::expect_error(
    qc:::.bpc_prepare_data(
      distribution = "t",
      x = NULL,
      mean = 1.5,
      sd = 2,
      N = 4L,
      allow_empty = TRUE
    ),
    regexp = "omit `mean`, `sd`, and `N`"
  )

  empty_t <- qc:::.bpc_prepare_data(
    distribution = "t",
    x = NULL,
    allow_empty = TRUE
  )

  testthat::expect_equal(empty_t$raw_data, numeric())
  testthat::expect_equal(empty_t$stan_data$N, 0L)
  testthat::expect_equal(empty_t$cached_state$n, 0L)
})

testthat::test_that("documented integration helpers are available through the package namespace", {
  expected_exports <- c(
    "compute_metric_value",
    "create_prior_generic",
    "create_prior_semi_mu",
    "create_prior_semi_sigma",
    "get_metric_constraints"
  )

  testthat::expect_true(all(expected_exports %in% getNamespaceExports("qc")))

  generic_prior <- qc::create_prior_generic(function(mu, sigma) rep(0, max(length(mu), length(sigma))))
  semi_mu_prior <- qc::create_prior_semi_mu(0, 1, function(sigma) rep(0, length(sigma)))
  semi_sigma_prior <- qc::create_prior_semi_sigma(1, 1, function(mu) rep(0, length(mu)))
  constraints <- qc::get_metric_constraints("Cp", c = 1, LSL = 0, USL = 10, target = 5)

  testthat::expect_s3_class(generic_prior, "PriorGeneric")
  testthat::expect_s3_class(semi_mu_prior, "PriorSemiConjugateMu")
  testthat::expect_s3_class(semi_sigma_prior, "PriorSemiConjugateSigma")
  testthat::expect_type(constraints$s_max_fn, "closure")
  testthat::expect_equal(constraints$s_max_fn(), 10 / 6)
  testthat::expect_equal(
    qc::compute_metric_value(mu = 5, sigma = 1, LSL = 0, USL = 10, target = 5, metric = "Cp"),
    10 / 6
  )
})

testthat::test_that("public integration prior constructors are accepted by bpc integration", {
  metric_names <- c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm")
  prior_cases <- list(
    generic = qc::create_prior_generic(function(mu, sigma) {
      ifelse(
        sigma <= 0,
        -Inf,
        stats::dnorm(mu, mean = 5, sd = 2, log = TRUE) +
          stats::dgamma(sigma, shape = 2, rate = 1, log = TRUE)
      )
    }),
    semi_mu = qc::create_prior_semi_mu(
      mu0 = 5,
      k0 = 1,
      log_dens_sigma = function(sigma) {
        stats::dgamma(sigma, shape = 2, rate = 1, log = TRUE)
      }
    ),
    semi_sigma = qc::create_prior_semi_sigma(
      alpha0 = 2,
      beta0 = 1,
      log_dens_mu = function(mu) {
        stats::dnorm(mu, mean = 5, sd = 2, log = TRUE)
      }
    )
  )

  for (case_name in names(prior_cases)) {
    prior_obj <- prior_cases[[case_name]]

    testthat::local_mocked_bindings(
      .bpc_fit_integration = function(distribution, data, LSL, USL, target, prior,
                                      sigma = 3, sample_priors = FALSE,
                                      cached_state = NULL) {
        testthat::expect_identical(distribution, "normal", info = case_name)
        testthat::expect_identical(prior, prior_obj, info = case_name)
        testthat::expect_identical(sample_priors, TRUE, info = case_name)

        metric_stub <- stats::setNames(rep(1, length(metric_names)), metric_names)
        list(
          metrics = as.list(metric_stub),
          coefficients = metric_stub,
          prior = prior,
          case = qc:::.integration_prior_case(prior),
          cached_state = cached_state,
          results = setNames(vector("list", length(metric_names)), metric_names)
        )
      },
      .package = "qc"
    )

    fit <- testthat::expect_no_error(
      qc::bpc(
        NULL,
        LSL = 2,
        target = 5,
        USL = 8,
        prior = prior_obj,
        method = "integration",
        sample_priors = TRUE
      )
    )

    testthat::expect_s3_class(fit, "bpc")
    testthat::expect_identical(fit$prior_resolved, prior_obj, info = case_name)
    testthat::expect_equal(unname(fit$coefficients), rep(1, length(metric_names)), info = case_name)
  }
})

testthat::test_that("capability metrics fall back to default E_abs_dev for extension distributions", {
  predictive_draws <- c(-1, 1, -2, 2)
  predictive_calls <- 0L

  base::registerS3method(
    "extract_samples",
    "qc_test_fit_default_e_abs_dev",
    function(fit, bootstrap) {
      structure(
        list(draw_id = seq_along(predictive_draws)),
        class = "qc_test_samples_default_e_abs_dev"
      )
    },
    envir = asNamespace("qc")
  )
  base::registerS3method(
    "samples_to_percentiles",
    "qc_test_samples_default_e_abs_dev",
    function(samples, sigma) {
      n <- length(samples[[1L]])
      list(
        LP = rep(-3, n),
        MP = rep(0, n),
        UP = rep(3, n)
      )
    },
    envir = asNamespace("qc")
  )
  base::registerS3method(
    "samples_to_posterior_predictives",
    "qc_test_samples_default_e_abs_dev",
    function(samples) {
      predictive_calls <<- predictive_calls + 1L
      predictive_draws
    },
    envir = asNamespace("qc")
  )

  fit <- structure(list(), class = "qc_test_fit_default_e_abs_dev")
  metrics <- qc:::.compute_capability_metrics(
    fit,
    LSL = -3,
    USL = 3,
    target = 0,
    sigma = 3
  )

  expected_cpc <- 1 / (sqrt(pi / 2) * mean(abs(predictive_draws)))

  testthat::expect_s3_class(metrics, "capability_metrics")
  testthat::expect_equal(as.numeric(metrics$Cp), rep(1, length(predictive_draws)))
  testthat::expect_equal(as.numeric(metrics$Cpk), rep(1, length(predictive_draws)))
  testthat::expect_equal(as.numeric(metrics$Cpm), rep(1, length(predictive_draws)))
  testthat::expect_equal(
    as.numeric(metrics$Cpc),
    rep(expected_cpc, length(predictive_draws))
  )
  testthat::expect_gt(predictive_calls, 0L)

  distributions <- attr(metrics, "distributions")
  testthat::expect_named(distributions, qc:::.qc_metric_names())
  testthat::expect_true(all(vapply(distributions, inherits, logical(1), "qc_metric_distribution")))
})

testthat::test_that("beta alpha extraction matches the BayesTools alpha/beta schema", {
  testthat::skip_if_not_installed("BayesTools")

  prior_beta <- BayesTools::prior("beta", list(alpha = 2, beta = 4))

  testthat::expect_named(prior_beta[["parameters"]], c("alpha", "beta"))
  testthat::expect_equal(
    qc:::.extract_alpha_parameter(prior_beta),
    prior_beta[["parameters"]][["alpha"]]
  )
})

testthat::test_that("beta prior log-density respects supported boundaries", {
  testthat::skip_if_not_installed("BayesTools")

  beta_uniform <- qc:::.make_prior_log_dens_fn(
    BayesTools::prior("beta", list(alpha = 1, beta = 1))
  )
  beta_left_singular <- qc:::.make_prior_log_dens_fn(
    BayesTools::prior("beta", list(alpha = 0.5, beta = 2))
  )
  beta_right_singular <- qc:::.make_prior_log_dens_fn(
    BayesTools::prior("beta", list(alpha = 2, beta = 0.5))
  )
  beta_zero_boundary <- qc:::.make_prior_log_dens_fn(
    BayesTools::prior("beta", list(alpha = 2, beta = 2))
  )

  testthat::expect_equal(unname(beta_uniform(c(0, 1))), c(0, 0))
  testthat::expect_true(is.infinite(unname(beta_left_singular(0))))
  testthat::expect_gt(unname(beta_left_singular(0)), 0)
  testthat::expect_true(is.infinite(unname(beta_right_singular(1))))
  testthat::expect_gt(unname(beta_right_singular(1)), 0)
  testthat::expect_equal(unname(beta_zero_boundary(c(0, 1))), c(-Inf, -Inf))
})
