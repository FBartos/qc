.make_stub_integration_result <- function(LSL,
                                          target,
                                          USL,
                                          sigma,
                                          prior,
                                          cached_state,
                                          case) {
  metric_names <- qc:::.qc_metric_names()
  metric_samples <- setNames(lapply(seq_along(metric_names), function(i) {
    seq(0.8 + 0.02 * i, 1.2 + 0.02 * i, length.out = 64)
  }), metric_names)

  distributions <- setNames(Map(function(metric, samples) {
    qc:::.new_qc_metric_distribution(metric = metric, samples = samples)
  }, metric_names, metric_samples), metric_names)

  divergence <- setNames(vector("list", length(metric_names)), metric_names)
  metrics <- qc:::.new_capability_metrics(
    metrics = metric_samples,
    LSL = LSL,
    USL = USL,
    target = target,
    sigma = sigma,
    method = "integration",
    distributions = distributions,
    prior = prior,
    cached_state = cached_state,
    divergence = divergence
  )

  results <- setNames(lapply(metric_names, function(metric) {
    samples <- metric_samples[[metric]]
    list(
      samples = samples,
      stats = qc:::.sample_metric_stats(samples)
    )
  }), metric_names)

  qc:::.new_integration_result(
    results = results,
    metrics = metrics,
    coefficients = stats::setNames(
      vapply(metric_samples, stats::median, numeric(1)),
      metric_names
    ),
    sigma = sigma,
    prior = prior,
    is_conjugate = identical(case, 1L),
    case = case,
    cached_state = cached_state,
    divergence = divergence
  )
}

testthat::test_that("integration summary recomputes from cached sufficient statistics after data leaves scope", {
  fit_model <- function() {
    set.seed(123)
    local_data <- rnorm(30, mean = 10, sd = 1)

    qc::bpc(
      local_data,
      LSL = 7,
      target = 10,
      USL = 13,
      method = "integration"
    )
  }

  fit <- fit_model()

  testthat::expect_false(exists("local_data", envir = environment(), inherits = FALSE))
  testthat::expect_true(all(c("n", "x_bar", "sse") %in% names(fit$integration_result$cached_state)))

  summ <- summary(fit, LSL = 7, target = 10, USL = 13)

  testthat::expect_s3_class(summ, "bpc_summary")
  testthat::expect_true(!is.null(summ$summary))

  cpk_idx <- which(summ$summary$metric == "Cpk")
  testthat::expect_equal(summ$summary$mean[cpk_idx], 1.0, tolerance = 0.5)
})

testthat::test_that("integration summary forwards cached generic state instead of re-reading out-of-scope data", {
  expected_mu_prior <- BayesTools::prior("normal", list(10, 5))
  expected_sigma_prior <- BayesTools::prior("gamma", list(1, 1))
  prior_info <- qc:::.bayestools_to_integration_prior(
    expected_mu_prior,
    expected_sigma_prior
  )
  prior_spec <- qc::prior_independent(mu = expected_mu_prior, sigma = expected_sigma_prior)

  expected_cached_state <- list(
    n = 30L,
    x_bar = 10,
    sse = 29,
    log_post = function(mu, sigma) mu + sigma,
    log_post_vec = function(mu, sigma) mu + sigma,
    h_max = -1,
    uni_s = 1,
    map_mu = 10,
    map_sig = 1,
    Z = 1,
    mu_scale = 1,
    int_2d = function(f) 0
  )

  initial_result <- .make_stub_integration_result(
    LSL = 7,
    target = 10,
    USL = 13,
    sigma = 3,
    prior = prior_info$prior,
    cached_state = expected_cached_state,
    case = 4L
  )

  fit <- qc:::.new_bpc_fit(
    call = quote(
      bpc(
        local_data,
        LSL = 7,
        target = 10,
        USL = 13,
        method = "integration",
        prior = prior_independent(
          mu = BayesTools::prior("normal", list(10, 5)),
          sigma = BayesTools::prior("gamma", list(1, 1))
        )
      )
    ),
    method = "integration",
    distribution = "normal",
    metrics = initial_result$metrics,
    coefficients = initial_result$coefficients,
    sigma = 3,
    prior = prior_spec,
    prior_resolved = prior_info$prior,
    prior_map = list(mu = expected_mu_prior, sigma = expected_sigma_prior),
    integration_result = initial_result
  )

  testthat::local_mocked_bindings(
    .bpc_fit_integration = function(distribution,
                                    data,
                                    LSL,
                                    USL,
                                    target,
                                    prior,
                                    sigma,
                                    sample_priors = FALSE,
                                    cached_state = NULL) {
      testthat::expect_identical(distribution, "normal")
      testthat::expect_type(data, "double")
      testthat::expect_length(data, 0)
      testthat::expect_identical(prior, prior_info$prior)
      testthat::expect_identical(sample_priors, FALSE)
      testthat::expect_identical(cached_state$n, expected_cached_state$n)
      testthat::expect_identical(cached_state$x_bar, expected_cached_state$x_bar)
      testthat::expect_identical(cached_state$sse, expected_cached_state$sse)
      testthat::expect_true(is.function(cached_state$log_post_vec))
      testthat::expect_true(is.function(cached_state$int_2d))

      .make_stub_integration_result(
        LSL = LSL,
        target = target,
        USL = USL,
        sigma = sigma,
        prior = prior_info$prior,
        cached_state = cached_state,
        case = 4L
      )
    },
    .package = "qc"
  )

  testthat::expect_false(exists("local_data", envir = environment(), inherits = FALSE))

  summ <- summary(fit, LSL = 6, target = 10, USL = 14)

  testthat::expect_s3_class(summ, "bpc_summary")
  testthat::expect_setequal(summ$summary$metric, qc:::.qc_metric_names())
})
