make_integration_adapter_stub <- function(adapter_calls_env) {
  metric_names <- qc:::.qc_metric_names()

  function(distribution, data, LSL, USL, target, prior,
           sigma = 3, sample_priors = FALSE, cached_state = NULL) {
    adapter_calls_env$calls <- adapter_calls_env$calls + 1L
    adapter_calls_env$distribution <- distribution

    list(
      metrics = setNames(
        as.list(rep(1, length(metric_names))),
        metric_names
      ),
      coefficients = setNames(rep(1, length(metric_names)), metric_names),
      prior = prior,
      cached_state = cached_state,
      results = setNames(vector("list", length(metric_names)), metric_names)
    )
  }
}

testthat::test_that("integration-capable distributions resolve through the central registry before the adapter", {
  adapter_calls <- new.env(parent = emptyenv())
  adapter_calls$calls <- 0L

  testthat::local_mocked_bindings(
    .bpc_fit_integration = make_integration_adapter_stub(adapter_calls),
    .package = "qc"
  )

  fit <- qc::bpc(
    c(1, 2, 3, 4),
    LSL = 0,
    target = 2,
    USL = 4,
    distribution = "normal",
    method = "integration"
  )

  testthat::expect_equal(qc:::.qc_distribution_names(method = "integration"), "normal")
  testthat::expect_equal(adapter_calls$calls, 1L)
  testthat::expect_identical(adapter_calls$distribution, "normal")
  testthat::expect_s3_class(fit, "bpc")
  testthat::expect_equal(fit$distribution, "normal")
  testthat::expect_equal(fit$method, "integration")
})

testthat::test_that("unsupported integration distributions fail cleanly through the registry gate", {
  adapter_calls <- new.env(parent = emptyenv())
  adapter_calls$calls <- 0L

  testthat::local_mocked_bindings(
    .bpc_fit_integration = make_integration_adapter_stub(adapter_calls),
    .package = "qc"
  )

  testthat::expect_error(
    qc::bpc(
      c(1, 2, 3, 4),
      LSL = 0,
      target = 2,
      USL = 4,
      distribution = "t",
      method = "integration"
    ),
    regexp = "The integration method currently only supports distribution = 'normal'"
  )

  testthat::expect_equal(adapter_calls$calls, 0L)
})
