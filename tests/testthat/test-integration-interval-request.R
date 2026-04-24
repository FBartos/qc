testthat::test_that("interval probability reconstruction preserves request-equivalent solver inputs", {
  LSL <- 2
  USL <- 10
  prior <- qc:::create_prior_conjugate()
  cached_state <- list(n = 6L, x_bar = 5, sse = 1.25)
  request <- qc:::.new_qc_integration_request(
    data = numeric(0),
    LSL = LSL,
    USL = USL,
    prior = prior,
    metric = "Cpk",
    target = NULL,
    cached_state = cached_state
  )

  captured_requests <- list()
  call_index <- 0L
  solve_probs <- c(
    "[-Inf,1]" = 0.2,
    "[1,2]" = 0.3,
    "[2,Inf]" = 0.5
  )

  testthat::local_mocked_bindings(
    compute_cpk_prob_integration = function(data = NULL, LSL, USL, bounds, prior,
                                            metric = "Cpk", target = NULL,
                                            cached_state = NULL,
                                            sigma_level = 3,
                                            request = NULL) {
      call_index <<- call_index + 1L
      captured_requests[[call_index]] <<- if (!is.null(request)) {
        request
      } else {
        qc:::.as_qc_integration_request(
          data = data,
          LSL = LSL,
          USL = USL,
          prior = prior,
          metric = metric,
          target = target,
          cached_state = cached_state,
          sigma_level = sigma_level
        )
      }
      solve_probs[[paste0("[", bounds[1], ",", bounds[2], "]")]]
    },
    .package = "qc"
  )

  loose_probs <- qc:::.integration_interval_probs_from_entry(
    entry = list(stats = c(Mean = 1, Median = 1, SD = 0.1)),
    interval_breaks = c(-Inf, 1, 2, Inf),
    metric = "Cpk",
    prior = prior,
    cached_state = cached_state,
    LSL = LSL,
    USL = USL,
    target = NULL,
    sigma_level = 3
  )

  request_probs <- qc:::.integration_interval_probs_from_entry(
    entry = list(stats = c(Mean = 1, Median = 1, SD = 0.1)),
    interval_breaks = c(-Inf, 1, 2, Inf),
    metric = "Cpk",
    prior = prior,
    cached_state = cached_state,
    LSL = LSL,
    USL = USL,
    target = NULL,
    sigma_level = 3,
    request = request
  )

  testthat::expect_equal(loose_probs, c(0.2, 0.3, 0.5))
  testthat::expect_equal(request_probs, loose_probs)
  testthat::expect_length(captured_requests, 6)
  for (captured in captured_requests) {
    testthat::expect_equal(captured, request)
  }
})

testthat::test_that("compute_cpk_prob_integration treats explicit request and loose args the same", {
  LSL <- 2
  USL <- 10
  prior <- qc:::create_prior_conjugate()
  cached_state <- list(n = 6L, x_bar = 5, sse = 1.25)
  request <- qc:::.new_qc_integration_request(
    data = numeric(0),
    LSL = LSL,
    USL = USL,
    prior = prior,
    metric = "Cpk",
    target = NULL,
    cached_state = cached_state
  )

  seen_requests <- list()
  call_index <- 0L

  testthat::local_mocked_bindings(
    .integration_backend_resolver = function(request = NULL, prefer_density = TRUE, ...) {
      call_index <<- call_index + 1L
      seen_requests[[call_index]] <<- request
      list(
        mode = "survival",
        pdf_fn = NULL,
        S = function(c) ifelse(c <= 0, 1, exp(-c)),
        can_use_density = FALSE,
        density_error = NULL
      )
    },
    .package = "qc"
  )

  loose <- qc:::compute_cpk_prob_integration(
    data = numeric(0),
    LSL = LSL,
    USL = USL,
    bounds = c(0.5, 1.5),
    prior = prior,
    metric = "Cpk",
    target = NULL,
    cached_state = cached_state,
    sigma_level = 3
  )

  explicit <- qc:::compute_cpk_prob_integration(
    data = numeric(0),
    LSL = LSL,
    USL = USL,
    bounds = c(0.5, 1.5),
    prior = prior,
    metric = "Cpk",
    target = NULL,
    cached_state = cached_state,
    sigma_level = 3,
    request = request
  )

  testthat::expect_equal(loose, explicit)
  testthat::expect_length(seen_requests, 2)
  for (seen in seen_requests) {
    testthat::expect_equal(seen, request)
  }
})
