testthat::test_that(".new_qc_integration_request resolves midpoint target for Cpc and Cpm", {
  LSL <- 2
  USL <- 10
  midpoint <- (LSL + USL) / 2

  for (metric in c("Cpc", "Cpm")) {
    request <- qc:::.new_qc_integration_request(
      data = numeric(0),
      LSL = LSL,
      USL = USL,
      prior = qc:::create_prior_conjugate(),
      metric = metric,
      target = NULL
    )

    testthat::expect_s3_class(request, "qc_integration_request")
    testthat::expect_identical(request$requested_target, NULL)
    testthat::expect_equal(request$target, midpoint)
    testthat::expect_identical(request$metric, metric)
  }
})

testthat::test_that(".as_qc_integration_request accepts target = NULL for metrics that do not use it", {
  request <- qc:::.new_qc_integration_request(
    data = c(4.8, 5.1, 5.4),
    LSL = 2,
    USL = 10,
    prior = qc:::create_prior_conjugate(),
    metric = "Cpk",
    target = NULL
  )

  testthat::expect_silent({
    validated <- qc:::.as_qc_integration_request(request = request)
    testthat::expect_identical(validated, request)
  })
})

testthat::test_that("integration requests validate metric names, limits, target, and sigma level", {
  prior <- qc:::create_prior_conjugate()

  testthat::expect_error(
    qc:::.new_qc_integration_request(
      data = numeric(0),
      LSL = 0,
      USL = 1,
      prior = prior,
      metric = "not-a-metric",
      target = NULL
    ),
    regexp = "not recognized by the 'metric' argument"
  )

  testthat::expect_error(
    qc:::.new_qc_integration_request(
      data = numeric(0),
      LSL = 1,
      USL = 0,
      prior = prior,
      metric = "Cp",
      target = NULL
    ),
    regexp = "LSL"
  )

  testthat::expect_error(
    qc:::.new_qc_integration_request(
      data = numeric(0),
      LSL = 0,
      USL = 1,
      prior = prior,
      metric = "Cpm",
      target = 2
    ),
    regexp = "`target` must lie within"
  )

  testthat::expect_error(
    qc:::.new_qc_integration_request(
      data = numeric(0),
      LSL = 0,
      USL = 1,
      prior = prior,
      metric = "Cp",
      target = NULL,
      sigma_level = 0
    ),
    regexp = "`sigma_level` must be greater than 0"
  )
})

testthat::test_that("analyze_capability_integration matches old-argument and request-based invocation", {
  x <- c(4.4, 4.9, 5.2, 5.6, 5.1, 4.7)
  LSL <- 2
  USL <- 10
  prior <- qc:::create_prior_conjugate()
  request <- qc:::.new_qc_integration_request(
    data = x,
    LSL = LSL,
    USL = USL,
    prior = prior,
    metric = "Cpm",
    target = NULL
  )

  direct <- qc:::analyze_capability_integration(
    x, LSL, USL, prior,
    metric = "Cpm",
    target = NULL,
    use_density_solver = TRUE
  )
  via_request <- qc:::analyze_capability_integration(
    x, LSL, USL, prior,
    metric = "Cpm",
    target = NULL,
    use_density_solver = TRUE,
    request = request
  )

  testthat::expect_equal(direct$stats, via_request$stats)
  testthat::expect_equal(direct$area, via_request$area)
  testthat::expect_equal(direct$grid, via_request$grid)
})
