test_that("integration interval helpers resolve modes from degenerate, stats, and density inputs", {
  expect_equal(
    qc:::.integration_result_mode(list(degenerate = list(type = "point", value = 2))),
    2
  )

  expect_equal(
    qc:::.integration_result_mode(c(Mean = Inf, Median = 1, SD = 1)),
    Inf
  )

  expect_equal(
    qc:::.integration_result_mode(c(Mean = 1, Median = 2, SD = 0)),
    2
  )

  expect_equal(
    qc:::.integration_result_mode(
      c(Mean = 1, Median = 1, SD = 0.5),
      x = c(0, 1, 2),
      density = c(1, 3, 2)
    ),
    1
  )
})

test_that("integration density helpers honor explicit cell geometry and tail masses", {
  entry <- list(
    grid = data.frame(
      x = c(0.3, 1.5),
      density = c(2.0, 0.6),
      x_left = c(0.2, 1.0),
      x_right = c(0.4, 2.0),
      width = c(0.2, 1.0)
    ),
    area = 1,
    grid_mass = 0.7,
    mass_nonpositive = 0.1,
    left_tail_mass = 0.15,
    right_tail_mass = 0.05,
    support_lower = 0
  )

  expect_equal(
    qc:::.integration_quantiles_from_density(
      x = entry$grid,
      density = entry$grid$density,
      probs = c(0.25, 0.5, 0.75),
      area = entry$area,
      entry = entry
    ),
    c(0.2, 0.3785714286, 1.5238095238),
    tolerance = 1e-10
  )
  expect_equal(
    qc:::.integration_hdi_from_density(
      x = entry$grid,
      density = entry$grid$density,
      ci_level = 0.5,
      area = entry$area,
      entry = entry
    ),
    c(0.2, 0.4)
  )

  breaks <- c(-Inf, 0, 0.5, 1.5, Inf)
  expect_equal(
    qc:::.integration_interval_probs_from_grid(entry, breaks),
    c(0.1, 0.43, 0.21, 0.26),
    tolerance = 1e-10
  )

  finite_tail_entry <- list(
    grid = data.frame(x = c(0, 1, 2), density = c(1, 2, 1)),
    area = 1,
    grid_mass = 0.7,
    mass_nonpositive = 0,
    left_tail_mass = 0,
    right_tail_mass = 0.3,
    support_lower = 0
  )
  finite_breaks <- c(-Inf, 0.5, 1.5, 2.5)
  expect_equal(
    qc:::.integration_interval_probs_from_grid(finite_tail_entry, finite_breaks),
    c(0.25, 0.5, 0.25),
    tolerance = 1e-10
  )
})

test_that("integration density helpers interpolate unrepresented tail mass", {
  entry <- list(
    grid = data.frame(
      x = c(0.3, 1.5),
      density = c(2.0, 0.6),
      x_left = c(0.2, 1.0),
      x_right = c(0.4, 2.0),
      width = c(0.2, 1.0)
    ),
    area = 1,
    grid_mass = 0.7,
    mass_nonpositive = 0.1,
    left_tail_mass = 0.15,
    right_tail_mass = 0.05,
    support_lower = 0
  )

  expect_equal(
    qc:::.integration_quantiles_from_density(
      x = entry$grid,
      density = entry$grid$density,
      probs = 0.975,
      area = entry$area,
      entry = entry
    ),
    2.0595238095,
    tolerance = 1e-10
  )

  expect_equal(
    qc:::.integration_hdi_from_density(
      x = entry$grid,
      density = entry$grid$density,
      ci_level = 0.20,
      area = entry$area,
      entry = entry
    ),
    c(0, 0.1333333333),
    tolerance = 1e-8
  )

  expect_equal(
    qc:::.integration_hdi_from_density(
      x = entry$grid,
      density = entry$grid$density,
      ci_level = 0.98,
      area = entry$area,
      entry = entry
    ),
    c(0.2, 2.0714285714),
    tolerance = 1e-10
  )
})

test_that("integration interval probability helper covers degenerate and sample-backed entries", {
  breaks <- c(-Inf, 1, 2, Inf)

  degenerate_probs <- qc:::.integration_interval_probs_from_entry(
    entry = list(degenerate = list(type = "point", value = 1.5)),
    interval_breaks = breaks,
    metric = "Cp",
    prior = NULL,
    cached_state = NULL,
    LSL = 0,
    USL = 1,
    target = NULL,
    sigma_level = 3
  )
  expect_equal(degenerate_probs, c(0, 1, 0))

  sample_probs <- qc:::.integration_interval_probs_from_entry(
    entry = list(samples = c(0.5, 1.5, 2.5, 1.5)),
    interval_breaks = breaks,
    metric = "Cp",
    prior = NULL,
    cached_state = NULL,
    LSL = 0,
    USL = 1,
    target = NULL,
    sigma_level = 3
  )
  expect_equal(sample_probs, c(0.25, 0.5, 0.25))
})

test_that("integration interval probability helper assigns point masses to the first matching finite bin", {
  breaks <- c(1, 2, 3)

  expect_equal(
    qc:::.integration_interval_probs_from_entry(
      entry = list(degenerate = list(type = "point", value = 1)),
      interval_breaks = breaks,
      metric = "Cp",
      prior = NULL,
      cached_state = NULL,
      LSL = 0,
      USL = 1,
      target = NULL,
      sigma_level = 3
    ),
    c(1, 0)
  )

  expect_equal(
    qc:::.integration_interval_probs_from_entry(
      entry = list(degenerate = list(type = "point", value = 2)),
      interval_breaks = breaks,
      metric = "Cp",
      prior = NULL,
      cached_state = NULL,
      LSL = 0,
      USL = 1,
      target = NULL,
      sigma_level = 3
    ),
    c(1, 0)
  )
})

test_that("integration interval probability helpers reject malformed break vectors", {
  point_entry <- list(degenerate = list(type = "point", value = 1))

  expect_error(
    qc:::.integration_interval_probs_from_entry(
      entry = point_entry,
      interval_breaks = c(0, Inf, 2),
      metric = "Cp",
      prior = NULL,
      cached_state = NULL,
      LSL = 0,
      USL = 1,
      target = NULL,
      sigma_level = 3
    ),
    "finite interior break points"
  )

  expect_error(
    qc:::.integration_interval_probs_from_entry(
      entry = point_entry,
      interval_breaks = c(0, 2, 1),
      metric = "Cp",
      prior = NULL,
      cached_state = NULL,
      LSL = 0,
      USL = 1,
      target = NULL,
      sigma_level = 3
    ),
    "strictly increasing"
  )
})

test_that("integration interval probability helper falls back to grid probabilities when solver output is unusable", {
  breaks <- c(-Inf, 0.5, 1.5, Inf)
  entry <- list(
    grid = data.frame(x = c(0, 1, 2), density = c(1, 2, 1)),
    area = 1
  )
  expected <- qc:::.integration_interval_probs_from_grid(entry, breaks)

  testthat::local_mocked_bindings(
    compute_cpk_prob_integration = function(...) NA_real_,
    .package = "qc"
  )

  actual <- qc:::.integration_interval_probs_from_entry(
    entry = entry,
    interval_breaks = breaks,
    metric = "Cp",
    prior = qc:::create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 1, beta0 = 1),
    cached_state = list(n = 10L, x_bar = 0, sse = 1),
    LSL = 0,
    USL = 1,
    target = NULL,
    sigma_level = 3
  )

  expect_equal(actual, expected)
})

test_that("integration interval probability helper uses solver probabilities when no grid is available", {
  breaks <- c(-Inf, 1, 2, Inf)
  values <- c(0.2, 0.3, 0.5)
  call_index <- 0L

  testthat::local_mocked_bindings(
    compute_cpk_prob_integration = function(...) {
      call_index <<- call_index + 1L
      values[[call_index]]
    },
    .package = "qc"
  )

  actual <- qc:::.integration_interval_probs_from_entry(
    entry = list(stats = c(Mean = 1, Median = 1, SD = 0.1)),
    interval_breaks = breaks,
    metric = "Cp",
    prior = qc:::create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 1, beta0 = 1),
    cached_state = list(n = 10L, x_bar = 0, sse = 1),
    LSL = 0,
    USL = 1,
    target = NULL,
    sigma_level = 3
  )

  expect_equal(actual, values)
})

test_that("integration interval extraction honors ci_level for sample and density-backed entries", {
  sample_entry <- list(samples = c(0, 1, 2, 3))
  expect_equal(
    qc:::.integration_interval_from_entry(sample_entry, ci = "central", ci_level = 0.50),
    c(0.75, 2.25)
  )

  density_entry <- list(
    grid = data.frame(
      x = c(0.3, 1.5),
      density = c(2.0, 0.6),
      x_left = c(0.2, 1.0),
      x_right = c(0.4, 2.0),
      width = c(0.2, 1.0)
    ),
    area = 1,
    grid_mass = 0.7,
    mass_nonpositive = 0.1,
    left_tail_mass = 0.15,
    right_tail_mass = 0.05,
    support_lower = 0
  )
  expect_equal(
    qc:::.integration_interval_from_entry(density_entry, ci = "central", ci_level = 0.50),
    c(0.2, 1.5238095238),
    tolerance = 1e-10
  )
  expect_equal(
    qc:::.integration_interval_from_entry(density_entry, ci = "HPD", ci_level = 0.50),
    c(0.2, 0.4)
  )
})

test_that("adaptive density analyses keep tail truncation separate from nonpositive mass", {
  prior <- structure(list(), class = "PriorConjugate")

  testthat::local_mocked_bindings(
    .compute_validated_conjugate_posterior = function(...) {
      list(post = list(mu_n = 0, k_n = 1), is_degenerate = FALSE)
    },
    compute_metric_moments = function(...) {
      list(mean = 1, sd = 0.5)
    },
    make_density_solver = function(...) {
      function(x) exp(-x)
    },
    .package = "qc"
  )

  result <- qc:::analyze_capability_integration(
    data = numeric(0),
    LSL = 0,
    USL = 1,
    prior = prior,
    metric = "Cp",
    use_density_solver = TRUE
  )

  expect_true(all(c("x_left", "x_right", "width") %in% names(result$grid)))
  expect_lt(result$grid_mass, 1)
  expect_equal(result$mass_nonpositive, 0)
  expect_gt(result$right_tail_mass, 0)

  probs <- qc:::.integration_interval_probs_from_entry(
    entry = result,
    interval_breaks = c(-Inf, 0, 1, Inf),
    metric = "Cp",
    prior = prior,
    cached_state = NULL,
    LSL = 0,
    USL = 1,
    target = NULL,
    sigma_level = 3
  )

  expect_equal(probs[1], 0)
  expect_gt(probs[3], 0)
})

test_that("survival-grid analyses keep uncovered right tails out of the last finite interval", {
  fake_survival <- function(c) {
    c <- as.numeric(c)
    ifelse(c <= 0, 1, exp(-c))
  }

  prior <- qc:::create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = 2, beta0 = 1)

  testthat::local_mocked_bindings(
    make_density_solver = function(...) stop("density boom", call. = FALSE),
    make_solver = function(...) fake_survival,
    compute_metric_moments = function(...) list(mean = 1, sd = 0.25),
    .package = "qc"
  )

  result <- qc:::analyze_capability_integration(
    data = numeric(0),
    LSL = 0,
    USL = 1,
    prior = prior,
    metric = "Cp",
    use_density_solver = TRUE
  )

  expect_lt(result$grid_mass, 1)
  expect_gt(result$right_tail_mass, 0)

  finite_probs <- qc:::.integration_interval_probs_from_grid(
    result,
    c(-Inf, 0.5, 1.5, 2.5)
  )
  open_probs <- qc:::.integration_interval_probs_from_grid(
    result,
    c(-Inf, 0.5, 1.5, Inf)
  )

  expect_lt(finite_probs[3], open_probs[3])
  expect_equal(sum(finite_probs), 1, tolerance = 1e-10)
  expect_equal(sum(open_probs), 1, tolerance = 1e-10)
})

test_that("integration interval extraction fails clearly for stats-only non-95 requests", {
  stats_entry <- list(stats = c(
    Mean = 1,
    Median = 1,
    SD = 0.5,
    Q2.5 = 0,
    Q97.5 = 2,
    HDI_Lo = 0.25,
    HDI_Hi = 1.75
  ))

  expect_equal(
    unname(qc:::.integration_interval_from_entry(stats_entry, ci = "central", ci_level = 0.95)),
    c(0, 2)
  )

  expect_equal(
    unname(qc:::.integration_interval_from_entry(stats_entry, ci = "HPD", ci_level = 0.95)),
    c(0.25, 1.75)
  )

  expect_error(
    qc:::.integration_interval_from_entry(stats_entry, ci = "central", ci_level = 0.50),
    "ci_level = 0.95"
  )

  expect_error(
    qc:::.integration_interval_from_entry(stats_entry, ci = "HPD", ci_level = 0.50),
    "ci_level = 0.95"
  )
})
