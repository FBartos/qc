testthat::test_that("extract_point_estimates keeps finite density-backed modes when stats diverge", {
  entry <- list(
    stats = c(Mean = Inf, Median = Inf, SD = 1),
    grid = data.frame(
      x = c(1, 2, 3),
      density = c(1, 5, 2)
    )
  )
  df_density <- tibble::tibble(
    x = entry$grid$x,
    density = entry$grid$density,
    metric = factor(rep("Cp", nrow(entry$grid)), levels = "Cp")
  )

  points <- qc::extract_point_estimates(
    obj = NULL,
    what = "Cp",
    point_estimate = "mode",
    dfDensity = df_density,
    stats_list = list(Cp = entry)
  )

  testthat::expect_equal(nrow(points), 1)
  testthat::expect_equal(points$x, 2)
  testthat::expect_equal(points$y, 5)
})

testthat::test_that("distribution interval probabilities keep point masses on the first finite boundary", {
  breaks <- c(1, 2, 3)
  dist <- qc:::.new_qc_metric_distribution(
    metric = "Cp",
    entry = list(degenerate = list(type = "point", value = 1))
  )

  testthat::expect_equal(
    qc:::.qc_metric_distribution_interval_probs(dist, breaks),
    c(1, 0)
  )
  testthat::expect_equal(
    qc:::.qc_metric_distribution_interval_probs(dist, breaks),
    qc:::.integration_interval_probs_from_entry(
      entry = dist$entry,
      interval_breaks = breaks,
      metric = "Cp",
      prior = NULL,
      cached_state = NULL,
      LSL = 0,
      USL = 1,
      target = NULL,
      sigma_level = 3
    )
  )
})

testthat::test_that("extract_ci_data drops non-finite integration intervals", {
  entry <- list(
    stats = c(
      Mean = Inf,
      Median = Inf,
      SD = Inf,
      Q2.5 = Inf,
      Q97.5 = Inf,
      HDI_Lo = Inf,
      HDI_Hi = Inf
    )
  )
  df_density <- tibble::tibble(
    x = c(0, 1),
    density = c(0, 0),
    metric = factor(c("Cp", "Cp"), levels = "Cp")
  )

  ci_data <- qc::extract_ci_data(
    obj = NULL,
    what = "Cp",
    ci = "HPD",
    ci_level = 0.95,
    dfDensity = df_density,
    stats_list = list(Cp = entry)
  )

  testthat::expect_equal(nrow(ci_data$dfCi), 0)
  testthat::expect_equal(nrow(ci_data$dfArea), 0)
})
