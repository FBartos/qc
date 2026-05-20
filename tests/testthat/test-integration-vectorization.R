testthat::test_that("adaptive density grids evaluate vectorized density batches", {
  call_lengths <- integer()
  pdf_fn <- function(x) {
    call_lengths <<- c(call_lengths, length(x))
    stats::dnorm(x)
  }

  result <- qc:::.adaptive_density_grid(pdf_fn, 0, 3, n_grid = 32L)

  testthat::expect_true(any(call_lengths > 1L))
  testthat::expect_lt(length(call_lengths), length(result$mid_x))
  testthat::expect_true(all(is.finite(result$pdf_vals)))
})

testthat::test_that("survival-derived densities evaluate vectorized survival batches", {
  call_lengths <- integer()
  S <- function(c) {
    call_lengths <<- c(call_lengths, length(c))
    ifelse(c <= 0, 1, exp(-c))
  }

  pdf_fn <- qc:::.density_from_survival_fn(S, support_lower = 0)
  vals <- pdf_fn(c(0, 0.5, 1))

  testthat::expect_true(any(call_lengths > 1L))
  testthat::expect_equal(length(vals), 3L)
  testthat::expect_true(all(is.finite(vals) & vals >= 0))
})

testthat::test_that("density backends used by adaptive grids accept vector inputs", {
  testthat::skip_if_not_installed("BayesTools")
  testthat::skip_if_not_installed("cubature")

  set.seed(11)
  x <- stats::rnorm(8)
  cases <- list(
    conjugate = list(
      prior = qc:::create_prior_conjugate(),
      metrics = c("Cp", "Cpm", "Cpc")
    ),
    semi_mu = list(
      prior = qc:::.bayestools_to_integration_prior(
        "Jeffreys_mu",
        BayesTools::prior("gamma", list(2, 1))
      )$prior,
      metrics = c("Cp", "Cpm", "Cpc")
    ),
    generic = list(
      prior = qc:::.bayestools_to_integration_prior(
        BayesTools::prior("normal", list(0, 1)),
        BayesTools::prior("gamma", list(2, 1))
      )$prior,
      metrics = c("Cp", "Cpm")
    )
  )

  points <- c(0.25, 0.5, 1)
  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    for (metric in case$metrics) {
      pdf_fn <- qc:::.integration_make_density_solver(
        data = x,
        LSL = -3,
        USL = 3,
        target = 0,
        prior = case$prior,
        metric = metric
      )
      vals <- pdf_fn(points)
      testthat::expect_equal(length(vals), length(points), info = paste(case_name, metric))
      testthat::expect_true(all(is.finite(vals) & vals >= 0), info = paste(case_name, metric))
    }
  }
})

testthat::test_that("conjugate smooth-metric survival solvers accept vector thresholds", {
  set.seed(12)
  x <- stats::rnorm(8)
  prior <- qc:::create_prior_conjugate()

  for (metric in c("Cpm", "Cpc")) {
    solver <- qc:::.integration_make_solver(
      data = x,
      LSL = -3,
      USL = 3,
      target = 0,
      prior = prior,
      metric = metric
    )
    vals <- solver(c(0.25, 0.5, 1))
    testthat::expect_equal(length(vals), 3L, info = metric)
    testthat::expect_true(all(is.finite(vals)), info = metric)
    testthat::expect_true(all(diff(vals) <= 1e-8), info = metric)
  }
})
