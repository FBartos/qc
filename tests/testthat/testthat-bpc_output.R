testthat::test_that("compute metrics", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(x, chains = 1, warmup = 100, iter = 200, silent = TRUE, seed = 1)

  metrics <- compute_capability_metrics(fit)

  summarize_capability_metrics(metrics)
  summarize_capability_metrics(fit)

  fit <- bpc(x, chains = 1, warmup = 1000, iter = 2000, silent = TRUE, seed = 1,
             prior_nu = prior("exp",  list(1)))

  qc:::extract_mu_and_sigma(fit)

  metrics <- compute_capability_metrics(fit)

  summarize_capability_metrics(metrics)
  summarize_capability_metrics(fit)

})

testthat::test_that("plot posterior density of metrics", {

  set.seed(1)
  x <- rnorm(100, 10, 2)

  fit <- bpc(x, chains = 1, warmup = 100, iter = 200, silent = TRUE, seed = 1)

  plot_density(fit)
  plot_density(fit, what = "Cp")

  fit_t <- bpc(x, chains = 1, warmup = 1000, iter = 2000, silent = TRUE, seed = 1,
             prior_nu = prior("exp",  list(1)))

  plot_density(fit_t)
  plot_density(fit_t, what = "Cp")

})

