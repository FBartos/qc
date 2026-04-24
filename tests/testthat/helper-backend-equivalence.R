get_backend_test_fixture <- local({
  cache <- new.env(parent = emptyenv())

  function(key, builder) {
    if (!exists(key, envir = cache, inherits = FALSE)) {
      assign(key, builder(), envir = cache)
    }

    get(key, envir = cache, inherits = FALSE)
  }
})

fit_backend_pair <- function(...,
                             mcmc_chains = 1,
                             mcmc_iter = 1500,
                             mcmc_warmup = 500,
                             seed = 1) {
  common_args <- list(...)
  if (is.null(common_args$prior)) {
    common_args$prior <- "Jeffreys"
  }

  fit_int <- do.call(
    qc::bpc,
    c(common_args, list(method = "integration"))
  )
  fit_mcmc <- do.call(
    qc::bpc,
    c(
      common_args,
      list(
        method = "mcmc",
        chains = mcmc_chains,
        iter = mcmc_iter,
        warmup = mcmc_warmup,
        cores = 1,
        silent = TRUE,
        seed = seed
      )
    )
  )

  list(integration = fit_int, mcmc = fit_mcmc)
}

summary_metric_rows <- function(summary_tbl, metrics) {
  rows <- summary_tbl[match(metrics, summary_tbl$metric), , drop = FALSE]
  rownames(rows) <- NULL
  rows
}

expect_named_numeric_close <- function(actual,
                                       expected,
                                       abs_tol,
                                       rel_tol = 0,
                                       context = "") {
  stopifnot(length(actual) == length(expected))

  for (i in seq_along(actual)) {
    tol <- max(abs_tol, rel_tol * max(abs(expected[[i]]), 1))
    ok <- abs(actual[[i]] - expected[[i]]) <= tol
    label <- names(actual)[i]
    if (is.na(label) || !nzchar(label)) {
      label <- paste0("[", i, "]")
    }

    testthat::expect_true(
      ok,
      info = sprintf(
        "%s%s actual=%.6f expected=%.6f tol=%.6f",
        if (nzchar(context)) paste0(context, " ") else "",
        label,
        actual[[i]],
        expected[[i]],
        tol
      )
    )
  }
}

expect_summary_columns_close <- function(summary_int,
                                         summary_mcmc,
                                         metrics,
                                         columns,
                                         abs_tol,
                                         rel_tol = 0) {
  rows_int <- summary_metric_rows(summary_int$summary, metrics)
  rows_mcmc <- summary_metric_rows(summary_mcmc$summary, metrics)

  for (column in columns) {
    actual <- stats::setNames(rows_int[[column]], rows_int$metric)
    expected <- stats::setNames(rows_mcmc[[column]], rows_mcmc$metric)

    expect_named_numeric_close(
      actual = actual,
      expected = expected,
      abs_tol = abs_tol,
      rel_tol = rel_tol,
      context = sprintf("summary column `%s`", column)
    )
  }
}

expect_interval_summary_close <- function(summary_int,
                                          summary_mcmc,
                                          metrics,
                                          abs_tol) {
  rows_int <- summary_metric_rows(summary_int$interval_summary, metrics)
  rows_mcmc <- summary_metric_rows(summary_mcmc$interval_summary, metrics)
  probability_columns <- setdiff(names(rows_int), "metric")

  for (column in probability_columns) {
    actual <- stats::setNames(rows_int[[column]], rows_int$metric)
    expected <- stats::setNames(rows_mcmc[[column]], rows_mcmc$metric)

    expect_named_numeric_close(
      actual = actual,
      expected = expected,
      abs_tol = abs_tol,
      context = sprintf("interval column `%s`", column)
    )
  }
}

expect_fit_summaries_close <- function(pair,
                                       metrics,
                                       columns = c("mean", "median", "lower", "upper"),
                                       abs_tol = 0.05,
                                       rel_tol = 0,
                                       ci_level = 0.95,
                                       interval_probability = c(1.00, 1.33, 1.50, 2.00),
                                       sigma = NULL) {
  summary_args <- list(
    ci.level = ci_level,
    interval_probability = interval_probability
  )
  if (!is.null(sigma)) {
    summary_args$sigma <- sigma
  }

  sum_int <- do.call(summary, c(list(object = pair$integration), summary_args))
  sum_mcmc <- do.call(summary, c(list(object = pair$mcmc), summary_args))

  expect_summary_columns_close(
    summary_int = sum_int,
    summary_mcmc = sum_mcmc,
    metrics = metrics,
    columns = columns,
    abs_tol = abs_tol,
    rel_tol = rel_tol
  )

  list(integration = sum_int, mcmc = sum_mcmc)
}

expect_predictive_samples_close <- function(samples_int,
                                            samples_mcmc,
                                            mean_sd_factor = 0.05,
                                            quantile_sd_factor = 0.10) {
  predictive_scale <- max(stats::sd(samples_int), stats::sd(samples_mcmc), .Machine$double.eps)

  testthat::expect_lte(abs(mean(samples_int) - mean(samples_mcmc)),
                       mean_sd_factor * predictive_scale)
  testthat::expect_lte(abs(stats::sd(samples_int) - stats::sd(samples_mcmc)),
                       mean_sd_factor * predictive_scale)

  for (prob in c(0.025, 0.5, 0.975)) {
    q_int <- unname(stats::quantile(samples_int, prob))
    q_mcmc <- unname(stats::quantile(samples_mcmc, prob))
    testthat::expect_true(
      abs(q_int - q_mcmc) <= quantile_sd_factor * predictive_scale,
      info = sprintf("predictive quantile %.3f", prob)
    )
  }
}
