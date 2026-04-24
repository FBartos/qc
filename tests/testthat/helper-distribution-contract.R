get_distribution_contract_fixture <- local({
  cache <- new.env(parent = emptyenv())

  function(key, builder) {
    if (!exists(key, envir = cache, inherits = FALSE)) {
      assign(key, builder(), envir = cache)
    }

    get(key, envir = cache, inherits = FALSE)
  }
})

distribution_contract_samples <- function(distribution) {
  switch(
    distribution,
    normal = 5 + 1.1 * stats::qnorm(seq(0.05, 0.95, length.out = 40)),
    t = 0.25 + 1.1 * stats::qt(seq(0.05, 0.95, length.out = 60), df = 8),
    stop(
      sprintf("No distribution-contract sample set is defined for `%s`.", distribution),
      call. = FALSE
    )
  )
}

distribution_contract_limits <- function(distribution) {
  switch(
    distribution,
    normal = list(LSL = 2, target = 5, USL = 8),
    t = list(LSL = -5, target = 0.25, USL = 6),
    stop(
      sprintf("No distribution-contract limits are defined for `%s`.", distribution),
      call. = FALSE
    )
  )
}

distribution_contract_summary_stats <- function(distribution) {
  x <- distribution_contract_samples(distribution)

  list(
    mean = mean(x),
    sd = stats::sd(x),
    N = length(x)
  )
}

distribution_contract_cached_state <- function(distribution) {
  x <- distribution_contract_samples(distribution)
  x_bar <- mean(x)

  list(
    n = length(x),
    x_bar = x_bar,
    sse = sum((x - x_bar)^2)
  )
}

distribution_contract_mcmc_args <- function(seed) {
  list(
    chains = 1,
    iter = 500,
    warmup = 200,
    cores = 1,
    silent = TRUE,
    seed = seed
  )
}

distribution_contract_nu_prior <- function() {
  qc::prior("exp", list(rate = 1 / 30), list(lower = 2, upper = Inf))
}

# Keep the supported matrix explicit so distribution/backend additions land here
# first and the test file stays mostly declarative.
distribution_contract_cases <- function(supported = NULL) {
  cases <- list(
    "pc-normal-raw" = list(
      label = "`pc()` supports `distribution = \"normal\"` raw fits",
      supported = TRUE,
      object_class = "pc",
      distribution = "normal",
      fit = function() {
        x <- distribution_contract_samples("normal")

        do.call(
          qc::pc,
          c(
            list(x = x),
            distribution_contract_limits("normal"),
            list(distribution = "normal", bootstrap = FALSE, seed = 1)
          )
        )
      },
      validate = function(fit) {
        testthat::expect_null(fit$metrics_boot)
      }
    ),
    "pc-t-raw" = list(
      label = "`pc()` supports `distribution = \"t\"` raw fits",
      supported = TRUE,
      object_class = "pc",
      distribution = "t",
      fit = function() {
        x <- distribution_contract_samples("t")

        suppressWarnings(
          do.call(
            qc::pc,
            c(
              list(x = x),
              distribution_contract_limits("t"),
              list(distribution = "t", bootstrap = FALSE, seed = 2)
            )
          )
        )
      },
      validate = function(fit) {
        testthat::expect_null(fit$metrics_boot)
        testthat::expect_true(all(is.finite(unlist(fit$fit$fit))))
        testthat::expect_gt(fit$fit$fit$scale, 0)
        testthat::expect_gt(fit$fit$fit$nu, 2)
      }
    ),
    "bpc-integration-normal-raw" = list(
      label = "`bpc(method = \"integration\")` supports `distribution = \"normal\"` raw fits",
      supported = TRUE,
      object_class = "bpc",
      distribution = "normal",
      method = "integration",
      fit = function() {
        x <- distribution_contract_samples("normal")

        do.call(
          qc::bpc,
          c(
            list(x = x),
            distribution_contract_limits("normal"),
            list(distribution = "normal", method = "integration")
          )
        )
      },
      validate = function(fit) {
        expected <- distribution_contract_cached_state("normal")
        actual <- fit$integration_result$cached_state

        testthat::expect_equal(actual$n, expected$n)
        testthat::expect_equal(actual$x_bar, expected$x_bar, tolerance = 1e-12)
        testthat::expect_equal(actual$sse, expected$sse, tolerance = 1e-12)
      }
    ),
    "bpc-integration-normal-summary" = list(
      label = "`bpc(method = \"integration\")` supports normal summary-statistics inputs",
      supported = TRUE,
      object_class = "bpc",
      distribution = "normal",
      method = "integration",
      fit = function() {
        summary_stats <- distribution_contract_summary_stats("normal")

        do.call(
          qc::bpc,
          c(
            list(x = NULL),
            distribution_contract_limits("normal"),
            summary_stats,
            list(distribution = "normal", method = "integration")
          )
        )
      },
      validate = function(fit) {
        expected <- distribution_contract_cached_state("normal")
        actual <- fit$integration_result$cached_state

        testthat::expect_equal(actual$n, expected$n)
        testthat::expect_equal(actual$x_bar, expected$x_bar, tolerance = 1e-12)
        testthat::expect_equal(actual$sse, expected$sse, tolerance = 1e-12)
      }
    ),
    "bpc-integration-t-raw" = list(
      label = "`bpc(method = \"integration\")` rejects `distribution = \"t\"`",
      supported = FALSE,
      error = "integration method currently only supports",
      fit = function() {
        x <- distribution_contract_samples("t")

        do.call(
          qc::bpc,
          c(
            list(x = x),
            distribution_contract_limits("t"),
            list(distribution = "t", method = "integration")
          )
        )
      }
    ),
    "bpc-mcmc-normal-raw" = list(
      label = "`bpc(method = \"mcmc\")` supports `distribution = \"normal\"` raw fits",
      supported = TRUE,
      object_class = "bpc",
      distribution = "normal",
      method = "mcmc",
      fit = function() {
        x <- distribution_contract_samples("normal")

        do.call(
          qc::bpc,
          c(
            list(x = x),
            distribution_contract_limits("normal"),
            list(distribution = "normal", method = "mcmc", prior = "Jeffreys"),
            distribution_contract_mcmc_args(seed = 3)
          )
        )
      },
      validate = function(fit) {
        testthat::expect_equal(fit$stan_data$is_ss, 0L)
      }
    ),
    "bpc-mcmc-normal-summary" = list(
      label = "`bpc(method = \"mcmc\")` supports normal summary-statistics inputs",
      supported = TRUE,
      object_class = "bpc",
      distribution = "normal",
      method = "mcmc",
      fit = function() {
        summary_stats <- distribution_contract_summary_stats("normal")

        do.call(
          qc::bpc,
          c(
            list(x = NULL),
            distribution_contract_limits("normal"),
            summary_stats,
            list(distribution = "normal", method = "mcmc", prior = "Jeffreys"),
            distribution_contract_mcmc_args(seed = 4)
          )
        )
      },
      validate = function(fit) {
        testthat::expect_equal(fit$stan_data$is_ss, 1L)
      }
    ),
    "bpc-mcmc-t-raw" = list(
      label = "`bpc(method = \"mcmc\")` supports `distribution = \"t\"` raw fits",
      supported = TRUE,
      object_class = "bpc",
      distribution = "t",
      method = "mcmc",
      fit = function() {
        x <- distribution_contract_samples("t")

        do.call(
          qc::bpc,
          c(
            list(x = x),
            distribution_contract_limits("t"),
            list(
              distribution = "t",
              method = "mcmc",
              prior = qc::prior_independent(nu = distribution_contract_nu_prior())
            ),
            distribution_contract_mcmc_args(seed = 5)
          )
        )
      },
      validate = function(fit) {
        testthat::expect_equal(fit$stan_data$is_ss, 0L)
      }
    ),
    "bpc-mcmc-t-summary" = list(
      label = "`bpc(method = \"mcmc\", distribution = \"t\")` rejects summary-statistics inputs",
      supported = FALSE,
      error = "Summary-statistics inputs .* `distribution = \"normal\"`",
      fit = function() {
        summary_stats <- distribution_contract_summary_stats("t")

        do.call(
          qc::bpc,
          c(
            list(x = NULL),
            distribution_contract_limits("t"),
            summary_stats,
            list(
              distribution = "t",
              method = "mcmc",
              prior = qc::prior_independent(nu = distribution_contract_nu_prior())
            ),
            distribution_contract_mcmc_args(seed = 6)
          )
        )
      }
    ),
    "bpc-mcmc-t-empty" = list(
      label = "`bpc(method = \"mcmc\", distribution = \"t\")` requires raw observations",
      supported = FALSE,
      error = "supply raw observations in `x`",
      fit = function() {
        do.call(
          qc::bpc,
          c(
            list(x = NULL),
            distribution_contract_limits("t"),
            list(
              distribution = "t",
              method = "mcmc",
              prior = qc::prior_independent(nu = distribution_contract_nu_prior())
            ),
            distribution_contract_mcmc_args(seed = 7)
          )
        )
      }
    )
  )

  if (!is.null(supported)) {
    cases <- cases[vapply(cases, function(case) identical(case$supported, supported), logical(1))]
  }

  ids <- names(cases)
  out <- Map(function(id, case) {
    case$id <- id
    case
  }, ids, unname(cases))
  names(out) <- ids
  out
}

distribution_contract_case <- function(id) {
  cases <- distribution_contract_cases()

  if (!id %in% names(cases)) {
    stop(sprintf("Unknown distribution-contract case `%s`.", id), call. = FALSE)
  }

  cases[[id]]
}

distribution_contract_fit <- function(case) {
  if (is.character(case)) {
    case <- distribution_contract_case(case)
  }

  get_distribution_contract_fixture(case$id, case$fit)
}

distribution_contract_pc_bootstrap_fit <- function() {
  get_distribution_contract_fixture("pc-normal-bootstrap", function() {
    x <- distribution_contract_samples("normal")

    do.call(
      qc::pc,
      c(
        list(x = x),
        distribution_contract_limits("normal"),
        list(distribution = "normal", samples = 20, seed = 8)
      )
    )
  })
}

expect_metric_distribution_cache <- function(metrics, info = NULL) {
  distributions <- attr(metrics, "distributions")

  testthat::expect_true(is.list(distributions), info = info)
  testthat::expect_equal(names(distributions), qc:::.qc_metric_names(), info = info)
  testthat::expect_true(
    all(vapply(distributions, inherits, logical(1), "qc_metric_distribution")),
    info = info
  )
}

expect_distribution_contract_support <- function(case) {
  if (is.character(case)) {
    case <- distribution_contract_case(case)
  }

  fit <- distribution_contract_fit(case)

  testthat::expect_s3_class(fit, case$object_class)
  testthat::expect_equal(fit$distribution, case$distribution, info = case$label)

  if (!is.null(case$method)) {
    testthat::expect_equal(fit$method, case$method, info = case$label)
  }

  testthat::expect_equal(names(fit$coefficients), qc:::.qc_metric_names(), info = case$label)
  testthat::expect_true(all(is.finite(unname(fit$coefficients))), info = case$label)
  expect_metric_distribution_cache(fit$metrics, info = case$label)

  if (!is.null(case$validate)) {
    case$validate(fit)
  }

  invisible(fit)
}

expect_distribution_contract_rejection <- function(case) {
  if (is.character(case)) {
    case <- distribution_contract_case(case)
  }

  testthat::expect_false(case$supported)
  testthat::expect_error(case$fit(), regexp = case$error, info = case$label)
}
