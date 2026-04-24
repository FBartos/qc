test_that("qc options round-trip and validate option names", {
  old_max_cores <- qc.get_option("max_cores")
  old_prior_mc_samples <- qc.get_option("prior_mc_samples")
  on.exit(qc.options(max_cores = old_max_cores, prior_mc_samples = old_prior_mc_samples), add = TRUE)

  current <- qc.options()
  expect_type(current, "list")
  expect_true("max_cores" %in% names(current))
  expect_true("prior_mc_samples" %in% names(current))
  expect_equal(current$max_cores, old_max_cores)

  updated_max_cores <- max(1L, as.integer(old_max_cores) - 1L)
  updated <- qc.options(max_cores = updated_max_cores)
  expect_equal(updated$max_cores, updated_max_cores)
  expect_equal(qc.get_option("max_cores"), updated_max_cores)
  expect_no_error(qc.options(prior_mc_samples = 1000L))

  expect_error(qc.options(unknown = 1), "Unmatched or ambiguous option")
  expect_error(qc.options(max_cores = 0L), "max_cores")
  expect_error(qc.options(prior_mc_samples = 999L), "prior_mc_samples")
  expect_error(qc.get_option(c("max_cores", "other")), "Only 1 option")
  expect_error(qc.get_option("unknown"), "Unmatched or ambiguous option")
})

test_that("Stan control helpers validate values and supply defaults", {
  expect_equal(
    set_control(),
    list(adapt_delta = 0.8, max_treedepth = 15L)
  )

  expect_error(set_control(adapt_delta = 1.2), "adapt_delta")
  expect_error(set_control(max_treedepth = 0), "max_treedepth")
})

test_that("Stan fit settings enforce iter greater than warmup and cap parallel cores", {
  old_max_cores <- qc.get_option("max_cores")
  on.exit(qc.options(max_cores = old_max_cores), add = TRUE)
  qc.options(max_cores = 2L)

  capped <- qc:::.stan_check_and_list_fit_settings(
    chains = 4L,
    warmup = 100L,
    iter = 250L,
    thin = 1L,
    parallel = TRUE,
    cores = 8L,
    silent = TRUE,
    seed = 11L,
    control = list(
      adapt_delta = NULL,
      max_treedepth = NULL
    )
  )

  expect_equal(capped$cores, 2L)
  expect_equal(capped$adapt_delta, 0.8)
  expect_equal(capped$max_treedepth, 15L)

  serial <- qc:::.stan_check_and_list_fit_settings(
    chains = 2L,
    warmup = 100L,
    iter = 250L,
    thin = 1L,
    parallel = FALSE,
    cores = 8L,
    silent = TRUE,
    seed = 11L,
    control = list(adapt_delta = 0.9)
  )

  expect_equal(serial$cores, 1L)
  expect_equal(serial$adapt_delta, 0.9)
  expect_equal(serial$max_treedepth, 15L)

  expect_error(
    qc:::.stan_check_and_list_fit_settings(
      chains = 1L,
      warmup = 100L,
      iter = 100L,
      thin = 1L,
      parallel = FALSE,
      cores = 1L,
      silent = TRUE,
      seed = 1L,
      control = list(
        adapt_delta = NULL,
        max_treedepth = NULL
      )
    ),
    "iter"
  )
})

test_that("optimization fit settings use qc max_cores for parallel fits", {
  old_max_cores <- qc.get_option("max_cores")
  on.exit(qc.options(max_cores = old_max_cores), add = TRUE)
  qc.options(max_cores = 3L)

  parallel_cfg <- qc:::.optim_check_and_list_fit_settings(
    bootstrap = TRUE,
    samples = 10L,
    parallel = TRUE,
    cores = NULL,
    seed = 7L,
    control = NULL
  )
  serial_cfg <- qc:::.optim_check_and_list_fit_settings(
    bootstrap = TRUE,
    samples = 10L,
    parallel = FALSE,
    cores = 4L,
    seed = 7L,
    control = NULL
  )

  expect_equal(parallel_cfg$cores, 3L)
  expect_equal(serial_cfg$cores, 1L)
})

test_that("Stan prior encoding handles Jeffreys, point, and truncated priors", {
  jeffreys_mu <- qc:::.stan_distribution("mu", "Jeffreys_mu", sample_priors = FALSE)
  expect_equal(jeffreys_mu$estimate_mu, 1)
  expect_equal(jeffreys_mu$prior_type_mu, 98)
  expect_equal(jeffreys_mu$bounds_mu, c(999, 999))
  expect_equal(jeffreys_mu$bounds_type_mu, c(0, 0))
  expect_equal(jeffreys_mu$prior_parameters_mu, c(999, 999, 999))

  expect_error(
    qc:::.stan_distribution("mu", "Jeffreys_mu", sample_priors = TRUE),
    "Improper prior distributions cannot be sampled from"
  )

  point_mu <- qc:::.stan_distribution(
    "mu",
    prior("point", list(location = 1.5)),
    sample_priors = FALSE
  )
  expect_equal(point_mu$estimate_mu, 0)
  expect_equal(drop(point_mu$fixed_mu), 1.5)
  expect_equal(point_mu$prior_type_mu, 0)
  expect_length(point_mu$bounds_mu, 0L)
  expect_length(point_mu$prior_parameters_mu, 0L)

  truncated_sigma <- qc:::.stan_distribution(
    "sigma",
    prior("gamma", list(shape = 2, rate = 0.5), list(lower = 0.25, upper = 4)),
    sample_priors = FALSE
  )
  expect_equal(truncated_sigma$estimate_sigma, 1)
  expect_equal(truncated_sigma$prior_type_sigma, 5)
  expect_equal(truncated_sigma$bounds_sigma, c(0.25, 4))
  expect_equal(truncated_sigma$bounds_type_sigma, c(1, 1))
  expect_equal(truncated_sigma$prior_parameters_sigma, c(2, 0.5, 999))
})

test_that("Stan prior parameter encoding matches the supported prior families", {
  prior_cases <- list(
    normal = list(prior = prior("normal", list(mean = 0, sd = 1)), expected = c(0, 1, 999)),
    lognormal = list(prior = prior("lognormal", list(meanlog = 0, sdlog = 1)), expected = c(0, 1, 999)),
    t = list(prior = prior("t", list(df = 5, location = 0, scale = 2)), expected = c(5, 0, 2)),
    gamma = list(prior = prior("gamma", list(shape = 2, rate = 1)), expected = c(2, 1, 999)),
    invgamma = list(prior = prior("invgamma", list(shape = 3, scale = 2)), expected = c(3, 2, 999)),
    uniform = list(prior = prior("uniform", list(a = 0, b = 10)), expected = c(0, 10, 999)),
    beta = list(prior = prior("beta", list(alpha = 2, beta = 4)), expected = c(2, 4, 999)),
    exp = list(prior = prior("exp", list(rate = 1.5)), expected = c(1.5, 999, 999))
  )

  for (case_name in names(prior_cases)) {
    case <- prior_cases[[case_name]]
    expect_equal(
      qc:::.stan_distribution_parameters(case$prior),
      case$expected,
      info = case_name
    )
  }
})

test_that("E_abs_dev helpers fall back to sampling on errors and non-finite outputs", {
  normal_samples <- structure(
    list(mu = c(0, 1), sigma = c(1, 2)),
    class = "normal"
  )
  t_samples <- structure(
    list(mu = c(0, 1), scale = c(1, 1.5), nu = c(6, 8)),
    class = "t"
  )

  testthat::local_mocked_bindings(
    pnorm = function(...) {
      args <- list(...)
      rep(NaN, length(args[[1L]]))
    },
    .package = "stats"
  )
  expect_silent(
    normal_fallback <- qc:::samples_to_E_abs_dev.normal(normal_samples, target = 0)
  )
  expect_length(normal_fallback, 2L)
  expect_true(all(is.finite(normal_fallback)))

  testthat::local_mocked_bindings(
    hyperg_2F1 = function(...) stop("boom"),
    .package = "gsl"
  )
  expect_silent(
    t_fallback <- qc:::samples_to_E_abs_dev.t(t_samples, target = 0)
  )
  expect_length(t_fallback, 2L)
  expect_true(all(is.finite(t_fallback)))
})
