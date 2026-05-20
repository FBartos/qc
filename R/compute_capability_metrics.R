.validate_LSL_USL_target <- function(LSL, USL, target) {
  BayesTools::check_real(LSL,    name = "LSL",    check_length = 1, allow_NA = FALSE, upper = USL)
  BayesTools::check_real(USL,    name = "USL",    check_length = 1, allow_NA = FALSE, lower = LSL)
  BayesTools::check_real(target, name = "target", check_length = 1, allow_NA = FALSE)

  if (target < LSL || target > USL) {
    stop(
      sprintf(
        "`target` must lie within [`LSL`, `USL`] so Cpm/Cpc are well-defined. Received target = %s with LSL = %s and USL = %s.",
        format(target, digits = 10),
        format(LSL, digits = 10),
        format(USL, digits = 10)
      ),
      call. = FALSE
    )
  }

  return(list(LSL = LSL, USL = USL, target = target))
}

.validate_sigma_level <- function(sigma_level, name = "sigma") {
  BayesTools::check_real(
    sigma_level,
    name = name,
    check_length = 1,
    lower = 0,
    allow_NA = FALSE
  )

  if (sigma_level <= 0) {
    stop("`", name, "` must be greater than 0.", call. = FALSE)
  }

  sigma_level
}

.validate_metric_name <- function(metric, name = "metric") {
  BayesTools::check_char(
    metric,
    name = name,
    check_length = 1,
    allow_values = .qc_metric_names()
  )

  metric
}

.validate_capability_request <- function(LSL, USL, target = NULL,
                                         sigma_level = 3,
                                         metric = NULL,
                                         target_required = TRUE,
                                         sigma_name = "sigma") {
  if (!is.null(metric)) {
    .validate_metric_name(metric)
  }

  if (target_required || !is.null(target)) {
    limits <- .validate_LSL_USL_target(LSL = LSL, USL = USL, target = target)
  } else {
    BayesTools::check_real(LSL, name = "LSL", check_length = 1, allow_NA = FALSE, upper = USL)
    BayesTools::check_real(USL, name = "USL", check_length = 1, allow_NA = FALSE, lower = LSL)
    limits <- list(LSL = LSL, USL = USL, target = target)
  }

  limits$sigma_level <- .validate_sigma_level(sigma_level, name = sigma_name)
  limits$metric <- metric
  limits
}

extract_samples     <- function(fit, bootstrap) {
  UseMethod("extract_samples")
}
#' @export
extract_samples.bpc <- function(fit, bootstrap) {

  # Check if this is an integration-based fit

  if (!is.null(fit$method) && fit$method == "integration") {
    stop("extract_samples() is not supported for integration method. ",
         "Use summary() to obtain posterior statistics, or refit with method = 'mcmc'.")
  }

  samples <- rstan::extract(fit[["stanfit"]], pars = .bpc_parameters(fit[["distribution"]]))
  class(samples) <- fit[["distribution"]]
  return(samples)
}
#' @export
extract_samples.pc  <- function(fit, bootstrap) {

  samples <- if (bootstrap) fit$fit[["boot_fit"]] else fit$fit[["fit"]]
  class(samples) <- fit[["distribution"]]
  return(samples)
}

#' Extract posterior (or prior) predictive samples from a fitted model
#'
#' Works for both \code{"mcmc"} and \code{"integration"} methods.  For MCMC
#' fits the existing Stan samples are reused.  For integration fits with a
#' conjugate NIG prior the posterior parameters are computed analytically and
#' \code{n_samples} draws are generated from the NIG predictive.
#'
#' @param fit A fitted object of class \code{bpc}.
#' @param n_samples Integer.  Number of predictive samples to draw.  Only
#'   used for \code{method = "integration"}; MCMC fits return one sample per
#'   posterior draw.
#' @param ... Currently unused.
#' @return A numeric vector of predictive samples.
#' @export
extract_predictive_samples <- function(fit, n_samples = 10000L, ...) {
  UseMethod("extract_predictive_samples")
}

#' @export
extract_predictive_samples.bpc <- function(fit, n_samples = 10000L, ...) {
  if (!is.null(fit$method) && fit$method == "integration") {
    .extract_predictive_samples_integration(fit, n_samples = n_samples)
  } else {
    raw_samples <- extract_samples(fit, bootstrap = FALSE)
    samples_to_posterior_predictives(raw_samples)
  }
}

.extract_predictive_samples_integration <- function(fit, n_samples) {
  ir    <- fit$integration_result
  distribution <- fit$distribution %||% ir$distribution %||% "normal"

  if (!identical(distribution, "normal")) {
    stop(
      "Predictive sampling from integration fits is currently only supported for `distribution = \"normal\"`. ",
      "Refit with `method = \"mcmc\"` or add a distribution-specific integration predictive sampler.",
      call. = FALSE
    )
  }

  prior <- ir$prior

  if (!inherits(prior, "PriorConjugate")) {
    stop(
      "Predictive sampling from integration fits is only supported for conjugate (NIG) priors. ",
      "Refit with 'method = \"mcmc\"' to obtain predictive samples with non-conjugate priors."
    )
  }

  posterior_info <- .compute_validated_conjugate_posterior(
    prior,
    data = numeric(0),
    cached_state = ir$cached_state,
    context = "The conjugate posterior predictive distribution"
  )
  post <- posterior_info$post

  if (posterior_info$is_degenerate) {
    return(rep(post$mu_n, n_samples))
  }

  sigma2  <- 1 / stats::rgamma(n_samples, shape = post$alpha_n, rate = post$beta_n)
  sigma   <- sqrt(sigma2)
  mu      <- stats::rnorm(n_samples, mean = post$mu_n, sd = sigma / sqrt(post$k_n))

  samples <- list(mu = mu, sigma = sigma)
  class(samples) <- "normal"
  samples_to_posterior_predictives(samples)
}

samples_to_percentiles          <- function(samples, sigma) {
  UseMethod("samples_to_percentiles")
}
#' @export
samples_to_percentiles.normal   <- function(samples, sigma) {

  return(list(
    LP = samples[["mu"]] - sigma * samples[["sigma"]],
    MP = samples[["mu"]],
    UP = samples[["mu"]] + sigma * samples[["sigma"]]
  ))
}
#' @export
samples_to_percentiles.t        <- function(samples, sigma) {

  return(list(
    LP = samples[["mu"]] + stats::qt(stats::pnorm(-sigma), df = samples[["nu"]]) * samples[["scale"]],
    MP = samples[["mu"]],
    UP = samples[["mu"]] - stats::qt(stats::pnorm(-sigma), df = samples[["nu"]]) * samples[["scale"]]
  ))
}

samples_to_E_abs_dev            <- function(samples, target, ...) {
  UseMethod("samples_to_E_abs_dev")
}

.should_fallback_E_abs_dev <- function(x) {
  inherits(x, "try-error") || any(!is.finite(x))
}

.metric_sample_draw_count <- function(samples) {
  lengths <- vapply(samples, length, integer(1))
  unique_lengths <- unique(lengths)
  if (length(unique_lengths) != 1L) {
    stop(
      "All posterior sample components must have the same length to compute capability metrics.",
      call. = FALSE
    )
  }

  unique_lengths[[1]]
}

.replicate_sample_draw <- function(samples, draw_index, n_inner) {
  replicated <- lapply(samples, function(component) {
    rep(component[[draw_index]], n_inner)
  })
  class(replicated) <- class(samples)
  replicated
}

#' @export
samples_to_E_abs_dev.default    <- function(samples, target, n_inner = 2048L, ...) {
  n_inner <- as.integer(n_inner)
  if (!is.finite(n_inner) || n_inner < 2L) {
    stop("`n_inner` must be a finite integer greater than or equal to 2.", call. = FALSE)
  }

  n_draws <- .metric_sample_draw_count(samples)
  E_abs_dev <- numeric(n_draws)

  for (i in seq_len(n_draws)) {
    post_pred_samples <- samples_to_posterior_predictives(
      .replicate_sample_draw(samples, draw_index = i, n_inner = n_inner)
    )
    E_abs_dev[i] <- mean(abs(target - post_pred_samples))
  }

  return(E_abs_dev)
}
#' @export
samples_to_E_abs_dev.normal     <- function(samples, target, ...) {

  E_abs_dev <- try(with(
    samples,
    .cpc_E_abs_dev_normal(mu, sigma, target)
  ), silent = TRUE)

  if (.should_fallback_E_abs_dev(E_abs_dev)) {
    return(samples_to_E_abs_dev.default(samples, target))
  }

  return(E_abs_dev)
}
#' @export
samples_to_E_abs_dev.t          <- function(samples, target, ...) {

  E_abs_dev <- try(with(
    samples,
    {
      # Equation 2.7 of https://arxiv.org/abs/1912.01607v3
      # note that in their notation (e.g., Equation 2.2) they specify \sigma / \nu * (t - \mu)^2
      # hence inv_scale_sq.
      inv_scale_sq <- 1 / (scale * scale)
      z <- -(mu - target)^2 * inv_scale_sq / nu
      gauss2F1 <- gsl::hyperg_2F1(-1 / 2, nu / 2 - 1 / 2, 1 / 2, z)
      # gamma((1 + 1) / 2) == 1, so dropped
      sqrt(nu / inv_scale_sq) * gamma(nu / 2 - 1 / 2) / (sqrt(pi) * gamma(nu / 2)) * gauss2F1
    }
  ), silent = TRUE)

  if (.should_fallback_E_abs_dev(E_abs_dev)) {
    return(samples_to_E_abs_dev.default(samples, target))
  }

  return(E_abs_dev)
}

samples_to_posterior_predictives        <- function(samples) {
  UseMethod("samples_to_posterior_predictives")
}
#' @export
samples_to_posterior_predictives.normal <- function(samples) {

  return(stats::rnorm(length(samples[["mu"]]), mean = samples[["mu"]], sd = samples[["sigma"]]))
}
#' @export
samples_to_posterior_predictives.t      <- function(samples) {

  return(samples[["mu"]] + stats::rt(length(samples[["mu"]]), df = samples[["nu"]]) * samples[["scale"]])
}

.capability_metrics_from_percentiles <- function(percentiles, E_abs_dev,
                                                 LSL, USL, target,
                                                 sigma_level = 3) {
  range <- USL - LSL
  two_sigma <- sigma_level + sigma_level

  Cp  <- range / (percentiles$UP - percentiles$LP)
  Cpu <- (USL - percentiles$MP) / (percentiles$UP - percentiles$MP)
  Cpl <- (percentiles$MP - LSL) / (percentiles$MP - percentiles$LP)
  Cpk <- pmin(Cpu, Cpl)

  Cpm <- min(USL - target, target - LSL) /
    (sigma_level * sqrt(((percentiles$UP - percentiles$LP) / two_sigma)^2 +
                          (percentiles$MP - target)^2))

  Cpc <- range / (two_sigma * sqrt(pi / 2) * E_abs_dev)

  list(
    Cp  = Cp,
    Cpu = Cpu,
    Cpl = Cpl,
    Cpk = Cpk,
    Cpc = Cpc,
    Cpm = Cpm
  )
}

.compute_capability_metrics <- function(fit, LSL, USL, target, sigma = 3, bootstrap = FALSE) {

  # validate input
  .validate_capability_request(
    LSL = LSL,
    USL = USL,
    target = target,
    sigma_level = sigma,
    sigma_name = "sigma"
  )

  # extract samples from the fitted objects
  raw_samples <- extract_samples(fit, bootstrap = bootstrap)
  method <- if (!is.null(fit$method)) {
    fit$method
  } else if (inherits(fit, "pc")) {
    "pc"
  } else {
    "mcmc"
  }

  # computes capability metrics using percentiles
  # (i.e., the q-version of the metric as a generalization to non-normal distributions)
  samples <- samples_to_percentiles(raw_samples, sigma = sigma)
  E_abs_dev <- samples_to_E_abs_dev(raw_samples, target)
  lst <- .capability_metrics_from_percentiles(
    percentiles = samples,
    E_abs_dev = E_abs_dev,
    LSL = LSL,
    USL = USL,
    target = target,
    sigma_level = sigma
  )

  .new_capability_metrics(
    lst,
    LSL = LSL,
    USL = USL,
    target = target,
    sigma = sigma,
    method = method,
    distributions = .sample_distributions_from_metrics(lst, what = names(lst))
  )
}
