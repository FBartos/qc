.validate_LSL_USL_target <- function(LSL, USL, target) {
  BayesTools::check_real(LSL,    name = "LSL",    check_length = 1, allow_NA = FALSE, upper = USL)
  BayesTools::check_real(USL,    name = "USL",    check_length = 1, allow_NA = FALSE, lower = LSL)
  BayesTools::check_real(target, name = "target", check_length = 1, allow_NA = FALSE)

  return(list(LSL = LSL, USL = USL, target = target))
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
    samples     <- samples_to_mu_and_sigma(raw_samples)
    samples_to_posterior_predictives(samples)
  }
}

.extract_predictive_samples_integration <- function(fit, n_samples) {
  ir    <- fit$integration_result
  prior <- ir$prior

  if (!inherits(prior, "PriorConjugate")) {
    stop(
      "Predictive sampling from integration fits is only supported for conjugate (NIG) priors. ",
      "Refit with 'method = \"mcmc\"' to obtain predictive samples with non-conjugate priors."
    )
  }

  cs   <- ir$cached_state
  post <- .nig_posterior(prior, n = cs$n, x_bar = cs$x_bar, SS = cs$sse)

  if (.is_improper_conjugate_posterior(post$k_n, post$alpha_n, post$beta_n)) {
    stop(
      "Predictive sampling is undefined for an improper conjugate prior/posterior. ",
      "Use a proper PriorConjugate for prior-only sampling or refit with data that yields a proper posterior.",
      call. = FALSE
    )
  }

  if (.is_degenerate_conjugate_posterior(post$beta_n)) {
    return(rep(post$mu_n, n_samples))
  }

  sigma2  <- 1 / stats::rgamma(n_samples, shape = post$alpha_n, rate = post$beta_n)
  sigma   <- sqrt(sigma2)
  mu      <- stats::rnorm(n_samples, mean = post$mu_n, sd = sigma / sqrt(post$k_n))

  samples <- list(mu = mu, sigma = sigma)
  class(samples) <- "normal"
  samples_to_posterior_predictives(samples)
}

samples_to_mu_and_sigma         <- function(samples) {
  UseMethod("samples_to_mu_and_sigma")
}
#' @export
samples_to_mu_and_sigma.normal  <- function(samples) {
  return(samples)
}
#' @export
samples_to_mu_and_sigma.t       <- function(samples) {
  return(list(
    mu    = samples[["mu"]],
    sigma = samples[["scale"]] * sqrt(samples[["nu"]] / (samples[["nu"]] - 2.0))
  ))
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

samples_to_E_abs_dev            <- function(samples, target) {
  UseMethod("samples_to_E_abs_dev")
}
#' @export
samples_to_E_abs_dev.default    <- function(samples, target) {

  # a slow but accurate way to compute E_abs_dev
  E_abs_dev <- numeric(length(samples[[1]]))
  for (i in seq_along(E_abs_dev)) {
    post_pred_samples <- samples_to_posterior_predictives(samples)
    E_abs_dev[i]      <- mean(abs(target - post_pred_samples))
  }

  return(E_abs_dev)
}
#' @export
samples_to_E_abs_dev.normal     <- function(samples, target) {

  E_abs_dev <- try(with(
    samples,
    {
      # Equation 17 of https://arxiv.org/abs/1209.4340
      # z <- -(mu - target)^2 / (2 * sigma^2)
      # sigma * sqrt(2 / pi) * gsl::hyperg_1F1(-1 / 2, 1 / 2, z)
      # simplification of hyperg_1F1
      z              <- (mu - target) / sigma
      abs_mu_minus_T <- abs(mu - target)
      sigma * sqrt(2 / pi) * exp(-0.5 * z^2) + abs_mu_minus_T * (1 - 2 * stats::pnorm(-abs(z)))
    }
  ))

  # TODO: Don (this is your previous comment from this function)
  # TODO: I'm confused about which of these we want, average across rows or columns?
  # E_abs_dev2 <- try(with(samples, {
  #   z      <- (mu - target) / sigma
  #   abs_mu_minus_T <- abs(mu - target)
  #
  #   sigma * sqrt(2 / pi) * exp(-0.5 * z^2) + abs_mu_minus_T * (1 - 2 * pnorm(-abs(z)))
  #   # z <- (mu - T) / sigma
  #   # a <- sigma * sqrt(2/pi) * exp(-0.5 * z^2) + abs(mu - T) * (1 - 2 * pnorm(-abs(z)))
  # }))
  #
  # E_abs_dev_ref  <- .bpc_compute_E_abs_dev.default(fit, target)
  # E_abs_dev_ref2 <- .bpc_compute_E_abs_dev.default(fit, target)
  #
  # par(mfrow = c(1, 3))
  # plot(density(E_abs_dev), main = "what we have")
  # lines(density(E_abs_dev_ref), col = "red")
  # plot(density(E_abs_dev2), col = "blue", main = "what we want")
  # lines(density(E_abs_dev_ref2), col = "green")
  # plot(density(E_abs_dev), main = "what we have")
  # lines(density(E_abs_dev2), col = "red")
  #
  # hh <- seq(.01, .99, .01)
  # plot(quantile(E_abs_dev, probs = hh), quantile(E_abs_dev_ref, probs = hh)); abline(0, 1)
  # plot(quantile(E_abs_dev2, probs = hh), quantile(E_abs_dev_ref, probs = hh)); abline(0, 1)
  # lm(E_abs_dev ~ E_abs_dev_ref)
  # plot(quantile((E_abs_dev - mean(E_abs_dev)) / (sd(E_abs_dev) * sqrt(length(E_abs_dev))) + mean(E_abs_dev), probs = hh), quantile(E_abs_dev_ref, probs = hh)); abline(0, 1)


  # E_abs_dev_mat <- matrix(nrow = with(fit$control, chains * (iter - warmup)), ncol = length(E_abs_dev))
  # for (i in seq_along(E_abs_dev)) {
  #   E_abs_dev_mat[, i] <- stats::rnorm(length(samples$mu), mean = samples$mu, sd = samples$sigma)
  # }
  # E_abs_dev_mat_abs    <- abs(target - E_abs_dev_mat)
  # E_abs_dev_mat_abs_rm <- colMeans(E_abs_dev_mat_abs)
  # E_abs_dev_mat_abs_cm <- rowMeans(E_abs_dev_mat_abs)
  #
  # this is what we've done analytically for the normal and the t...
  # plot(E_abs_dev, E_abs_dev_mat_abs_cm)
  # this is a Rao-blackwellized estimate that agrees in mean, but I'm not so sure about any other statistics...
  # but we want/ need E_abs_dev_mat_abs_rm!
  # is there a trick we can use?

  # hh <- seq(.01, .99, .01)
  # plot(quantile(E_abs_dev, probs = hh), quantile(E_abs_dev_mat_abs_cm, probs = hh)); abline(0, 1)
  # lm(E_abs_dev ~ E_abs_dev_ref)
  # plot(quantile((E_abs_dev - mean(E_abs_dev)) / (sd(E_abs_dev) * sqrt(length(E_abs_dev))) + mean(E_abs_dev), probs = hh), quantile(E_abs_dev_ref, probs = hh)); abline(0, 1)

  # E_abs_dev_ref <- .bpc_compute_E_abs_dev.default(fit, target)
  # par(mfrow = c(1, 2))
  # plot(density(E_abs_dev), main = "what we have")
  # lines(density(E_abs_dev_mat_abs_cm), col = "red")
  # plot(density(E_abs_dev_ref), col = "blue", main = "what we want")
  # lines(density(E_abs_dev_mat_abs_rm), col = "green")

  if (inherits(E_abs_dev, "try-error")) {
    # if the hypergeometric function fails, we fall back to the slow method
    warning("Failed to compute E_abs_dev using hypergeometric function, using the slow sampling based method as a fallback.")
    return(samples_to_E_abs_dev.default(fit, target))
  }

  return(E_abs_dev)
}
#' @export
samples_to_E_abs_dev.t          <- function(samples, target) {

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
  ))

  if (inherits(E_abs_dev, "try-error")) {
    # if the hypergeometric function fails, we fall back to the slow method
    warning("Failed to compute E_abs_dev using hypergeometric function, using the slow sampling based method as a fallback.")
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

.compute_capability_metrics <- function(fit, LSL, USL, target, sigma = 3, force_normal = FALSE, bootstrap = FALSE) {

  # validate input
  .validate_LSL_USL_target(LSL, USL, target)

  # precomputed settings
  range     <- USL - LSL
  one_sigma <- sigma
  two_sigma <- sigma + sigma

  # extract samples from the fitted objects
  raw_samples <- extract_samples(fit, bootstrap = bootstrap)

  # compute mean and standard deviation if normal distributions calculation is required
  if (force_normal) {

    # computes standard capability metrics assuming normal distribution
    samples <- samples_to_mu_and_sigma(raw_samples)

    three_sigma <- one_sigma * samples$sigma
    six_sigma   <- two_sigma * samples$sigma

    Cp  <- range / six_sigma
    Cpu <- (USL - samples$mu) / three_sigma
    Cpl <- (samples$mu - LSL) / three_sigma
    Cpk <- pmin(Cpu, Cpl)

    # Eq. 8.14 of Montgomery, 8th edition
    xi <- (samples$mu - target) / (samples$sigma)
    # Eq. 8.13 of Montgomery, 8th edition
    Cpm <- Cp / sqrt(1 + xi^2)

    delta <- (target - samples$mu) / samples$sigma
    E_abs_dev_samples <- samples$sigma * (sqrt(2 / pi) * dnorm(delta) + abs(delta) * (1 - 2 * pnorm(-abs(delta))))
    E_abs_dev <- mean(E_abs_dev_samples)

  } else {

    # computes capability metrics using percentiles
    # (i.e., (q) version of the metric = generalization to non-normal distributions)
    samples <- samples_to_percentiles(raw_samples, sigma = one_sigma)

    Cp  <- range / (samples$UP - samples$LP)
    Cpu <- (USL - samples$MP) / (samples$UP - samples$MP)
    Cpl <- (samples$MP - LSL) / (samples$MP - samples$LP)
    Cpk <- pmin(Cpu, Cpl)

    Cpm <- min(USL - target, target - LSL) / (one_sigma * sqrt( ( (samples$UP - samples$LP) / two_sigma)^2  + (samples$MP - target)^2) )

    E_abs_dev <- samples_to_E_abs_dev(raw_samples, target)

  }

  Cpc <- range / (two_sigma * sqrt(pi /  2) * E_abs_dev)

  lst <- list(
    Cp  = Cp,
    Cpu = Cpu,
    Cpl = Cpl,
    Cpk = Cpk,
    Cpc = Cpc,
    Cpm = Cpm
  )

  class(lst) <- "capability_metrics"
  attr(lst, "LSL")    <- LSL
  attr(lst, "USL")    <- USL
  attr(lst, "target") <- target
  attr(lst, "method") <- "mcmc"
  return(lst)
}
