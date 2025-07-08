.validate_LSL_USL_target <- function(LSL, USL, target) {
  BayesTools::check_real(LSL,    name = "LSL",    check_length = 1, allow_NA = FALSE, upper = USL)
  BayesTools::check_real(USL,    name = "USL",    check_length = 1, allow_NA = FALSE, lower = LSL)
  BayesTools::check_real(target, name = "target", check_length = 1, allow_NA = FALSE)

  return(list(LSL = LSL, USL = USL, target = target))
}

extract_mu_and_sigma            <- function(fit) {
  UseMethod("extract_mu_and_sigma")
}

#' @export
extract_mu_and_sigma.bpc_normal <- function(fit) {
  rstan::extract(fit$stanfit, pars = c("mu", "sigma"))
}

#' @export
extract_mu_and_sigma.bpc_t      <- function(fit) {

  all <- rstan::extract(fit$stanfit, pars = c("mu", "scale", "nu"))
  return(with(all, {
    sigma <- scale * sqrt(nu / (nu - 2.0))
    return(list(mu = mu, sigma = sigma))
  }))
}

extract_percentiles            <- function(fit, sigma) {
  UseMethod("extract_percentiles")
}

#' @export
extract_percentiles.bpc_normal <- function(fit, sigma) {

  samples <- rstan::extract(fit$stanfit, pars = c("mu", "sigma"))

  # for normal we can compute the percentiles directly
  return(list(
    LP = samples$mu - sigma * samples$sigma,
    MP = samples$mu,
    UP = samples$mu + sigma * samples$sigma
  ))
}

#' @export
extract_percentiles.bpc_t      <- function(fit, sigma) {

  samples <- rstan::extract(fit$stanfit, pars = c("mu", "scale", "nu"))

  return(list(
    LP = samples$mu + stats::qt(stats::pnorm(-sigma), df = samples$nu) * samples$scale,
    MP = samples$mu,
    UP = samples$mu - stats::qt(stats::pnorm(-sigma), df = samples$nu) * samples$scale
  ))
}

.bpc_posterior_predictive <- function(fit) {
  UseMethod(".bpc_posterior_predictive")
}

#' @export
.bpc_posterior_predictive.bpc_normal <- function(fit) {
  samples <- rstan::extract(fit$stanfit, pars = c("mu", "sigma"))
  return(stats::rnorm(length(samples$mu), mean = samples$mu, sd = samples$sigma))
}

#' @export
.bpc_posterior_predictive.bpc_t <- function(fit) {
  samples <- rstan::extract(fit$stanfit, pars = c("mu", "scale", "nu"))
  return(samples$mu + stats::rt(length(samples$mu), df = samples$nu) * samples$scale)
}

.bpc_posterior_predictive <- function(fit) {
  UseMethod(".bpc_posterior_predictive")
}

#' @export
.bpc_posterior_predictive.bpc_normal <- function(fit) {
  samples <- rstan::extract(fit$stanfit, pars = c("mu", "sigma"))
  return(stats::rnorm(length(samples$mu), mean = samples$mu, sd = samples$sigma))
}

#' @export
.bpc_posterior_predictive.bpc_t <- function(fit) {
  samples <- rstan::extract(fit$stanfit, pars = c("mu", "scale", "nu"))
  return(samples$mu + stats::rt(length(samples$mu), df = samples$nu) * samples$scale)
}

.bpc_compute_E_abs_dev <- function(fit, target) {
  UseMethod(".bpc_compute_E_abs_dev")
}

#' @export
.bpc_compute_E_abs_dev.default <- function(fit, target) {
  # a slow but accurate way to compute E_abs_dev
  E_abs_dev <- numeric(with(fit$control, chains * (iter - warmup)))
  for (i in seq_along(E_abs_dev)) {
    post_pred_samples <- .bpc_posterior_predictive(fit)
    E_abs_dev[i] <- mean(abs(target - post_pred_samples))
  }
  return(E_abs_dev)
}

.bpc_compute_E_abs_dev.bpc_normal <- function(fit, target) {

  samples <- rstan::extract(fit$stanfit, pars = c("mu", "sigma"))
  E_abs_dev <- try(with(
    samples,
    {
      # Equation 17 of https://arxiv.org/abs/1209.4340
      # z <- -(mu - target)^2 / (2 * sigma^2)
      # sigma * sqrt(2 / pi) * gsl::hyperg_1F1(-1 / 2, 1 / 2, z)
      # simplification of hyperg_1F1
      z      <- (mu - target) / sigma
      abs_mu_minus_T <- abs(mu - target)
      sigma * sqrt(2 / pi) * exp(-0.5 * z^2) + abs_mu_minus_T * (1 - 2 * pnorm(-abs(z)))
    }
  ))

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
    return(.bpc_compute_E_abs_dev.default(fit, target))
  }

  return(E_abs_dev)
}

#' @export
.bpc_compute_E_abs_dev.bpc_t <- function(fit, target) {

  samples <- rstan::extract(fit$stanfit, pars = c("mu", "scale", "nu"))
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
    return(.bpc_compute_E_abs_dev.default(fit, target))
  }

  return(E_abs_dev)
}


.bpc_compute_capability_metrics <- function(fit, LSL, USL, target, sigma = 3, force_normal = FALSE) {

  # validate input
  .validate_LSL_USL_target(LSL, USL, target)

  range <- USL - LSL

  one_sigma <- sigma
  two_sigma <- sigma + sigma

  # compute mean and standard deviation if normal distributions calculation is required
  if (force_normal) {

    # computes standard capability metrics assuming normal distribution
    samples <- extract_mu_and_sigma(fit)

    three_sigma <- one_sigma * samples$sigma
    six_sigma   <- two_sigma * samples$sigma

    Cp  <- range / six_sigma
    CpU <- (USL - samples$mu) / three_sigma
    CpL <- (samples$mu - LSL) / three_sigma
    Cpk <- pmin(CpU, CpL)

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
    samples <- extract_percentiles(fit, sigma = one_sigma)

    Cp  <- range / (samples$UP - samples$LP)
    CpU <- (USL - samples$MP) / (samples$UP - samples$MP)
    CpL <- (samples$MP - LSL) / (samples$MP - samples$LP)
    Cpk <- pmin(CpU, CpL)

    Cpm <- min(USL - target, target - LSL) / (one_sigma * sqrt( ( (samples$UP - samples$LP) / two_sigma)^2  + (samples$MP - target)^2) )

    E_abs_dev <- .bpc_compute_E_abs_dev(fit, target)

  }

  Cpc <- range / (two_sigma * sqrt(pi /  2) * E_abs_dev)

  lst <- list(
    Cp  = Cp,
    CpU = CpU,
    CpL = CpL,
    Cpk = Cpk,
    Cpc = Cpc,
    Cpm = Cpm
  )

  class(lst) <- "bpc_capability_metrics"
  attr(lst, "LSL")    <- LSL
  attr(lst, "USL")    <- USL
  attr(lst, "target") <- target
  return(lst)
}
