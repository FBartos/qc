.validate_LSL_USL_target <- function(LSL, USL, target) {
  BayesTools::check_real(LSL,    name = "LSL",    check_length = 1, allow_NA = FALSE, upper = USL)
  BayesTools::check_real(USL,    name = "USL",    check_length = 1, allow_NA = FALSE, lower = LSL)
  BayesTools::check_real(target, name = "target", check_length = 1, allow_NA = FALSE)

  return(list(LSL = LSL, USL = USL, target = target))
}

extract_mu_and_sigma            <- function(fit) {
  UseMethod("extract_mu_and_sigma")
}
extract_mu_and_sigma.bpc_normal <- function(fit) {
  rstan::extract(fit$stanfit, pars = c("mu", "sigma"))
}
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
extract_percentiles.bpc_normal <- function(fit, sigma) {

  samples <- rstan::extract(fit$stanfit, pars = c("mu", "sigma"))

  # for normal we can compute the percentiles directly
  return(list(
    LP = samples$mu - sigma * samples$sigma,
    MP = samples$mu,
    UP = samples$mu + sigma * samples$sigma
  ))
}
extract_percentiles.bpc_t      <- function(fit, sigma) {

  samples <- rstan::extract(fit$stanfit, pars = c("mu", "scale", "nu"))

  return(list(
    LP = samples$mu + stats::qt(stats::pnorm(-sigma), df = samples$nu) * samples$scale,
    MP = samples$mu,
    UP = samples$mu - stats::qt(stats::pnorm(-sigma), df = samples$nu) * samples$scale
  ))
}

.bpc_compute_capability_metrics <- function(fit, LSL, USL, target, sigma = 3, force_normal = FALSE) {

  # validate input
  .validate_LSL_USL_target(LSL, USL, target)

  # compute mean and standard deviation if normal distributions calculation is required
  if (force_normal) {

    # computes standard capability metrics assuming normal distribution
    samples <- extract_mu_and_sigma(fit)

    range <- USL - LSL

    three_sigma <- sigma * samples$sigma
    six_sigma   <- (2*sigma) * samples$sigma

    Cp  <- range / six_sigma
    CpU <- (USL - samples$mu) / three_sigma
    CpL <- (samples$mu - LSL) / three_sigma
    Cpk <- pmin(CpU, CpL)

    # TODO: I replaced 6 with '(2*sigma)' here to generalize to different sigmas, check it's correct
    Cpc <- range / ((2*sigma) * sqrt(pi /  2) * samples$mu - target)

    # Eq. 8.14 of Montgomery, 8th edition
    xi <- (samples$mu - target) / (samples$sigma)
    # Eq. 8.13 of Montgomery, 8th edition
    Cpm <- Cp / sqrt(1 + xi^2)

    # TODO: add Cpmk? (that also adjusts for the wrong location of the target?)

  } else {

    # computes capability metrics using percentiles
    # (i.e., (q) version of the metric = generalization to non-normal distributions)
    samples <- extract_percentiles(fit, sigma = sigma)

    range <- USL - LSL

    Cp  <- range / (samples$UP - samples$LP)
    CpU <- (USL - samples$MP) / (samples$UP - samples$MP)
    CpL <- (samples$MP - LSL) / (samples$MP - samples$LP)
    Cpk <- pmin(CpU, CpL)

    # TODO: I didn't find this one in the manuscript
    Cpc <- rep(NA, length(samples$LP))

    Cpm <- min(USL - target, target - LSL) / (sigma * sqrt( ( (samples$UP - samples$LP) / (2 * sigma) )^2  + (samples$MP - target)^2) )
  }

  lst <- list(
    Cp  = Cp,
    CpU = CpU,
    CpL = CpL,
    Cpk = Cpk,
    Cpc = Cpc,
    Cpm = Cpm
  )

  class(lst) <- "capability_metrics"
  attr(lst, "LSL")    <- LSL
  attr(lst, "USL")    <- USL
  attr(lst, "target") <- target
  return(lst)
}
