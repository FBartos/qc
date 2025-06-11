validate_LSL_USL_target <- function(LSL, USL, target) {
  BayesTools::check_real(LSL,    name = "LSL",    check_length = 1, allow_NA = FALSE, upper = USL)
  BayesTools::check_real(USL,    name = "USL",    check_length = 1, allow_NA = FALSE, lower = LSL)
  BayesTools::check_real(target, name = "target", check_length = 1, allow_NA = FALSE)
}

#'@export
extract_mu_and_sigma <- function(fit) {
  UseMethod("extract_mu_and_sigma")
}

#'@export
extract_mu_and_sigma.bpc_normal <- function(fit) {
  rstan::extract(fit$stanfit, pars = c("mu", "sigma"))
}

#'@export
extract_mu_and_sigma.bpc_t <- function(fit) {

  all <- rstan::extract(fit$stanfit, pars = c("mu", "scale", "nu"))
  return(with(all, {
    sigma <- scale * sqrt(nu / (nu - 2.0))
    return(list(mu = mu, sigma = sigma))
  }))
}

#' Title
#'
#' @param fit
#' @param LSL
#' @param USL
#' @param target
#'
#' @returns
#' @export
#'
#' @examples
compute_capability_metrics <- function(fit, LSL = -1, USL = 1, target = 0) {

  validate_LSL_USL_target(LSL, USL, target)
  samples <- extract_mu_and_sigma(fit)

  range <- USL - LSL

  three_sigma <- 3 * samples$sigma
  six_sigma   <- 6 * samples$sigma

  Cp  <- range / six_sigma
  CpU <- (USL - samples$mu) / three_sigma
  CpL <- (samples$mu - LSL) / three_sigma
  Cpk <- pmin(CpU, CpL)

  Cpc <- range / (6 * sqrt(pi /  2) * samples$mu - target)

  # Eq. 8.14 of Montgomery, 8th edition
  xi <- (samples$mu - target) / (samples$sigma)
  # Eq. 8.13 of Montgomery, 8th edition
  Cpm <- Cp / sqrt(1 + xi^2)

  lst <- list(
    Cp  = Cp,
    CpU = CpU,
    CpL = CpL,
    Cpk = Cpk,
    Cpc = Cpc,
    Cpm = Cpm
  )

  class(lst) <- "capability_metrics"
  return(lst)

}



#' Title
#'
#' @param fit
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
summarize_capability_metrics <- function(fit, ...) {
  UseMethod("summarize_capability_metrics")
}

#' @export
summarize_capability_metrics.bpc <- function(fit, LSL = -1, USL = 1, target = 0, ...) {

  metrics <- compute_capability_metrics(fit, LSL, USL, target)
  return(summarize_capability_metrics(metrics))

}

#' @export
summarize_capability_metrics.capability_metrics <- function(fit, ...) {

  out <- t(vapply(fit, function(x) {
    quantiles <- unname(stats::quantile(x, probs = c(0.025, .5, 0.975)))
    c(mean = mean(x), median = quantiles[2], sd = stats::sd(x), lower = quantiles[1], upper = quantiles[3])
  }, numeric(5L)))

  return(tibble::as_tibble(out, rownames = "metric"))

}
