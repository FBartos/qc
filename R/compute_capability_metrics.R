compute_capability_metrics <- function(fit, LSL = -1, USL = 1, target = 0) {

  samples <- rstan::extract(fit$stanfit, pars = c("mu", "sigma"))

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

  ls <- list(
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
