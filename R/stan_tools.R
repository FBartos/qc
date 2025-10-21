### this file contains common stan functions ###
### stan control functions ----

#' @title Convergence checks of the fitting process
#'
#' @description Set values for the convergence checks of the fitting process.
#'
#' @param max_Rhat maximum value of the R-hat diagnostic.
#' Defaults to \code{1.05}.
#' @param min_ESS minimum estimated sample size.
#' Defaults to \code{500}.
#' @param adapt_delta tuning parameter of HMC.
#' Defaults to \code{0.80}.
#' @param max_treedepth tuning parameter of HMC.
#' Defaults to \code{15}.
#' @param bridge_max_iter maximum number of iterations for the
#' \link[bridgesampling]{bridge_sampler} function. Defaults to \code{10000}
#'
#'
#' @return \code{set_control} returns a list of control settings
#' and \code{set_convergence_checks} returns a list of convergence checks settings.
#'
#' @export set_control
#' @export set_convergence_checks
#' @name stan_control
#' @aliases set_control, set_convergence_checks
NULL

#' @rdname stan_control
set_convergence_checks  <- function(max_Rhat = 1.05, min_ESS = 500){

  BayesTools::check_real(max_Rhat, "max_Rhat", lower = 1)
  BayesTools::check_real(min_ESS,  "min_ESS",  lower = 0)

  convergence_checks <- list(
    max_Rhat            = max_Rhat,
    min_ESS             = min_ESS
  )

  return(convergence_checks)
}
#' @rdname stan_control
set_control             <- function(adapt_delta = 0.80, max_treedepth = 15, bridge_max_iter = 1000){

  BayesTools::check_real(adapt_delta, "adapt_delta", lower = 0, upper = 1)
  BayesTools::check_int(max_treedepth, "max_treedepth", lower = 1)
  BayesTools::check_int(bridge_max_iter, "bridge_max_iter", lower = 1)

  control <- list(
    adapt_delta     = adapt_delta,
    max_treedepth   = max_treedepth,
    bridge_max_iter = bridge_max_iter
  )

  return(control)

}

.stan_check_and_list_fit_settings  <- function(chains, warmup, iter, thin, parallel, cores, silent, seed, control, check_mins = list(chains = 1, warmup = 50, iter = 50, thin = 1), call = ""){

  BayesTools::check_int(chains, "chains",  lower = check_mins[["chains"]],  call = call)
  BayesTools::check_int(warmup, "warmup",  lower = check_mins[["warmup"]],  call = call)
  BayesTools::check_int(iter,   "iter",    lower = min(check_mins[["iter"]], warmup + 1), call = call)
  BayesTools::check_int(thin,   "thin",    lower = check_mins[["thin"]],    call = call)
  BayesTools::check_list(control, "control", check_names = c("adapt_delta", "max_treedepth", "bridge_max_iter"))

  BayesTools::check_bool(parallel, "parallel",                call = call)
  BayesTools::check_int(cores,     "cores", lower = 1,        call = call)
  BayesTools::check_bool(silent,   "silent",                  call = call)
  BayesTools::check_int(seed,      "seed", allow_NULL = TRUE, call = call)

  if(!parallel){
    cores <- 1
  }else if(cores > qc.get_option("max_cores")){
    cores <- qc.get_option("max_cores")
  }

  if(is.null(control[["adapt_delta"]])){
    control[["adapt_delta"]] <- 0.80
  }else{
    BayesTools::check_real(control[["adapt_delta"]], "adapt_delta", lower = 0, upper = 1)
  }
  if(is.null(control[["max_treedepth"]])){
    control[["max_treedepth"]] <- 15
  }else{
    BayesTools::check_int(control[["max_treedepth"]], "max_treedepth", lower = 1)
  }
  if(is.null(control[["bridge_max_iter"]])){
    control[["bridge_max_iter"]] <- 1000
  }else{
    BayesTools::check_int(control[["bridge_max_iter"]], "bridge_max_iter", lower = 1)
  }

  return(invisible(list(
    chains   = chains,
    warmup   = warmup,
    iter     = iter,
    thin     = thin,
    parallel = parallel,
    cores    = cores,
    silent   = silent,
    adapt_delta     = control[["adapt_delta"]],
    max_treedepth   = control[["max_treedepth"]],
    bridge_max_iter = control[["bridge_max_iter"]],
    seed     = seed
  )))
}

### stan prior functions ----
# transforms BayesTools priors into pre-compiled stan code
.stan_distribution            <- function(parameter, prior, sample_priors){

  out <- list()

  # special handling of the Jeffreys priors pseudo-distributions
  if(parameter %in% c("mu", "sigma", "nu") && is.character(prior) && prior %in% c("Jeffreys_mu", "Jeffreys_sigma", "uniform_nu")){

    if (sample_priors)
      stop("Improper prior distributions cannot be sampled from with `sample_priors = TRUE`")

    out[[paste0("estimate_", parameter)]]  <- 1
    out[[paste0("fixed_",    parameter)]]  <- numeric()

    out[[paste0("prior_type_", parameter)]] <- switch(
      prior,
      "Jeffreys_mu"     = 98,
      "Jeffreys_sigma"  = 99,
      "uniform_nu"      = 98
    )

    out[[paste0("bounds_", parameter)]]      <- switch(
      prior,
      "Jeffreys_mu"     = c(999, 999),
      "Jeffreys_sigma"  = c(0,   999),
      "uniform_nu"      = c(0,   999)
    )
    out[[paste0("bounds_type_", parameter)]] <- switch(
      prior,
      "Jeffreys_mu"     = c(0, 0),
      "Jeffreys_sigma"  = c(1, 0),
      "uniform_nu"      = c(1, 0)
    )

    out[[paste0("prior_parameters_", parameter)]] <- c(999, 999, 999)

    return(out)
  }

  out[[paste0("prior_type_", parameter)]] <- switch(
    prior[["distribution"]],
    "point"           = 0,
    "normal"          = 1,
    "lognormal"       = 2,
    "t"               = 4,
    "gamma"           = 5,
    "invgamma"        = 6,
    "uniform"         = 7,
    "beta"            = 8,
    "exp"             = 9
  )

  if(is.prior.point(prior)){
    # point priors
    out[[paste0("estimate_", parameter)]]    <- 0
    out[[paste0("fixed_", parameter)]]       <- as.array(prior$parameters[["location"]])

    out[[paste0("bounds_", parameter)]]      <- numeric()
    out[[paste0("bounds_type_", parameter)]] <- numeric()

    out[[paste0("prior_parameters_", parameter)]] <- numeric()

  }else if(is.prior.simple(prior)){
    # non-point priors
    out[[paste0("estimate_", parameter)]]  <- 1
    out[[paste0("fixed_",    parameter)]]  <- numeric()

    out[[paste0("bounds_", parameter)]] <- c(
      if(is.infinite(prior$truncation[["lower"]])) 999 else prior$truncation[["lower"]],
      if(is.infinite(prior$truncation[["upper"]])) 999 else prior$truncation[["upper"]]
    )

    out[[paste0("bounds_type_", parameter)]] <- c(
      if(is.infinite(prior$truncation[["lower"]])) 0 else 1,
      if(is.infinite(prior$truncation[["upper"]])) 0 else 1
    )

    out[[paste0("prior_parameters_", parameter)]] <- .stan_distribution_parameters(prior)

  }else{
    stop("Other prior distributions are not implemented for stan.")
  }

  return(out)
}
.stan_distribution_parameters <- function(prior){

  # a vector of length three always needs to be passed - filling the redundant values with 999
  prior_parameters <- rep(999, 3)

  if(prior[["distribution"]] == "normal"){
    prior_parameters[1] <- prior$parameters[["mean"]]
    prior_parameters[2] <- prior$parameters[["sd"]]
  }else if(prior[["distribution"]] == "lognormal"){
    prior_parameters[1] <- prior$parameters[["meanlog"]]
    prior_parameters[2] <- prior$parameters[["sdlog"]]
  }else if(prior[["distribution"]] == "t"){
    prior_parameters[1] <- prior$parameters[["df"]]
    prior_parameters[2] <- prior$parameters[["location"]]
    prior_parameters[3] <- prior$parameters[["scale"]]
  }else if(prior[["distribution"]] == "gamma"){
    prior_parameters[1] <- prior$parameters[["shape"]]
    prior_parameters[2] <- prior$parameters[["rate"]]
  }else if(prior[["distribution"]] == "invgamma"){
    prior_parameters[1] <- prior$parameters[["shape"]]
    prior_parameters[2] <- prior$parameters[["scale"]]
  }else if(prior[["distribution"]] == "uniform"){
    prior_parameters[1] <- prior$parameters[["a"]]
    prior_parameters[2] <- prior$parameters[["b"]]
  }else if(prior[["distribution"]] == "beta"){
    prior_parameters[1] <- prior$parameters[["alpha"]]
    prior_parameters[2] <- prior$parameters[["beta"]]
  }else if(prior[["distribution"]] == "exp"){
    prior_parameters[1] <- prior$parameters[["rate"]]
  }

  return(prior_parameters)
}
