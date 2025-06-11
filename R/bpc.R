
#' @export
bpc <- function(
    x,

    # prior settings
    prior_mu    = "Jeffreys_mu",
    prior_sigma = "Jeffreys_sigma",
    prior_nu    = NULL,

    # stan control settings
    chains  = 4, iter = 10000, warmup = 5000, thin = 1, parallel = FALSE,
    control = set_control(), convergence_checks = set_convergence_checks(),
    seed = NULL, silent = TRUE,
    ...) {

  dots         <- list(...)
  object       <- list()
  object$call  <- match.call()

  # prepare data
  object$stan_data   <- .bpc_data(x = x, mean = dots$mean, sd = dots$sd, N = dots$N)

  # prepare priors
  object$stan_priors <- .bpc_priors(prior_mu = prior_mu, prior_sigma = prior_sigma, prior_nu = prior_nu)

  # collect control settings
  object$control <- .stan_check_and_list_fit_settings(
    chains = chains, warmup = warmup, iter = iter, thin = thin,
    parallel = parallel, cores = chains, silent = silent, seed = seed, control = control
  )

  # fit stan model
  object$stanfit   <- .bpc_fit(data = object$stan_data, priors = object$stan_priors, control = object$control)

  class(object) <- "bpc"
  return(object)
}

.bpc_data   <- function(x = NULL, mean = NULL, sd = NULL, N = NULL) {

  # use raw data if supplied, otherwise use summary statistics
  if (!is.null(x)) {

    # check input
    BayesTools::check_real(x, name = "x", check_length = 0)

    # remove NAs
    x <- na.omit(x)

    # return stan formatted data
    return(list(
      x = as.array(x),
      N = length(x),

      is_ss   = 0,
      ss_mean = numeric(),
      ss_sd   = numeric()
    ))

  } else {

    # check input
    BayesTools::check_real(mean, name = "mean", check_length = 1, allow_NA = FALSE)
    BayesTools::check_real(sd,   name = "sd",   check_length = 1, allow_NA = FALSE)
    BayesTools::check_integer(N, name = "N",    check_length = 1, allow_NA = FALSE)

    # return stan formatted data
    return(list(
      x = numeric(),
      N = N,

      is_ss   = 0,
      ss_mean = ss_mean,
      ss_sd   = ss_sd
    ))
  }
}
.bpc_priors <- function(prior_mu = NULL, prior_sigma = NULL, prior_nu = NULL) {

  # transform priors into stan format
  out <- c(
    .stan_distribution("mu",    prior_mu),
    .stan_distribution("sigma", prior_sigma),
    if (!is.null(prior_nu)) .stan_distribution("nu",    prior_nu)
  )

  return(out)
}
.bpc_fit    <- function(data, priors, control) {

  model_call <- list(
    object    = stanmodels[[if (is.null(priors[["prior_type_nu"]])) "normal" else "t"]],
    data      = c(data, priors),
    pars      = c("mu", "sigma", if (!is.null(priors[["prior_type_nu"]])) "nu"),
    chains    = control[["chains"]],
    warmup    = control[["warmup"]],
    iter      = control[["iter"]],
    thin      = control[["thin"]],
    cores     = control[["cores"]],
    control   = list(
      adapt_delta   = control[["adapt_delta"]],
      max_treedepth = control[["max_treedepth"]]
    )
  )

  if(control[["silent"]]){
    model_call$refresh <- -1
    model_call$open_progress <- FALSE
    model_call$show_messages <- FALSE
  }

  if(!is.null(control[["seed"]])){
    set.seed(control[["seed"]])
    model_call$seed <- control[["seed"]]
  }

  fit <- tryCatch(suppressWarnings(do.call(rstan::sampling, model_call)), error = function(e)e)

  return(fit)
}


#' @export
print.bpc <- function(x, ...) {

  cat("Bayesian Process Capability\n")
  print(x$call)

  invisible(x)
}
