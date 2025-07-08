
#' Bayesian Process Capability
#'
#' @param x the data
#' @param LSL Lower Specification Limit
#' @param target Target value
#' @param USL Upper Specification Limit
#' @param prior_mu prior for the process mean, can be a string or a prior object
#' @param prior_sigma prior for the process standard deviation, can be a string or a prior object
#' @param prior_nu prior for the degrees of freedom, can be a string or a prior object, only used if the t-distribution is selected
#' @param chains number of chains to run, defaults to 4
#' @param iter number of iterations per chain, defaults to 10000
#' @param warmup number of warmup iterations per chain, defaults to 5000
#' @param thin thinning factor, defaults to 1
#' @param parallel whether to run chains in parallel, defaults to FALSE
#' @param control a list of control settings for the Stan model, see \code{\link[rstan]{sampling}} for details, defaults to \code{set_control()}
#' @param convergence_checks a list of convergence checks for the Stan model, see \code{\link[rstan]{check_convergence}} for details, defaults to \code{set_convergence_checks()}
#' @param seed a random seed for reproducibility, defaults to NULL
#' @param silent whether to suppress output during fitting, defaults to TRUE
#' @param sigma the number of standard deviations to use for the capability metrics, defaults to 3
#' @param force_normal whether to force the calculation of capability metrics assuming normal distribution, defaults to FALSE
#' @param ...
#'
#' @export
bpc <- function(
    x,
    LSL, target, USL,

    # prior settings
    prior_mu    = "Jeffreys_mu",
    prior_sigma = "Jeffreys_sigma",
    prior_nu    = NULL,

    # stan control settings
    chains  = 4, iter = 10000, warmup = 5000, thin = 1, parallel = FALSE,
    control = set_control(), convergence_checks = set_convergence_checks(),
    seed = NULL, silent = TRUE,

    # capability metrics settings
    sigma = 3, force_normal = FALSE,
    ...) {

  dots         <- list(...)
  object       <- list()
  object$call  <- match.call()

  # check input for capability metrics calculation to fail fast before estimation
  .validate_LSL_USL_target(LSL = LSL, USL = USL, target = target)
  BayesTools::check_real(sigma, name = "sigma", check_length = 1, lower = 0, allow_NA = FALSE)
  BayesTools::check_bool(force_normal, name = "force_normal", check_length = 1, allow_NA = FALSE)

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
  object$stanfit <- .bpc_fit(data = object$stan_data, priors = object$stan_priors, control = object$control)

  # add class to the object because it's required for dispatching in compute_capability_metrics
  class(object) <- c("bpc", paste0("bpc_", .bpc_distribution(object$stan_priors)))

  # compute capability metrics
  object$metrics <- .bpc_compute_capability_metrics(object, LSL = LSL, USL = USL, target = target, sigma = sigma, force_normal = force_normal)

  # add coefficients
  object$coefficients <- sapply(object$metrics, mean)

  return(object)
}

### internal functions ----
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
    object    = stanmodels[[.bpc_distribution(priors)]],
    data      = c(data, priors),
    pars      = c("mu", if (.bpc_distribution(priors) == "t") c("scale", "nu") else "sigma"),
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
  attr(fit, "distribution") <- .bpc_distribution(priors)

  return(fit)
}

# helpers
.bpc_distribution <- function(priors) {
  if (!is.null(priors[["prior_type_nu"]])) {
    return("t")
  } else {
    return("normal")
  }
}

### print and summary functions ----
#' @export
print.bpc <- function(x, ...) {

  cat("\nCall:\n")
  print(x$call)
  cat("\nBayesian Process Capability:\n")
  print(round(coef(x), 4))

  invisible(x)
}

#' @export
summary.bpc <- function(object, LSL, target, USL, sigma = 3, force_normal = FALSE, ci.level = 0.95, interval_probability = c(1.00, 1.33, 1.50, 2.00), ...) {

  BayesTools::check_real(sigma, name = "sigma", check_length = 1, lower = 0, allow_NA = FALSE)
  BayesTools::check_real(ci.level, name = "ci.level", check_length = 1, allow_NA = FALSE, lower = 0, upper = 1)
  BayesTools::check_real(interval_probability, name = "interval_probability", check_length = 0, allow_NA = FALSE)
  BayesTools::check_bool(force_normal, name = "force_normal", check_length = 1, allow_NA = FALSE)

  # recompute capability metrics if specification limits are set
  if (!missing(LSL) && !missing(target) && !missing(USL)) {
    metrics <- .bpc_compute_capability_metrics(object, LSL = LSL, USL = USL, target = target, sigma = sigma, force_normal = force_normal)
  } else {
    metrics <- object$metrics
  }

  ### compute mean, median, and credible intervals for the metrics
  h     <- (1 - ci.level) / 2
  probs <- c(h, .5, 1 - h)

  summary <- t(vapply(metrics, function(x) {
    quantiles <- unname(stats::quantile(x, probs = probs, na.rm = TRUE))
    c(mean = mean(x), median = quantiles[2], sd = stats::sd(x), lower = quantiles[1], upper = quantiles[3])
  }, numeric(5L)))
  summary <- tibble::as_tibble(summary, rownames = "metric")

  ### compute interval summaries
  interval_summary <- t(vapply(metrics, FUN = function(x) {
    table(cut(x, breaks = c(-Inf, interval_probability, Inf), include.lowest = TRUE)) / length(x)
  }, FUN.VALUE = numeric(length(interval_probability) + 1L)))
  interval_summary <- tibble::as_tibble(interval_summary, rownames = "metric")

  out <- list(
    call             = object$call,
    summary          = summary,
    interval_summary = interval_summary,
    metrics          = metrics
  )
  class(out) <- "bpc_summary"

  return(out)
}

#' @export
print.bpc_summary <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)

  cat("\nBayesian Process Capability:\n")
  print(as.data.frame(round(x$summary[,-1], 4)), quote = FALSE, right = TRUE, row.names = unlist(x$summary[,1]))

  cat("\nInterval Probability:\n")
  print(as.data.frame(round(x$interval_summary, 4)), quote = FALSE, right = TRUE, row.names = rownames(x$interval_summary))

  invisible(x$summary)
}
