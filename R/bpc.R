#' Bayesian Process Capability
#'
#' @param x the data
#' @param LSL lower specification limit
#' @param target target value
#' @param USL upper specification limit
#' @param distribution distribution for fitting the data
#' @param method computation method: "mcmc" for Stan-based MCMC sampling, or "integration" for numerical integration
#' @param prior_mu prior for the process mean, can be a string or a prior object
#' @param prior_sigma prior for the process standard deviation, can be a string or a prior object
#' @param prior_nu prior for the degrees of freedom, can be a string or a prior object
#' @param chains number of chains to run, defaults to 4
#' @param iter number of iterations per chain, defaults to 10000
#' @param warmup number of warmup iterations per chain, defaults to 5000
#' @param thin thinning factor, defaults to 1
#' @param cores number of cores for parallel processing
#' @param parallel whether to run chains in parallel, defaults to FALSE
#' @param control a list of control settings for the Stan model, see \code{\link[rstan]{sampling}} for details, defaults to \code{set_control()}
#' @param convergence_checks a list of convergence checks for the Stan model, see \code{\link[rstan]{check_convergence}} for details, defaults to \code{set_convergence_checks()}
#' @param seed a random seed for reproducibility, defaults to NULL
#' @param silent whether to suppress output during fitting, defaults to TRUE
#' @param sigma the number of standard deviations to use for the capability metrics, defaults to 3
#' @param force_normal whether to force the calculation of capability metrics assuming normal distribution, defaults to FALSE
#' @param sample_priors whether samples should be obtained from the prior distribution only (cannot be combined with improper prior distributions)
#' @param ...
#'
#' @export
bpc <- function(
    x,
    LSL, target, USL,
    distribution = "normal",
    method = "mcmc",

    # prior settings
    prior_mu    = "Jeffreys_mu",
    prior_sigma = "Jeffreys_sigma",
    prior_nu    = NULL, # TODO: set up default prior distribution?

    # stan control settings
    chains  = 4, iter = 10000, warmup = 5000, thin = 1, cores = chains, parallel = FALSE,
    control = set_control(), convergence_checks = set_convergence_checks(),
    seed = NULL, silent = TRUE,

    sample_priors = FALSE,

    # capability metrics settings
    sigma = 3, force_normal = FALSE,
    ...) {

  ### Instructions for adding a new distribution
  # the following functions need to be created:
  # samples_to_mu_and_sigma.<distribution> (S3 Class)
  # samples_to_percentiles.<distribution> (S3 Class)
  # samples_to_posterior_predictives.<distribution> (S3 Class)
  # samples_to_E_abs_dev.<distribution> (S3 Class)
  # stanfit model definition that are compiled to stanmodels (stored in inst/stan/<distribution>.stan)
  # the following functions need to be extended:
  # .bpc_parameters (Dispatch)
  # .bpc_priors (possibly including definition of new priors in the function calls)

  dots         <- list(...)
  object       <- list()
  object$call  <- match.call()

  # check input for capability metrics calculation to fail fast before estimation
  .validate_LSL_USL_target(LSL = LSL, USL = USL, target = target)
  BayesTools::check_char(distribution, name = "distribution", check_length = 1, allow_values = c("normal", "t"))
  BayesTools::check_char(method, name = "method", check_length = 1, allow_values = c("mcmc", "integration"))
  BayesTools::check_real(sigma, name = "sigma", check_length = 1, lower = 0, allow_NA = FALSE)
  BayesTools::check_bool(force_normal, name = "force_normal", check_length = 1, allow_NA = FALSE)
  BayesTools::check_bool(sample_priors, name = "sample_priors", check_length = 1, allow_NA = FALSE)

  # Check method-distribution compatibility

  if (method == "integration" && distribution != "normal") {
    stop("The integration method currently only supports distribution = 'normal'. ",
         "Use method = 'mcmc' for t-distribution.")
  }

  object$method <- method
  object$distribution <- distribution

  # Dispatch based on method
  if (method == "integration") {

    # Integration method: bypass Stan, use numerical integration
    # Remove NAs from data
    x <- na.omit(x)

    # Fit using integration
    int_result <- .bpc_fit_integration(
      data = x, LSL = LSL, USL = USL, target = target,
      prior_mu = prior_mu, prior_sigma = prior_sigma, sigma = sigma
    )

    object$integration_result <- int_result
    object$metrics <- int_result$metrics
    object$coefficients <- int_result$coefficients

    class(object) <- "bpc"
    return(object)

  } else {
    # MCMC method: use Stan

    # prepare data
    object$stan_data   <- .bpc_data(x = x, mean = dots$mean, sd = dots$sd, N = dots$N)

    # prepare priors
    object$stan_priors <- .bpc_priors(distribution = distribution, prior_mu = prior_mu, prior_sigma = prior_sigma, prior_nu = prior_nu, sample_priors = sample_priors)

    # collect control settings
    object$control <- .stan_check_and_list_fit_settings(
      chains = chains, warmup = warmup, iter = iter, thin = thin,
      parallel = parallel, cores = cores, silent = silent, seed = seed, control = control
    )

    # fit stan model
    object$stanfit      <- .bpc_fit(distribution = distribution, data = object$stan_data, priors = object$stan_priors, control = object$control)

    # add class to the object because it's required for dispatching in compute_capability_metrics
    class(object) <- "bpc"

    # compute capability metrics
    object$metrics <- .compute_capability_metrics(object, LSL = LSL, USL = USL, target = target, sigma = sigma, force_normal = force_normal)

    # add coefficients
    object$coefficients <- sapply(object$metrics, mean)

    return(object)
  }
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
.bpc_priors <- function(distribution, prior_mu = NULL, prior_sigma = NULL, prior_nu = NULL, sample_priors = FALSE) {

  # transform priors into stan format
  out <- c(
    if (distribution %in% c("normal", "t")) .stan_distribution("mu",    prior_mu,    sample_priors),
    if (distribution %in% c("normal", "t")) .stan_distribution("sigma", prior_sigma, sample_priors),
    if (distribution %in% c("t"))           .stan_distribution("nu",    prior_nu,    sample_priors)
  )

  out[["sample_priors"]] <- ifelse(sample_priors, 1, 0)

  return(out)
}
.bpc_fit    <- function(distribution, data, priors, control) {

  model_call <- list(
    object    = stanmodels[[distribution]],
    data      = c(data, priors),
    pars      = .bpc_parameters(distribution),
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
  attr(fit, "distribution") <- distribution

  return(fit)
}

# helpers
.bpc_parameters <- function(distribution) {

  if (distribution == c("normal"))
    return(c("mu", "sigma"))
  if (distribution == c("t"))
    return(c("mu", "scale", "nu"))
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

  # Handle integration method separately
  if (!is.null(object$method) && object$method == "integration") {

    # For integration, recompute if new limits provided
    if (!missing(LSL) && !missing(target) && !missing(USL)) {

      # Validate new specification limits
      .validate_LSL_USL_target(LSL = LSL, USL = USL, target = target)

      # Get original data and priors from the call
      data <- eval(object$call$x, envir = parent.frame())
      prior_mu <- eval(object$call$prior_mu, envir = parent.frame()) %||% "Jeffreys_mu"
      prior_sigma <- eval(object$call$prior_sigma, envir = parent.frame()) %||% "Jeffreys_sigma"

      # Re-fit with new specification limits
      int_result <- .bpc_fit_integration(
        data = data, LSL = LSL, USL = USL, target = target,
        prior_mu = prior_mu, prior_sigma = prior_sigma, sigma = sigma
      )

    } else {

      int_result <- object$integration_result
    }

    # Build summary from pre-computed stats
    metric_names <- names(int_result$results)
    h <- (1 - ci.level) / 2

    summary <- do.call(rbind, lapply(metric_names, function(m) {
      stats <- int_result$results[[m]]$stats
      data.frame(
        metric = m,
        mean = stats["Mean"],
        median = stats["Median"],
        sd = stats["SD"],
        lower = stats["Q2.5"],
        upper = stats["Q97.5"],
        row.names = NULL
      )
    }))
    summary <- tibble::as_tibble(summary)

    # Compute interval probabilities analytically
    prior_info <- .bayestools_to_integration_prior(
      attr(object$call, "prior_mu") %||% "Jeffreys_mu",
      attr(object$call, "prior_sigma") %||% "Jeffreys_sigma"
    )
    prior <- prior_info$prior
    cached_state <- int_result$cached_state

    # Get data from the original call
    data <- eval(object$call$x, envir = parent.frame())
    orig_LSL <- eval(object$call$LSL, envir = parent.frame())
    orig_USL <- eval(object$call$USL, envir = parent.frame())
    orig_target <- eval(object$call$target, envir = parent.frame())

    # Compute interval probabilities for each metric
    interval_breaks <- c(-Inf, interval_probability, Inf)
    n_intervals <- length(interval_probability) + 1

    interval_summary <- do.call(rbind, lapply(metric_names, function(m) {
      probs <- numeric(n_intervals)
      for (i in seq_len(n_intervals)) {
        bounds <- c(interval_breaks[i], interval_breaks[i + 1])
        probs[i] <- compute_cpk_prob_integration(
          data, orig_LSL, orig_USL, bounds, prior, metric = m,
          target = orig_target, cached_state = cached_state
        )
      }
      df <- as.data.frame(t(probs))
      names(df) <- levels(cut(0, breaks = interval_breaks, include.lowest = TRUE))
      df$metric <- m
      df[, c("metric", names(df)[names(df) != "metric"])]
    }))
    interval_summary <- tibble::as_tibble(interval_summary)

    out <- list(
      call               = object$call,
      summary            = summary,
      interval_summary   = interval_summary,
      metrics            = object$metrics,
      integration_result = int_result
    )
    class(out) <- "bpc_summary"
    return(out)
  }

  # MCMC method: original behavior
  # recompute capability metrics if specification limits are set
  if (!missing(LSL) && !missing(target) && !missing(USL)) {
    metrics <- .compute_capability_metrics(object, LSL = LSL, USL = USL, target = target, sigma = sigma, force_normal = force_normal)
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
  print(as.data.frame(round(x$interval_summary[,-1], 4)), quote = FALSE, right = TRUE, row.names = unlist(x$interval_summary[,1]))

  invisible(x$summary)
}
