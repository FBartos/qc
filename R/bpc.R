#' Bayesian Process Capability
#'
#' @param x the data
#' @param LSL lower specification limit
#' @param target target value
#' @param USL upper specification limit
#' @param distribution distribution for fitting the data
#' @param method computation method: "mcmc" for Stan-based MCMC sampling,
#'   "integration" for numerical integration, or \code{NULL} to use the
#'   distribution-specific default.
#' @param prior unified prior specification. When omitted, the normal
#'   integration path uses \code{"DCSI"}, while MCMC and Student-t fits use
#'   \code{"Jeffreys"}. Use \code{prior_independent()} for parameter-wise
#'   priors.
#' @param chains number of chains to run, defaults to 4
#' @param iter number of iterations per chain, defaults to 10000
#' @param warmup number of warmup iterations per chain, defaults to 5000
#' @param thin thinning factor, defaults to 1
#' @param cores number of cores for parallel processing
#' @param parallel whether to run chains in parallel, defaults to FALSE
#' @param control a list of control settings for the Stan model, see \code{\link[rstan]{sampling}} for details, defaults to \code{set_control()}
#' @param seed a random seed for reproducibility, defaults to NULL
#' @param silent whether to suppress output during fitting, defaults to TRUE
#' @param sigma the number of standard deviations to use for the capability metrics, defaults to 3
#' @param sample_priors whether samples should be obtained from the prior distribution only (cannot be combined with improper prior distributions)
#' @param ... additional arguments. When \code{x} is \code{NULL}, \code{mean},
#' \code{sd}, and \code{N} can be supplied instead of raw observations for
#' \code{distribution = "normal"}. No other named arguments are accepted. The
#' Student-t model requires raw observations in \code{x}.
#'
#' @details
#' Priors are specified through the single \code{prior} argument. When
#' \code{method} is omitted, the normal likelihood defaults to
#' \code{method = "integration"} and \code{prior = "DCSI"}, which resolves to a
#' conjugate Normal-Inverse-Gamma prior centered at the midpoint of the
#' specification limits. Explicit MCMC fits with omitted \code{prior} use
#' \code{prior = "Jeffreys"} because MCMC currently supports parameter-wise
#' priors. Use \code{prior_independent()} for custom parameter-wise priors, or
#' \code{prior_conjugate()}, \code{prior_joint()}, and
#' \code{prior_semi_conjugate()} for normal-likelihood integration priors.
#'
#' @examples
#' data("pistonrings", package = "qc")
#' diameter <- subset(pistonrings, trial)[["diameter"]]
#'
#' # Default: normal likelihood, numerical integration, and DCSI prior.
#' fit <- bpc(diameter, LSL = 73.98, target = 74.00, USL = 74.02)
#' coef(fit)
#' summary(fit)
#'
#' # Switch to the parameter-wise Jeffreys prior.
#' fit_jeffreys <- bpc(
#'   diameter,
#'   LSL = 73.98, target = 74.00, USL = 74.02,
#'   prior = "Jeffreys"
#' )
#' coef(fit_jeffreys)
#'
#' # Use a custom conjugate normal-likelihood prior with integration.
#' fit_conjugate <- bpc(
#'   diameter,
#'   LSL = 73.98, target = 74.00, USL = 74.02,
#'   prior = prior_conjugate(mu0 = 74, k0 = 2, alpha0 = 3, beta0 = 8.5e-5)
#' )
#' summary(fit_conjugate, LSL = 73.975, target = 74.00, USL = 74.025)
#'
#' \dontrun{
#' data("pistonrings", package = "qc")
#' diameter <- subset(pistonrings, trial)[["diameter"]]
#'
#' # MCMC uses parameter-wise priors.
#' fit_mcmc <- bpc(
#'   diameter,
#'   LSL = 73.98, target = 74.00, USL = 74.02,
#'   method = "mcmc",
#'   prior = prior_independent(
#'     mu = prior("normal", list(74, 0.02)),
#'     sigma = prior("gamma", list(2, 200))
#'   ),
#'   chains = 2, iter = 2000, warmup = 1000, seed = 1
#' )
#' summary(fit_mcmc)
#' }
#'
#' @export
bpc <- function(
    x = NULL,
    LSL, target, USL,
    distribution = "normal",
    method = NULL,

    # prior settings
    prior = "DCSI",

    # stan control settings
    chains  = 4, iter = 10000, warmup = 5000, thin = 1, cores = chains, parallel = FALSE,
    control = set_control(),
    seed = NULL, silent = TRUE,

    sample_priors = FALSE,

    # capability metrics settings
    sigma = 3,
    ...) {

  dots <- list(...)
  call <- match.call()
  prior_missing <- missing(prior)

  # Process method argument first (handles default value)
  .qc_distribution_spec(distribution)
  method <- .resolve_bpc_method(method, distribution)
  .check_bpc_dots(dots)

  # check input for capability metrics calculation to fail fast before estimation
  .validate_capability_request(
    LSL = LSL,
    USL = USL,
    target = target,
    sigma_level = sigma,
    sigma_name = "sigma"
  )
  BayesTools::check_bool(sample_priors, name = "sample_priors", check_length = 1, allow_NA = FALSE)

  .qc_distribution_spec(distribution, method = method)

  prepared_data <- .bpc_prepare_data(
    distribution = distribution,
    x = x,
    mean = dots$mean,
    sd = dots$sd,
    N = dots$N,
    allow_empty = sample_priors
  )

  prior_info <- .resolve_bpc_prior(
    prior = prior,
    prior_missing = prior_missing,
    distribution = distribution,
    method = method,
    LSL = LSL,
    USL = USL,
    cached_state = prepared_data$cached_state
  )

  .validate_bpc_prior_configuration(
    method = method,
    distribution = distribution,
    prior_info = prior_info,
    sample_priors = sample_priors
  )

  # Dispatch based on method
  if (method == "integration") {

    int_result <- .bpc_fit_integration(
      distribution = distribution,
      data = prepared_data$raw_data, LSL = LSL, USL = USL, target = target,
      prior = prior_info$prior, sigma = sigma,
      sample_priors = sample_priors,
      cached_state = prepared_data$cached_state
    )

    return(.new_bpc_fit(
      call = call,
      method = method,
      distribution = distribution,
      metrics = int_result$metrics,
      coefficients = int_result$coefficients,
      sigma = sigma,
      prior = prior_info$specification,
      prior_resolved = prior_info$prior,
      prior_map = prior_info$prior_map,
      integration_result = int_result
    ))

  } else {
    # MCMC method: use Stan

    # prepare priors
    stan_data <- prepared_data$stan_data
    stan_priors <- .bpc_priors(
      distribution = distribution,
      prior_map = prior_info$prior_map,
      sample_priors = sample_priors
    )

    # collect control settings
    control <- .stan_check_and_list_fit_settings(
      chains = chains, warmup = warmup, iter = iter, thin = thin,
      parallel = parallel, cores = cores, silent = silent, seed = seed, control = control
    )

    # fit stan model
    stanfit <- .bpc_fit(
      distribution = distribution,
      data = stan_data,
      priors = stan_priors,
      control = control
    )

    object <- .new_bpc_fit(
      call = call,
      method = method,
      distribution = distribution,
      sigma = sigma,
      prior = prior_info$specification,
      prior_resolved = prior_info$prior,
      prior_map = prior_info$prior_map,
      stan_data = stan_data,
      stan_priors = stan_priors,
      control = control,
      stanfit = stanfit
    )

    # compute capability metrics
    object$metrics <- .compute_capability_metrics(object, LSL = LSL, USL = USL, target = target, sigma = sigma)

    # add coefficients
    object$coefficients <- sapply(object$metrics, mean)

    return(object)
  }
}

### internal functions ----
.resolve_bpc_method <- function(method, distribution) {
  if (is.null(method)) {
    if (identical(distribution, "normal")) {
      return("integration")
    }
    return("mcmc")
  }

  match.arg(method, choices = c("mcmc", "integration"))
}

.check_bpc_dots <- function(dots) {
  dot_names <- names(dots) %||% character()
  if (length(dots) == 0L) {
    return(invisible(NULL))
  }

  if (any(!nzchar(dot_names))) {
    stop("All arguments in `...` must be named.", call. = FALSE)
  }

  unknown_args <- setdiff(dot_names, c("mean", "sd", "N"))
  if (length(unknown_args) > 0L) {
    stop(
      sprintf(
        "Unsupported argument%s in `...`: %s.",
        if (length(unknown_args) == 1L) "" else "s",
        paste(sprintf("`%s`", unknown_args), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  invisible(NULL)
}

.is_proper_prior_conjugate <- function(prior) {
  inherits(prior, "PriorConjugate") &&
    is.numeric(prior$k0) && length(prior$k0) == 1L && is.finite(prior$k0) &&
    is.numeric(prior$alpha0) && length(prior$alpha0) == 1L && is.finite(prior$alpha0) &&
    is.numeric(prior$beta0) && length(prior$beta0) == 1L && is.finite(prior$beta0) &&
    prior$k0 > 0 && prior$alpha0 > 0 && prior$beta0 > 0
}

.uses_unbounded_uniform_prior <- function(prior) {
  if (!inherits(prior, "prior") || !identical(prior[["distribution"]], "uniform")) {
    return(FALSE)
  }

  params <- prior[["parameters"]]
  truncation <- prior[["truncation"]]
  lower <- max(params[["a"]] %||% -Inf, truncation[["lower"]] %||% -Inf)
  upper <- min(params[["b"]] %||% Inf, truncation[["upper"]] %||% Inf)

  !is.finite(lower) || !is.finite(upper)
}

.is_proper_prior_for_sampling <- function(prior) {
  if (is.null(prior)) {
    return(TRUE)
  }

  if (is.character(prior)) {
    return(!prior %in% c("Jeffreys_mu", "Jeffreys_sigma", "uniform_nu"))
  }

  if (inherits(prior, "PriorConjugate")) {
    return(.is_proper_prior_conjugate(prior))
  }

  if (.uses_unbounded_uniform_prior(prior)) {
    return(FALSE)
  }

  TRUE
}

.format_prior_sampling_issue <- function(parameter_name, prior) {
  if (inherits(prior, "PriorConjugate")) {
    return(sprintf(
      "`%s` uses an improper PriorConjugate (requires k0 > 0, alpha0 > 0, beta0 > 0)",
      parameter_name
    ))
  }

  if (.uses_unbounded_uniform_prior(prior)) {
    return(sprintf("`%s` uses an unbounded uniform prior", parameter_name))
  }

  if (is.character(prior)) {
    return(sprintf("`%s = \"%s\"` is improper", parameter_name, prior))
  }

  sprintf("`%s` is improper", parameter_name)
}

.validate_bpc_prior_configuration <- function(method, distribution,
                                              prior_info,
                                              sample_priors,
                                              prior_map = NULL) {
  issues <- character()
  prior_map <- prior_map %||% prior_info$prior_map

  if (sample_priors) {
    improper <- character()
    if (!is.null(prior_map)) {
      for (parameter in names(prior_map)) {
        prior <- prior_map[[parameter]]
        if (!.is_proper_prior_for_sampling(prior)) {
          improper <- c(
            improper,
            .format_prior_sampling_issue(sprintf("prior$%s", parameter), prior)
          )
        }
      }
    } else if (!.is_proper_prior_for_sampling(prior_info$prior)) {
      improper <- c(
        improper,
        .format_prior_sampling_issue("prior", prior_info$prior)
      )
    }

    if (identical(prior_info$kind, "joint") && inherits(prior_info$prior, "PriorGeneric")) {
      # Properness of arbitrary joint priors is checked by the integration backend.
      improper <- setdiff(improper, .format_prior_sampling_issue("prior", prior_info$prior))
    }

    if (length(improper) > 0L) {
      issues <- c(
      issues,
      paste0(
          "Improper prior distributions cannot be sampled from with `sample_priors = TRUE`:",
          " ",
          paste(improper, collapse = "; "),
          "."
        )
      )
    }
  }

  if (length(issues) > 0L) {
    stop(paste(issues, collapse = " "), call. = FALSE)
  }

  invisible(NULL)
}

.bpc_priors <- function(distribution, prior_map, sample_priors = FALSE) {
  if (is.null(prior_map)) {
    stop(
      sprintf(
        "`method = \"mcmc\"` requires parameter-wise priors for `distribution = \"%s\"`.",
        distribution
      ),
      call. = FALSE
    )
  }
  prior_parameters <- .qc_distribution_parameter_names(distribution, type = "prior")
  prior_map <- prior_map[prior_parameters]

  # transform priors into stan format
  out <- unlist(
    lapply(prior_parameters, function(parameter) {
      .stan_distribution(parameter, prior_map[[parameter]], sample_priors)
    }),
    recursive = FALSE
  )

  out[["sample_priors"]] <- ifelse(sample_priors, 1, 0)

  return(out)
}

.summary_statistics_unsupported_message <- function(distribution, allow_empty) {
  model_name <- .qc_distribution_display_name(distribution)
  supported <- paste(sprintf('"%s"', .qc_distribution_summary_stat_names()), collapse = ", ")

  paste(
    sprintf("For `distribution = \"%s\"`, supply raw observations in `x`.", distribution),
    sprintf(
      "Summary-statistics inputs (`mean`, `sd`, and `N`) are only supported for `distribution = %s`.",
      supported
    ),
    if (allow_empty) {
      sprintf("For prior-only sampling with `distribution = \"%s\"`, omit `mean`, `sd`, and `N`.", distribution)
    } else {
      sprintf("The %s likelihood cannot be reconstructed from `mean`, `sd`, and `N` alone.", model_name)
    }
  )
}

.bpc_prepare_data <- function(distribution = "normal", x = NULL, mean = NULL, sd = NULL, N = NULL, allow_empty = FALSE) {
  supports_summary_statistics <- .qc_distribution_supports_summary_statistics(distribution)

  if (!is.null(x)) {
    BayesTools::check_real(x, name = "x", check_length = 0)
    x <- stats::na.omit(x)
    suff_state <- .as_qc_suff_stats_state(data = x)

    return(list(
      raw_data = as.numeric(x),
      stan_data = list(
        x = as.array(x),
        N = suff_state$n,
        is_ss = 0L,
        ss_mean = numeric(),
        ss_sd = numeric()
      ),
      cached_state = suff_state
    ))
  }

  has_summary <- !is.null(mean) || !is.null(sd) || !is.null(N)
  if (!has_summary) {
    if (!allow_empty) {
      if (supports_summary_statistics) {
        stop("When 'x' is NULL, supply all of 'mean', 'sd', and 'N'.", call. = FALSE)
      }

      stop(.summary_statistics_unsupported_message(distribution, allow_empty = FALSE), call. = FALSE)
    }

    return(list(
      raw_data = numeric(),
      stan_data = list(
        x = numeric(),
        N = 0L,
        is_ss = 0L,
        ss_mean = numeric(),
        ss_sd = numeric()
      ),
      cached_state = .as_qc_suff_stats_state(data = numeric())
    ))
  }

  if (!supports_summary_statistics) {
    stop(.summary_statistics_unsupported_message(distribution, allow_empty), call. = FALSE)
  }

  if (is.null(mean) || is.null(sd) || is.null(N)) {
    stop("When 'x' is NULL, supply all of 'mean', 'sd', and 'N'.", call. = FALSE)
  }

  BayesTools::check_real(mean, name = "mean", check_length = 1, allow_NA = FALSE)
  BayesTools::check_real(sd,   name = "sd",   check_length = 1, lower = 0, allow_NA = FALSE)
  BayesTools::check_int(N,     name = "N",    check_length = 1, lower = 1, allow_NA = FALSE)

  N <- as.integer(N)
  if (N == 1L && sd != 0) {
    stop("When 'N' is 1, 'sd' must be 0.", call. = FALSE)
  }

  suff_state <- .as_qc_suff_stats_state(
    cached_state = list(
      n = N,
      x_bar = mean,
      sse = if (N > 1L) (N - 1L) * sd^2 else 0
    )
  )

  list(
    raw_data = numeric(),
    stan_data = list(
      x = numeric(),
      N = N,
      is_ss = 1L,
      ss_mean = as.array(c(mean)),
      ss_sd = as.array(c(sd))
    ),
    cached_state = suff_state
  )
}

.bpc_data   <- function(distribution = "normal", x = NULL, mean = NULL, sd = NULL, N = NULL, allow_empty = FALSE) {
  .bpc_prepare_data(distribution = distribution, x = x, mean = mean, sd = sd, N = N, allow_empty = allow_empty)$stan_data
}
.bpc_fit    <- function(distribution, data, priors, control) {
  stan_model_name <- .qc_distribution_stan_model(distribution)
  stan_model <- stanmodels[[stan_model_name]]
  if (is.null(stan_model)) {
    stop(
      sprintf(
        "No Stan model is registered for `distribution = \"%s\"`.",
        distribution
      ),
      call. = FALSE
    )
  }

  model_call <- list(
    object    = stan_model,
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
    model_call$refresh <- 0
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
  .qc_distribution_parameter_names(distribution, type = "sample")
}

.qc_requested_limits <- function(LSL, target, USL,
                                 LSL_missing, target_missing, USL_missing) {
  provided <- c(!LSL_missing, !target_missing, !USL_missing)

  if (!any(provided)) {
    return(NULL)
  }

  if (!all(provided)) {
    stop("If any of LSL, USL, or target are provided, all three must be specified.", call. = FALSE)
  }

  .validate_LSL_USL_target(LSL = LSL, USL = USL, target = target)
}

.bpc_requested_limits <- .qc_requested_limits

.bpc_query_metrics <- function(object,
                               limits = NULL,
                               sigma = object$sigma %||% object$integration_result$sigma %||% 3) {
  is_integration <- identical(object$method, "integration")

  if (is.null(limits)) {
    return(list(
      metrics = object$metrics,
      integration_result = if (is_integration) object$integration_result else NULL
    ))
  }

  if (is_integration) {
    prior <- .bpc_query_prior(
      object = object,
      limits = limits
    )

    int_result <- .bpc_fit_integration(
      distribution = object$distribution %||% "normal",
      data = numeric(0),
      LSL = limits$LSL,
      USL = limits$USL,
      target = limits$target,
      prior = prior,
      sigma = sigma,
      cached_state = object$integration_result$cached_state
    )

    return(list(
      metrics = int_result$metrics,
      integration_result = int_result
    ))
  }

  list(
    metrics = .compute_capability_metrics(
      object,
      LSL = limits$LSL,
      USL = limits$USL,
      target = limits$target,
      sigma = sigma
    ),
    integration_result = NULL
  )
}

.bpc_query_prior <- function(object, limits) {
  prior_specification <- object$prior %||% object$prior_resolved
  if (is.null(prior_specification)) {
    return(NULL)
  }

  prior_info <- .resolve_bpc_prior(
    prior = prior_specification,
    prior_missing = FALSE,
    distribution = object$distribution %||% "normal",
    method = "integration",
    LSL = limits$LSL,
    USL = limits$USL,
    cached_state = object$integration_result$cached_state
  )

  prior_info$prior
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
summary.bpc <- function(object, LSL, target, USL, sigma = 3, ci.level = 0.95, interval_probability = c(1.00, 1.33, 1.50, 2.00), ...) {

  if (missing(sigma)) {
    sigma <- object$sigma %||% 3
  }
  .validate_sigma_level(sigma, name = "sigma")
  BayesTools::check_real(ci.level, name = "ci.level", check_length = 1, allow_NA = FALSE, lower = 0, upper = 1)
  BayesTools::check_real(interval_probability, name = "interval_probability", check_length = 0, allow_NA = FALSE)

  limits <- .qc_requested_limits(
    LSL = LSL,
    target = target,
    USL = USL,
    LSL_missing = missing(LSL),
    target_missing = missing(target),
    USL_missing = missing(USL)
  )

  is_integration <- identical(object$method, "integration")
  query <- .bpc_query_metrics(
    object,
    limits = limits,
    sigma = sigma
  )
  metrics <- query$metrics
  int_result <- query$integration_result

  distributions <- .as_qc_metric_distributions(metrics, what = .qc_metric_names())
  tables <- .qc_metric_distributions_summary_bundle(
    distributions,
    ci_level = ci.level,
    interval_probability = interval_probability
  )

  .new_bpc_summary(
    call = object$call,
    metrics = metrics,
    summary = tables$summary,
    interval_summary = tables$interval_summary,
    integration_result = if (is_integration) int_result else NULL,
    divergence_diagnostics = tables$divergence_diagnostics
  )
}

#' @export
print.bpc_summary <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)

  cat("\nBayesian Process Capability:\n")
  print(as.data.frame(round(x$summary[,-1], 4)), quote = FALSE, right = TRUE, row.names = unlist(x$summary[,1]))

  cat("\nInterval Probability:\n")
  print(as.data.frame(round(x$interval_summary[,-1], 4)), quote = FALSE, right = TRUE, row.names = unlist(x$interval_summary[,1]))

  if (isTRUE(attr(x, "has_divergent_moments"))) {
    cat("\nNote: Some moments are analytically infinite due to prior specification.\n")
    for (i in seq_len(nrow(x$divergence_diagnostics))) {
      cat("  ", x$divergence_diagnostics$reason[i], "\n")
    }
  }

  invisible(x$summary)
}

# ==============================================================================
# Divergence Query API
# ==============================================================================

#' Check whether a statistic was analytically determined to be infinite
#'
#' For certain prior-metric combinations, the posterior mean or standard deviation
#' is analytically known to diverge. This function queries the divergence
#' diagnostics attached to a bpc_summary object.
#'
#' @param x A bpc_summary object (from \code{summary(bpc(...))})
#' @param metric One of "Cp", "Cpu", "Cpl", "Cpk", "Cpm", "Cpc"
#' @param statistic One of "mean", "sd"
#' @return Logical: TRUE if the statistic was analytically determined to be infinite
#' @export
is_analytic <- function(x, metric, statistic = c("mean", "sd")) {
  UseMethod("is_analytic")
}

#' @export
is_analytic.bpc_summary <- function(x, metric, statistic = c("mean", "sd")) {
  statistic <- match.arg(statistic)
  diag <- x$divergence_diagnostics
  if (is.null(diag) || nrow(diag) == 0)
    return(FALSE)
  row <- diag[diag$metric == metric, , drop = FALSE]
  if (nrow(row) == 0)
    return(FALSE)
  switch(statistic,
    "mean" = row$mean_divergent[1],
    "sd"   = row$sd_divergent[1]
  )
}

#' Retrieve the full divergence diagnostics table
#'
#' Returns a tibble indicating which metric-statistic combinations were
#' analytically determined to be infinite, along with the effective shape
#' parameter alpha and a human-readable explanation.
#'
#' @param x A bpc_summary object
#' @return A tibble with columns: metric, mean_divergent, sd_divergent, alpha, reason.
#'   NULL if no divergent moments were detected.
#' @export
get_analytic_flags <- function(x) {
  UseMethod("get_analytic_flags")
}

#' @export
get_analytic_flags.bpc_summary <- function(x) {
  x$divergence_diagnostics
}
