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
#' @param ... additional arguments. When \code{x} is \code{NULL}, \code{mean},
#' \code{sd}, and \code{N} can be supplied instead of raw observations for
#' \code{distribution = "normal"}. The Student-t model requires raw
#' observations in \code{x}.
#'
#' @export
bpc <- function(
    x,
    LSL, target, USL,
    distribution = "normal",
    method = c("mcmc", "integration"),

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

  # Process method argument first (handles default value)
  method <- match.arg(method)

  # check input for capability metrics calculation to fail fast before estimation
  .validate_LSL_USL_target(LSL = LSL, USL = USL, target = target)
  BayesTools::check_char(distribution, name = "distribution", check_length = 1, allow_values = c("normal", "t"))
  BayesTools::check_char(method, name = "method", check_length = 1, allow_values = c("mcmc", "integration"))
  BayesTools::check_real(sigma, name = "sigma", check_length = 1, lower = 0, allow_NA = FALSE)
  BayesTools::check_bool(force_normal, name = "force_normal", check_length = 1, allow_NA = FALSE)
  BayesTools::check_bool(sample_priors, name = "sample_priors", check_length = 1, allow_NA = FALSE)

  .validate_bpc_prior_configuration(
    method = method,
    distribution = distribution,
    prior_mu = prior_mu,
    prior_sigma = prior_sigma,
    prior_nu = prior_nu,
    sample_priors = sample_priors
  )


  # Check method-distribution compatibility
  if (method == "integration" && distribution != "normal") {
    stop("The integration method currently only supports distribution = 'normal'. ",
         "Use method = 'mcmc' for t-distribution.")
  }

  object$method <- method
  object$distribution <- distribution
  object$prior_mu <- prior_mu
  object$prior_sigma <- prior_sigma
  object$sigma <- sigma

  prepared_data <- .bpc_prepare_data(
    distribution = distribution,
    x = x,
    mean = dots$mean,
    sd = dots$sd,
    N = dots$N,
    allow_empty = sample_priors
  )

  # Dispatch based on method
  if (method == "integration") {

    int_result <- .bpc_fit_integration(
      data = prepared_data$raw_data, LSL = LSL, USL = USL, target = target,
      prior_mu = prior_mu, prior_sigma = prior_sigma, sigma = sigma,
      sample_priors = sample_priors,
      cached_state = prepared_data$cached_state
    )

    object$integration_result <- int_result
    object$metrics <- int_result$metrics
    object$coefficients <- int_result$coefficients

    class(object) <- "bpc"
    return(object)

  } else {
    # MCMC method: use Stan

    object$stan_data <- prepared_data$stan_data

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

.is_default_prior_sigma_placeholder <- function(prior_sigma) {
  is.null(prior_sigma) || identical(prior_sigma, "Jeffreys_sigma")
}

.validate_bpc_prior_configuration <- function(method, distribution,
                                              prior_mu, prior_sigma, prior_nu,
                                              sample_priors) {
  issues <- character()

  conjugate_args <- c(
    if (inherits(prior_mu, "PriorConjugate")) "`prior_mu`",
    if (inherits(prior_sigma, "PriorConjugate")) "`prior_sigma`"
  )

  if (method == "mcmc" && length(conjugate_args) > 0L) {
    issues <- c(
      issues,
      sprintf(
        "%s %s only supported with `method = \"integration\"`; use BayesTools priors for `method = \"mcmc\"`.",
        paste(conjugate_args, collapse = " and "),
        if (length(conjugate_args) == 1L) "is" else "are"
      )
    )
  }

  if (inherits(prior_mu, "PriorConjugate") &&
      !.is_default_prior_sigma_placeholder(prior_sigma)) {
    issues <- c(
      issues,
      "`prior_mu` is a full PriorConjugate and already contains the sigma prior; leave `prior_sigma` at its default placeholder."
    )
  }

  if (sample_priors) {
    improper <- c(
      if (!.is_proper_prior_for_sampling(prior_mu)) .format_prior_sampling_issue("prior_mu", prior_mu),
      if (!inherits(prior_mu, "PriorConjugate") &&
          !.is_proper_prior_for_sampling(prior_sigma)) .format_prior_sampling_issue("prior_sigma", prior_sigma),
      if (distribution == "t" &&
          !.is_proper_prior_for_sampling(prior_nu)) .format_prior_sampling_issue("prior_nu", prior_nu)
    )

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

.summary_statistics_unsupported_message <- function(distribution, allow_empty) {
  paste(
    sprintf("For `distribution = \"%s\"`, supply raw observations in `x`.", distribution),
    "Summary-statistics inputs (`mean`, `sd`, and `N`) are only supported for `distribution = \"normal\"`.",
    if (allow_empty) {
      sprintf("For prior-only sampling with `distribution = \"%s\"`, omit `mean`, `sd`, and `N`.", distribution)
    } else {
      sprintf("The %s likelihood cannot be reconstructed from `mean`, `sd`, and `N` alone.", if (distribution == "t") "Student-t" else distribution)
    }
  )
}

.bpc_prepare_data <- function(distribution = "normal", x = NULL, mean = NULL, sd = NULL, N = NULL, allow_empty = FALSE) {

  if (!is.null(x)) {
    BayesTools::check_real(x, name = "x", check_length = 0)
    x <- stats::na.omit(x)
    n <- length(x)
    x_bar <- if (n > 0L) mean(x) else 0
    sse <- if (n > 0L) sum((x - x_bar)^2) else 0

    return(list(
      raw_data = as.numeric(x),
      stan_data = list(
        x = as.array(x),
        N = n,
        is_ss = 0L,
        ss_mean = numeric(),
        ss_sd = numeric()
      ),
      cached_state = list(
        n = n,
        x_bar = x_bar,
        sse = sse
      )
    ))
  }

  has_summary <- !is.null(mean) || !is.null(sd) || !is.null(N)
  if (!has_summary) {
    if (!allow_empty) {
      if (distribution == "normal") {
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
      cached_state = list(
        n = 0L,
        x_bar = 0,
        sse = 0
      )
    ))
  }

  if (distribution != "normal") {
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

  list(
    raw_data = numeric(),
    stan_data = list(
      x = numeric(),
      N = N,
      is_ss = 1L,
      ss_mean = as.array(c(mean)),
      ss_sd = as.array(c(sd))
    ),
    cached_state = list(
      n = N,
      x_bar = mean,
      sse = if (N > 1L) (N - 1L) * sd^2 else 0
    )
  )
}

.bpc_data   <- function(distribution = "normal", x = NULL, mean = NULL, sd = NULL, N = NULL, allow_empty = FALSE) {
  .bpc_prepare_data(distribution = distribution, x = x, mean = mean, sd = sd, N = N, allow_empty = allow_empty)$stan_data
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

  if (missing(sigma)) {
    sigma <- object$sigma %||% 3
  }

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

      # Re-fit with new specification limits, reusing the cached posterior state
      int_result <- .bpc_fit_integration(
        data = numeric(0), LSL = LSL, USL = USL, target = target,
        prior_mu = object$prior_mu %||% "Jeffreys_mu",
        prior_sigma = object$prior_sigma %||% "Jeffreys_sigma",
        sigma = sigma,
        cached_state = object$integration_result$cached_state
      )

    } else {

      int_result <- object$integration_result
    }

    metrics <- int_result$metrics

    # Build summary from pre-computed stats
    metric_names <- names(int_result$results)

    summary <- do.call(rbind, lapply(metric_names, function(m) {
      stats <- int_result$results[[m]]$stats
      interval <- .integration_interval_from_entry(
        int_result$results[[m]],
        ci = "central",
        ci_level = ci.level
      )
      data.frame(
        metric = m,
        mean = stats["Mean"],
        median = stats["Median"],
        sd = stats["SD"],
        lower = interval[1],
        upper = interval[2],
        row.names = NULL
      )
    }))
    summary <- tibble::as_tibble(summary)

    # Use the prior and cached state from the integration result
    # (these were already computed during bpc() and should be reused)
    prior <- int_result$prior
    cached_state <- int_result$cached_state
    metric_sigma <- int_result$sigma %||% object$sigma %||% sigma

    orig_LSL    <- attr(metrics, "LSL")
    orig_USL    <- attr(metrics, "USL")
    orig_target <- attr(metrics, "target")

    # Compute interval probabilities for each metric
    interval_breaks <- c(-Inf, interval_probability, Inf)

    interval_summary <- do.call(rbind, lapply(metric_names, function(m) {
      r <- int_result$results[[m]]
      probs <- .integration_interval_probs_from_entry(
        entry = r,
        interval_breaks = interval_breaks,
        metric = m,
        prior = prior,
        cached_state = cached_state,
        LSL = orig_LSL,
        USL = orig_USL,
        target = orig_target,
        sigma_level = metric_sigma
      )

      df <- as.data.frame(t(probs))
      names(df) <- levels(cut(0, breaks = interval_breaks, include.lowest = TRUE))
      df$metric <- m
      df[, c("metric", names(df)[names(df) != "metric"])]
    }))
    interval_summary <- tibble::as_tibble(interval_summary)

    # Build divergence diagnostics table from integration results
    divergence_diagnostics <- NULL
    if (!is.null(int_result$divergence)) {
      div_rows <- lapply(metric_names, function(m) {
        d <- int_result$divergence[[m]]
        if (is.null(d)) d <- int_result$results[[m]]$divergence_info
        if (is.null(d)) return(NULL)
        if (!d$mean_divergent && !d$sd_divergent) return(NULL)
        data.frame(
          metric = m,
          mean_divergent = d$mean_divergent,
          sd_divergent   = d$sd_divergent,
          alpha          = d$alpha,
          reason         = d$reason %||% "",
          stringsAsFactors = FALSE
        )
      })
      div_rows <- Filter(Negate(is.null), div_rows)
      if (length(div_rows) > 0)
        divergence_diagnostics <- tibble::as_tibble(do.call(rbind, div_rows))
    }

    out <- list(
      call               = object$call,
      summary            = summary,
      interval_summary   = interval_summary,
      metrics            = metrics,
      integration_result = int_result,
      divergence_diagnostics = divergence_diagnostics
    )
    attr(out, "has_divergent_moments") <- !is.null(divergence_diagnostics) && nrow(divergence_diagnostics) > 0
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
