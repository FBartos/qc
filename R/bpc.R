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

  # Identification of improper priors
  # 1. String "Jeffreys..."
  # 2. Uniform unbounded? (BayesTools checks this usually)
  is_improper_mu <- identical(prior_mu, "Jeffreys_mu")
  is_improper_sigma <- identical(prior_sigma, "Jeffreys_sigma")

  # When prior_mu is a PriorConjugate, the sigma prior is already encoded in
  # the joint NIG structure, so the Jeffreys check on prior_sigma is not relevant.
  if (sample_priors && !inherits(prior_mu, "PriorConjugate")) {
    if (is_improper_mu) {
      stop("Improper prior for mu (Jeffreys) cannot be used without data (or with sample_priors = TRUE).")
    }
    if (is_improper_sigma) {
      stop("Improper prior for sigma (Jeffreys) cannot be used without data (or with sample_priors = TRUE).")
    }
  }


  # Check method-distribution compatibility
  if (method == "integration" && distribution != "normal") {
    stop("The integration method currently only supports distribution = 'normal'. ",
         "Use method = 'mcmc' for t-distribution.")
  }

  object$method <- method
  object$distribution <- distribution
  object$prior_mu <- prior_mu
  object$prior_sigma <- prior_sigma

  # Dispatch based on method
  if (method == "integration") {

    # Integration method: bypass Stan, use numerical integration
    # Remove NAs from data
    x <- na.omit(x)

    # Fit using integration
    int_result <- .bpc_fit_integration(
      data = x, LSL = LSL, USL = USL, target = target,
      prior_mu = prior_mu, prior_sigma = prior_sigma, sigma = sigma,
      sample_priors = sample_priors
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
      ss_mean = mean,
      ss_sd   = sd
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

    orig_LSL    <- attr(metrics, "LSL")
    orig_USL    <- attr(metrics, "USL")
    orig_target <- attr(metrics, "target")

    # Compute interval probabilities for each metric
    interval_breaks <- c(-Inf, interval_probability, Inf)
    n_intervals <- length(interval_probability) + 1

    interval_summary <- do.call(rbind, lapply(metric_names, function(m) {
      r <- int_result$results[[m]]

      # If MC samples are available (prior-only mode), use direct binning
      if (!is.null(r$samples)) {
        bin_counts <- table(cut(r$samples,
                                breaks = interval_breaks, include.lowest = TRUE))
        probs <- as.numeric(bin_counts) / sum(bin_counts)
      } else {
        # Trapezoidal-rule CDF from the pre-computed density grid (fallback)
        grid_probs <- NULL
        if (!is.null(r$grid)) {
          g <- r$grid
          area <- r$area
          if (is.null(area)) area <- 1
          n_g <- nrow(g)
          dxx <- diff(g$x)
          avg_dens <- (g$density[-n_g] + g$density[-1]) / 2
          cdf_cond <- c(0, cumsum(avg_dens * dxx))

          grid_probs <- numeric(n_intervals)
          for (j in seq_len(n_intervals)) {
            lo <- interval_breaks[j]
            hi <- interval_breaks[j + 1]
            lo_cdf <- if (is.finite(lo) && lo >= min(g$x))
              (1 - area) + area * stats::approx(g$x, cdf_cond, xout = lo, rule = 2)$y
            else 0
            hi_cdf <- if (is.finite(hi) && hi <= max(g$x))
              (1 - area) + area * stats::approx(g$x, cdf_cond, xout = hi, rule = 2)$y
            else 1
            grid_probs[j] <- max(0, hi_cdf - lo_cdf)
          }
          gt <- sum(grid_probs)
          if (gt > 0) grid_probs <- grid_probs / gt
        }

        # Try proper numerical integration on the density solver
        probs <- numeric(n_intervals)
        use_integration <- TRUE
        for (j in seq_len(n_intervals)) {
          lo <- interval_breaks[j]
          hi <- interval_breaks[j + 1]
          probs[j] <- tryCatch(
            compute_cpk_prob_integration(
              numeric(0), orig_LSL, orig_USL, c(lo, hi), prior,
              metric = m, target = orig_target,
              cached_state = cached_state
            ),
            error = function(e) NA_real_
          )
        }

        # Validate: the density solver for Cpk can fail structurally for some
        # priors (returning ~0 where the grid shows non-trivial mass)
        if (anyNA(probs)) {
          use_integration <- FALSE
        } else if (!is.null(grid_probs)) {
          suspect <- any(grid_probs > 0.05 & probs < 0.001)
          if (suspect) use_integration <- FALSE
        }

        if (!use_integration) {
          probs <- if (!is.null(grid_probs)) grid_probs else {
            probs[is.na(probs)] <- 0
            probs
          }
        }

        total <- sum(probs)
        if (total > 0) probs <- probs / total
      }

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
