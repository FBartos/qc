# ==============================================================================
# Integration Entry Points and Public Prior Constructors
# ==============================================================================

# ==============================================================================
# Prior Classes
# ==============================================================================

#' Create a conjugate prior for Normal-InverseGamma model
#' @param mu0 Prior mean location
#' @param k0 Prior precision multiplier (0 = noninformative)
#' @param alpha0 Prior shape for sigma^2
#' @param beta0 Prior rate for sigma^2
#' @return PriorConjugate object
#' @export
create_prior_conjugate <- function(mu0 = 0, k0 = 0, alpha0 = -0.5, beta0 = 0) {
  structure(list(mu0 = mu0, k0 = k0, alpha0 = alpha0, beta0 = beta0),
            class = "PriorConjugate")
}

#' Create a generic prior with custom log-density function
#' @param log_dens_fn Function(mu, sigma) returning log prior density
#' @param bayestools_priors Optional list with original BayesTools prior objects (mu and sigma)
#' @return PriorGeneric object
#' @export
create_prior_generic <- function(log_dens_fn, bayestools_priors = NULL) {
  structure(list(log_dens = log_dens_fn, bayestools_priors = bayestools_priors),
            class = "PriorGeneric")
}

#' Create a semi-conjugate prior (conjugate mu, non-conjugate sigma)
#' @param mu0 Prior mean location
#' @param k0 Prior precision multiplier (0 = noninformative)
#' @param log_dens_sigma Function(sigma) returning log prior density for sigma
#' @param bayestools_priors Optional list with original BayesTools prior objects
#' @return PriorSemiConjugateMu object
#' @export
create_prior_semi_mu <- function(mu0, k0, log_dens_sigma, bayestools_priors = NULL) {
  structure(list(mu0 = mu0, k0 = k0, log_dens_sigma = log_dens_sigma,
                 bayestools_priors = bayestools_priors),
            class = "PriorSemiConjugateMu")
}

#' Create a semi-conjugate prior (non-conjugate mu, conjugate sigma)
#' @param alpha0 Prior shape for sigma^2
#' @param beta0 Prior rate for sigma^2
#' @param log_dens_mu Function(mu) returning log prior density for mu
#' @param bayestools_priors Optional list with original BayesTools prior objects
#' @return PriorSemiConjugateSigma object
#' @export
create_prior_semi_sigma <- function(alpha0, beta0, log_dens_mu, bayestools_priors = NULL) {
  structure(list(alpha0 = alpha0, beta0 = beta0, log_dens_mu = log_dens_mu,
                 bayestools_priors = bayestools_priors),
            class = "PriorSemiConjugateSigma")
}

#' Create a unit information prior for a normal model
#'
#' Constructs a Normal-Inverse-Gamma (NIG) conjugate prior that carries
#' approximately one unit of Fisher information.  The prior parameters are
#' derived from the observed data:
#' \itemize{
#'   \item \eqn{\mu \mid \sigma^2 \sim \mathrm{Normal}(\bar{x},\, \sigma^2 / 1)}
#'   \item \eqn{\sigma^2 \sim \mathrm{Inv-Gamma}(1/2,\, s^2/2)}
#' }
#' which corresponds to \code{create_prior_conjugate(mu0 = xbar, k0 = 1,
#' alpha0 = 0.5, beta0 = s2 / 2)}.
#'
#' The returned \code{PriorConjugate} object can be passed as \code{prior}
#' to \code{\link{bpc}} when \code{method = "integration"}.
#'
#' @param x Numeric vector of observations.  \code{NA} values are removed.
#'   At least 2 finite observations are required.
#' @return A \code{PriorConjugate} object.
#' @export
create_prior_unit_information <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2L)
    stop("At least 2 finite observations are required to compute a unit information prior.")
  xbar <- mean(x)
  s2   <- stats::var(x)
  create_prior_conjugate(mu0 = xbar, k0 = 1, alpha0 = 0.5, beta0 = s2 / 2)
}

.is_integration_prior <- function(prior) {
  inherits(
    prior,
    c("PriorConjugate", "PriorGeneric", "PriorSemiConjugateMu", "PriorSemiConjugateSigma")
  )
}

.integration_prior_case <- function(prior) {
  if (inherits(prior, "PriorConjugate")) {
    return(1L)
  }
  if (inherits(prior, "PriorSemiConjugateMu")) {
    return(2L)
  }
  if (inherits(prior, "PriorSemiConjugateSigma")) {
    return(3L)
  }
  if (inherits(prior, "PriorGeneric")) {
    return(4L)
  }

  NULL
}

.conjugate_sigma_divergence_signature <- function(alpha0, beta0) {
  # A positive beta0 contributes exp(-beta0 / sigma^2), which suppresses the
  # sigma -> 0 singularity faster than any power law.
  if (is.numeric(beta0) && length(beta0) == 1L && is.finite(beta0) && beta0 > 0) {
    return(list(
      alpha = Inf,
      label = sprintf(
        "NIG(alpha0 = %s, beta0 = %s)",
        format(alpha0, digits = 3),
        format(beta0, digits = 3)
      )
    ))
  }

  list(
    alpha = alpha0,
    label = sprintf("NIG(alpha0 = %s)", format(alpha0, digits = 3))
  )
}

.integration_sigma_divergence_signature <- function(prior, sigma_prior = NULL) {
  if (inherits(prior, "PriorConjugate")) {
    return(.conjugate_sigma_divergence_signature(prior$alpha0, prior$beta0))
  }

  if (inherits(prior, "PriorSemiConjugateSigma")) {
    return(.conjugate_sigma_divergence_signature(prior$alpha0, prior$beta0))
  }

  bt_priors <- prior$bayestools_priors %||% NULL
  sigma_source <- bt_priors$sigma %||% sigma_prior

  if (inherits(sigma_source, "PriorConjugate")) {
    return(.conjugate_sigma_divergence_signature(
      sigma_source$alpha0,
      sigma_source$beta0
    ))
  }

  if (identical(sigma_source, "Jeffreys_sigma") || inherits(sigma_source, "prior")) {
    return(list(
      alpha = .extract_alpha_parameter(sigma_source),
      label = .sigma_prior_label(sigma_source)
    ))
  }

  list(alpha = Inf, label = "embedded sigma prior")
}

# ==============================================================================
# Integration Fit Function (called from bpc)
# ==============================================================================

#' Fit using numerical integration method
#' @param data Numeric vector of observations
#' @param LSL Lower specification limit
#' @param USL Upper specification limit
#' @param target Target value
#' @param prior Resolved integration prior object
#' @param sigma Number of standard deviations for capability metrics
#' @return List with metrics and integration results
#' @keywords internal
.bpc_fit_integration <- function(distribution = "normal",
                                 data, LSL, USL, target, prior,
                                 sigma = 3, sample_priors = FALSE,
                                 cached_state = NULL) {
  .qc_distribution_fit_integration(
    distribution = distribution,
    data = data,
    LSL = LSL,
    USL = USL,
    target = target,
    prior = prior,
    sigma = sigma,
    sample_priors = sample_priors,
    cached_state = cached_state
  )
}

.bpc_fit_integration_normal <- function(data, LSL, USL, target, prior,
                                        sigma = 3, sample_priors = FALSE,
                                        cached_state = NULL) {
  distribution <- "normal"

  # If sampling from priors, ignore data (likelihood becomes flat/unity effectively)
  # Ideally we should pass sample_priors down, but clearing data works if the
  # functions handle empty data correctly.
  # However, the conjugate update uses n and sufficiency stats.
  # If we set data to empty, n=0, and the posterior parameters equal prior parameters.
  if (sample_priors) {
    data <- numeric(0)
    cached_state <- NULL
  }

  if (!.is_integration_prior(prior)) {
    stop(
      "`method = \"integration\"` requires a resolved normal likelihood prior object.",
      call. = FALSE
    )
  }

  is_conjugate <- inherits(prior, "PriorConjugate")
  case <- .integration_prior_case(prior)

  # Pre-compute state (sufficient statistics) unless already provided
  if (is.null(cached_state)) {
    suff_state <- .integration_require_suff_state(data = data)

    if (case == 4L) {
      cached_state <- precompute_generic_state(data, prior)
    } else {
      cached_state <- suff_state
    }
  } else if (case == 4L && !inherits(cached_state, "qc_generic_cached_state")) {
    cached_state <- precompute_generic_state(data, prior, cached_state = cached_state)
  }

  n <- cached_state$n %||% 0L

  if (is_conjugate) {
    .compute_validated_conjugate_posterior(
      prior, data, cached_state,
      context = "The conjugate posterior for integration"
    )
  }

  # Analyze all metrics
  # Match order of MCMC results for consistency (Cp, Cpu, Cpl, Cpk, Cpc, Cpm)
  metrics <- c("Cp", "Cpu", "Cpl", "Cpk", "Cpc", "Cpm")

  # Pre-compute analytic divergence info for each metric.
  # With data (n > 0), the likelihood provides superexponential decay at sigma = 0,
  # so all moments are finite regardless of the prior.
  is_prior_only <- n == 0
  sigma_divergence <- .integration_sigma_divergence_signature(
    prior = prior
  )
  alpha_sigma <- sigma_divergence$alpha
  sigma_label <- sigma_divergence$label
  divergence_map <- if (is_prior_only) {
    setNames(lapply(metrics, function(m)
      .check_moment_divergence(m, alpha_sigma, sigma_label)
    ), metrics)
  } else {
    no_div <- list(mean_divergent = FALSE, sd_divergent = FALSE,
                   alpha = alpha_sigma, reason = NULL)
    setNames(rep(list(no_div), length(metrics)), metrics)
  }

  # For prior-only mode with BayesTools priors, pre-generate shared samples
  # so all metrics use the same (mu, sigma) draws for consistency.
  mc_samples <- NULL
  has_bt_priors <- inherits(prior, "PriorGeneric") &&
                   !is.null(prior$bayestools_priors) &&
                   inherits(prior$bayestools_priors$sigma, "prior")
  if (is_prior_only && has_bt_priors) {
    bt <- prior$bayestools_priors
    n_mc <- 2000000L
    mc_samples <- list(
      mu  = BayesTools::rng(bt$mu, n_mc),
      sig = BayesTools::rng(bt$sigma, n_mc)
    )
  }

  results <- lapply(metrics, function(m) {
    analyze_capability_integration(
      mc_samples = mc_samples,
      divergence_info = divergence_map[[m]],
      request = .new_qc_integration_request(
        data = data,
        LSL = LSL,
        USL = USL,
        prior = prior,
        metric = m,
        target = target,
        cached_state = cached_state,
        sigma_level = sigma
      )
    )
  })
  names(results) <- metrics

  # Extract coefficients (posterior means)
  coefficients <- sapply(results, function(r) r$stats["Mean"])
  names(coefficients) <- metrics

  # Build metrics list compatible with bpc structure
  # For integration, we store stats rather than samples
  metrics_list <- lapply(metrics, function(m) {
    results[[m]]$stats
  })
  names(metrics_list) <- metrics
  metric_context <- .new_qc_metric_context(
    prior = prior,
    cached_state = cached_state,
    LSL = LSL,
    USL = USL,
    target = target,
    sigma_level = sigma,
    divergence = divergence_map
  )
  metrics_list <- .new_capability_metrics(
    metrics_list,
    LSL = LSL,
    USL = USL,
    target = target,
    sigma = sigma,
    distribution = distribution,
    method = "integration",
    distributions = .entry_distributions_from_results(
      results,
      what = metrics,
      context = metric_context
    ),
    results = results,
    prior = prior,
    cached_state = cached_state,
    divergence = divergence_map
  )

  .new_integration_result(
    results = results,
    metrics = metrics_list,
    coefficients = coefficients,
    sigma = sigma,
    distribution = distribution,
    prior = prior,
    is_conjugate = is_conjugate,
    case = case,
    cached_state = cached_state,
    divergence = divergence_map
  )
}

