.new_bpc_fit <- function(call,
                         method,
                         distribution,
                         metrics = NULL,
                         coefficients = NULL,
                         sigma = NULL,
                         prior = NULL,
                         prior_resolved = NULL,
                         prior_map = NULL,
                         stan_data = NULL,
                         stan_priors = NULL,
                         control = NULL,
                         stanfit = NULL,
                         integration_result = NULL) {
  out <- list(
    call = call,
    method = method,
    distribution = distribution,
    metrics = metrics,
    coefficients = coefficients
  )

  if (!is.null(sigma)) {
    out$sigma <- sigma
  }
  if (!is.null(prior)) {
    out$prior <- prior
  }
  if (!is.null(prior_resolved)) {
    out$prior_resolved <- prior_resolved
  }
  if (!is.null(prior_map)) {
    out$prior_map <- prior_map
  }
  if (!is.null(stan_data)) {
    out$stan_data <- stan_data
  }
  if (!is.null(stan_priors)) {
    out$stan_priors <- stan_priors
  }
  if (!is.null(control)) {
    out$control <- control
  }
  if (!is.null(stanfit)) {
    out$stanfit <- stanfit
  }
  if (!is.null(integration_result)) {
    out$integration_result <- integration_result
  }

  class(out) <- "bpc"
  out
}

.new_pc_fit <- function(call,
                        distribution,
                        data,
                        control,
                        fit,
                        metrics = NULL,
                        metrics_boot = NULL,
                        coefficients = NULL,
                        sigma = NULL) {
  out <- list(
    call = call,
    distribution = distribution,
    data = data,
    control = control,
    fit = fit,
    metrics = metrics,
    metrics_boot = metrics_boot,
    coefficients = coefficients
  )

  if (!is.null(sigma)) {
    out$sigma <- sigma
  }

  class(out) <- "pc"
  out
}

.new_integration_result <- function(results,
                                    metrics,
                                    coefficients,
                                    sigma,
                                    distribution = NULL,
                                    prior,
                                    is_conjugate,
                                    case,
                                    cached_state,
                                    divergence) {
  structure(
    list(
      results = results,
      metrics = metrics,
      coefficients = coefficients,
      sigma = sigma,
      distribution = distribution,
      prior = prior,
      is_conjugate = is_conjugate,
      case = case,
      cached_state = cached_state,
      divergence = divergence
    ),
    class = "qc_integration_result"
  )
}

.new_capability_metrics <- function(metrics,
                                    LSL,
                                    USL,
                                    target,
                                    sigma = NULL,
                                    distribution = NULL,
                                    method = NULL,
                                    distributions = NULL,
                                    results = NULL,
                                    prior = NULL,
                                    cached_state = NULL,
                                    divergence = NULL) {
  class(metrics) <- "capability_metrics"
  attr(metrics, "LSL") <- LSL
  attr(metrics, "USL") <- USL
  attr(metrics, "target") <- target

  if (!is.null(sigma)) {
    attr(metrics, "sigma") <- sigma
  }
  if (!is.null(distribution)) {
    attr(metrics, "distribution") <- distribution
  }
  if (!is.null(method)) {
    attr(metrics, "method") <- method
  }
  if (!is.null(distributions)) {
    attr(metrics, "distributions") <- distributions
  }
  if (!is.null(results)) {
    attr(metrics, "results") <- results
  }
  if (!is.null(prior)) {
    attr(metrics, "prior") <- prior
  }
  if (!is.null(cached_state)) {
    attr(metrics, "cached_state") <- cached_state
  }
  if (!is.null(divergence)) {
    attr(metrics, "divergence") <- divergence
  }

  metrics
}

.new_bpc_summary <- function(call,
                             metrics,
                             summary,
                             interval_summary,
                             integration_result = NULL,
                             divergence_diagnostics = NULL) {
  out <- list(
    call = call,
    summary = summary,
    interval_summary = interval_summary,
    metrics = metrics
  )

  if (!is.null(integration_result)) {
    out$integration_result <- integration_result
  }

  out$divergence_diagnostics <- divergence_diagnostics
  attr(out, "has_divergent_moments") <- !is.null(divergence_diagnostics) &&
    nrow(divergence_diagnostics) > 0
  class(out) <- "bpc_summary"
  out
}

.new_pc_summary <- function(call,
                            metrics,
                            metrics_boot,
                            summary,
                            interval_summary) {
  out <- list(
    call = call,
    summary = summary,
    interval_summary = interval_summary,
    metrics = metrics,
    metrics_boot = metrics_boot
  )

  attr(out, "has_bootstrap") <- !is.null(metrics_boot)
  class(out) <- "pc_summary"
  out
}
