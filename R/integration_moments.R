# ==============================================================================
# Moment Computation for Integration Backends
# ==============================================================================

#' Compute posterior mean and SD of a capability metric
#' @param data Numeric vector of observations
#' @param LSL Lower specification limit
#' @param USL Upper specification limit
#' @param prior Prior object (PriorConjugate or PriorGeneric)
#' @param metric Capability index name
#' @param target Target value for Cpm
#' @param use_analytic Use analytic formulas (TRUE) or numerical fallback (FALSE)
#' @param cached_state Pre-computed state for PriorGeneric
#' @param sigma_level Number of process standard deviations used to define the
#'   capability metric scale.
#' @return List with mean and sd of the posterior distribution of the metric
#' @keywords internal
compute_metric_moments <- function(data, LSL, USL, prior, metric = "Cpk",
                                   target = NULL, use_analytic = TRUE,
                                   cached_state = NULL,
                                   sigma_level = 3) {
  .validate_capability_request(
    LSL = LSL,
    USL = USL,
    target = target,
    sigma_level = sigma_level,
    metric = metric,
    target_required = FALSE,
    sigma_name = "sigma_level"
  )

  UseMethod("compute_metric_moments", prior)
}

#' Compute grid bounds for a metric using sigma quantiles (prior-only case)
#'
#' When sampling from priors only (no data), we can compute grid bounds directly
#' from the sigma prior quantiles using BayesTools::quant, avoiding expensive
#' 2D numerical integration for moment computation.
#'
#' @param sigma_prior BayesTools prior object for sigma
#' @param LSL Lower specification limit
#' @param USL Upper specification limit
#' @param metric Capability index name
#' @param target Target value for Cpm
#' @param p_low Lower quantile probability (default 0.001)
#' @param p_high Upper quantile probability (default 0.999)
#' @return List with x_start and x_end for the grid
#' @keywords internal
.compute_metric_grid_bounds_from_quantiles <- function(sigma_prior, LSL, USL,
                                                        metric, target = NULL,
                                                        p_low = 0.001,
                                                        p_high = 0.999,
                                                        sigma_level = 3) {
  tol <- USL - LSL
  mid <- (LSL + USL) / 2
  if (is.null(target)) target <- mid

  # Get sigma quantiles using BayesTools::quant
  sigma_low <- BayesTools::quant(sigma_prior, p_low)
  sigma_high <- BayesTools::quant(sigma_prior, p_high)

  # For metrics inversely related to sigma (Cp, Cpk, Cpu, Cpl, Cpm, Cpc),

# metric_high corresponds to sigma_low and vice versa
  cpm_dist <- .cpm_spec_distance(LSL, USL, target)
  metric_at_sigma_low <- switch(metric,
    "Cp" = tol / ((2 * sigma_level) * sigma_low),
    "Cpk" = tol / ((2 * sigma_level) * sigma_low),  # Upper bound (assumes mu = mid)
    "Cpu" = (USL - mid) / (sigma_level * sigma_low),
    "Cpl" = (mid - LSL) / (sigma_level * sigma_low),
    "Cpm" = cpm_dist / (sigma_level * sigma_low),  # Upper bound (assumes mu = target)
    "Cpc" = tol / ((2 * sigma_level) * sigma_low),  # Upper bound approximation
    tol / ((2 * sigma_level) * sigma_low)  # Default
  )

  metric_at_sigma_high <- switch(metric,
    "Cp" = tol / ((2 * sigma_level) * sigma_high),
    "Cpk" = tol / ((2 * sigma_level) * sigma_high),
    "Cpu" = (USL - mid) / (sigma_level * sigma_high),
    "Cpl" = (mid - LSL) / (sigma_level * sigma_high),
    "Cpm" = cpm_dist / (sigma_level * sigma_high),
    "Cpc" = tol / ((2 * sigma_level) * sigma_high),
    tol / ((2 * sigma_level) * sigma_high)
  )

  x_start <- max(0, metric_at_sigma_high)
  x_end <- metric_at_sigma_low

  # Ensure valid range
  if (x_end <= x_start || !is.finite(x_end)) {
    x_end <- max(x_start + 3, 10)
  }

  list(x_start = x_start, x_end = x_end)
}

#' @export
compute_metric_moments.PriorConjugate <- function(data, LSL, USL, prior,
                                                   metric = "Cpk", target = NULL,
                                                   use_analytic = TRUE,
                                                   cached_state = NULL,
                                                   sigma_level = 3) {
  request <- .new_qc_integration_request(
    data = data,
    LSL = LSL,
    USL = USL,
    prior = prior,
    metric = metric,
    target = target,
    cached_state = cached_state,
    sigma_level = sigma_level
  )
  target <- request$target

  posterior_info <- .compute_validated_conjugate_posterior(
    prior, data, cached_state,
    context = "The conjugate posterior for metric computation"
  )
  ss <- posterior_info$ss
  post <- posterior_info$post
  n <- ss$n; k_n <- post$k_n; mu_n <- post$mu_n
  alpha_n <- post$alpha_n; beta_n <- post$beta_n

  tol <- USL - LSL
  mid <- (LSL + USL) / 2
  if (is.null(target)) target <- mid

  if (posterior_info$is_degenerate) {
    dist <- .degenerate_conjugate_metric_distribution(
      mu_n, k_n, LSL, USL, target, metric,
      sigma_level = sigma_level
    )
    if (!is.null(dist)) {
      return(.degenerate_metric_moments(dist))
    }
  }

  if (metric %in% c("Cpm", "Cpc")) {
    # Smooth conjugate metrics do not admit a stable closed-form shortcut.
    # Compute moments from the exact survival function instead of truncating
    # the diffuse joint posterior to a finite box.
    return(.positive_moments_from_survival(
      .integration_make_solver(request = request)
    ))
  }

  if (!use_analytic) {
    return(.compute_moments_numerical_conjugate(
      mu_n, k_n, alpha_n, beta_n, LSL, USL, target, metric,
      sigma_level = sigma_level
    ))
  }

  # Analytic computation
  # E[1/sigma] for Inverse-Gamma(alpha, beta): sqrt(beta) * Gamma(alpha-0.5) / Gamma(alpha) / sqrt(2)
  # Using the chi-square parameterization: sigma^2 = 2*beta/y where y ~ chi^2(2*alpha)
  # E[1/sigma] = E[sqrt(y/(2*beta))] = E[sqrt(y)] / sqrt(2*beta)
  # E[sqrt(y)] for y ~ chi^2(df) = sqrt(2) * Gamma((df+1)/2) / Gamma(df/2)
  df_p <- 2 * alpha_n
  E_inv_sigma <- sqrt(2) * exp(lgamma((df_p + 1) / 2) - lgamma(df_p / 2)) / sqrt(2 * beta_n)
  E_inv_sigma2 <- df_p / (2 * beta_n)  # E[y/(2*beta)] = E[y]/(2*beta) = df/(2*beta)

  if (metric == "Cp") {
    # Cp = tol / (2 * sigma_level * sigma)
    E1 <- (tol / (2 * sigma_level)) * E_inv_sigma
    E2 <- (tol / (2 * sigma_level))^2 * E_inv_sigma2
    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }

  if (metric == "Cpu") {
    # Cpu = (USL - mu) / (sigma_level * sigma)
    # E[Cpu] = E[(USL - mu)/sigma] / sigma_level
    # mu|sigma ~ N(mu_n, sigma^2/k_n), so E[mu|sigma] = mu_n
    # E[(USL - mu)/sigma] = (USL - mu_n) * E[1/sigma]
    E1 <- (USL - mu_n) / sigma_level * E_inv_sigma

    # E[Cpu^2] = E[(USL - mu)^2 / sigma^2] / sigma_level^2
    # (USL - mu)^2 = (USL - mu_n)^2 - 2*(USL - mu_n)*(mu - mu_n) + (mu - mu_n)^2
    # E[(mu - mu_n)^2 | sigma] = sigma^2/k_n
    # E[(USL - mu)^2 / sigma^2] = (USL - mu_n)^2 * E[1/sigma^2] + E[1/k_n] = (USL - mu_n)^2 * E[1/sigma^2] + 1/k_n
    E2 <- ((USL - mu_n)^2 * E_inv_sigma2 + 1 / k_n) / sigma_level^2
    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }

  if (metric == "Cpl") {
    # Cpl = (mu - LSL) / (sigma_level * sigma)
    E1 <- (mu_n - LSL) / sigma_level * E_inv_sigma
    E2 <- ((mu_n - LSL)^2 * E_inv_sigma2 + 1 / k_n) / sigma_level^2
    return(list(mean = E1, sd = sqrt(max(0, E2 - E1^2))))
  }

  if (metric == "Cpk") {
    return(.compute_moments_numerical_conjugate(
      mu_n, k_n, alpha_n, beta_n, LSL, USL, target, metric,
      sigma_level = sigma_level
    ))
  }

  stop("Unknown metric: ", metric)
}

#' Compute the first two moments of a non-negative metric from its survival function
#' @keywords internal
.positive_moments_from_survival <- function(S,
                                            rel.tol = 1e-5,
                                            subdivisions = 400) {
  S_vec <- function(c) {
    vals <- S(c)
    if (length(vals) != length(c)) {
      stop("The integration survival solver must return one value per input point.")
    }
    vals[!is.finite(vals) | vals < 0] <- 0
    pmin(vals, 1)
  }

  E1 <- stats::integrate(
    function(c) S_vec(c),
    lower = 0,
    upper = Inf,
    rel.tol = rel.tol,
    subdivisions = subdivisions
  )$value

  E2 <- stats::integrate(
    function(c) {
      vals <- 2 * c * S_vec(c)
      vals[!is.finite(vals)] <- 0
      vals
    },
    lower = 0,
    upper = Inf,
    rel.tol = rel.tol,
    subdivisions = subdivisions
  )$value

  list(mean = E1, sd = sqrt(max(0, E2 - E1^2)))
}

#' Conditional metric moments under mu | sigma ~ Normal(...)
#' @keywords internal
.conjugate_metric_conditional_moments <- function(sigma, mu_n, k_n,
                                                  LSL, USL, target, metric,
                                                  sigma_level = 3) {
  tol <- USL - LSL
  mid <- (LSL + USL) / 2
  if (is.null(target)) target <- mid

  sd_mu <- sigma / sqrt(k_n)
  a <- mu_n - LSL
  cc <- USL - mu_n
  inv_sigma_scale <- 1 / (sigma_level * sigma)

  switch(metric,
    "Cp" = list(
      E1 = tol * inv_sigma_scale / 2,
      E2 = (tol * inv_sigma_scale / 2)^2
    ),
    "Cpu" = list(
      E1 = cc * inv_sigma_scale,
      E2 = (cc^2 + sd_mu^2) * inv_sigma_scale^2
    ),
    "Cpl" = list(
      E1 = a * inv_sigma_scale,
      E2 = (a^2 + sd_mu^2) * inv_sigma_scale^2
    ),
    "Cpk" = {
      zs <- (mid - mu_n) / sd_mu
      Phi <- stats::pnorm(zs)
      phi <- stats::dnorm(zs)
      list(
        E1 = (a * Phi + cc * (1 - Phi) - 2 * sd_mu * phi) * inv_sigma_scale,
        E2 = ((a^2 + sd_mu^2) * Phi + (cc^2 + sd_mu^2) * (1 - Phi) -
                2 * sd_mu * (a + cc) * phi) * inv_sigma_scale^2
      )
    },
    stop("Unsupported conjugate conditional metric: ", metric)
  )
}

#' Numerical fallback for conjugate prior moments (2D integration)
#' @keywords internal
.compute_moments_numerical_conjugate <- function(mu_n, k_n, alpha_n, beta_n,
                                                  LSL, USL, target, metric,
                                                  sigma_level = 3) {
  if (.is_improper_conjugate_posterior(k_n, alpha_n, beta_n)) {
    stop(
      .conjugate_posterior_error_message(
        k_n, alpha_n, beta_n,
        context = "The conjugate posterior for numerical moment computation"
      ),
      call. = FALSE
    )
  }

  if (.is_degenerate_conjugate_posterior(beta_n)) {
    dist <- .degenerate_conjugate_metric_distribution(
      mu_n, k_n, LSL, USL, target, metric,
      sigma_level = sigma_level
    )
    if (!is.null(dist)) {
      return(.degenerate_metric_moments(dist))
    }
  }

  if (!metric %in% c("Cp", "Cpu", "Cpl", "Cpk")) {
    stop("Unsupported numerical conjugate metric: ", metric)
  }

  log_sigma_kernel <- function(t) {
    log(2) + alpha_n * log(beta_n) - lgamma(alpha_n) -
      2 * alpha_n * t - beta_n * exp(-2 * t)
  }
  t_mode <- 0.5 * (log(beta_n) - log(alpha_n))
  h_max <- log_sigma_kernel(t_mode)

  tail_prob <- 1e-10
  y_lower <- stats::qchisq(tail_prob, df = 2 * alpha_n)
  y_upper <- stats::qchisq(1 - tail_prob, df = 2 * alpha_n)
  if (!is.finite(y_lower) || y_lower <= 0 || !is.finite(y_upper) || y_upper <= y_lower) {
    t_lower <- -Inf
    t_upper <- Inf
  } else {
    t_lower <- 0.5 * log(2 * beta_n / y_upper)
    t_upper <- 0.5 * log(2 * beta_n / y_lower)
  }

  scaled_weight <- function(t) {
    vals <- exp(log_sigma_kernel(t) - h_max)
    vals[!is.finite(vals)] <- 0
    vals
  }

  moment_integrand <- function(t, power) {
    sigma <- exp(t)
    cond <- .conjugate_metric_conditional_moments(
      sigma, mu_n, k_n, LSL, USL, target, metric,
      sigma_level = sigma_level
    )
    base <- if (power == 1) cond$E1 else cond$E2
    vals <- base * scaled_weight(t)
    vals[!is.finite(vals)] <- 0
    vals
  }

  Z <- stats::integrate(
    function(t) scaled_weight(t),
    lower = t_lower,
    upper = t_upper,
    rel.tol = 1e-5,
    subdivisions = 400
  )$value

  E1 <- stats::integrate(
    function(t) moment_integrand(t, power = 1),
    lower = t_lower,
    upper = t_upper,
    rel.tol = 1e-5,
    subdivisions = 400
  )$value / Z

  E2 <- stats::integrate(
    function(t) moment_integrand(t, power = 2),
    lower = t_lower,
    upper = t_upper,
    rel.tol = 1e-5,
    subdivisions = 400
  )$value / Z

  list(mean = E1, sd = sqrt(max(0, E2 - E1^2)))
}

.semi_mu_as_generic_prior <- function(prior) {
  bt_priors <- prior$bayestools_priors %||% list()
  bt_sigma <- bt_priors$sigma %||% NULL

  log_dens <- function(mu, sigma) {
    result <- rep(-Inf, length(mu))
    valid <- is.finite(mu) & is.finite(sigma) & sigma > 0
    if (!any(valid)) {
      return(result)
    }

    sigma_v <- sigma[valid]
    log_sigma <- prior$log_dens_sigma(sigma_v)
    if (prior$k0 > 0) {
      log_mu <- stats::dnorm(
        mu[valid],
        mean = prior$mu0,
        sd = sigma_v / sqrt(prior$k0),
        log = TRUE
      )
      result[valid] <- log_sigma + log_mu
    } else {
      result[valid] <- log_sigma
    }

    result
  }

  create_prior_generic(
    log_dens_fn = log_dens,
    bayestools_priors = list(mu = bt_priors$mu %||% NULL, sigma = bt_sigma)
  )
}

.semi_sigma_as_generic_prior <- function(prior) {
  bt_priors <- prior$bayestools_priors %||% list()
  sigma_prior <- bt_priors$sigma %||%
    create_prior_conjugate(mu0 = 0, k0 = 1, alpha0 = prior$alpha0, beta0 = prior$beta0)
  log_dens_sigma <- .integration_make_sigma_log_dens_fn(sigma_prior)

  create_prior_generic(
    function(mu, sigma) {
      prior$log_dens_mu(mu) + log_dens_sigma(sigma)
    },
    bayestools_priors = c(bt_priors, list(sigma = sigma_prior))
  )
}


#' @export
compute_metric_moments.PriorSemiConjugateMu <- function(data, LSL, USL, prior,
                                                         metric = "Cpk", target = NULL,
                                                         use_analytic = TRUE,
                                                         cached_state = NULL,
                                                         sigma_level = 3) {
  request <- .new_qc_integration_request(
    data = data,
    LSL = LSL,
    USL = USL,
    prior = prior,
    metric = metric,
    target = target,
    cached_state = cached_state,
    sigma_level = sigma_level
  )
  target <- request$target

  state <- .prepare_semi_mu_state(data, prior, cached_state)
  k_n <- state$k_n
  mu_n <- state$mu_n

  spec_tol <- USL - LSL
  M <- (LSL + USL) / 2
  if (is.null(target)) target <- M

  if (!use_analytic || metric %in% c("Cpm", "Cpc")) {
    generic_prior <- .semi_mu_as_generic_prior(prior)
    return(.integration_compute_metric_moments(
      use_analytic = FALSE,
      request = .integration_request_update(
        request,
        prior = generic_prior
      )
    ))
  }

  log_sw <- function(sigma) {
    state$log_p_sigma(sigma)
  }

  t_lo <- state$t_lower
  t_hi <- state$t_upper

  base_w <- function(t) {
    sigma <- exp(t)
    vals <- exp(log_sw(sigma) + t - state$log_Z)
    vals[!is.finite(vals)] <- 0
    vals
  }

  if (!metric %in% c("Cp", "Cpu", "Cpl", "Cpk")) {
    stop("Unknown metric: ", metric)
  }

  # Smooth metrics route earlier through the generic fallback, so only the
  # closed-form family reaches this block.
  inner_moments <- function(sigma) {
    sd_mu <- sigma / sqrt(k_n)
    a <- mu_n - LSL
    cc <- USL - mu_n
    b <- sd_mu
    inv_sigma_scale <- 1 / (sigma_level * sigma)
    switch(metric,
      "Cp" = list(
        E1 = spec_tol * inv_sigma_scale / 2,
        E2 = (spec_tol * inv_sigma_scale / 2)^2
      ),
      "Cpu" = list(
        E1 = cc * inv_sigma_scale,
        E2 = (cc^2 + b^2) * inv_sigma_scale^2
      ),
      "Cpl" = list(
        E1 = a * inv_sigma_scale,
        E2 = (a^2 + b^2) * inv_sigma_scale^2
      ),
      "Cpk" = {
        zs <- (M - mu_n) / sd_mu
        Phi <- stats::pnorm(zs)
        phi <- stats::dnorm(zs)
        list(
          E1 = (a * Phi + cc * (1 - Phi) - 2 * b * phi) * inv_sigma_scale,
          E2 = ((a^2 + b^2) * Phi + (cc^2 + b^2) * (1 - Phi) -
                  2 * b * (a + cc) * phi) * inv_sigma_scale^2
        )
      }
    )
  }

  E1 <- stats::integrate(function(t) {
    m <- inner_moments(exp(t))
    base_w(t) * m$E1
  }, t_lo, t_hi, rel.tol = 1e-5)$value
  E2 <- stats::integrate(function(t) {
    m <- inner_moments(exp(t))
    base_w(t) * m$E2
  }, t_lo, t_hi, rel.tol = 1e-5)$value

  list(mean = E1, sd = sqrt(max(0, E2 - E1^2)))
}

.safe_prior_quantiles <- function(bt_prior, probs) {
  if (is.null(bt_prior) || !inherits(bt_prior, "prior")) {
    return(rep(NA_real_, length(probs)))
  }

  tryCatch(
    as.numeric(BayesTools::quant(bt_prior, probs)),
    error = function(e) rep(NA_real_, length(probs))
  )
}

.semi_sigma_mu_domain <- function(prior, n, x_bar, sse, beta0,
                                  quantile_probs = c(1e-4, 1 - 1e-4),
                                  scale_multiplier = 8) {
  mu_prior <- prior$bayestools_priors$mu %||% NULL
  support <- if (!is.null(mu_prior) && inherits(mu_prior, "prior")) {
    .extract_prior_bounds(mu_prior)
  } else {
    list(lower = -Inf, upper = Inf)
  }

  # Without prior metadata the mu prior's scale is unknown; do not invent one.
  init <- if (!is.null(mu_prior)) {
    .extract_prior_init(mu_prior)
  } else {
    list(value = x_bar, scale = NA_real_)
  }

  qs <- .safe_prior_quantiles(mu_prior, quantile_probs)
  qs_finite <- qs[is.finite(qs)]

  data_scale <- if (n > 1 && sse > 0) sqrt(sse / (n - 1)) else 0
  beta_scale <- if (beta0 > 0) {
    if (n > 0) sqrt((sse + 2 * beta0) / n) else sqrt(2 * beta0)
  } else {
    0
  }
  quantile_scale <- if (length(qs_finite) == 2L && diff(qs_finite) > 0) {
    diff(qs_finite) / 6
  } else {
    0
  }

  # The scale must follow the data's units: an absolute floor (e.g. 1) makes
  # the probe window far too coarse for measurements with sd ~ 0.01. Fall
  # back to 1 only when nothing informs the scale.
  scales <- c(init$scale, data_scale, beta_scale, quantile_scale)
  scales <- scales[is.finite(scales) & scales > 0]
  ref_scale <- if (length(scales) > 0L) max(scales) else 1
  data_radius <- scale_multiplier * max(data_scale, beta_scale, 0)
  prior_radius <- scale_multiplier * ref_scale

  finite_lower <- c(qs[1], x_bar - data_radius)
  finite_upper <- c(qs[2], x_bar + data_radius)
  if (!all(is.finite(qs))) {
    finite_lower <- c(finite_lower, init$value - prior_radius)
    finite_upper <- c(finite_upper, init$value + prior_radius)
  }
  finite_lower <- finite_lower[is.finite(finite_lower)]
  finite_upper <- finite_upper[is.finite(finite_upper)]

  lower <- if (is.finite(support$lower)) support$lower else min(finite_lower)
  upper <- if (is.finite(support$upper)) support$upper else max(finite_upper)

  if (!is.finite(lower)) lower <- x_bar - prior_radius
  if (!is.finite(upper)) upper <- x_bar + prior_radius

  lower <- max(lower, support$lower)
  upper <- min(upper, support$upper)

  if (!is.finite(lower) || !is.finite(upper) || upper <= lower) {
    center_candidates <- c(x_bar, init$value, qs_finite)
    center_candidates <- center_candidates[is.finite(center_candidates)]
    center <- if (length(center_candidates) > 0L) stats::median(center_candidates) else 0
    radius <- ref_scale
    lower <- max(center - radius, support$lower)
    upper <- min(center + radius, support$upper)
  }

  list(lower = lower, upper = upper, scale = ref_scale)
}

.semi_sigma_probe_interval <- function(mu_domain, lower = -Inf, upper = Inf) {
  probe_lower <- mu_domain$lower
  probe_upper <- mu_domain$upper

  if (!is.finite(probe_lower) || !is.finite(probe_upper) || probe_upper <= probe_lower) {
    radius <- max(mu_domain$scale, 1, na.rm = TRUE)
    probe_lower <- -radius
    probe_upper <- radius
  }

  width <- probe_upper - probe_lower
  if (!is.finite(width) || width <= 0) {
    width <- max(mu_domain$scale, 1, na.rm = TRUE)
  }

  if (is.finite(lower) && probe_upper < lower) {
    probe_lower <- lower
    probe_upper <- if (is.finite(upper)) min(upper, lower + width) else lower + width
  } else if (is.finite(upper) && probe_lower > upper) {
    probe_upper <- upper
    probe_lower <- if (is.finite(lower)) max(lower, upper - width) else upper - width
  } else {
    if (is.finite(lower)) probe_lower <- max(probe_lower, lower)
    if (is.finite(upper)) probe_upper <- min(probe_upper, upper)
  }

  if (!is.finite(probe_lower) || !is.finite(probe_upper) || probe_upper <= probe_lower) {
    if (is.finite(lower) && is.finite(upper) && upper > lower) {
      probe_lower <- lower
      probe_upper <- upper
    } else {
      radius <- max(mu_domain$scale, 1, na.rm = TRUE)
      center <- if (is.finite(lower) && is.finite(upper)) {
        (lower + upper) / 2
      } else if (is.finite(lower)) {
        lower + radius
      } else if (is.finite(upper)) {
        upper - radius
      } else {
        0
      }
      probe_lower <- if (is.finite(lower)) lower else center - radius
      probe_upper <- if (is.finite(upper)) upper else center + radius
    }
  }

  c(lower = probe_lower, upper = probe_upper)
}

.semi_sigma_find_log_peak <- function(log_integrand, mu_domain,
                                      lower = -Inf, upper = Inf,
                                      probe_n = 1025L, max_expansions = 8L) {
  bounds <- .semi_sigma_probe_interval(mu_domain, lower, upper)

  for (iter in seq_len(max_expansions + 1L)) {
    probe_mu <- seq(bounds[["lower"]], bounds[["upper"]], length.out = probe_n)
    probe_log <- log_integrand(probe_mu)
    finite_probe <- is.finite(probe_log)

    # The probe window can still be widened on a side it has not yet pushed
    # to the integration bound.
    can_left <- bounds[["lower"]] > lower
    can_right <- bounds[["upper"]] < upper

    if (any(finite_probe)) {
      finite_idx <- which(finite_probe)
      best_idx <- finite_idx[which.max(probe_log[finite_probe])]
      h_max <- probe_log[best_idx]
      mode <- probe_mu[best_idx]

      local_lower <- probe_mu[max(best_idx - 1L, 1L)]
      local_upper <- probe_mu[min(best_idx + 1L, length(probe_mu))]
      if (local_upper > local_lower) {
        opt <- tryCatch(
          stats::optimize(
            function(mu) {
              val <- log_integrand(mu)
              if (is.finite(val)) val else -.Machine$double.xmax
            },
            interval = c(local_lower, local_upper),
            maximum = TRUE
          ),
          error = function(e) NULL
        )
        if (!is.null(opt) && is.finite(opt$objective) && opt$objective > h_max) {
          h_max <- opt$objective
          mode <- opt$maximum
        }
      }

      at_left_edge <- best_idx <= 2L && can_left
      at_right_edge <- best_idx >= length(probe_mu) - 1L && can_right
      if (!at_left_edge && !at_right_edge) {
        return(list(h_max = h_max, mode = mode))
      }
    } else {
      if (!can_left && !can_right) break
      at_left_edge <- can_left
      at_right_edge <- can_right
    }

    width <- bounds[["upper"]] - bounds[["lower"]]
    if (!is.finite(width) || width <= 0) {
      width <- max(mu_domain$scale, 1, na.rm = TRUE)
    }
    if (at_left_edge) bounds[["lower"]] <- max(lower, bounds[["lower"]] - width)
    if (at_right_edge) bounds[["upper"]] <- min(upper, bounds[["upper"]] + width)
  }

  NULL
}

# Locate the mass of the mu marginal of the semi-conjugate sigma posterior,
#   log w(mu) = log p(mu) - alpha_n * log(beta0 + (sse + n * (mu - x_bar)^2) / 2),
# in the units of the data. With n observations w is roughly a Student-t with
# scale sqrt((beta0 + sse / 2) / (alpha_n * n)), i.e. about s / sqrt(n), which
# can be orders of magnitude narrower than the mu prior. Returns the mode,
# a per-side local scale (within a factor of two of the half-width at which
# log w drops by 1/2), and the range outside which log w has dropped by
# `tail_log_drop` (infinite when a heavy tail never does).
.semi_sigma_mu_frame <- function(log_weight, mu_domain, n, x_bar, sse, beta0,
                                 alpha_n_eff, tail_log_drop = 40,
                                 max_steps = 64L) {
  eval_lw <- function(mu) {
    vals <- as.numeric(log_weight(mu))
    vals[is.na(vals)] <- -Inf
    vals
  }

  beta_floor <- beta0 + sse / 2
  lik_scale <- if (n > 0 && beta_floor > 0) {
    sqrt(beta_floor / (alpha_n_eff * n))
  } else {
    NA_real_
  }

  domain <- .semi_sigma_probe_interval(mu_domain)
  probe <- seq(domain[["lower"]], domain[["upper"]], length.out = 1025L)
  if (is.finite(lik_scale) && lik_scale > 0) {
    probe <- c(probe, x_bar + lik_scale * seq(-64, 64, length.out = 513L))
  }
  probe <- sort(unique(probe[is.finite(probe)]))
  probe_lw <- eval_lw(probe)
  if (!any(is.finite(probe_lw))) {
    return(NULL)
  }

  best <- which.max(probe_lw)
  h_max <- probe_lw[best]
  bracket <- probe[c(max(best - 1L, 1L), min(best + 1L, length(probe)))]

  # A best probe on the window edge means the kernel is still rising there.
  # March outward until it turns down so the mode is bracketed.
  if (best == 1L || best == length(probe)) {
    dir <- if (best == 1L) -1 else 1
    step <- if (is.finite(lik_scale) && lik_scale > 0) {
      lik_scale
    } else {
      diff(range(probe)) / 64
    }
    prev <- probe[best]
    prev_lw <- h_max
    for (i in seq_len(max_steps)) {
      cand <- prev + dir * step
      cand_lw <- eval_lw(cand)
      if (!is.finite(cand_lw) || cand_lw <= prev_lw) {
        bracket <- sort(c(prev - dir * step / 2, cand))
        break
      }
      bracket <- sort(c(prev, cand))
      prev <- cand
      prev_lw <- cand_lw
      step <- 2 * step
    }
    h_max <- prev_lw
    best_mu <- prev
  } else {
    best_mu <- probe[best]
  }

  mode <- best_mu
  if (bracket[2] > bracket[1]) {
    opt <- tryCatch(
      stats::optimize(
        function(mu) {
          val <- eval_lw(mu)
          if (is.finite(val)) val else -.Machine$double.xmax
        },
        interval = bracket,
        maximum = TRUE,
        tol = 1e-10 * max(abs(bracket), 1e-300)
      ),
      error = function(e) NULL
    )
    if (!is.null(opt) && is.finite(opt$objective) && opt$objective > h_max) {
      h_max <- opt$objective
      mode <- opt$maximum
    }
  }

  drop_at <- function(dir, d) {
    val <- eval_lw(mode + dir * d)
    if (is.finite(val)) h_max - val else Inf
  }

  # Find d with drop(d / 2) < 1/2 <= drop(d), starting from the data scale.
  side_scale <- function(dir) {
    d <- if (is.finite(lik_scale) && lik_scale > 0) {
      lik_scale
    } else {
      diff(range(probe)) / 1024
    }
    if (drop_at(dir, d) >= 0.5) {
      for (i in seq_len(60L)) {
        if (drop_at(dir, d / 2) < 0.5) break
        d <- d / 2
      }
    } else {
      for (i in seq_len(max_steps)) {
        d <- 2 * d
        if (drop_at(dir, d) >= 0.5) break
      }
    }
    d
  }
  scale <- c(lower = side_scale(-1), upper = side_scale(1))

  tail_bound <- function(dir, d) {
    for (i in seq_len(max_steps)) {
      if (drop_at(dir, d) >= tail_log_drop) {
        return(mode + dir * d)
      }
      d <- 2 * d
    }
    dir * Inf
  }
  lower <- tail_bound(-1, scale[["lower"]])
  upper <- tail_bound(1, scale[["upper"]])

  # Keep any separate mode the probe saw above the tail threshold.
  heavy <- which(probe_lw > h_max - tail_log_drop)
  if (length(heavy) > 0L) {
    lower <- min(lower, probe[max(min(heavy) - 1L, 1L)])
    upper <- max(upper, probe[min(max(heavy) + 1L, length(probe))])
  }

  list(
    mode = mode,
    h_max = h_max,
    scale = scale,
    lower = lower,
    upper = upper
  )
}

# Quadrature breakpoints anchored on the frame's mode and local scale, so a
# narrow peak is never left between Gauss-Kronrod nodes of a wide interval.
.semi_sigma_mu_breaks <- function(frame, lower = -Inf, upper = Inf) {
  a <- max(lower, frame$lower)
  b <- min(upper, frame$upper)
  if (is.na(a) || is.na(b) || !(b > a)) {
    return(NULL)
  }

  offsets <- c(1, 4, 16, 64, 256)
  pts <- c(
    frame$mode,
    frame$mode - offsets * frame$scale[["lower"]],
    frame$mode + offsets * frame$scale[["upper"]]
  )
  pts <- pts[is.finite(pts) & pts > a & pts < b]
  sort(unique(c(a, pts, b)))
}

.semi_sigma_integrate_breaks <- function(fn, breaks,
                                         rel.tol = 1e-8,
                                         subdivisions = 400L) {
  if (is.null(breaks)) {
    return(0)
  }

  integrate_piece <- function(lower, upper, tol) {
    tryCatch(
      stats::integrate(
        fn, lower, upper,
        rel.tol = tol,
        subdivisions = subdivisions,
        stop.on.error = FALSE
      ),
      error = function(e) NULL
    )
  }

  total <- 0
  for (k in seq_len(length(breaks) - 1L)) {
    res <- integrate_piece(breaks[k], breaks[k + 1L], rel.tol)
    # Integrands built from interpolated geometry (e.g. Cpc) are not smooth
    # enough for a tight tolerance; QUADPACK then flags roundoff. Retry looser.
    if (is.null(res) || !identical(res$message, "OK")) {
      res <- integrate_piece(breaks[k], breaks[k + 1L], max(rel.tol, 1e-5))
    }
    ok <- !is.null(res) && is.finite(res$value) &&
      !grepl("divergent|non-finite", res$message)
    piece <- if (ok) res$value else NA_real_
    total <- total + piece
  }
  total
}

# Integrate exp(log_integrand(mu)) over [lower, upper] within the frame's
# effective range. The integrand is rescaled by its own maximum over probes
# placed between the breakpoints, so a factor that shifts mass away from the
# frame's mode (e.g. a sigma-interval probability) does not underflow.
.semi_sigma_integrate_log_mu <- function(log_integrand, frame,
                                         lower = -Inf, upper = Inf,
                                         rel.tol = 1e-8,
                                         subdivisions = 400L) {
  breaks <- .semi_sigma_mu_breaks(frame, lower, upper)
  if (is.null(breaks)) {
    return(list(value = 0, log_value = -Inf))
  }

  finite_breaks <- breaks[is.finite(breaks)]
  probe <- finite_breaks
  if (length(finite_breaks) > 1L) {
    probe <- unlist(lapply(seq_len(length(finite_breaks) - 1L), function(k) {
      seq(finite_breaks[k], finite_breaks[k + 1L], length.out = 17L)
    }))
  }
  probe_vals <- as.numeric(log_integrand(probe))
  probe_vals <- probe_vals[is.finite(probe_vals)]
  if (length(probe_vals) == 0L) {
    return(list(value = 0, log_value = -Inf))
  }
  h <- max(probe_vals)

  scaled <- .semi_sigma_integrate_breaks(
    function(mu) {
      vals <- exp(log_integrand(mu) - h)
      vals[!is.finite(vals)] <- 0
      vals
    },
    breaks,
    rel.tol = rel.tol,
    subdivisions = subdivisions
  )

  if (!is.finite(scaled)) {
    return(list(value = NA_real_, log_value = NA_real_))
  }
  if (scaled <= 0) {
    return(list(value = 0, log_value = -Inf))
  }

  log_value <- log(scaled) + h
  list(value = exp(log_value), log_value = log_value)
}

# E[g(mu, sigma)^k | mu] for k = 1, 2 under sigma^2 | mu ~ InvGamma(alpha,
# beta_n(mu)), on a fixed log(sigma) grid around each conditional mode. The
# log(sigma) density is smooth and decays double-exponentially on the left
# and at rate 2 * alpha on the right, so a self-normalized trapezoid rule on
# that grid is accurate.
.semi_sigma_conditional_metric_moments <- function(mu, beta_n, alpha, metric_fn,
                                                   n_points = 401L) {
  sd_t <- 1 / (2 * sqrt(alpha))
  offsets <- seq(-12 * sd_t, max(12 * sd_t, 50 / (2 * alpha)),
                 length.out = n_points)
  t <- outer(0.5 * log(beta_n / alpha), offsets, "+")
  log_p <- -2 * alpha * t - beta_n * exp(-2 * t)
  w <- exp(log_p - apply(log_p, 1L, max))
  w[, c(1L, n_points)] <- w[, c(1L, n_points)] / 2
  w <- w / rowSums(w)

  vals <- metric_fn(matrix(mu, nrow = length(mu), ncol = n_points), exp(t))
  list(E1 = rowSums(w * vals), E2 = rowSums(w * vals^2))
}

.semi_sigma_support_contains_point <- function(log_dens_mu, point) {
  probe_point <- function(x) {
    val <- tryCatch(log_dens_mu(x), error = function(e) NA_real_)
    val <- as.numeric(val)[1]
    !is.na(val) && val > -Inf
  }

  if (probe_point(point)) {
    return(TRUE)
  }

  step_scale <- max(1, abs(point))
  probe_steps <- step_scale * sqrt(.Machine$double.eps) * 2^seq(0, 8)
  for (step in probe_steps) {
    if (probe_point(point - step) || probe_point(point + step)) {
      return(TRUE)
    }
  }

  FALSE
}

.validate_semi_sigma_posterior <- function(n, x_bar, sse, alpha0, beta0,
                                           log_dens_mu, jeffreys_adj = 0,
                                           context = "The semi-conjugate sigma posterior") {
  alpha_n_eff <- alpha0 + n / 2 + jeffreys_adj
  beta_floor <- beta0 + sse / 2

  if (!is.finite(alpha_n_eff) || alpha_n_eff <= 0) {
    stop(
      sprintf(
        "%s is improper for the supplied data and prior (alpha_n = %s).",
        context,
        .format_conjugate_posterior_parameter(alpha_n_eff)
      ),
      call. = FALSE
    )
  }

  if (!is.finite(beta_floor) || beta_floor < 0) {
    stop(
      sprintf(
        "%s is improper for the supplied data and prior (beta_n = %s).",
        context,
        .format_conjugate_posterior_parameter(beta_floor)
      ),
      call. = FALSE
    )
  }

  zero_singularity_in_support <- beta_floor == 0 &&
    (n == 0 || .semi_sigma_support_contains_point(log_dens_mu, x_bar))
  if (zero_singularity_in_support && alpha_n_eff >= 0.5) {
    stop(
      sprintf(
        "%s is improper because beta_n reaches zero inside the mu support (alpha_n = %s).",
        context,
        .format_conjugate_posterior_parameter(alpha_n_eff)
      ),
      call. = FALSE
    )
  }

  invisible(alpha_n_eff)
}

.prepare_semi_sigma_state <- function(data, prior, cached_state = NULL,
                                      context = "The semi-conjugate sigma posterior") {
  ss <- .extract_suff_stats(data, cached_state)
  n <- ss$n
  x_bar <- ss$x_bar
  sse <- ss$SS
  alpha0 <- prior$alpha0
  beta0 <- prior$beta0
  jeffreys_adj <- if (alpha0 == -0.5 && beta0 == 0) 0.5 else 0
  alpha_n <- alpha0 + n / 2
  alpha_n_eff <- .validate_semi_sigma_posterior(
    n, x_bar, sse, alpha0, beta0, prior$log_dens_mu,
    jeffreys_adj = jeffreys_adj,
    context = context
  )
  mu_domain <- .semi_sigma_mu_domain(prior, n, x_bar, sse, beta0)

  # Unnormalized log marginal of mu after integrating sigma analytically.
  log_weight <- function(mu) {
    beta_n <- beta0 + (sse + n * (mu - x_bar)^2) / 2
    vals <- prior$log_dens_mu(mu) - alpha_n_eff * log(beta_n)
    vals[!is.finite(vals)] <- -Inf
    vals
  }

  mu_frame <- .semi_sigma_mu_frame(
    log_weight, mu_domain, n, x_bar, sse, beta0, alpha_n_eff
  )
  if (is.null(mu_frame) || !is.finite(mu_frame$h_max)) {
    stop(
      sprintf("%s is improper for the supplied data and prior.", context),
      call. = FALSE
    )
  }
  # The frame doubles its tail search out to ~2^64 local scales; a marginal
  # that has still not decayed by then (e.g. a flat mu prior without data)
  # cannot be normalized.
  if (!is.finite(mu_frame$lower) || !is.finite(mu_frame$upper)) {
    stop(
      sprintf(
        "%s is improper for the supplied data and prior (the mu marginal does not decay).",
        context
      ),
      call. = FALSE
    )
  }

  log_Z <- .semi_sigma_integrate_log_mu(log_weight, mu_frame)$log_value
  if (!is.finite(log_Z)) {
    stop(
      sprintf("%s could not be normalized numerically.", context),
      call. = FALSE
    )
  }

  # Conditional sigma scale at the mu mode: sigma^2 | mu ~ InvGamma(alpha_n,
  # beta_n(mu)) has mode beta_n / (alpha_n + 1); sqrt(beta_n / alpha_n) is a
  # matching reference in the data's units.
  beta_at_mode <- beta0 + (sse + n * (mu_frame$mode - x_bar)^2) / 2
  sigma_scale <- sqrt(beta_at_mode / alpha_n_eff)
  if (!is.finite(sigma_scale) || sigma_scale <= 0) {
    sigma_scale <- mu_domain$scale
  }

  list(
    n = n,
    x_bar = x_bar,
    sse = sse,
    alpha0 = alpha0,
    beta0 = beta0,
    alpha_n = alpha_n,
    alpha_n_eff = alpha_n_eff,
    jeffreys_adj = jeffreys_adj,
    alpha_0_times_logbeta0 = if (beta0 == 0) 0 else alpha0 * log(beta0),
    mu_domain = mu_domain,
    log_weight = log_weight,
    mu_frame = mu_frame,
    log_Z = log_Z,
    sigma_scale = sigma_scale
  )
}

#' @export
compute_metric_moments.PriorSemiConjugateSigma <- function(data, LSL, USL, prior,
                                                            metric = "Cpk", target = NULL,
                                                            use_analytic = TRUE,
                                                            cached_state = NULL,
                                                            sigma_level = 3) {
  request <- .new_qc_integration_request(
    data = data,
    LSL = LSL,
    USL = USL,
    prior = prior,
    metric = metric,
    target = target,
    cached_state = cached_state,
    sigma_level = sigma_level
  )
  target <- request$target

  # Case 3: non-conjugate mu, conjugate sigma (InvGamma/Jeffreys).
  # use_analytic = TRUE keeps the semi-analytic reduction over mu and uses
  # exact conditional sigma moments where available. use_analytic = FALSE
  # falls back to the generic 2D reference path.

  if (!use_analytic) {
    return(.integration_compute_metric_moments(
      use_analytic = FALSE,
      request = .integration_request_update(
        request,
        prior = .semi_sigma_as_generic_prior(prior)
      )
    ))
  }

  state <- .prepare_semi_sigma_state(
    data, prior, cached_state,
    context = "The semi-conjugate sigma posterior for metric computation"
  )
  n <- state$n
  x_bar <- state$x_bar
  sse <- state$sse
  beta0 <- state$beta0
  alpha_n_eff <- state$alpha_n_eff
  log_weight <- state$log_weight
  frame <- state$mu_frame
  h_max <- frame$h_max

  if (!metric %in% c("Cp", "Cpu", "Cpl", "Cpk", "Cpm", "Cpc")) {
    stop("Unknown metric: ", metric)
  }

  # E[metric^k | mu] after integrating sigma^2 | mu ~ InvGamma(alpha_n, beta_n).
  metric_conditional_moments <- function(mu, beta_n) {
    if (metric %in% c("Cpm", "Cpc")) {
      cpm_dist <- .cpm_spec_distance(LSL, USL, target)
      tol <- USL - LSL
      metric_fn <- if (metric == "Cpm") {
        function(mu, sigma) cpm_dist / (sigma_level * sqrt(sigma^2 + (mu - target)^2))
      } else {
        function(mu, sigma) {
          tol / ((2 * sigma_level) * sqrt(pi / 2) *
                   .cpc_E_abs_dev_normal(mu, sigma, target))
        }
      }
      return(.semi_sigma_conditional_metric_moments(
        mu, beta_n, alpha_n_eff, metric_fn
      ))
    }

    E_inv_sigma <- exp(
      lgamma(alpha_n_eff + 0.5) - lgamma(alpha_n_eff) - 0.5 * log(beta_n)
    )
    E_inv_sigma2 <- alpha_n_eff / beta_n

    numerator <- switch(metric,
      "Cp" = rep((USL - LSL) / (2 * sigma_level), length(mu)),
      "Cpu" = (USL - mu) / sigma_level,
      "Cpl" = (mu - LSL) / sigma_level,
      "Cpk" = pmin(USL - mu, mu - LSL) / sigma_level
    )

    list(
      E1 = numerator * E_inv_sigma,
      E2 = numerator^2 * E_inv_sigma2
    )
  }

  integrand_Ek <- function(mu, power) {
    beta_n <- beta0 + (sse + n * (mu - x_bar)^2) / 2
    scaled_weight <- exp(log_weight(mu) - h_max)
    result <- numeric(length(mu))
    ok <- is.finite(scaled_weight) & scaled_weight > 0 &
      is.finite(beta_n) & beta_n > 0
    if (!any(ok)) {
      return(result)
    }

    if (power == 0) {
      result[ok] <- scaled_weight[ok]
      return(result)
    }

    inner <- metric_conditional_moments(mu[ok], beta_n[ok])
    metric_moment <- if (power == 1) inner$E1 else inner$E2
    result[ok] <- scaled_weight[ok] * metric_moment
    result[!is.finite(result)] <- 0
    result
  }

  # Normalize with the same breakpoints so quadrature error cancels in the
  # ratios.
  breaks <- .semi_sigma_mu_breaks(frame)
  moment_integral <- function(power) {
    .semi_sigma_integrate_breaks(
      function(mu) integrand_Ek(mu, power),
      breaks,
      rel.tol = 1e-8
    )
  }
  Z_scaled <- moment_integral(0)
  if (!is.finite(Z_scaled) || Z_scaled <= 0) {
    stop(
      "The semi-conjugate sigma posterior for metric computation could not be normalized numerically.",
      call. = FALSE
    )
  }
  E1 <- moment_integral(1) / Z_scaled
  E2 <- moment_integral(2) / Z_scaled

  list(mean = E1, sd = sqrt(max(0, E2 - E1^2)))
}

#' @export
compute_metric_moments.PriorGeneric <- function(data, LSL, USL, prior,
                                                 metric = "Cpk", target = NULL,
                                                 use_analytic = TRUE,
                                                 cached_state = NULL,
                                                 sigma_level = 3) {
  # For generic priors, always use numerical 2D integration
  # (use_analytic is ignored - no closed form available)

  cached_state <- precompute_generic_state(data, prior, cached_state = cached_state)

  Z <- cached_state$Z

  metric_power <- function(power) {
    function(mu, sigma) {
      compute_metric_value(
        mu, sigma, LSL, USL, target, metric,
        sigma_level = sigma_level
      )^power
    }
  }

  # The cached integrator (MAP-anchored nested quadrature with data, a finite
  # box for prior-only state) is the one that defines Z, so the h_max
  # normalization cancels.
  E1 <- cached_state$int_2d(Inf, NULL, metric_power(1)) / Z
  E2 <- cached_state$int_2d(Inf, NULL, metric_power(2)) / Z

  list(mean = E1, sd = sqrt(max(0, E2 - E1^2)))
}


# ==============================================================================
# Solver Factory: Returns function P(Index > c)
# ==============================================================================


