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

  init <- if (!is.null(mu_prior)) {
    .extract_prior_init(mu_prior)
  } else {
    list(value = x_bar, scale = 1)
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

  ref_scale <- max(c(init$scale, data_scale, beta_scale, quantile_scale, 1), na.rm = TRUE)
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
    radius <- max(ref_scale, 1)
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

    if (any(finite_probe)) {
      finite_idx <- which(finite_probe)
      best_idx <- finite_idx[which.max(probe_log[finite_probe])]
      h_max <- probe_log[best_idx]

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
        if (!is.null(opt) && is.finite(opt$objective)) {
          h_max <- max(h_max, opt$objective)
        }
      }

      at_left_edge <- best_idx <= 2L && !is.finite(lower)
      at_right_edge <- best_idx >= length(probe_mu) - 1L && !is.finite(upper)
      if (!at_left_edge && !at_right_edge) {
        return(list(h_max = h_max))
      }

      width <- bounds[["upper"]] - bounds[["lower"]]
      if (!is.finite(width) || width <= 0) {
        width <- max(mu_domain$scale, 1, na.rm = TRUE)
      }
      if (at_left_edge) bounds[["lower"]] <- bounds[["lower"]] - width
      if (at_right_edge) bounds[["upper"]] <- bounds[["upper"]] + width
      next
    }

    width <- bounds[["upper"]] - bounds[["lower"]]
    if (!is.finite(width) || width <= 0) {
      width <- max(mu_domain$scale, 1, na.rm = TRUE)
    }
    if (!is.finite(lower)) bounds[["lower"]] <- bounds[["lower"]] - width
    if (!is.finite(upper)) bounds[["upper"]] <- bounds[["upper"]] + width
  }

  NULL
}

.semi_sigma_integrate_log_kernel <- function(log_integrand, mu_domain,
                                             lower = -Inf, upper = Inf,
                                             rel.tol = 1e-6,
                                             subdivisions = 400,
                                             context = NULL) {
  peak <- .semi_sigma_find_log_peak(
    log_integrand, mu_domain,
    lower = lower, upper = upper
  )

  if (is.null(peak) || !is.finite(peak$h_max)) {
    if (is.null(context)) {
      return(list(value = 0, log_value = -Inf, h_max = -Inf))
    }
    stop(
      sprintf("%s is improper for the supplied data and prior.", context),
      call. = FALSE
    )
  }

  scaled <- tryCatch(
    stats::integrate(
      function(mu) {
        vals <- exp(log_integrand(mu) - peak$h_max)
        vals[!is.finite(vals)] <- 0
        vals
      },
      lower = lower,
      upper = upper,
      rel.tol = rel.tol,
      subdivisions = subdivisions
    )$value,
    error = function(e) NA_real_
  )

  if (!is.finite(scaled) || scaled <= 0) {
    if (is.null(context)) {
      return(list(value = 0, log_value = -Inf, h_max = peak$h_max))
    }
    stop(
      sprintf("%s could not be normalized numerically.", context),
      call. = FALSE
    )
  }

  log_value <- log(scaled) + peak$h_max
  list(value = exp(log_value), log_value = log_value, h_max = peak$h_max)
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

  sigma_scale_num <- sse + 2 * max(beta0, 0)
  if (n > 0) {
    sigma_scale_num <- sigma_scale_num + n * mu_domain$scale^2
  }
  sigma_scale <- sqrt(max(sigma_scale_num, 0) / max(n + 1, 1))
  if (!is.finite(sigma_scale) || sigma_scale <= 0) {
    sigma_scale <- max(mu_domain$scale, 1)
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

  if (metric %in% c("Cpm", "Cpc")) {
    return(.integration_compute_metric_moments(
      use_analytic = FALSE,
      request = .integration_request_update(
        request,
        prior = .semi_sigma_as_generic_prior(prior)
      )
    ))
  }

  if (!metric %in% c("Cp", "Cpu", "Cpl", "Cpk")) {
    stop("Unknown metric: ", metric)
  }

  log_kernel <- function(mu) {
    sse_mu <- sse + n * (mu - x_bar)^2
    beta_n <- beta0 + sse_mu / 2
    vals <- prior$log_dens_mu(mu) - alpha_n_eff * log(beta_n)
    vals[!is.finite(vals)] <- -Inf
    vals
  }

  weight_info <- .semi_sigma_integrate_log_kernel(
    log_kernel,
    state$mu_domain,
    lower = -Inf,
    upper = Inf,
    rel.tol = 1e-5,
    subdivisions = 400,
    context = "The semi-conjugate sigma posterior for metric computation"
  )
  log_Z <- weight_info$log_value
  h_max <- weight_info$h_max

  metric_conditional_moments <- function(mu, beta_n) {
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
    sse_mu <- sse + n * (mu - x_bar)^2
    beta_n <- beta0 + sse_mu / 2
    scaled_weight <- exp(log_kernel(mu) - h_max)
    result <- numeric(length(mu))
    ok <- is.finite(scaled_weight) & scaled_weight > 0 &
      is.finite(beta_n) & beta_n > 0
    if (!any(ok)) {
      return(result)
    }

    inner <- metric_conditional_moments(mu[ok], beta_n[ok])
    metric_moment <- if (power == 1) inner$E1 else inner$E2
    result[ok] <- scaled_weight[ok] * metric_moment
    result[!is.finite(result)] <- 0
    result
  }

  scale_back <- exp(h_max - log_Z)
  E1 <- scale_back * stats::integrate(
    function(mu) integrand_Ek(mu, 1),
    lower = -Inf,
    upper = Inf,
    rel.tol = 1e-5,
    subdivisions = 400
  )$value
  E2 <- scale_back * stats::integrate(
    function(mu) integrand_Ek(mu, 2),
    lower = -Inf,
    upper = Inf,
    rel.tol = 1e-5,
    subdivisions = 400
  )$value

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

  log_post_vec <- cached_state$log_post_vec
  h_max <- cached_state$h_max
  Z <- cached_state$Z
  mu_lower <- cached_state$mu_lower
  mu_upper <- cached_state$mu_upper
  sigma_lower <- max(cached_state$sigma_lower, 1e-10)
  sigma_upper <- max(cached_state$sigma_upper, sigma_lower * 1.1)

  # Fully vectorized integrand
  integrand <- function(x, power = 1) {
    mu <- x[1, ]
    sigma <- x[2, ]

    # Vectorized metric and log posterior computation
    metric_val <- compute_metric_value(
      mu, sigma, LSL, USL, target, metric,
      sigma_level = sigma_level
    )
    log_p <- log_post_vec(mu, sigma)

    result <- exp(log_p - h_max) * (metric_val^power)
    result[!is.finite(result)] <- 0
    matrix(result, nrow = 1)
  }

  # Both the integral and Z use the same h_max normalization, so they cancel out
  E1 <- cubature::pcubature(
    function(x) integrand(x, power = 1),
    lowerLimit = c(mu_lower, sigma_lower),
    upperLimit = c(mu_upper, sigma_upper),
    tol = 1e-4, vectorInterface = TRUE
  )$integral / Z

  E2 <- cubature::pcubature(
    function(x) integrand(x, power = 2),
    lowerLimit = c(mu_lower, sigma_lower),
    upperLimit = c(mu_upper, sigma_upper),
    tol = 1e-4, vectorInterface = TRUE
  )$integral / Z

  list(mean = E1, sd = sqrt(max(0, E2 - E1^2)))
}


# ==============================================================================
# Solver Factory: Returns function P(Index > c)
# ==============================================================================


