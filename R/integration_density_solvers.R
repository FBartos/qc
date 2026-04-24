# ==============================================================================
# Density Solvers
# ==============================================================================

#' Create density solver for direct PDF computation via contour integration
#'
#' Instead of computing P(metric > c) and differentiating, this computes the
#' marginal PDF p(c) directly by integrating along the contour lines where
#' metric = c. This is faster and avoids finite-difference noise.
#'
#' @param data Numeric vector of observations
#' @param LSL Lower specification limit
#' @param USL Upper specification limit
#' @param prior Prior object (PriorConjugate)
#' @param metric Capability metric (currently "Cpk", "Cpu", "Cpl", "Cp" supported)
#' @param target Target value for Cpm/Cpc
#' @param ... Additional arguments
#' @return Function pdf(c) returning the marginal posterior density at c
#' @keywords internal
make_density_solver <- function(data, LSL, USL, prior, metric = "Cpk",
                                target = NULL, sigma_level = 3, ...) {
  .validate_capability_request(
    LSL = LSL,
    USL = USL,
    target = target,
    sigma_level = sigma_level,
    metric = metric,
    target_required = FALSE,
    sigma_name = "sigma_level"
  )

  UseMethod("make_density_solver", prior)
}

.density_from_survival_fn <- function(S,
                                      rel_step = 1e-4,
                                      abs_step = 1e-5,
                                      support_lower = 0) {
  force(S)
  force(support_lower)

  function(c) {
    c <- as.numeric(c)
    result <- numeric(length(c))
    valid <- is.finite(c)
    if (!any(valid)) {
      return(result)
    }

    c_valid <- c[valid]
    h <- pmax(abs(c_valid) * rel_step, abs_step)
    lower <- c_valid - h
    upper <- c_valid + h
    density <- numeric(length(c_valid))

    outside_support <- is.finite(support_lower) & c_valid < support_lower
    use_forward <- !outside_support &
      is.finite(support_lower) &
      lower < support_lower
    use_central <- !outside_support & !use_forward

    if (any(use_forward)) {
      base <- pmax(c_valid[use_forward], support_lower)
      density[use_forward] <- (
        vapply(base, S, numeric(1L)) -
          vapply(upper[use_forward], S, numeric(1L))
      ) / (upper[use_forward] - base)
    }

    if (any(use_central)) {
      density[use_central] <- (
        vapply(lower[use_central], S, numeric(1L)) -
          vapply(upper[use_central], S, numeric(1L))
      ) / (upper[use_central] - lower[use_central])
    }

    density[!is.finite(density) | density < 0] <- 0
    result[valid] <- density
    result
  }
}

#' @export
#' @export
make_density_solver.PriorConjugate <- function(data, LSL, USL, prior,
                                                 metric = "Cpk", target = NULL,
                                                 cached_state = NULL,
                                                sigma_level = 3, ...) {
  posterior_info <- .compute_validated_conjugate_posterior(
    prior, data, cached_state,
    context = "The conjugate posterior for density computation"
  )
  ss <- posterior_info$ss
  post <- posterior_info$post
  n <- ss$n; k_n <- post$k_n; mu_n <- post$mu_n
  alpha_n <- post$alpha_n; beta_n <- post$beta_n

  M <- (LSL + USL) / 2
  if (is.null(target)) target <- M
  tol <- USL - LSL
  cpm_dist <- .cpm_spec_distance(LSL, USL, target)
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

  if (posterior_info$is_degenerate) {
    dist <- .degenerate_conjugate_metric_distribution(
      mu_n, k_n, LSL, USL, target, metric,
      sigma_level = sigma_level
    )
    if (!is.null(dist)) {
      return(.degenerate_metric_density(dist))
    }
  }

  log_Z_sigma <- lgamma(alpha_n) - log(2) - alpha_n * log(beta_n)

  # Cp: no integration needed (contour is a single sigma point)
  if (metric == "Cp") {
    return(function(c) {
      if (c <= 0) return(0)
      sigma_c <- tol / ((2 * sigma_level) * c)
      if (sigma_c <= 0) return(0)
      log_p_sigma <- -(2 * alpha_n + 1) * log(sigma_c) - beta_n / sigma_c^2
      log_jacobian <- log(tol / (2 * sigma_level)) - 2 * log(c)
      exp(log_p_sigma + log_jacobian - log_Z_sigma)
    })
  }

  if (metric %in% c("Cpu", "Cpl", "Cpk")) {
    return(.density_from_survival_fn(
      .integration_make_solver(request = request),
      support_lower = -Inf
    ))
  }

  if (metric == "Cpm") {
    n_theta <- 1024L
    theta_grid <- seq(1e-8, pi - 1e-8, length.out = n_theta)
    d_theta <- theta_grid[2] - theta_grid[1]

    return(function(c) {
      if (c <= 0) return(0)
      R <- cpm_dist / (sigma_level * c)
      if (R <= 0) return(0)
      log_jacobian <- 2 * log(R) - log(c)

      sigma_v <- R * sin(theta_grid)
      valid <- sigma_v > 1e-10
      if (!any(valid)) return(0)
      sv <- sigma_v[valid]
      mu_v <- target + R * cos(theta_grid[valid])
      sd_mu <- sv / sqrt(k_n)

      log_p_sigma <- -(2 * alpha_n + 1) * log(sv) - beta_n / sv^2
      log_p_mu <- stats::dnorm(mu_v, mu_n, sd_mu, log = TRUE)
      vals <- exp(log_p_sigma + log_p_mu + log_jacobian - log_Z_sigma)
      vals[!is.finite(vals)] <- 0
      sum(vals) * d_theta
    })
  }

  if (metric == "Cpc") {
    if (is.null(target)) target <- M

    return(function(c) {
      if (c <= 0) return(0)

      integrand <- function(z) {
        contour <- .cpc_contour_from_z(
          z, c, tol, target,
          sigma_level = sigma_level
        )
        sd_mu <- contour$sigma / sqrt(k_n)
        log_p_sigma <- -(2 * alpha_n + 1) * log(contour$sigma) -
          beta_n / contour$sigma^2
        log_p_mu_L <- stats::dnorm(contour$mu_lower, mu_n, sd_mu, log = TRUE)
        log_p_mu_U <- stats::dnorm(contour$mu_upper, mu_n, sd_mu, log = TRUE)

        vals <- exp(log_p_sigma + log_p_mu_L + contour$log_jacobian - log_Z_sigma) +
          exp(log_p_sigma + log_p_mu_U + contour$log_jacobian - log_Z_sigma)
        vals[!is.finite(vals)] <- 0
        vals
      }

      stats::integrate(
        integrand, 0, 20,
        rel.tol = 1e-5,
        subdivisions = 400
      )$value
    })
  }

  stop("Unsupported metric for density solver: ", metric)
}

#' @export
make_density_solver.PriorSemiConjugateMu <- function(data, LSL, USL, prior,
                                                      metric = "Cpk", target = NULL,
                                                      cached_state = NULL,
                                                      sigma_level = 3, ...) {
  # Case 2: Contour integration with semi-analytical mu integration
  # Uses the same contour approach as PriorConjugate but with
  # log p(sigma | data) = log_lik(sigma) + log_sigma_prior(sigma)

  state <- .prepare_semi_mu_state(data, prior, cached_state)
  k_n <- state$k_n
  mu_n <- state$mu_n
  sse_n <- state$sse_n

  M <- (LSL + USL) / 2
  if (is.null(target)) target <- M
  tol <- USL - LSL
  cpm_dist <- .cpm_spec_distance(LSL, USL, target)
  log_Z <- state$log_Z
  log_p_sigma <- state$log_p_sigma
  sigma_support <- state$sigma_support
  sigma_domain <- state$sigma_domain
  sigma_lower <- max(sigma_support$lower, sigma_domain$lower)
  sigma_upper <- min(sigma_support$upper, sigma_domain$upper)

  if (metric %in% c("Cpk", "Cpu", "Cpl")) {
    return(.density_from_survival_fn(
      make_solver(
        data, LSL, USL, prior,
        metric = metric,
        target = target,
        cached_state = cached_state,
        sigma_level = sigma_level
      ),
      support_lower = -Inf
    ))
  }

  # Cp: sigma is fixed on the contour, no integration needed
  if (metric == "Cp") {
    return(function(c) {
      if (c <= 0) return(0)

      sigma_c <- tol / ((2 * sigma_level) * c)
      if (sigma_c <= 0 ||
          sigma_c < sigma_lower ||
          sigma_c > sigma_upper) {
        return(0)
      }

      lps <- log_p_sigma(sigma_c)
      log_jacobian <- log(tol / (2 * sigma_level)) - 2 * log(c)
      exp(lps + log_jacobian - log_Z)
    })
  }

  if (metric == "Cpm") {
    n_theta <- 1024L
    theta_grid <- seq(1e-8, pi - 1e-8, length.out = n_theta)
    d_theta <- theta_grid[2] - theta_grid[1]

    return(function(c) {
      if (c <= 0) return(0)

      R <- cpm_dist / (sigma_level * c)
      if (R <= 0) return(0)

      log_jacobian <- 2 * log(R) - log(c)

      sigma_v <- R * sin(theta_grid)
      mu_v <- target + R * cos(theta_grid)
      valid <- sigma_v >= sigma_lower & sigma_v <= sigma_upper & sigma_v > 1e-10
      if (!any(valid)) return(0)
      sigma_v <- sigma_v[valid]
      mu_v <- mu_v[valid]
      sd_mu <- sigma_v / sqrt(k_n)

      lps <- log_p_sigma(sigma_v)
      log_p_mu <- stats::dnorm(mu_v, mu_n, sd_mu, log = TRUE)

      vals <- exp(lps + log_p_mu + log_jacobian - log_Z)
      vals[!is.finite(vals)] <- 0
      sum(vals) * d_theta
    })
  }

  if (metric == "Cpc") {
    return(function(c) {
      if (c <= 0) return(0)

      integrand <- function(z) {
        contour <- .cpc_contour_from_z(
          z, c, tol, target,
          sigma_level = sigma_level
        )
        valid <- contour$sigma >= sigma_lower &
          contour$sigma <= sigma_upper &
          contour$sigma > 0
        vals <- numeric(length(contour$sigma))
        if (!any(valid)) return(vals)

        sigma_v <- contour$sigma[valid]
        sd_mu <- sigma_v / sqrt(k_n)
        lps <- log_p_sigma(sigma_v)
        log_p_mu_L <- stats::dnorm(contour$mu_lower[valid], mu_n, sd_mu, log = TRUE)
        log_p_mu_U <- stats::dnorm(contour$mu_upper[valid], mu_n, sd_mu, log = TRUE)

        vals[valid] <- exp(lps + log_p_mu_L + contour$log_jacobian[valid] - log_Z) +
          exp(lps + log_p_mu_U + contour$log_jacobian[valid] - log_Z)
        vals[!is.finite(vals)] <- 0
        vals
      }

      stats::integrate(
        integrand, 0, 20,
        rel.tol = 1e-5,
        subdivisions = 400
      )$value
    })
  }

  stop("Unsupported metric for density solver: ", metric)
}

#' Helper to compute log normalization constant for semi-conjugate mu prior
#' @keywords internal
.semi_mu_sigma_support <- function(prior) {
  sigma_prior <- prior$bayestools_priors$sigma %||% NULL
  sigma_support <- if (!is.null(sigma_prior) && inherits(sigma_prior, "prior")) {
    .extract_prior_bounds(sigma_prior)
  } else {
    list(lower = 0, upper = Inf)
  }
  sigma_support$lower <- max(sigma_support$lower, 0)
  list(prior = sigma_prior, support = sigma_support)
}

.semi_mu_sigma_log_domain <- function(n_eff, sse_n, log_dens_sigma,
                                      sigma_support = list(lower = 0, upper = Inf),
                                      sigma_prior = NULL,
                                      tail_log_drop = 36,
                                      context = "The semi-conjugate mu posterior") {
  sigma_floor <- if (is.finite(sigma_support$lower) && sigma_support$lower > 0) {
    sigma_support$lower
  } else {
    1e-12
  }
  upper_bound <- if (is.finite(sigma_support$upper) && sigma_support$upper > 0) {
    log(sigma_support$upper)
  } else {
    Inf
  }

  init <- .extract_prior_init(sigma_prior)
  qs <- .safe_prior_quantiles(sigma_prior, c(1e-6, 0.5, 1 - 1e-6))
  sigma_ref <- sqrt(max(sse_n, 0) / max(n_eff + 1, 1))
  if (!is.finite(sigma_ref) || sigma_ref <= 0) sigma_ref <- NA_real_

  sigma_candidates <- c(sigma_ref, init$value, qs)
  if (is.finite(sigma_support$lower) && sigma_support$lower > 0) {
    sigma_candidates <- c(sigma_candidates, sigma_support$lower)
  }
  if (is.finite(sigma_support$upper) && sigma_support$upper > 0) {
    sigma_candidates <- c(sigma_candidates, sigma_support$upper)
  }
  sigma_candidates <- sort(unique(sigma_candidates[is.finite(sigma_candidates) & sigma_candidates > 0]))
  if (length(sigma_candidates) == 0L) sigma_candidates <- 1

  lower_bound <- log(sigma_floor)
  t_lo_init <- max(lower_bound, min(log(sigma_candidates)) - 4)
  t_hi_init <- max(t_lo_init + 1e-6, max(log(sigma_candidates)) + 4)
  if (is.finite(upper_bound)) {
    t_hi_init <- min(t_hi_init, upper_bound)
  }
  if (t_hi_init <= t_lo_init) {
    t_hi_init <- if (is.finite(upper_bound) && upper_bound > t_lo_init) {
      upper_bound
    } else {
      t_lo_init + 8
    }
  }

  log_base_t <- function(t) {
    sigma <- exp(t)
    vals <- -n_eff * t - sse_n / (2 * sigma^2) + log_dens_sigma(sigma) + t
    vals[!is.finite(vals)] <- -Inf
    vals
  }

  opt <- stats::optimize(function(t) {
    val <- log_base_t(t)
    if (!is.finite(val)) Inf else -val
  }, c(t_lo_init, t_hi_init))
  t_mode <- opt$minimum
  h_max <- log_base_t(t_mode)
  if (!is.finite(h_max)) {
    stop(
      sprintf(
        "%s is improper for the supplied data and prior (no finite sigma mode within support).",
        context
      ),
      call. = FALSE
    )
  }

  t_lower <- if (is.finite(sigma_support$lower) && sigma_support$lower > 0) {
    lower_bound
  } else {
    max(lower_bound, min(t_lo_init, t_mode - 8))
  }
  lower_tail_captured <- TRUE
  if (!is.finite(sigma_support$lower) || sigma_support$lower == 0) {
    lower_tail_captured <- FALSE
    repeat {
      val <- log_base_t(t_lower)
      if (!is.finite(val) || val <= h_max - tail_log_drop) {
        lower_tail_captured <- TRUE
        break
      }
      if (t_lower <= lower_bound + 1e-8) {
        break
      }
      t_lower <- max(lower_bound, t_lower - 2)
    }
  }

  t_upper <- if (is.finite(upper_bound)) {
    upper_bound
  } else {
    max(t_hi_init, t_mode + 8)
  }
  upper_tail_captured <- TRUE
  if (!is.finite(upper_bound)) {
    upper_tail_captured <- FALSE
    for (i in seq_len(60)) {
      val <- log_base_t(t_upper)
      if (!is.finite(val) || val <= h_max - tail_log_drop) {
        upper_tail_captured <- TRUE
        break
      }
      t_upper <- t_upper + 2
    }
  }

  if (t_upper <= t_lower) {
    t_upper <- t_lower + 1e-6
  }

  if (!lower_tail_captured || !upper_tail_captured) {
    stop(
      sprintf(
        "%s is improper for the supplied data and prior (sigma tail mass does not decay fast enough to normalize).",
        context
      ),
      call. = FALSE
    )
  }

  list(lower = t_lower, upper = t_upper, mode = t_mode, h_max = h_max)
}

.prepare_semi_mu_state <- function(data, prior, cached_state = NULL) {
  ss <- .extract_suff_stats(data, cached_state)
  n <- ss$n
  x_bar <- ss$x_bar
  sse <- ss$SS
  smp <- .semi_mu_posterior(prior, n, x_bar, sse)
  k_n <- smp$k_n
  mu_n <- smp$mu_n
  sse_n <- smp$sse_n

  # The sigma exponent depends on whether the mu prior contributes sigma^{-1}.
  n_eff <- if (prior$k0 == 0) n - 1 else n
  sigma_info <- .semi_mu_sigma_support(prior)
  log_domain <- .semi_mu_sigma_log_domain(
    n_eff, sse_n, prior$log_dens_sigma,
    sigma_support = sigma_info$support,
    sigma_prior = sigma_info$prior,
    context = "The semi-conjugate mu posterior"
  )

  log_p_sigma <- function(sigma_v) {
    vals <- -n_eff * log(sigma_v) - sse_n / (2 * sigma_v^2) + prior$log_dens_sigma(sigma_v)
    vals[!is.finite(vals)] <- -Inf
    vals
  }

  log_Z <- .compute_semi_mu_log_Z(
    n_eff, sse_n, prior$log_dens_sigma,
    sigma_support = sigma_info$support,
    sigma_prior = sigma_info$prior,
    log_domain = log_domain
  )

  list(
    n = n,
    x_bar = x_bar,
    sse = sse,
    k_n = k_n,
    mu_n = mu_n,
    sse_n = sse_n,
    n_eff = n_eff,
    sigma_prior = sigma_info$prior,
    sigma_support = sigma_info$support,
    sigma_domain = list(lower = exp(log_domain$lower), upper = exp(log_domain$upper)),
    t_lower = log_domain$lower,
    t_upper = log_domain$upper,
    base_log_h_max = log_domain$h_max,
    log_p_sigma = log_p_sigma,
    log_Z = log_Z
  )
}

.semi_mu_t_grid <- function(t_lower, t_upper, n_points = 2048L) {
  if (!is.finite(t_lower) || !is.finite(t_upper) || t_upper <= t_lower) {
    return(NULL)
  }

  t_grid <- seq(t_lower, t_upper, length.out = n_points)
  dt <- (t_upper - t_lower) / max(n_points - 1, 1)
  weights <- rep(dt, n_points)
  weights[c(1, n_points)] <- dt / 2

  list(t = t_grid, sigma = exp(t_grid), weights = weights)
}

.integrate_log_grid <- function(log_vals, weights, log_norm = 0) {
  valid <- is.finite(log_vals) & is.finite(weights) & weights > 0
  if (!any(valid)) return(0)

  h_max <- max(log_vals[valid])
  total <- sum(exp(log_vals[valid] - h_max) * weights[valid])
  if (!is.finite(total) || total <= 0) return(0)

  exp(h_max - log_norm) * total
}

.compute_semi_mu_log_Z <- function(n_eff, sse_n, log_dens_sigma,
                                   sigma_support = list(lower = 0, upper = Inf),
                                   sigma_prior = NULL,
                                   log_domain = NULL) {
  if (is.null(log_domain)) {
    log_domain <- .semi_mu_sigma_log_domain(
      n_eff, sse_n, log_dens_sigma,
      sigma_support = sigma_support,
      sigma_prior = sigma_prior
    )
  }

  integrand <- function(t) {
    sigma <- exp(t)
    vals <- -n_eff * log(sigma) - sse_n / (2 * sigma^2) + log_dens_sigma(sigma) + t
    vals <- exp(vals - log_domain$h_max)
    vals[!is.finite(vals)] <- 0
    vals
  }

  Z_scaled <- stats::integrate(
    integrand,
    log_domain$lower,
    log_domain$upper,
    rel.tol = 1e-5,
    subdivisions = 400
  )$value

  log(max(Z_scaled, 1e-300)) + log_domain$h_max
}

#' @export
make_density_solver.PriorSemiConjugateSigma <- function(data, LSL, USL, prior,
                                                         metric = "Cpk", target = NULL,
                                                         cached_state = NULL,
                                                         sigma_level = 3, ...) {
  # Case 3: Contour integration with semi-analytical sigma integration
  # For each mu on the contour, p(sigma|mu,data) is Inverse-Gamma

  state <- .prepare_semi_sigma_state(
    data, prior, cached_state,
    context = "The semi-conjugate sigma posterior for density computation"
  )
  n <- state$n
  x_bar <- state$x_bar
  sse <- state$sse
  alpha0 <- state$alpha0
  beta0 <- state$beta0
  alpha_n_eff <- state$alpha_n_eff

  M <- (LSL + USL) / 2
  if (is.null(target)) target <- M
  tol <- USL - LSL

  # Compute normalization constant Z
  log_Z <- .compute_semi_sigma_log_Z_survival(
    n, x_bar, sse, alpha0, beta0,
    prior$log_dens_mu, state$jeffreys_adj,
    mu_domain = state$mu_domain,
    context = "The semi-conjugate sigma posterior for normalization"
  )

  # Constant for the unnormalized sigma kernel
  log_C <- log(2) - lgamma(alpha_n_eff)

  survival_density <- function(metric_name) {
    .density_from_survival_fn(
      make_solver(
        data, LSL, USL, prior,
        metric = metric_name,
        target = target,
        cached_state = cached_state,
        sigma_level = sigma_level
      ),
      support_lower = if (.metric_can_be_negative(metric_name)) -Inf else 0
    )
  }

  if (metric %in% c("Cpu", "Cpl", "Cpk")) {
    return(survival_density(metric))
  }

  if (metric == "Cp") {
    return(function(c) {
      if (c <= 0) return(0)
      sigma_c <- tol / ((2 * sigma_level) * c)
      if (sigma_c <= 0) return(0)

      log_jacobian <- log(tol / (2 * sigma_level)) - 2 * log(c)

      log_integrand <- function(mu) {
        sse_mu <- sse + n * (mu - x_bar)^2
        beta_n <- beta0 + sse_mu / 2
        log_sigma_kernel <- log_C - (2 * alpha_n_eff + 1) * log(sigma_c) - beta_n / sigma_c^2
        vals <- prior$log_dens_mu(mu) + log_sigma_kernel + log_jacobian - log_Z
        vals[!is.finite(vals)] <- -Inf
        vals
      }
      eff_bounds <- .semi_sigma_effective_bounds(
        log_integrand,
        state$mu_domain,
        lower = -Inf,
        upper = Inf
      )

      .semi_sigma_integrate_log_kernel(
        log_integrand,
        state$mu_domain,
        lower = eff_bounds[["lower"]],
        upper = eff_bounds[["upper"]],
        rel.tol = 1e-5,
        subdivisions = 400
      )$value
    })
  }

  if (metric == "Cpm") {
    return(survival_density("Cpm"))
  }

  if (metric == "Cpc") {
    return(survival_density("Cpc"))
  }

  stop("Unsupported metric for density solver: ", metric)
}


