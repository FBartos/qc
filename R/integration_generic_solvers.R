# ==============================================================================
# Generic Prior Solvers and Cached State
# ==============================================================================

.generic_axis_bounds <- function(log_fn, center,
                                 step = 1,
                                 tail_log_drop = 30,
                                 max_steps = 80L) {
  if (!is.finite(center) || !is.finite(step) || step <= 0) {
    return(list(
      lower = NA_real_,
      upper = NA_real_,
      lower_tail_captured = FALSE,
      upper_tail_captured = FALSE
    ))
  }

  h_center <- log_fn(center)
  if (!is.finite(h_center)) {
    return(list(
      lower = NA_real_,
      upper = NA_real_,
      lower_tail_captured = FALSE,
      upper_tail_captured = FALSE
    ))
  }

  lower <- center
  upper <- center
  lower_tail_captured <- FALSE
  upper_tail_captured <- FALSE

  for (i in seq_len(max_steps)) {
    candidate <- center - i * step
    h_candidate <- log_fn(candidate)
    lower <- candidate
    if (!is.finite(h_candidate) || h_candidate <= h_center - tail_log_drop) {
      lower_tail_captured <- TRUE
      break
    }
  }

  for (i in seq_len(max_steps)) {
    candidate <- center + i * step
    h_candidate <- log_fn(candidate)
    upper <- candidate
    if (!is.finite(h_candidate) || h_candidate <= h_center - tail_log_drop) {
      upper_tail_captured <- TRUE
      break
    }
  }

  list(
    lower = lower,
    upper = upper,
    lower_tail_captured = lower_tail_captured,
    upper_tail_captured = upper_tail_captured
  )
}

.generic_positive_axis_bounds <- function(log_fn, center,
                                          log_step = 0.1,
                                          tail_log_drop = 30,
                                          max_steps = 80L) {
  if (!is.finite(center) || center <= 0 || !is.finite(log_step) || log_step <= 0) {
    return(list(
      lower = NA_real_,
      upper = NA_real_,
      lower_tail_captured = FALSE,
      upper_tail_captured = FALSE
    ))
  }

  center_t <- log(center)
  bounds_t <- .generic_axis_bounds(
    function(t) log_fn(exp(t)),
    center = center_t,
    step = log_step,
    tail_log_drop = tail_log_drop,
    max_steps = max_steps
  )

  list(
    lower = exp(bounds_t$lower),
    upper = exp(bounds_t$upper),
    lower_tail_captured = bounds_t$lower_tail_captured,
    upper_tail_captured = bounds_t$upper_tail_captured
  )
}

.generic_prior_improper_error <- function(context = "The generic posterior for integration") {
  stop(
    sprintf(
      "%s is improper for the supplied data and prior (normalization depends on an arbitrary finite box).",
      context
    ),
    call. = FALSE
  )
}

#' @export
make_density_solver.PriorGeneric <- function(data, LSL, USL, prior,
                                              metric = "Cpk", target = NULL,
                                              cached_state = NULL,
                                              sigma_level = 3, ...) {
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

  # For generic priors, we still use 1D integration along contours,
  # but evaluate the joint posterior numerically using vectorized operations

  cached_state <- precompute_generic_state(data, prior, cached_state = cached_state)

  n <- cached_state$n
  x_bar <- cached_state$x_bar
  h_max <- cached_state$h_max
  Z <- cached_state$Z
  uni_s <- cached_state$uni_s
  map_mu <- cached_state$map_mu
  mu_scale <- cached_state$mu_scale
  mu_lower <- cached_state$mu_lower
  mu_upper <- cached_state$mu_upper
  sigma_lower <- cached_state$sigma_lower
  sigma_upper <- cached_state$sigma_upper
  log_post_vec <- cached_state$log_post_vec
  sigma_floor <- max(sigma_lower, 1e-10)

  contour_density <- function(mu, sigma, log_jacobian) {
    result <- numeric(length(mu))
    if (length(result) == 0L) return(result)

    if (length(log_jacobian) == 1L) {
      log_jacobian <- rep(log_jacobian, length(mu))
    }

    valid <- is.finite(mu) & is.finite(sigma) &
      sigma >= sigma_floor & sigma <= sigma_upper &
      mu >= mu_lower & mu <= mu_upper
    if (!any(valid)) return(result)

    log_p <- log_post_vec(mu[valid], sigma[valid])
    contrib <- exp(log_p + log_jacobian[valid] - h_max) / Z
    contrib[!is.finite(contrib)] <- 0

    result[valid] <- contrib
    result
  }

  survival_density <- function(metric_name) {
    .density_from_survival_fn(
      .integration_make_solver(
        request = .integration_request_update(
          request,
          metric = metric_name,
          cached_state = cached_state
        )
      ),
      support_lower = if (.metric_can_be_negative(metric_name)) -Inf else 0
    )
  }

  M <- (LSL + USL) / 2
  target <- .integration_resolve_target(metric, target, LSL, USL)
  tol <- USL - LSL
  cpm_dist <- .cpm_spec_distance(LSL, USL, target)

  if (metric %in% c("Cpk", "Cpu", "Cpl")) {
    # These metrics can have negative support, so differentiating the survival
    # function is the reliable way to cover the full line.
    return(survival_density(metric))
  }

  # Cp: single point
  if (metric == "Cp") {
    return(function(c) {
      if (c <= 0) return(0)

      sigma_c <- tol / ((2 * sigma_level) * c)
      if (sigma_c <= 0) return(0)
      if (sigma_c < sigma_floor || sigma_c > sigma_upper) return(0)
      if (!is.finite(mu_lower) || !is.finite(mu_upper) || mu_upper <= mu_lower) return(0)

      # At sigma_c, mu can be anything inside the cached normalization box.
      # p(Cp = c) = ∫ p(mu, sigma_c) |dsigma/dc| dmu
      # |dsigma/dc| = tol / (2 * sigma_level * c^2)
      log_jacobian <- log(tol / (2 * sigma_level)) - 2 * log(c)

      # Integrate over mu
      mu_integrand <- function(mu) {
        sigma_vec <- rep(sigma_c, length(mu))
        log_p <- log_post_vec(mu, sigma_vec)
        exp(log_p - h_max) / Z
      }

      tryCatch({
        mu_integral <- stats::integrate(mu_integrand, mu_lower, mu_upper, rel.tol = 1e-5)$value
        mu_integral * exp(log_jacobian)
      }, error = function(e) {
        0
      })
    })
  }

  if (metric == "Cpm") {
    return(function(c) {
      if (c <= 0) return(0)

      # Contour radius uses the nearest target-to-spec distance.
      R <- cpm_dist / (sigma_level * c)
      if (R <= 0) return(0)

      # Jacobian determinant for polar coordinates wrt c
      # |J| = R^2 / c
      log_jacobian <- 2 * log(R) - log(c)

      integrand <- function(theta) {
        sigma_v <- R * sin(theta)
        mu_v <- target + R * cos(theta)
        contour_density(mu_v, sigma_v, log_jacobian)
      }

      tryCatch({
        stats::integrate(integrand, 0, pi, rel.tol = 1e-5, subdivisions = 200)$value
      }, error = function(e) {
        0
      })
    })
  }

  if (metric == "Cpc") {
    if (is.null(target)) target <- M
    return(survival_density("Cpc"))
  }

  stop("Unsupported metric for density solver: ", metric)
}


#' @export
make_solver.PriorGeneric <- function(data, LSL, USL, prior, metric = "Cpk",
                                      target = NULL, cached_state = NULL,
                                      sigma_level = 3, ...) {
  cached_state <- precompute_generic_state(data, prior, cached_state = cached_state)
  Z <- cached_state$Z
  int_2d <- cached_state$int_2d

  # Return closure
  function(c) {
    if (!.metric_can_be_negative(metric) && c <= 0) return(1.0)
    constr <- get_metric_constraints(metric, c, LSL, USL, target,
                                     sigma_level = sigma_level)
    num <- int_2d(constr$s_max_fn(), constr$mu_b_fn_vec)
    return(num / Z)
  }
}


# ==============================================================================
# Pre-computation for Generic Priors
# ==============================================================================

#' Precompute expensive posterior state for Generic priors
#' @param data Numeric vector of observations
#' @param prior PriorGeneric object
#' @return cached_state object to pass to make_solver
#' @keywords internal
precompute_generic_state <- function(data, prior, cached_state = NULL) {
  if (!is.null(cached_state) && inherits(cached_state, "qc_generic_cached_state")) {
    .validate_qc_generic_cached_state(cached_state)
    return(cached_state)
  }

  suff_state <- .as_qc_suff_stats_state(data = data, cached_state = cached_state)
  n <- suff_state$n
  x_bar <- suff_state$x_bar
  sse <- suff_state$sse

  if (n > 0) {
    init_sd <- max(sqrt(sse / max(n - 1, 1)), 1e-6)
    init_mu <- x_bar
  } else {
    bt <- prior$bayestools_priors
    if (!is.null(bt)) {
      mu_init  <- .extract_prior_init(bt$mu)
      sig_init <- .extract_prior_init(bt$sigma)
      init_mu <- mu_init$value
      init_sd <- sig_init$value
    } else {
      init_mu <- 0
      init_sd <- 1
    }
  }

  # Scalar log posterior (used by optim)
  log_post <- function(mu, sigma) {
    if (sigma <= 0) return(-Inf)
    log_lik <- if (n > 0) {
      -n * log(sigma) - (sse + n * (mu - x_bar)^2) / (2 * sigma^2)
    } else {
      0
    }
    log_lik + prior$log_dens(mu, sigma)
  }

  init <- .find_feasible_generic_init(init_mu, init_sd, prior, log_post)
  init_mu <- init$mu
  init_sd <- init$sigma

  # Vectorized log posterior (used by cubature)
  log_post_vec <- function(mu, sigma) {
    log_lik <- rep(0, length(mu))
    if (n > 0) {
      log_lik[] <- -Inf
      valid <- sigma > 0
      if (any(valid)) {
        log_lik[valid] <- -n * log(sigma[valid]) -
                          (sse + n * (mu[valid] - x_bar)^2) / (2 * sigma[valid]^2)
      }
    }
    log_prior <- prior$log_dens(mu, sigma)
    log_lik + log_prior
  }

  # Find MAP using constrained optimization (sigma > 0)
  # Use L-BFGS-B with reasonable bounds
  opt <- tryCatch({
    stats::optim(c(init_mu, init_sd), function(p) -log_post(p[1], p[2]),
                 method = "L-BFGS-B",
                 lower = c(-1e6, 1e-6),
                 upper = c(1e6, 1e6))
  }, error = function(e) {
    # Fallback to Nelder-Mead if L-BFGS-B fails
    stats::optim(c(init_mu, max(init_sd, 0.1)), function(p) {
      if (p[2] <= 0) return(1e10)
      -log_post(p[1], p[2])
    })
  })
  map_mu <- opt$par[1]
  map_sig <- max(opt$par[2], 1e-3)  # Ensure positive
  h_max <- -opt$value

  # Compute integration bounds based on posterior concentration
  # For Normal likelihood with data:
  #   - Posterior of sigma is concentrated around map_sig
  #   - Effective SD of sigma is roughly map_sig / sqrt(2*n)
  #   - Posterior of mu|sigma is N(x_bar, sigma^2/n)
  # Use bounds that capture ~6 SDs of the posterior to get >99.99% of mass
  mu_bounds_stable <- TRUE
  sigma_bounds_stable <- TRUE
  if (n > 0) {
    # Sigma bounds: use multiplier based on posterior concentration
    # sigma_sd ~ map_sig / sqrt(2n), so 6 sigma_sd ~ 6 * map_sig / sqrt(2n)
    # For n=30, this is about 0.77 * map_sig, so bounds [0.5 * map_sig, 3 * map_sig] should work
    sigma_lower <- max(0.2 * map_sig, 1e-6)
    sigma_upper <- 4 * map_sig

    bt <- prior$bayestools_priors
    if (!is.null(bt) && inherits(bt$sigma, "prior")) {
      sig_bounds <- .extract_prior_bounds(bt$sigma)
      if (is.finite(sig_bounds$lower)) {
        sigma_lower <- max(sigma_lower, sig_bounds$lower + 1e-6)
      }
      if (is.finite(sig_bounds$upper)) {
        sigma_upper <- min(sigma_upper, sig_bounds$upper)
      }
    }
    sigma_upper <- max(sigma_upper, sigma_lower * 1.1)

    # Mu bounds: at sigma_upper, SD of mu|sigma is sigma_upper/sqrt(n)
    # Use 6 SDs for safety
    mu_sd_at_sigma_upper <- sigma_upper / sqrt(n)
    mu_lower <- map_mu - 6 * mu_sd_at_sigma_upper
    mu_upper <- map_mu + 6 * mu_sd_at_sigma_upper

    # For other functions (like make_solver), uni_s is used as a reference scale
    uni_s <- sigma_upper
  } else {
    # No data (prior-only): derive bounds from prior quantiles
    bt <- prior$bayestools_priors
    if (!is.null(bt) && inherits(bt$sigma, "prior")) {
      sig_q_lo <- BayesTools::quant(bt$sigma, 0.0001)
      sig_q_hi <- BayesTools::quant(bt$sigma, 0.9999)
      sigma_lower <- max(1e-6, sig_q_lo)
      sigma_upper <- sig_q_hi
    } else if (is.finite(map_sig) && map_sig > 0) {
      sigma_bounds <- .generic_positive_axis_bounds(
        function(sigma) log_post(map_mu, sigma),
        center = map_sig
      )
      sigma_lower <- max(1e-6, sigma_bounds[["lower"]])
      sigma_upper <- max(sigma_bounds[["upper"]], sigma_lower * 1.1)
      sigma_bounds_stable <- isTRUE(sigma_bounds$lower_tail_captured) &&
        isTRUE(sigma_bounds$upper_tail_captured)
    } else {
      sigma_lower <- 0.01
      sigma_upper <- 100
      sigma_bounds_stable <- FALSE
    }
    uni_s <- sigma_upper
    if (!is.null(bt) && inherits(bt$mu, "prior")) {
      mu_init_info <- .extract_prior_init(bt$mu)
      mu_support <- .extract_prior_bounds(bt$mu)
      mu_scale <- mu_init_info$scale
      mu_lower <- map_mu - 6 * mu_scale
      mu_upper <- map_mu + 6 * mu_scale
      if (is.finite(mu_support$lower)) {
        mu_lower <- max(mu_lower, mu_support$lower)
      }
      if (is.finite(mu_support$upper)) {
        mu_upper <- min(mu_upper, mu_support$upper)
      }
    } else {
      mu_scale <- NULL
      mu_bounds <- .generic_axis_bounds(
        function(mu) log_post(mu, max(map_sig, 1e-6)),
        center = map_mu,
        step = max(1, abs(map_mu) * 0.05)
      )
      mu_lower <- mu_bounds[["lower"]]
      mu_upper <- mu_bounds[["upper"]]
      mu_bounds_stable <- isTRUE(mu_bounds$lower_tail_captured) &&
        isTRUE(mu_bounds$upper_tail_captured)
    }
  }

  integrate_box <- function(mu_lower_box, mu_upper_box,
                            sigma_lower_box, sigma_upper_box,
                            s_lim, m_fn) {
    s_top <- if (is.infinite(s_lim)) sigma_upper_box else min(s_lim, sigma_upper_box)
    sigma_floor <- max(sigma_lower_box, 1e-10)
    if (!is.finite(s_top) || s_top <= sigma_floor || mu_upper_box <= mu_lower_box) {
      return(0)
    }

    integrand <- function(x) {
      mu <- x[1, ]
      sigma <- x[2, ]

      log_p <- log_post_vec(mu, sigma)
      mb <- m_fn(sigma)
      valid <- (mu >= mb$lower) & (mu <= mb$upper)

      result <- rep(0, length(mu))
      result[valid] <- exp(log_p[valid] - h_max)
      result[!is.finite(result)] <- 0
      matrix(result, nrow = 1)
    }

    tryCatch(
      cubature::pcubature(
        integrand,
        lowerLimit = c(mu_lower_box, sigma_floor),
        upperLimit = c(mu_upper_box, s_top),
        tol = 1e-4,
        vectorInterface = TRUE
      )$integral,
      error = function(e) NA_real_
    )
  }

  # cubature 2D integration with fully vectorized interface
  int_2d <- function(s_lim, m_fn) {
    integrate_box(mu_lower, mu_upper, sigma_lower, sigma_upper, s_lim, m_fn)
  }

  unconstrained_mu_box <- function(s) {
    list(lower = rep(-Inf, length(s)), upper = rep(Inf, length(s)))
  }

  Z <- int_2d(Inf, unconstrained_mu_box)
  if ((!is.finite(Z) || Z <= 0) &&
      n == 0 &&
      (is.null(bt) || !inherits(bt$sigma, "prior")) &&
      is.finite(map_sig) &&
      map_sig > 0) {
    sigma_bounds <- .generic_positive_axis_bounds(
      function(sigma) log_post(map_mu, sigma),
      center = map_sig,
      log_step = 0.05,
      tail_log_drop = 40
    )
    sigma_lower <- max(1e-6, sigma_bounds[["lower"]])
    sigma_upper <- max(sigma_bounds[["upper"]], sigma_lower * 1.1)
    sigma_bounds_stable <- isTRUE(sigma_bounds$lower_tail_captured) &&
      isTRUE(sigma_bounds$upper_tail_captured)
    uni_s <- sigma_upper
    mu_bounds <- .generic_axis_bounds(
      function(mu) log_post(mu, max(map_sig, 1e-6)),
      center = map_mu,
      step = max(0.5, abs(map_mu) * 0.02),
      tail_log_drop = 40
    )
    mu_lower <- mu_bounds[["lower"]]
    mu_upper <- mu_bounds[["upper"]]
    mu_bounds_stable <- isTRUE(mu_bounds$lower_tail_captured) &&
      isTRUE(mu_bounds$upper_tail_captured)
    Z <- int_2d(Inf, unconstrained_mu_box)
  }

  if (n == 0) {
    # Optional metadata helps choose finite bounds, but it does not guarantee
    # that the embedded custom log-density is proper. Validate prior-only
    # normalization by checking that it stays stable under box expansion.
    if (!is.finite(Z) || Z <= 0) {
      .generic_prior_improper_error()
    }

    mu_mid <- (mu_lower + mu_upper) / 2
    mu_half_width <- max((mu_upper - mu_lower) / 2, 1)
    expanded_mu_lower <- mu_mid - 2 * mu_half_width
    expanded_mu_upper <- mu_mid + 2 * mu_half_width
    expanded_sigma_lower <- max(1e-10, sigma_lower / 10)
    expanded_sigma_upper <- max(sigma_upper * 10, expanded_sigma_lower * 1.1)

    Z_expanded <- integrate_box(
      expanded_mu_lower,
      expanded_mu_upper,
      expanded_sigma_lower,
      expanded_sigma_upper,
      Inf,
      unconstrained_mu_box
    )

    box_stable <- is.finite(Z_expanded) && Z_expanded > 0 && (Z_expanded / Z) <= 1.1
    if (!(mu_bounds_stable && sigma_bounds_stable) && !box_stable) {
      .generic_prior_improper_error()
    }
  }

  # mu_scale: prior SD for mu (used by density solver for bounds when n=0)
  if (!exists("mu_scale")) mu_scale <- NULL

  state <- .new_generic_backend_state(
    log_post = log_post,
    log_post_vec = log_post_vec,
    h_max = h_max,
    uni_s = uni_s,
    map_mu = map_mu,
    mu_scale = mu_scale,
    int_2d = int_2d,
    Z = Z,
    n = n,
    x_bar = x_bar,
    sse = sse,
    mu_lower = mu_lower,
    mu_upper = mu_upper,
    sigma_lower = sigma_lower,
    sigma_upper = sigma_upper,
    map_sig = map_sig,
    bayestools_priors = prior$bayestools_priors %||% NULL
  )
  .validate_qc_generic_cached_state(state)
  state
}


# ==============================================================================
# Main Analysis Functions
# ==============================================================================

