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

  survival_density <- function(metric_name) {
    .density_from_survival_fn(
      .integration_make_solver(
        request = .integration_request_update(
          request,
          metric = metric_name,
          cached_state = cached_state
        )
      ),
      support_lower = 0
    )
  }

  M <- (LSL + USL) / 2
  target <- .integration_resolve_target(metric, target, LSL, USL)
  tol <- USL - LSL
  cpm_dist <- .cpm_spec_distance(LSL, USL, target)

  if (metric %in% c("Cpk", "Cpu", "Cpl")) {
    return(.linear_metric_density(
      function(mu, sigma) log_post_vec(mu, sigma) - h_max - log(Z),
      LSL, USL, metric, sigma_level,
      sigma_ref = cached_state$map_sig,
      sigma_lower = sigma_floor, sigma_upper = sigma_upper,
      mu_lower = mu_lower, mu_upper = mu_upper
    ))
  }

  if (metric == "Cp") {
    return(function(c) {
      c <- as.numeric(c)
      result <- numeric(length(c))
      valid_c <- is.finite(c) & c > 0
      if (!any(valid_c)) return(result)
      if (!is.finite(mu_lower) || !is.finite(mu_upper) || mu_upper <= mu_lower) {
        return(result)
      }

      cv <- c[valid_c]
      sigma_c <- tol / ((2 * sigma_level) * cv)
      valid_sigma <- sigma_c > 0 &
        sigma_c >= sigma_floor &
        sigma_c <= sigma_upper
      vals_out <- numeric(length(cv))
      if (!any(valid_sigma)) {
        result[valid_c] <- vals_out
        return(result)
      }

      sv <- sigma_c[valid_sigma]
      cv_valid <- cv[valid_sigma]
      log_jacobian <- log(tol / (2 * sigma_level)) - 2 * log(cv_valid)

      if (is.function(cached_state$int_mu)) {
        vals_out[valid_sigma] <- cached_state$int_mu(sv) * exp(log_jacobian) / Z
        result[valid_c] <- vals_out
        return(result)
      }

      integrand <- function(mu) {
        mu <- if (is.matrix(mu)) as.numeric(mu[1, ]) else as.numeric(mu)
        mu_mat <- matrix(rep(mu, each = length(sv)), nrow = length(sv))
        sigma_mat <- matrix(rep(sv, times = length(mu)), nrow = length(sv))
        log_p <- matrix(
          log_post_vec(as.vector(mu_mat), as.vector(sigma_mat)),
          nrow = length(sv)
        )
        vals <- exp(log_p - h_max + log_jacobian) / Z
        vals[!is.finite(vals)] <- 0
        vals
      }

      vals <- tryCatch(
        cubature::pcubature(
          integrand,
          lowerLimit = mu_lower,
          upperLimit = mu_upper,
          tol = 1e-5,
          fDim = length(sv),
          vectorInterface = TRUE
        )$integral,
        error = function(e) rep(0, length(sv))
      )
      vals_out[valid_sigma] <- vals
      result[valid_c] <- vals_out
      result
    })
  }

  if (metric == "Cpm") {
    return(function(c) {
      c <- as.numeric(c)
      result <- numeric(length(c))
      valid_c <- is.finite(c) & c > 0
      if (!any(valid_c)) return(result)

      cv <- c[valid_c]
      R <- cpm_dist / (sigma_level * cv)
      valid_R <- is.finite(R) & R > 0
      vals_out <- numeric(length(cv))
      if (!any(valid_R)) {
        result[valid_c] <- vals_out
        return(result)
      }

      Rv <- R[valid_R]
      cv_valid <- cv[valid_R]
      log_jacobian <- 2 * log(Rv) - log(cv_valid)

      integrand <- function(theta) {
        theta <- if (is.matrix(theta)) as.numeric(theta[1, ]) else as.numeric(theta)
        sigma_mat <- outer(Rv, sin(theta), `*`)
        mu_mat <- target + outer(Rv, cos(theta), `*`)
        log_jacobian_mat <- matrix(
          rep(log_jacobian, times = length(theta)),
          nrow = length(Rv)
        )
        valid <- is.finite(mu_mat) & is.finite(sigma_mat) &
          sigma_mat >= sigma_floor & sigma_mat <= sigma_upper &
          mu_mat >= mu_lower & mu_mat <= mu_upper
        vals <- matrix(0, nrow = length(Rv), ncol = length(theta))
        if (any(valid)) {
          log_p <- log_post_vec(mu_mat[valid], sigma_mat[valid])
          vals[valid] <- exp(log_p + log_jacobian_mat[valid] - h_max) / Z
        }
        vals[!is.finite(vals)] <- 0
        vals
      }

      vals <- tryCatch(
        cubature::pcubature(
          integrand,
          lowerLimit = 0,
          upperLimit = pi,
          tol = 1e-5,
          fDim = length(Rv),
          vectorInterface = TRUE
        )$integral,
        error = function(e) rep(0, length(Rv))
      )
      vals_out[valid_R] <- vals
      result[valid_c] <- vals_out
      result
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
# Local (MAP-anchored) quadrature for Generic Posteriors with Data
# ==============================================================================

.gauss_legendre_rule <- function(k) {
  i <- seq_len(k - 1L)
  off_diag <- i / sqrt(4 * i^2 - 1)
  jacobi <- matrix(0, k, k)
  jacobi[cbind(i, i + 1L)] <- off_diag
  jacobi[cbind(i + 1L, i)] <- off_diag
  eig <- eigen(jacobi, symmetric = TRUE)
  ord <- order(eig$values)
  list(nodes = eig$values[ord], weights = 2 * eig$vectors[1, ord]^2)
}

# Composite Gauss-Legendre rule on [0, 1] with equal-width panels.
.composite_gauss_legendre_unit <- function(panels, k = 10L) {
  rule <- .gauss_legendre_rule(k)
  edges <- seq(0, 1, length.out = panels + 1L)
  half <- diff(edges) / 2
  mid <- edges[-1L] - half
  list(
    nodes = as.vector(outer(rule$nodes, half) + rep(mid, each = k)),
    weights = as.vector(outer(rule$weights, half))
  )
}

# Builds the integration machinery for a generic posterior with informative data
# (n >= 2, sse > 0). Everything is anchored on the posterior mode in
# (mu, log sigma) and scaled by local posterior widths, so a narrow posterior is
# resolved regardless of the measurement units.
#
# The posterior is integrated as an outer Gauss-Legendre rule over
# t = log(sigma) and, for every sigma node, an inner Gauss-Legendre rule over a
# mu window centred on the conditional mode and sized by the conditional
# curvature. Metric constraints enter as limits of the inner window rather than
# as indicators, keeping survival functions smooth in the metric value.
#
# Returns NULL when the local structure cannot be established, in which case
# the caller falls back to the box-based integration.
.generic_local_posterior <- function(log_post_vec, n, x_bar, sse,
                                     init_mu, init_sd,
                                     mu_window = 12,
                                     tail_log_drop = 45,
                                     n_profile = 129L,
                                     t_panels = 16L,
                                     mu_panels = 8L,
                                     gl_order = 10L) {
  s_ref <- sqrt(sse / (n - 1))
  t_ref <- log(s_ref)
  sentinel <- -.Machine$double.xmax

  lp_t <- function(mu, t) {
    vals <- log_post_vec(mu, exp(t)) + t
    vals[is.na(vals)] <- -Inf
    vals
  }

  # Mode in standardized coordinates z = ((mu - x_bar) / s_ref, t - t_ref)
  objective <- function(z) {
    val <- lp_t(x_bar + s_ref * z[1], t_ref + z[2])
    if (is.finite(val)) -val else .Machine$double.xmax
  }
  starts <- list(c(0, 0))
  if (is.finite(init_mu) && is.finite(init_sd) && init_sd > 0) {
    starts[[2L]] <- c((init_mu - x_bar) / s_ref, log(init_sd) - t_ref)
  }
  start_vals <- vapply(starts, objective, numeric(1))
  if (!any(start_vals < .Machine$double.xmax)) return(NULL)
  opt <- stats::optim(
    starts[[which.min(start_vals)]], objective,
    method = "Nelder-Mead",
    control = list(reltol = 1e-12, maxit = 2000L)
  )
  opt_bfgs <- tryCatch(
    stats::optim(opt$par, objective, method = "BFGS",
                 control = list(reltol = 1e-14, maxit = 500L)),
    error = function(e) NULL
  )
  if (!is.null(opt_bfgs) && is.finite(opt_bfgs$value) &&
      opt_bfgs$value <= opt$value) {
    opt <- opt_bfgs
  }
  if (!is.finite(opt$value) || opt$value >= .Machine$double.xmax) return(NULL)
  z_star <- opt$par
  mu_star <- x_bar + s_ref * z_star[1]
  t_star <- t_ref + z_star[2]

  # Laplace covariance in standardized coordinates
  hess <- tryCatch(stats::optimHess(z_star, objective), error = function(e) NULL)
  cov_z <- NULL
  if (!is.null(hess) && all(is.finite(hess))) {
    eig <- eigen(hess, symmetric = TRUE, only.values = TRUE)$values
    if (all(eig > 0)) cov_z <- solve(hess)
  }
  if (is.null(cov_z)) cov_z <- diag(c(1 / n, 1 / (2 * n)))
  sd_t <- sqrt(cov_z[2, 2])
  slope_mu_t <- s_ref * cov_z[1, 2] / cov_z[2, 2]

  # Conditional mode and curvature scale of mu given t
  lik_scale <- function(t) exp(t) / sqrt(n)
  cond_peak <- function(t) {
    h0 <- lik_scale(t)
    guess <- mu_star + slope_mu_t * (t - t_star)
    lo <- min(guess, x_bar) - 30 * h0
    hi <- max(guess, x_bar) + 30 * h0
    probe <- seq(lo, hi, length.out = 121L)
    heights <- lp_t(probe, rep(t, length(probe)))
    if (!any(is.finite(heights))) {
      return(list(center = guess, scale = h0, height = -Inf))
    }
    i <- which.max(heights)
    f <- function(mu) {
      val <- lp_t(mu, t)
      if (is.finite(val)) val else sentinel
    }
    opt_mu <- stats::optimize(
      f,
      c(probe[max(1L, i - 1L)], probe[min(length(probe), i + 1L)]),
      maximum = TRUE,
      tol = h0 * 1e-8
    )
    center <- probe[i]
    height <- heights[i]
    if (opt_mu$objective > height) {
      center <- opt_mu$maximum
      height <- opt_mu$objective
    }
    delta <- 0.05 * h0
    curv <- (f(center + delta) - 2 * height + f(center - delta)) / delta^2
    scale <- if (is.finite(curv) && curv < 0) 1 / sqrt(-curv) else h0
    scale <- min(max(scale, 1e-3 * h0), 3 * h0)
    list(center = center, scale = scale, height = height)
  }

  p_star <- cond_peak(t_star)$height
  if (!is.finite(p_star)) return(NULL)
  h_max <- p_star - t_star

  # Extend log(sigma) until the profile has dropped by tail_log_drop or the
  # prior support ends (located by bisection).
  find_edge <- function(direction) {
    step <- sd_t / 2
    inside <- t_star
    for (i in seq_len(400L)) {
      t <- t_star + direction * i * step
      p <- cond_peak(t)$height
      if (!is.finite(p) || p <= sentinel) {
        outside <- t
        for (j in seq_len(60L)) {
          mid <- (inside + outside) / 2
          if (is.finite(cond_peak(mid)$height)) inside <- mid else outside <- mid
        }
        return(inside)
      }
      if (p <= p_star - tail_log_drop) return(t)
      inside <- t
    }
    NA_real_
  }
  t_lo <- find_edge(-1)
  t_hi <- find_edge(1)
  if (!is.finite(t_lo) || !is.finite(t_hi) || t_hi <= t_lo) return(NULL)

  # Conditional mu windows on a log(sigma) grid, clipped to the prior support
  t_grid <- seq(t_lo, t_hi, length.out = n_profile)
  peaks <- lapply(t_grid, cond_peak)
  centers <- vapply(peaks, `[[`, numeric(1), "center")
  scales <- vapply(peaks, `[[`, numeric(1), "scale")
  win_lo <- centers - mu_window * scales
  win_hi <- centers + mu_window * scales
  support_edge <- function(t, inside, outside) {
    if (is.finite(lp_t(outside, t))) return(outside)
    for (j in seq_len(60L)) {
      mid <- (inside + outside) / 2
      if (is.finite(lp_t(mid, t))) inside <- mid else outside <- mid
    }
    inside
  }
  for (i in seq_along(t_grid)) {
    if (!is.finite(lp_t(centers[i], t_grid[i]))) next
    win_lo[i] <- support_edge(t_grid[i], centers[i], win_lo[i])
    win_hi[i] <- support_edge(t_grid[i], centers[i], win_hi[i])
  }
  win_lo_fn <- stats::approxfun(t_grid, win_lo, rule = 2)
  win_hi_fn <- stats::approxfun(t_grid, win_hi, rule = 2)

  rule_t <- .composite_gauss_legendre_unit(t_panels, gl_order)
  rule_mu <- .composite_gauss_legendre_unit(mu_panels, gl_order)
  m_nodes <- length(rule_mu$nodes)

  # Integral over mu of exp(log_post - h_max) * g at each sigma
  # (sigma-scale density, no Jacobian), restricted by the optional constraint.
  int_mu <- function(sigma, m_fn = NULL, g = NULL) {
    sigma <- as.numeric(sigma)
    out <- numeric(length(sigma))
    t <- suppressWarnings(log(sigma))
    ok <- is.finite(t) & t >= t_lo & t <= t_hi
    if (!any(ok)) return(out)
    a <- win_lo_fn(t[ok])
    b <- win_hi_fn(t[ok])
    if (!is.null(m_fn)) {
      mb <- m_fn(sigma[ok])
      a <- pmax(a, mb$lower)
      b <- pmin(b, mb$upper)
    }
    keep <- is.finite(a) & is.finite(b) & b > a
    if (!any(keep)) return(out)
    width <- b[keep] - a[keep]
    mu <- as.vector(outer(rule_mu$nodes, width) + rep(a[keep], each = m_nodes))
    sig <- rep(sigma[ok][keep], each = m_nodes)
    vals <- exp(log_post_vec(mu, sig) - h_max)
    if (!is.null(g)) vals <- vals * g(mu, sig)
    vals[!is.finite(vals)] <- 0
    idx <- which(ok)[keep]
    out[idx] <- colSums(matrix(vals * rule_mu$weights, nrow = m_nodes)) * width
    out
  }

  int_2d <- function(s_lim, m_fn, g = NULL) {
    t_top <- if (is.finite(s_lim)) {
      if (s_lim > 0) min(t_hi, log(s_lim)) else -Inf
    } else {
      t_hi
    }
    if (!(t_top > t_lo)) return(0)
    span <- t_top - t_lo
    t_nodes <- t_lo + span * rule_t$nodes
    sigma_nodes <- exp(t_nodes)
    sum(span * rule_t$weights * sigma_nodes * int_mu(sigma_nodes, m_fn, g))
  }

  list(
    h_max = h_max,
    map_mu = mu_star,
    map_sig = exp(t_star),
    mu_lower = min(win_lo),
    mu_upper = max(win_hi),
    sigma_lower = exp(t_lo),
    sigma_upper = exp(t_hi),
    int_mu = int_mu,
    int_2d = int_2d
  )
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

  # With informative data, integrate around the posterior mode on local
  # scales; a fixed (mu, sigma) box misses narrow posteriors.
  if (n >= 2 && is.finite(sse) && sse > 0) {
    local <- .generic_local_posterior(
      log_post_vec, n = n, x_bar = x_bar, sse = sse,
      init_mu = init_mu, init_sd = init_sd
    )
    if (!is.null(local)) {
      Z <- local$int_2d(Inf, NULL)
      if (is.finite(Z) && Z > 0) {
        state <- .new_generic_backend_state(
          log_post = log_post,
          log_post_vec = log_post_vec,
          h_max = local$h_max,
          uni_s = local$sigma_upper,
          map_mu = local$map_mu,
          mu_scale = NULL,
          int_2d = local$int_2d,
          int_mu = local$int_mu,
          Z = Z,
          n = n,
          x_bar = x_bar,
          sse = sse,
          mu_lower = local$mu_lower,
          mu_upper = local$mu_upper,
          sigma_lower = local$sigma_lower,
          sigma_upper = local$sigma_upper,
          map_sig = local$map_sig,
          bayestools_priors = prior$bayestools_priors %||% NULL
        )
        .validate_qc_generic_cached_state(state)
        return(state)
      }
    }
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
                            s_lim, m_fn, g = NULL) {
    s_top <- if (is.infinite(s_lim)) sigma_upper_box else min(s_lim, sigma_upper_box)
    sigma_floor <- max(sigma_lower_box, 1e-10)
    if (!is.finite(s_top) || s_top <= sigma_floor || mu_upper_box <= mu_lower_box) {
      return(0)
    }

    integrand <- function(x) {
      mu <- x[1, ]
      sigma <- x[2, ]

      log_p <- log_post_vec(mu, sigma)
      valid <- rep(TRUE, length(mu))
      if (!is.null(m_fn)) {
        mb <- m_fn(sigma)
        valid <- (mu >= mb$lower) & (mu <= mb$upper)
      }

      result <- rep(0, length(mu))
      result[valid] <- exp(log_p[valid] - h_max)
      if (!is.null(g)) result[valid] <- result[valid] * g(mu[valid], sigma[valid])
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
  int_2d <- function(s_lim, m_fn, g = NULL) {
    integrate_box(mu_lower, mu_upper, sigma_lower, sigma_upper, s_lim, m_fn, g)
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

