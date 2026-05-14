# ==============================================================================
# Survival-Function Solvers
# ==============================================================================

# ==============================================================================
# Solver Factory: Returns function P(Index > c)
# ==============================================================================

#' Create solver function for P(Index > c)
#' @param data Numeric vector of observations
#' @param LSL Lower specification limit
#' @param USL Upper specification limit
#' @param prior Prior object (PriorConjugate or PriorGeneric)
#' @param metric Capability index name
#' @param target Target value for Cpm
#' @param cached_state Pre-computed state for PriorGeneric (optional)
#' @param sigma_level Number of process standard deviations used to define the
#'   capability metric scale.
#' @param ... Additional arguments.
#' @return Function that takes threshold c and returns P(Index > c)
#' @keywords internal
make_solver <- function(data, LSL, USL, prior, metric = "Cpk", target = NULL,
                        cached_state = NULL, sigma_level = 3, ...) {
  .validate_capability_request(
    LSL = LSL,
    USL = USL,
    target = target,
    sigma_level = sigma_level,
    metric = metric,
    target_required = FALSE,
    sigma_name = "sigma_level"
  )

  UseMethod("make_solver", prior)
}

#' @export
make_solver.PriorConjugate <- function(data, LSL, USL, prior, metric = "Cpk",
                                        target = NULL, cached_state = NULL,
                                        sigma_level = 3, ...) {
  posterior_info <- .compute_validated_conjugate_posterior(
    prior, data, cached_state,
    context = "The conjugate posterior for survival-function computation"
  )
  ss <- posterior_info$ss
  post <- posterior_info$post
  n <- ss$n; k_n <- post$k_n; mu_n <- post$mu_n
  alpha_n <- post$alpha_n; beta_n <- post$beta_n
  df_p <- 2 * alpha_n
  tol <- USL - LSL

  if (posterior_info$is_degenerate) {
    dist <- .degenerate_conjugate_metric_distribution(
      mu_n, k_n, LSL, USL, target, metric,
      sigma_level = sigma_level
    )
    if (!is.null(dist)) {
      return(.degenerate_metric_solver(dist))
    }
  }

  # Analytic fast paths for Cp, Cpu, Cpl
  if (metric == "Cp") {
    return(function(c) {
      if (c <= 0) return(1.0)
      stats::pchisq(8 * beta_n * sigma_level^2 * c^2 / tol^2,
                    df = df_p, lower.tail = FALSE)
    })
  }
  if (metric == "Cpu") {
    x_U <- (USL - mu_n) * sqrt(k_n * alpha_n / beta_n)
    sqrt_k_n <- sqrt(k_n)
    return(function(c) {
      suppressWarnings(stats::pt(x_U, df = df_p,
                                 ncp = sigma_level * c * sqrt_k_n))
    })
  }
  if (metric == "Cpl") {
    x_L <- (mu_n - LSL) * sqrt(k_n * alpha_n / beta_n)
    sqrt_k_n <- sqrt(k_n)
    return(function(c) {
      suppressWarnings(stats::pt(x_L, df = df_p,
                                 ncp = sigma_level * c * sqrt_k_n))
    })
  }

  if (metric %in% c("Cpm", "Cpc")) {
    if (is.null(target)) target <- (LSL + USL) / 2
    cpm_dist <- .cpm_spec_distance(LSL, USL, target)

    cpc_width_matrix <- function(sigma, c_vec) {
      A <- tol / ((2 * sigma_level) * sqrt(pi / 2))
      sigma_mat <- matrix(rep(sigma, each = length(c_vec)), nrow = length(c_vec))
      c_mat <- matrix(rep(c_vec, times = length(sigma)), nrow = length(c_vec))
      K <- A / (c_mat * sigma_mat)

      width <- matrix(NA_real_, nrow = length(c_vec), ncol = length(sigma))
      feasible <- is.finite(K) & K >= .cpc_lookup$g0
      if (!any(feasible)) {
        return(width)
      }

      Kf <- K[feasible]
      z <- numeric(length(Kf))
      use_interp <- Kf <= .cpc_lookup$g_max
      if (any(use_interp)) {
        z[use_interp] <- stats::approx(
          x = .cpc_lookup$g,
          y = .cpc_lookup$z,
          xout = Kf[use_interp],
          ties = "ordered",
          rule = 2
        )$y
      }
      if (any(!use_interp)) {
        z[!use_interp] <- Kf[!use_interp]
      }

      width[feasible] <- sigma_mat[feasible] * z
      width
    }

    return(function(c) {
      c <- as.numeric(c)
      result <- numeric(length(c))
      finite <- is.finite(c)
      result[finite & c <= 0] <- 1

      use <- finite & c > 0
      if (!any(use)) {
        return(result)
      }

      c_vec <- c[use]
      if (length(c_vec) == 1L) {
        ci <- c_vec
        constr <- get_metric_constraints(metric, ci, LSL, USL, target,
                                         sigma_level = sigma_level)
        s_max <- constr$s_max_fn()
        vals <- if (!is.infinite(s_max) && s_max <= 0) {
          0
        } else {
          y_min <- if (is.infinite(s_max)) 0 else (2 * beta_n) / (s_max^2)
          h_max <- stats::dchisq(max(df_p - 2, 1e-6), df_p, log = TRUE)
          log_int <- function(y) {
            sigma <- sqrt((2 * beta_n) / y)
            sd_mu <- sigma / sqrt(k_n)
            mb <- constr$mu_b_fn_vec(sigma)
            z_U <- ifelse(is.infinite(mb$upper), Inf, (mb$upper - mu_n) / sd_mu)
            z_L <- ifelse(is.infinite(mb$lower), -Inf, (mb$lower - mu_n) / sd_mu)
            log_prob <- log_diff_exp(
              stats::pnorm(z_U, log.p = TRUE),
              stats::pnorm(z_L, log.p = TRUE)
            )
            invalid <- mb$lower >= mb$upper
            log_prob[invalid] <- -Inf
            log_prob + stats::dchisq(y, df_p, log = TRUE)
          }
          scaled <- stats::integrate(
            function(y) {
              vals <- exp(log_int(y) - h_max)
              vals[!is.finite(vals)] <- 0
              vals
            },
            y_min,
            Inf
          )$value
          if (scaled <= 0) 0 else exp(h_max + log(scaled))
        }
        result[use] <- pmin(pmax(vals, 0), 1)
        return(result)
      }

      eps <- 1e-10
      u_grid <- sort(unique(c(
        seq(eps, 1 - eps, length.out = 4096L),
        1 - 10^seq(-10, -3, length.out = 512L)
      )))
      u_edges <- c(0, (u_grid[-1] + u_grid[-length(u_grid)]) / 2, 1)
      u_weights <- diff(u_edges)

      survival_matrix <- function(u) {
        y <- stats::qchisq(u, df = df_p)
        sigma <- sqrt((2 * beta_n) / y)
        sd_mu <- sigma / sqrt(k_n)

        if (metric == "Cpm") {
          R <- cpm_dist / (sigma_level * c_vec)
          width_sq <- outer(R^2, sigma^2, "-")
          valid_width <- is.finite(width_sq) & width_sq > 0
          width <- sqrt(pmax(width_sq, 0))
        } else {
          width <- cpc_width_matrix(sigma, c_vec)
          valid_width <- is.finite(width)
        }

        lower <- target - width
        upper <- target + width
        lower[!valid_width] <- 0
        upper[!valid_width] <- -1

        z_U <- sweep(upper - mu_n, 2, sd_mu, "/")
        z_L <- sweep(lower - mu_n, 2, sd_mu, "/")
        log_prob <- log_diff_exp(
          stats::pnorm(z_U, log.p = TRUE),
          stats::pnorm(z_L, log.p = TRUE)
        )
        prob <- exp(log_prob)
        prob[!is.finite(prob) | z_L >= z_U] <- 0
        prob
      }

      vals <- tryCatch({
        prob <- survival_matrix(u_grid)
        as.numeric(prob %*% u_weights)
      }, error = function(e) rep(NA_real_, length(c_vec)))
      vals[!is.finite(vals)] <- 0
      result[use] <- pmin(pmax(vals, 0), 1)
      result
    })
  }

  # Pre-compute global h_max (chi-square mode density for numerical stability)
  y_mode <- max(df_p - 2, 1e-6)
  h_max_global <- stats::dchisq(y_mode, df_p, log = TRUE)

  function(c) {
    if (!.metric_can_be_negative(metric) && c <= 0) return(1.0)

    constr <- get_metric_constraints(metric, c, LSL, USL, target,
                                     sigma_level = sigma_level)
    s_max <- constr$s_max_fn()
    if (!is.infinite(s_max) && s_max <= 0) return(0.0)

    y_min <- if (is.infinite(s_max)) 0 else (2 * beta_n) / (s_max^2)

    # Vectorized log integrand
    log_int <- function(y) {
      # All operations vectorized
      sigma <- sqrt((2 * beta_n) / y)
      sd_mu <- sigma / sqrt(k_n)

      # Get bounds for all sigma values at once
      mb <- constr$mu_b_fn_vec(sigma)
      mb_L <- mb$lower
      mb_U <- mb$upper

      # Identify valid intervals
      valid <- mb_L < mb_U

      # Initialize result with -Inf

      result <- rep(-Inf, length(y))

      if (!any(valid)) return(result)

      # Compute z-scores (vectorized)
      z_U <- ifelse(is.infinite(mb_U), Inf, (mb_U - mu_n) / sd_mu)
      z_L <- ifelse(is.infinite(mb_L), -Inf, (mb_L - mu_n) / sd_mu)

      # log_diff_exp(stats::pnorm(z_U, log.p=TRUE), stats::pnorm(z_L, log.p=TRUE)) + stats::dchisq(y, df_p, log=TRUE)
      log_prob <- log_diff_exp(stats::pnorm(z_U, log.p = TRUE), stats::pnorm(z_L, log.p = TRUE))
      result[valid] <- log_prob[valid] + stats::dchisq(y[valid], df_p, log = TRUE)
      result
    }

    h_max <- h_max_global

    safe_integrand <- function(y) {
      vals <- exp(log_int(y) - h_max)
      vals[!is.finite(vals)] <- 0
      return(vals)
    }

    res <- stats::integrate(safe_integrand, y_min, Inf)$value
    if (res <= 0) return(0.0)
    return(exp(h_max + log(res)))
  }
}

#' @export
make_solver.PriorSemiConjugateMu <- function(data, LSL, USL, prior,
                                              metric = "Cpk", target = NULL,
                                              cached_state = NULL,
                                              sigma_level = 3, ...) {
  # Case 2: Conjugate mu (Normal/Jeffreys), non-conjugate sigma
  # Integrate mu analytically using pnorm, then 1D numerical over sigma

  state <- .prepare_semi_mu_state(data, prior, cached_state)
  k_n <- state$k_n
  mu_n <- state$mu_n
  log_Z <- state$log_Z
  log_p_sigma <- state$log_p_sigma
  sigma_support <- state$sigma_support
  sigma_domain <- state$sigma_domain
  base_log_h_max <- state$base_log_h_max

  function(c) {
    if (!.metric_can_be_negative(metric) && c <= 0) return(1.0)
    constr <- get_metric_constraints(metric, c, LSL, USL, target,
                                     sigma_level = sigma_level)
    s_max <- constr$s_max_fn()
    if (!is.infinite(s_max) && s_max <= 0) return(0.0)

    sigma_lower <- max(sigma_support$lower, sigma_domain$lower)
    sigma_upper <- min(
      if (is.infinite(s_max)) sigma_domain$upper else max(s_max, 0),
      sigma_support$upper,
      sigma_domain$upper
    )
    if (!is.finite(sigma_upper) || sigma_upper <= sigma_lower) return(0.0)

    integrand <- function(t) {
      sigma <- exp(t)
      # P(mu in constraint region | sigma, data)
      sd_mu <- sigma / sqrt(k_n)
      mb <- constr$mu_b_fn_vec(sigma)
      z_U <- (mb$upper - mu_n) / sd_mu
      z_L <- (mb$lower - mu_n) / sd_mu

      # Safety check for impossible regions
      log_prob_mu <- ifelse(z_L >= z_U, -Inf,
                             log_diff_exp(stats::pnorm(z_U, log.p = TRUE),
                                          stats::pnorm(z_L, log.p = TRUE)))

      vals <- exp(log_p_sigma(sigma) + log_prob_mu + t - base_log_h_max)
      vals[!is.finite(vals)] <- 0
      vals
    }

    prob_scaled <- stats::integrate(
      integrand,
      log(sigma_lower),
      log(sigma_upper),
      rel.tol = 1e-5,
      subdivisions = 400
    )$value
    prob <- exp(base_log_h_max - log_Z) * prob_scaled
    min(max(prob, 0), 1)
  }
}

.sigma_interval_log_prob_inv_gamma <- function(lower, upper, shape, rate) {
  log_prob <- rep(-Inf, length(lower))

  full_region <- lower <= 0 & is.infinite(upper) & upper > 0
  log_prob[full_region] <- 0

  upper_only <- lower <= 0 & is.finite(upper) & upper > 0
  if (any(upper_only)) {
    log_prob[upper_only] <- stats::pgamma(
      1 / upper[upper_only]^2,
      shape = shape,
      rate = rate[upper_only],
      lower.tail = FALSE,
      log.p = TRUE
    )
  }

  lower_only <- is.finite(lower) & lower > 0 & is.infinite(upper) & upper > 0
  if (any(lower_only)) {
    log_prob[lower_only] <- stats::pgamma(
      1 / lower[lower_only]^2,
      shape = shape,
      rate = rate[lower_only],
      lower.tail = TRUE,
      log.p = TRUE
    )
  }

  finite_band <- is.finite(lower) & lower > 0 &
    is.finite(upper) & upper > lower
  if (any(finite_band)) {
    lower_band <- lower[finite_band]
    upper_band <- upper[finite_band]
    rate_band <- rate[finite_band]

    log_cdf_lower <- stats::pgamma(
      1 / lower_band^2,
      shape = shape,
      rate = rate_band,
      lower.tail = TRUE,
      log.p = TRUE
    )
    log_cdf_upper <- stats::pgamma(
      1 / upper_band^2,
      shape = shape,
      rate = rate_band,
      lower.tail = TRUE,
      log.p = TRUE
    )
    log_surv_upper <- stats::pgamma(
      1 / upper_band^2,
      shape = shape,
      rate = rate_band,
      lower.tail = FALSE,
      log.p = TRUE
    )
    log_surv_lower <- stats::pgamma(
      1 / lower_band^2,
      shape = shape,
      rate = rate_band,
      lower.tail = FALSE,
      log.p = TRUE
    )

    log_mass_cdf <- log_diff_exp(log_cdf_lower, log_cdf_upper)
    log_mass_surv <- log_diff_exp(log_surv_upper, log_surv_lower)

    log_band <- log_mass_cdf
    use_surv <- !is.finite(log_band) | (is.finite(log_mass_surv) & log_mass_surv > log_band)
    log_band[use_surv] <- log_mass_surv[use_surv]
    log_prob[finite_band] <- log_band
  }

  log_prob[!is.finite(log_prob) & !full_region] <- -Inf
  log_prob
}

.compute_semi_sigma_log_Z_survival <- function(n, x_bar, sse, alpha0, beta0,
                                               log_dens_mu, jeffreys_adj = 0,
                                               mu_domain,
                                               context = "The semi-conjugate sigma posterior") {
  alpha_n_eff <- alpha0 + n / 2 + jeffreys_adj

  .validate_semi_sigma_posterior(
    n, x_bar, sse, alpha0, beta0, log_dens_mu,
    jeffreys_adj = jeffreys_adj,
    context = context
  )

  log_integrand <- function(mu) {
    sse_mu <- sse + n * (mu - x_bar)^2
    beta_n <- beta0 + sse_mu / 2
    vals <- log_dens_mu(mu) - alpha_n_eff * log(beta_n)
    vals[!is.finite(vals)] <- -Inf
    vals
  }

  eff_bounds <- .semi_sigma_effective_bounds(
    log_integrand,
    mu_domain,
    lower = -Inf,
    upper = Inf
  )

  .semi_sigma_integrate_log_kernel(
    log_integrand,
    mu_domain,
    lower = eff_bounds[["lower"]],
    upper = eff_bounds[["upper"]],
    rel.tol = 1e-6,
    subdivisions = 400,
    context = context
  )$log_value
}

.semi_sigma_effective_bounds <- function(log_integrand, mu_domain,
                                         lower = -Inf, upper = Inf,
                                         probe_n = 1025L,
                                         tail_log_drop = 36,
                                         max_expansions = 32L) {
  probe_bounds <- .semi_sigma_probe_interval(mu_domain, lower, upper)
  lower_limit <- probe_bounds[["lower"]]
  upper_limit <- probe_bounds[["upper"]]
  if (!is.finite(lower_limit) || !is.finite(upper_limit) || upper_limit <= lower_limit) {
    return(c(lower = lower, upper = upper))
  }

  probe_mu <- seq(lower_limit, upper_limit, length.out = probe_n)
  probe_log <- log_integrand(probe_mu)
  finite_idx <- which(is.finite(probe_log))
  if (!length(finite_idx)) {
    return(c(lower = lower_limit, upper = upper_limit))
  }

  best_idx <- finite_idx[which.max(probe_log[finite_idx])]
  center <- probe_mu[best_idx]
  h_max <- probe_log[best_idx]

  local_lower <- probe_mu[max(best_idx - 1L, 1L)]
  local_upper <- probe_mu[min(best_idx + 1L, probe_n)]
  if (local_upper > local_lower) {
    opt <- tryCatch(
      stats::optimize(
        function(mu) {
          val <- log_integrand(mu)
          if (is.finite(val)) -val else .Machine$double.xmax
        },
        interval = c(local_lower, local_upper)
      ),
      error = function(e) NULL
    )
    if (!is.null(opt)) {
      center <- opt$minimum
      h_max <- max(h_max, -opt$objective)
    }
  }

  step <- max((upper_limit - lower_limit) / max(probe_n - 1L, 1L), 1) * 8
  eff_lower <- center
  eff_upper <- center

  for (i in seq_len(max_expansions)) {
    if (is.finite(lower) && eff_lower <= lower) break
    candidate <- eff_lower - step
    if (is.finite(lower)) {
      candidate <- max(lower, candidate)
    }
    val <- log_integrand(candidate)
    eff_lower <- candidate
    if (!is.finite(val) || val <= h_max - tail_log_drop) break
  }

  for (i in seq_len(max_expansions)) {
    if (is.finite(upper) && eff_upper >= upper) break
    candidate <- eff_upper + step
    if (is.finite(upper)) {
      candidate <- min(upper, candidate)
    }
    val <- log_integrand(candidate)
    eff_upper <- candidate
    if (!is.finite(val) || val <= h_max - tail_log_drop) break
  }

  c(lower = eff_lower, upper = eff_upper)
}

#' @export
make_solver.PriorSemiConjugateSigma <- function(data, LSL, USL, prior,
                                                 metric = "Cpk", target = NULL,
                                                 cached_state = NULL,
                                                 sigma_level = 3, ...) {
  # Case 3: Non-conjugate mu, conjugate sigma (InvGamma/Jeffreys)
  # Integrate sigma analytically using Gamma functions, then 1D over mu

  state <- .prepare_semi_sigma_state(
    data, prior, cached_state,
    context = "The semi-conjugate sigma posterior for survival-function computation"
  )
  n <- state$n
  x_bar <- state$x_bar
  sse <- state$sse
  beta0 <- state$beta0
  jeffreys_adj <- state$jeffreys_adj
  alpha_n_eff <- state$alpha_n_eff
  if (is.null(target)) target <- (LSL + USL) / 2
  # Normalize to the same posterior mass used by the density solver.
  log_Z <- .compute_semi_sigma_log_Z_survival(
    n, x_bar, sse, state$alpha0, beta0,
    prior$log_dens_mu, jeffreys_adj,
    mu_domain = state$mu_domain,
    context = "The semi-conjugate sigma posterior for normalization"
  )

  function(c) {
    if (!.metric_can_be_negative(metric) && c <= 0) return(1.0)

    mu_low <- -Inf
    mu_high <- Inf
    if (metric == "Cpm") {
      mu_radius <- .cpm_spec_distance(LSL, USL, target) / (sigma_level * c)
      mu_low <- target - mu_radius
      mu_high <- target + mu_radius
    } else if (metric == "Cpc") {
      mu_radius <- (USL - LSL) / ((2 * sigma_level) * sqrt(pi / 2) * c)
      mu_low <- target - mu_radius
      mu_high <- target + mu_radius
    }

    if (is.na(mu_low) || is.na(mu_high) ||
        (is.finite(mu_low) && is.finite(mu_high) && mu_high <= mu_low)) {
      return(0)
    }

    log_integrand <- function(mu) {
      sse_mu <- sse + n * (mu - x_bar)^2
      beta_n <- beta0 + sse_mu / 2
      log_weight_mu <- prior$log_dens_mu(mu) - alpha_n_eff * log(beta_n)

      sigma_region <- .metric_sigma_region_from_mu(
        metric, mu, c, LSL, USL, target,
        sigma_level = sigma_level
      )
      log_prob_sigma <- .sigma_interval_log_prob_inv_gamma(
        sigma_region$lower, sigma_region$upper,
        shape = alpha_n_eff, rate = beta_n
      )

      vals <- log_weight_mu + log_prob_sigma - log_Z
      vals[!is.finite(vals)] <- -Inf
      vals
    }
    eff_bounds <- .semi_sigma_effective_bounds(
      log_integrand,
      state$mu_domain,
      lower = mu_low,
      upper = mu_high
    )

    prob <- .semi_sigma_integrate_log_kernel(
      log_integrand,
      state$mu_domain,
      lower = eff_bounds[["lower"]],
      upper = eff_bounds[["upper"]],
      rel.tol = 1e-5,
      subdivisions = 400
    )$value
    if (!is.finite(prob) || prob <= 0) return(0)

    min(max(prob, 0), 1)
  }
}

# ==============================================================================
# Density Solver Factory: Direct PDF via Contour Integration
# ==============================================================================

