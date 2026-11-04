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
#' @param cached_state Optional precomputed integration state.
#' @param sigma_level Number of process standard deviations used to define the
#'   capability metric scale.
#' @param ... Additional arguments
#' @return Function pdf(c) returning the marginal posterior density at c
#' @keywords internal
make_density_solver <- function(data, LSL, USL, prior, metric = "Cpk",
                                target = NULL, cached_state = NULL,
                                sigma_level = 3, ...) {
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
        .integration_eval(
          S, base,
          context = "The integration survival solver"
        ) -
          .integration_eval(
            S, upper[use_forward],
            context = "The integration survival solver"
          )
      ) / (upper[use_forward] - base)
    }

    if (any(use_central)) {
      density[use_central] <- (
        .integration_eval(
          S, lower[use_central],
          context = "The integration survival solver"
        ) -
          .integration_eval(
            S, upper[use_central],
            context = "The integration survival solver"
          )
      ) / (upper[use_central] - lower[use_central])
    }

    density[!is.finite(density) | density < 0] <- 0
    result[valid] <- density
    result
  }
}

.contour_rules <- function() {
  list(
    rule = .composite_gauss_legendre_unit(3L, 10L),
    check_rule = .composite_gauss_legendre_unit(4L, 10L)
  )
}

# Integrates exp(log_kernel) along many contours at once. log_kernel(t, index)
# receives a matrix with one column of contour positions per element of index.
# Every step works on such a matrix, so the kernel is evaluated in a handful of
# batches regardless of the number of contours.
.contour_mass <- function(log_kernel, rules, center, t_lower, t_upper,
                          tail_log_drop = 40) {
  probe_unit <- seq(0, 1, length.out = 65L)
  zoom_unit <- seq(0, 1, length.out = 17L)
  reach_offsets <- 2^seq(-20, 10)
  n_contours <- length(center)

  column_grid <- function(unit, from, to) {
    outer(unit, to - from) + rep(from, each = length(unit))
  }
  best_row <- function(vals) {
    max.col(t(pmax(vals, -.Machine$double.xmax)), ties.method = "first")
  }

  # Locate the contour peak so that the quadrature does not miss a narrow
  # posterior, widening the probe while the peak sits on its edge.
  center <- pmin(pmax(center, t_lower), t_upper)
  probe_lower <- pmax(t_lower, center - 8)
  probe_upper <- pmin(t_upper, center + 8)
  mode <- center
  height <- rep(-Inf, n_contours)
  active <- seq_len(n_contours)
  for (iter in seq_len(17L)) {
    t <- column_grid(probe_unit, probe_lower[active], probe_upper[active])
    vals <- log_kernel(t, active)
    best <- best_row(vals)
    peak <- cbind(best, seq_along(active))
    mode[active] <- t[peak]
    height[active] <- vals[peak]

    missed <- vals[peak] == -Inf
    widen_left <- probe_lower[active] > t_lower[active] &
      (best <= 2L | missed)
    widen_right <- probe_upper[active] < t_upper[active] &
      (best >= length(probe_unit) - 1L | missed)
    if (!any(widen_left | widen_right)) break

    width <- probe_upper[active] - probe_lower[active]
    probe_lower[active[widen_left]] <- pmax(
      t_lower[active[widen_left]],
      probe_lower[active[widen_left]] - width[widen_left]
    )
    probe_upper[active[widen_right]] <- pmin(
      t_upper[active[widen_right]],
      probe_upper[active[widen_right]] + width[widen_right]
    )
    active <- active[widen_left | widen_right]
  }

  # No finite kernel value anywhere on a contour: it misses the support.
  density <- numeric(n_contours)
  on_support <- which(height > -Inf)
  if (length(on_support) == 0L) return(density)
  mode <- mode[on_support]
  height <- height[on_support]
  t_lower <- t_lower[on_support]
  t_upper <- t_upper[on_support]
  step <- ((probe_upper - probe_lower) / (length(probe_unit) - 1L))[on_support]

  for (pass in seq_len(5L)) {
    from <- pmax(t_lower, mode - step)
    to <- pmin(t_upper, mode + step)
    t <- column_grid(zoom_unit, from, to)
    vals <- log_kernel(t, on_support)
    peak <- cbind(best_row(vals), seq_along(on_support))
    mode <- t[peak]
    height <- vals[peak]
    step <- (to - from) / (length(zoom_unit) - 1L)
  }

  scan_offsets <- function(direction) {
    limit <- if (direction > 0) t_upper else t_lower
    clamp <- if (direction > 0) pmin else pmax
    t_reach <- clamp(
      rep(mode, each = length(reach_offsets)) +
        direction * reach_offsets,
      rep(limit, each = length(reach_offsets))
    )
    dim(t_reach) <- c(length(reach_offsets), length(on_support))
    log_kernel(t_reach, on_support)
  }
  scan_below <- scan_offsets(-1)
  scan_above <- scan_offsets(1)

  # Distance from the peak at which the kernel has dropped by log_drop,
  # keeping all of the support when it never does.
  reach_of_drop <- function(direction, scan, log_drop, passes) {
    limit <- if (direction > 0) t_upper else t_lower
    clamp <- if (direction > 0) pmin else pmax
    above <- scan > rep(height - log_drop, each = length(reach_offsets))
    last <- max.col(t(above * seq_along(reach_offsets)), ties.method = "first")
    inside <- ifelse(above[cbind(last, seq_along(last))], reach_offsets[last], 0)
    outside <- reach_offsets[pmin(last + 1L, length(reach_offsets))]
    # A smooth drop needs only a rough position, but a hard support edge such
    # as a truncated prior has to be pinned down.
    active <- seq_along(on_support)
    for (pass in seq_len(passes)) {
      if (pass == 9L) {
        t_outside <- matrix(clamp(mode + direction * outside, limit), nrow = 1L)
        active <- which(is.infinite(log_kernel(t_outside, on_support)))
        if (length(active) == 0L) break
      }
      middle <- (inside[active] + outside[active]) / 2
      t_middle <- matrix(
        clamp(mode[active] + direction * middle, limit[active]),
        nrow = 1L
      )
      reached <- as.vector(log_kernel(t_middle, on_support[active])) >
        height[active] - log_drop
      inside[active[reached]] <- middle[reached]
      outside[active[!reached]] <- middle[!reached]
    }
    abs(clamp(mode + direction * inside, limit) - mode)
  }

  reach_below <- reach_of_drop(-1, scan_below, tail_log_drop, 30L)
  reach_above <- reach_of_drop(1, scan_above, tail_log_drop, 30L)
  scale_below <- reach_of_drop(-1, scan_below, 1, 4L)
  scale_above <- reach_of_drop(1, scan_above, 1, 4L)

  side_mass <- function(direction, rule, reach, scale) {
    side <- .geometric_side_rule(rule, reach, scale)
    t_nodes <- rep(mode, each = nrow(side$offsets)) + direction * side$offsets
    scaled <- exp(log_kernel(t_nodes, on_support) -
                    rep(height, each = nrow(side$offsets)))
    colSums(scaled * side$weights)
  }
  peak_mass <- function(rule) {
    side_mass(-1, rule, reach_below, scale_below) +
      side_mass(1, rule, reach_above, scale_above)
  }

  # The fixed rules can miss a second narrow mode on the contour, which shows
  # as disagreement between them; such columns are integrated adaptively
  # unless they are negligible next to the others.
  scaled <- peak_mass(rules$check_rule)
  log_density <- height + log(scaled)
  negligible <- is.finite(log_density) &
    log_density < max(log_density[is.finite(log_density)], -Inf) - 30
  unresolved <- which(
    !negligible &
      (!is.finite(scaled) | abs(peak_mass(rules$rule) - scaled) > 1e-6 * scaled)
  )
  for (i in unresolved) {
    integrand <- function(u) {
      exp(as.vector(log_kernel(matrix(u + mode[i]), on_support[i])) - height[i])
    }
    integrate_part <- function(lower, upper) {
      if (lower >= upper) return(0)
      stats::integrate(integrand, lower, upper, rel.tol = 1e-7,
                       abs.tol = 1e-10, subdivisions = 400L)$value
    }
    # A failed quadrature is reported as NA (not 0) so callers can fall back
    # to the survival solver instead of silently dropping density.
    scaled[i] <- tryCatch(
      integrate_part(-reach_below[i], 0) + integrate_part(0, reach_above[i]),
      error = function(e) NA_real_
    )
  }

  density[on_support] <- ifelse(
    is.finite(scaled),
    exp(height + log(scaled)),
    NA_real_
  )
  density
}

# For Cpu/Cpl, mu is USL - sigma_level*c*sigma or
# LSL + sigma_level*c*sigma. The change-of-variables Jacobian is
# sigma_level*sigma. Cpk adds both branches up to their meeting at the
# specification midpoint (sigma <= (USL-LSL)/(2*sigma_level*c) when c > 0).
.linear_metric_density <- function(log_joint, LSL, USL, metric, sigma_level,
                                   sigma_ref, sigma_lower = 0, sigma_upper = Inf,
                                   mu_lower = -Inf, mu_upper = Inf) {
  force(log_joint)
  branches <- switch(metric, Cpu = "upper", Cpl = "lower",
                     Cpk = c("lower", "upper"))
  rules <- .contour_rules()

  branch_density <- function(c, branch) {
    intercept <- if (branch == "upper") USL else LSL
    slope <- if (branch == "upper") -sigma_level * c else sigma_level * c
    lower <- rep(sigma_lower, length(c))
    upper <- rep(sigma_upper, length(c))
    if (metric == "Cpk") {
      positive <- c > 0
      upper[positive] <- pmin(
        upper[positive],
        (USL - LSL) / (2 * sigma_level * c[positive])
      )
    }
    # Generic posterior normalization uses a finite (mu, sigma) box.
    # Restrict the contour to that same box, including for negative c.
    flat <- slope == 0
    if (intercept < mu_lower || intercept > mu_upper) upper[flat] <- 0
    if (is.finite(mu_lower) || is.finite(mu_upper)) {
      at_mu_lower <- (mu_lower - intercept) / slope[!flat]
      at_mu_upper <- (mu_upper - intercept) / slope[!flat]
      lower[!flat] <- pmax(lower[!flat], pmin(at_mu_lower, at_mu_upper))
      upper[!flat] <- pmin(upper[!flat], pmax(at_mu_lower, at_mu_upper))
    }

    density <- numeric(length(c))
    contour <- which(upper > lower & upper > 0)
    if (length(contour) == 0L) return(density)
    contour_slope <- slope[contour]

    log_kernel <- function(t, index) {
      sigma <- exp(t)
      mu <- intercept + rep(contour_slope[index], each = nrow(t)) * sigma
      vals <- matrix(-Inf, nrow(t), ncol(t))
      valid <- is.finite(mu) & is.finite(sigma) & sigma > 0
      if (any(valid)) {
        # The second sigma factor is d sigma / d log(sigma).
        vals[valid] <- log_joint(mu[valid], sigma[valid]) +
          log(sigma_level) + 2 * t[valid]
      }
      vals[!is.finite(vals)] <- -Inf
      vals
    }

    density[contour] <- .contour_mass(
      log_kernel, rules,
      center = rep(log(sigma_ref), length(contour)),
      t_lower = log(pmax(lower[contour], 0)),
      t_upper = log(upper[contour])
    )
    density
  }

  function(c) {
    c <- as.numeric(c)
    result <- numeric(length(c))
    valid <- is.finite(c)
    if (!any(valid)) return(result)
    result[valid] <- Reduce(`+`, lapply(branches, function(branch) {
      branch_density(c[valid], branch)
    }))
    result
  }
}

# Cpm and Cpc contours are arcs around (target, 0). With z = |mu - target| /
# sigma they are sigma(z, c), one branch on either side of the target, and
# both metrics have the Jacobian sigma^2 / c in z. The peak of a contour can
# lie at any z, which grows without bound as the posterior narrows away from
# the target, so the contours are integrated over log(z).
.radial_metric_density <- function(log_joint, LSL, USL, target, metric,
                                   sigma_level, sigma_ref, mu_ref,
                                   sigma_lower = 0, sigma_upper = Inf) {
  force(log_joint)
  rules <- .contour_rules()
  tol <- USL - LSL
  cpm_dist <- .cpm_spec_distance(LSL, USL, target)

  sigma_of_z <- function(z, c) {
    if (metric == "Cpm") {
      cpm_dist / (sigma_level * c * sqrt(1 + z^2))
    } else {
      tol / ((2 * sigma_level) * sqrt(pi / 2) * c * .cpc_g(z))
    }
  }
  z_of_sigma <- function(sigma, c) {
    z <- if (metric == "Cpm") {
      sqrt(pmax((cpm_dist / (sigma_level * c * sigma))^2 - 1, 0))
    } else {
      .cpc_mu_width_from_sigma(sigma, c, tol, sigma_level) / sigma
    }
    z[is.na(z)] <- 0
    z[sigma <= 0] <- Inf
    z
  }

  branch_density <- function(c, direction, t_lower, t_upper) {
    log_kernel <- function(t, index) {
      z <- exp(t)
      c_mat <- rep(c[index], each = nrow(t))
      sigma <- sigma_of_z(z, c_mat)
      vals <- matrix(-Inf, nrow(t), ncol(t))
      valid <- is.finite(z) & is.finite(sigma) & sigma > 0
      if (any(valid)) {
        mu <- target + direction * sigma[valid] * z[valid]
        vals[valid] <- log_joint(mu, sigma[valid]) +
          2 * log(sigma[valid]) - log(c_mat[valid]) + t[valid]
      }
      vals[!is.finite(vals)] <- -Inf
      vals
    }

    z_ref <- max(abs(mu_ref - target) / sigma_ref, 0.1)
    .contour_mass(
      log_kernel, rules,
      center = rep(log(z_ref), length(c)),
      t_lower = t_lower,
      t_upper = t_upper
    )
  }

  function(c) {
    c <- as.numeric(c)
    result <- numeric(length(c))
    valid <- is.finite(c) & c > 0
    if (!any(valid)) return(result)

    cv <- c[valid]
    # Padded so that the interpolated Cpc geometry never clips the contour.
    t_lower <- suppressWarnings(log(z_of_sigma(rep(sigma_upper, length(cv)), cv))) - 1e-6
    t_upper <- suppressWarnings(log(z_of_sigma(rep(sigma_lower, length(cv)), cv))) + 1e-6
    has_arc <- t_upper > t_lower
    if (!any(has_arc)) return(result)

    vals <- numeric(length(cv))
    vals[has_arc] <- Reduce(`+`, lapply(c(-1, 1), function(direction) {
      branch_density(cv[has_arc], direction, t_lower[has_arc], t_upper[has_arc])
    }))
    result[valid] <- vals
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
  conjugate_log_joint <- function(mu, sigma) {
    -(2 * alpha_n + 1) * log(sigma) - beta_n / sigma^2 -
      log_Z_sigma + stats::dnorm(mu, mu_n, sigma / sqrt(k_n), log = TRUE)
  }

  # Cp: no integration needed (contour is a single sigma point)
  if (metric == "Cp") {
    return(function(c) {
      c <- as.numeric(c)
      result <- numeric(length(c))
      valid <- is.finite(c) & c > 0
      if (!any(valid)) return(result)

      sigma_c <- tol / ((2 * sigma_level) * c[valid])
      valid_sigma <- is.finite(sigma_c) & sigma_c > 0
      vals <- numeric(length(sigma_c))
      if (any(valid_sigma)) {
        sv <- sigma_c[valid_sigma]
        cv <- c[valid][valid_sigma]
        log_p_sigma <- -(2 * alpha_n + 1) * log(sv) - beta_n / sv^2
        log_jacobian <- log(tol / (2 * sigma_level)) - 2 * log(cv)
        vals[valid_sigma] <- exp(log_p_sigma + log_jacobian - log_Z_sigma)
      }
      result[valid] <- vals
      result
    })
  }

  if (metric %in% c("Cpu", "Cpl", "Cpk")) {
    return(.linear_metric_density(
      conjugate_log_joint,
      LSL, USL, metric, sigma_level,
      sigma_ref = sqrt(beta_n / alpha_n)
    ))
  }

  if (metric %in% c("Cpm", "Cpc")) {
    return(.radial_metric_density(
      conjugate_log_joint,
      LSL, USL, target, metric, sigma_level,
      sigma_ref = sqrt(beta_n / alpha_n),
      mu_ref = mu_n
    ))
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
  log_Z <- state$log_Z
  log_p_sigma <- state$log_p_sigma
  sigma_support <- state$sigma_support
  sigma_domain <- state$sigma_domain
  sigma_lower <- max(sigma_support$lower, sigma_domain$lower)
  sigma_upper <- min(sigma_support$upper, sigma_domain$upper)
  semi_mu_log_joint <- function(mu, sigma) {
    log_p_sigma(sigma) - log_Z +
      stats::dnorm(mu, mu_n, sigma / sqrt(k_n), log = TRUE)
  }

  if (metric %in% c("Cpk", "Cpu", "Cpl")) {
    return(.linear_metric_density(
      semi_mu_log_joint,
      LSL, USL, metric, sigma_level,
      sigma_ref = exp((state$t_lower + state$t_upper) / 2),
      sigma_lower = sigma_lower, sigma_upper = sigma_upper
    ))
  }

  # Cp: sigma is fixed on the contour, no integration needed
  if (metric == "Cp") {
    return(function(c) {
      c <- as.numeric(c)
      result <- numeric(length(c))
      valid <- is.finite(c) & c > 0
      if (!any(valid)) return(result)

      sigma_c <- tol / ((2 * sigma_level) * c[valid])
      valid_sigma <- sigma_c > 0 &
        sigma_c >= sigma_lower &
        sigma_c <= sigma_upper
      vals <- numeric(length(sigma_c))
      if (any(valid_sigma)) {
        sv <- sigma_c[valid_sigma]
        cv <- c[valid][valid_sigma]
        lps <- log_p_sigma(sv)
        log_jacobian <- log(tol / (2 * sigma_level)) - 2 * log(cv)
        vals[valid_sigma] <- exp(lps + log_jacobian - log_Z)
      }
      result[valid] <- vals
      result
    })
  }

  if (metric %in% c("Cpm", "Cpc")) {
    return(.radial_metric_density(
      semi_mu_log_joint,
      LSL, USL, target, metric, sigma_level,
      sigma_ref = exp((state$t_lower + state$t_upper) / 2),
      mu_ref = mu_n,
      sigma_lower = sigma_lower, sigma_upper = sigma_upper
    ))
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

  # The search above steps in wide units, which leaves a narrow posterior as a
  # tiny fraction of the domain that adaptive quadrature cannot find.
  tail_edge <- function(outer) {
    inner <- t_mode
    probe <- seq(outer, inner, length.out = 257L)
    keep <- which(log_base_t(probe) > h_max - tail_log_drop)
    if (length(keep) == 0L || keep[1] == 1L) return(outer)
    outer <- probe[keep[1] - 1L]
    inner <- probe[keep[1]]
    for (i in seq_len(30L)) {
      middle <- (outer + inner) / 2
      if (log_base_t(middle) > h_max - tail_log_drop) {
        inner <- middle
      } else {
        outer <- middle
      }
    }
    outer
  }
  if (t_lower < t_mode && t_mode < t_upper) {
    t_lower <- tail_edge(t_lower)
    t_upper <- tail_edge(t_upper)
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

  # Normalization constant of the mu marginal (sigma integrated analytically).
  log_Z <- state$log_Z

  # Constant for the unnormalized sigma kernel
  log_C <- log(2) - lgamma(alpha_n_eff)

  # Used only for Cpm and Cpc, which are supported on [0, Inf).
  survival_density <- function(metric_name) {
    .density_from_survival_fn(
      make_solver(
        data, LSL, USL, prior,
        metric = metric_name,
        target = target,
        cached_state = cached_state,
        sigma_level = sigma_level
      ),
      support_lower = 0
    )
  }

  if (metric %in% c("Cpu", "Cpl", "Cpk")) {
    mu_support <- .extract_prior_bounds(prior$bayestools_priors$mu)
    return(.linear_metric_density(
      function(mu, sigma) {
        beta_n <- beta0 + (sse + n * (mu - x_bar)^2) / 2
        prior$log_dens_mu(mu) + log_C -
          (2 * alpha_n_eff + 1) * log(sigma) - beta_n / sigma^2 - log_Z
      },
      LSL, USL, metric, sigma_level,
      sigma_ref = state$sigma_scale,
      mu_lower = mu_support$lower, mu_upper = mu_support$upper
    ))
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
      .semi_sigma_integrate_log_mu(log_integrand, state$mu_frame)$value
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


