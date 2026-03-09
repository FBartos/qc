rm(list = ls())
# ==============================================================================
# 1. SETUP: HELPER FUNCTIONS & PRIOR CLASSES
# ==============================================================================

log_diff_exp <- function(x, y) ifelse(x <= y, -Inf, x + log1p(-exp(y - x)))

create_prior_conjugate <- function(mu0=0, k0=0, alpha0=-0.5, beta0=0) {
  structure(list(mu0=mu0, k0=k0, alpha0=alpha0, beta0=beta0), class = "PriorConjugate")
}

create_prior_generic <- function(log_dens_fn) {
  structure(list(log_dens = log_dens_fn), class = "PriorGeneric")
}

get_metric_constraints <- function(metric, c, LSL, USL, target) {
  tol <- USL - LSL; mid <- (LSL + USL)/2; if(is.null(target)) target <- mid

  list(
    s_max_fn = function() {
      if (metric %in% c("Cp", "Cpk", "Cpm", "Cpc")) return(tol / (6 * c))
      return(Inf)
    },
    mu_b_fn = function(s) {
      if (metric=="Cp") return(c(-Inf, Inf))
      if (metric=="Cpk") return(c(LSL+3*c*s, USL-3*c*s))
      if (metric=="Cpu"||metric=="Cpu") return(c(-Inf, USL-3*c*s))
      if (metric=="Cpl"||metric=="Cpl") return(c(LSL+3*c*s, Inf))
      if (metric=="Cpm"||metric=="Cpc") {
        T_val <- if(metric=="Cpc") mid else target
        R <- tol/(6*c); if(s>=R) return(c(0,-1))
        w <- sqrt(R^2-s^2); return(c(T_val-w, T_val+w))
      }
      stop("Unknown Metric")
    }
  )
}

# ==============================================================================
# 2. THE SOLVER FACTORY (Returns a function P(Index > c))
# ==============================================================================

make_solver <- function(data, LSL, USL, prior, metric="Cpk", target=NULL, ...) {
  UseMethod("make_solver", prior)
}

make_solver.PriorConjugate <- function(data, LSL, USL, prior, metric="Cpk", target=NULL, ...) {
  n <- length(data); x_bar <- mean(data); SS <- sum((data - x_bar)^2)
  k_n <- prior$k0 + n; mu_n <- (prior$k0*prior$mu0 + n*x_bar)/k_n
  alpha_n <- prior$alpha0 + n/2
  beta_n <- prior$beta0 + 0.5*SS + (prior$k0*n*(x_bar-prior$mu0)^2)/(2*k_n)
  df_p <- 2 * alpha_n

  # Pre-compute global h_max_global ONCE (chi-square mode density, no truncation)
  # This is a valid upper bound for all c values since the log_diff_exp term <= 0
  y_mode <- max(df_p - 2, 1e-6)  # Mode of chi-square(df) is df-2 for df>2
  h_max_global <- dchisq(y_mode, df_p, log=TRUE)

  function(c) {
    if (c <= 0) return(1.0)
    constr <- get_metric_constraints(metric, c, LSL, USL, target)
    s_max <- constr$s_max_fn(); if (!is.infinite(s_max) && s_max <= 0) return(0.0)
    y_min <- if (is.infinite(s_max)) 0 else (2*beta_n)/(s_max^2)

    log_int <- function(y) {
      sapply(y, function(y_v) {
        sigma <- sqrt((2*beta_n)/y_v); sd_mu <- sigma/sqrt(k_n)
        mb <- constr$mu_b_fn(sigma); if (mb[1] >= mb[2]) return(-Inf)
        z_U <- if(is.infinite(mb[2])) Inf else (mb[2]-mu_n)/sd_mu
        z_L <- if(is.infinite(mb[1])) -Inf else (mb[1]-mu_n)/sd_mu
        log_diff_exp(pnorm(z_U,log.p=TRUE), pnorm(z_L,log.p=TRUE)) + dchisq(y_v, df_p, log=TRUE)
      })
    }

    # Use cached global h_max (valid upper bound for stable numerics)
    h_max <- h_max_global

    # Safe integrand: convert non-finite values to 0
    safe_integrand <- function(y) {
      vals <- exp(log_int(y) - h_max)
      vals[!is.finite(vals)] <- 0
      return(vals)
    }

    res <- integrate(safe_integrand, y_min, Inf)$value
    if (res <= 0) return(0.0)
    return(exp(h_max + log(res)))
  }
}

make_solver.PriorGeneric <- function(data, LSL, USL, prior, metric="Cpk", target=NULL,
                                     cached_state=NULL) {

  # Use cached state if available (massive speedup for multiple metrics)
  if (!is.null(cached_state)) {
    log_post <- cached_state$log_post
    h_max <- cached_state$h_max
    uni_s <- cached_state$uni_s
    uni_m <- cached_state$uni_m
    Z <- cached_state$Z
    int_2d <- cached_state$int_2d
  } else {
    # Compute from scratch
    n <- length(data); x_bar <- mean(data); sse <- sum((data - x_bar)^2)

    log_post <- function(mu, sigma) {
      if (sigma <= 0) return(-Inf)
      -n*log(sigma) - (sse+n*(mu-x_bar)^2)/(2*sigma^2) + prior$log_dens(mu, sigma)
    }

    # 1. MAP and Normalization (Computed ONCE)
    init_sd <- sqrt(sse/(n-1))
    opt <- optim(c(x_bar, init_sd), function(p) -log_post(p[1], p[2]))
    map_mu <- opt$par[1]; map_sig <- opt$par[2]; h_max <- -opt$value

    uni_s <- map_sig*10; uni_m <- function(s) c(map_mu-20*s, map_mu+20*s)

    int_2d <- function(s_lim, m_fn) {
      s_top <- if(is.infinite(s_lim)) uni_s else min(s_lim, uni_s)
      inner <- function(s_vec) {
        sapply(s_vec, function(s) {
          rb <- m_fn(s); ub <- uni_m(s); l <- max(rb[1],ub[1]); u <- min(rb[2],ub[2])
          if (l >= u) return(0)
          integrate(function(m) exp(log_post(m,s)-h_max), l, u, rel.tol=1e-4)$value
        })
      }
      integrate(inner, 0, s_top, rel.tol=1e-4)$value
    }

    Z <- int_2d(Inf, function(s) c(-Inf, Inf))
  }

  # 2. Return Closure
  function(c) {
    if (c <= 0) return(1.0)
    constr <- get_metric_constraints(metric, c, LSL, USL, target)
    num <- int_2d(constr$s_max_fn(), constr$mu_b_fn)
    return(num / Z)
  }
}

# ==============================================================================
# 3. ANALYSIS TOOL: DENSITY, STATS, PLOTS
# ==============================================================================

#' Compute Posterior Probability that a Capability Index is in a Region
#' @param data Numeric vector of observations
#' @param LSL Lower specification limit
#' @param USL Upper specification limit
#' @param bounds Vector c(lower, upper) defining the interval
#' @param prior Prior object (PriorConjugate or PriorGeneric)
#' @param metric Capability index: "Cp", "Cpk", "Cpm", "Cpc", "Cpu", "Cpl"
#' @param target Target value for Cpm (defaults to midpoint)
#' @param cached_state Pre-computed state from precompute_generic_state (for PriorGeneric only)
#' @return Probability that metric is in (bounds[1], bounds[2])
compute_cpk_prob <- function(data, LSL, USL, bounds, prior, metric="Cpk", target=NULL,
                              cached_state=NULL) {
  # Create solver function P(Index > c)
  if (inherits(prior, "PriorGeneric") && !is.null(cached_state)) {
    S <- make_solver(data, LSL, USL, prior, metric, target, cached_state=cached_state)
  } else {
    S <- make_solver(data, LSL, USL, prior, metric, target)
  }

  # P(lower < Index < upper) = P(Index > lower) - P(Index > upper)
  p_lower <- S(min(bounds))
  p_upper <- S(max(bounds))

  return(p_lower - p_upper)
}

#' Precompute expensive posterior state for Generic priors (speedup for multiple metrics)
#' @return cached_state object to pass to make_solver
precompute_generic_state <- function(data, prior) {
  n <- length(data); x_bar <- mean(data); sse <- sum((data - x_bar)^2)

  log_post <- function(mu, sigma) {
    if (sigma <= 0) return(-Inf)
    -n*log(sigma) - (sse+n*(mu-x_bar)^2)/(2*sigma^2) + prior$log_dens(mu, sigma)
  }

  init_sd <- sqrt(sse/(n-1))
  opt <- optim(c(x_bar, init_sd), function(p) -log_post(p[1], p[2]))
  map_mu <- opt$par[1]; map_sig <- opt$par[2]; h_max <- -opt$value

  uni_s <- map_sig*10
  uni_m <- function(s) c(map_mu-20*s, map_mu+20*s)

  int_2d <- function(s_lim, m_fn) {
    s_top <- if(is.infinite(s_lim)) uni_s else min(s_lim, uni_s)
    inner <- function(s_vec) {
      sapply(s_vec, function(s) {
        rb <- m_fn(s); ub <- uni_m(s); l <- max(rb[1],ub[1]); u <- min(rb[2],ub[2])
        if (l >= u) return(0)
        integrate(function(m) exp(log_post(m,s)-h_max), l, u, rel.tol=1e-4)$value
      })
    }
    integrate(inner, 0, s_top, rel.tol=1e-4)$value
  }

  Z <- int_2d(Inf, function(s) c(-Inf, Inf))

  list(log_post=log_post, h_max=h_max, uni_s=uni_s, uni_m=uni_m, int_2d=int_2d, Z=Z)
}

#' Analyze Capability Posterior
#' @param range vector c(min, max) to scan for density
#' @param n_grid number of points to evaluate
#' @param cached_state Pre-computed state from precompute_generic_state (for PriorGeneric only)
# analyze_capability <- function(data, LSL, USL, prior, metric="Cpk", target=NULL,
#                                range=c(0, 3), n_grid=200, cached_state=NULL) {

#   # 1. Create Solver
#   if (inherits(prior, "PriorGeneric")) {
#     solve_survival <- make_solver(data, LSL, USL, prior, metric, target, cached_state)
#   } else {
#     solve_survival <- make_solver(data, LSL, USL, prior, metric, target)
#   }

#   # 2. Evaluate Grid (Survival Function)
#   grid_x <- seq(range[1], range[2], length.out = n_grid)
#   # Computes P(Index > x)
#   S_vals <- sapply(grid_x, solve_survival)

#   # 3. Compute PDF via Finite Differences
#   # PDF(x) = - d/dx S(x)
#   pdf_vals <- -diff(S_vals) / diff(grid_x)
#   mid_x    <- (grid_x[-1] + grid_x[-n_grid]) / 2

#   # Normalize numerically to ensure area = 1 (removes minor discretization error)
#   area <- sum(pdf_vals * diff(grid_x))
#   pdf_vals <- pdf_vals / area

#   # 4. Compute Statistics from Grid PDF
#   post_mean <- sum(mid_x * pdf_vals * diff(grid_x))

#   # Standard Deviation: sqrt( E[x^2] - E[x]^2 )
#   post_var  <- sum((mid_x^2) * pdf_vals * diff(grid_x)) - post_mean^2
#   post_sd   <- sqrt(max(0, post_var))

#   # Quantiles (Inverse CDF)
#   # Reconstruct CDF from grid
#   cdf_vals <- cumsum(pdf_vals * diff(grid_x))
#   get_quant <- function(q) mid_x[which.min(abs(cdf_vals - q))]

#   q_low  <- get_quant(0.025)
#   q_high <- get_quant(0.975)

#   # HDI (Highest Density Interval)
#   # Simple grid search for 95% mass with shortest width
#   sorted_idx <- order(pdf_vals, decreasing = TRUE)
#   sorted_mass <- pdf_vals[sorted_idx] * diff(grid_x)[1] # assume uniform grid step
#   cum_mass <- cumsum(sorted_mass)
#   cutoff_idx <- which(cum_mass >= 0.95)[1]

#   # The indices belonging to HDI
#   hdi_indices <- sorted_idx[1:cutoff_idx]
#   hdi_x <- mid_x[hdi_indices]

#   return(list(
#     metric = metric,
#     grid = data.frame(x = mid_x, density = pdf_vals),
#     stats = c(Mean = post_mean, SD = post_sd,
#               Q2.5 = q_low, Q97.5 = q_high,
#               HDI_Lo = min(hdi_x), HDI_Hi = max(hdi_x))
#   ))
# }

# ==============================================================================
# 3. ANALYSIS TOOL: AUTO-RANGING DENSITY & STATS
# ==============================================================================

#' Analyze Capability with Automatic Grid Detection
#' @param alpha_tail Probability mass to leave in the tails for grid detection (default 0.0001)
#' @param cached_state Pre-computed state from precompute_generic_state (for PriorGeneric only)
analyze_capability <- function(data, LSL, USL, prior, metric="Cpk", target=NULL,
                               n_grid=128, alpha_tail=0.0001, cached_state=NULL) {

  # 1. Create Solver Function: P(Index > c)
  # This function is monotonic decreasing, which makes it perfect for root finding.
  if (inherits(prior, "PriorGeneric") && !is.null(cached_state)) {
    S <- make_solver(data, LSL, USL, prior, metric, target, cached_state=cached_state)
  } else {
    S <- make_solver(data, LSL, USL, prior, metric, target)
  }

  # 2. Heuristic Search for Bracket Interval [0, max_c]
  # We need to find a rough upper limit where Prob effectively drops to 0
  # to give 'uniroot' a valid bracket.
  max_c <- 3.0
  while (S(max_c) > alpha_tail) {
    max_c <- max_c * 2
    if (max_c > 100) break # Safety brake for degenerate priors
  }

  # 3. Find Grid Bounds using Root Finding (Inverse CDF)
  # We want the grid to cover the range from quantile(alpha) to quantile(1-alpha)

  # Find Lower Bound: c where P(Index > c) approx 1 - alpha
  # i.e., S(c) - (1 - alpha_tail) = 0
  get_quantile <- function(target_prob) {
    tryCatch({
      uniroot(function(c) S(c) - target_prob,
              interval = c(0, max_c),
              extendInt = "downX", # allow searching near 0
              tol = 1e-4)$root
    }, error = function(e) NA)
  }

  # Note: S(c) is "Prob Greater", so:
  # Lower limit of grid (0.1% quantile) corresponds to S(c) = 0.999
  # Upper limit of grid (99.9% quantile) corresponds to S(c) = 0.001

  x_start <- get_quantile(1 - alpha_tail)
  x_end   <- get_quantile(alpha_tail)

  if (is.na(x_start)) x_start <- 0
  if (is.na(x_end))   x_end <- max_c

  # 4. Evaluate Grid
  grid_x <- seq(x_start, x_end, length.out = n_grid)
  S_vals <- sapply(grid_x, S)

  # 5. Compute PDF via Finite Differences
  # PDF(x) = - d/dx S(x)
  # We use centered differencing for the interior points for slightly better accuracy
  pdf_vals <- -diff(S_vals) / diff(grid_x)
  mid_x    <- (grid_x[-1] + grid_x[-n_grid]) / 2

  # Normalize area to 1.0
  area <- sum(pdf_vals * diff(grid_x))
  pdf_vals <- pdf_vals / area

  # 6. Compute Statistics
  # Mean
  post_mean <- sum(mid_x * pdf_vals * diff(grid_x))

  # Variance / SD
  post_var <- sum((mid_x^2) * pdf_vals * diff(grid_x)) - post_mean^2
  post_sd  <- sqrt(max(0, post_var))

  # Quantiles (Reconstruct CDF from grid)
  cdf_vals <- cumsum(pdf_vals * diff(grid_x))
  get_q <- function(q) mid_x[which.min(abs(cdf_vals - q))]

  q2.5  <- get_q(0.025)
  q97.5 <- get_q(0.975)
  median_val <- get_q(0.5)

  # HDI (Highest Density Interval)
  sorted_idx <- order(pdf_vals, decreasing = TRUE)
  # Mass of each bar is height * width
  sorted_mass <- pdf_vals[sorted_idx] * diff(grid_x)[1]
  cum_mass <- cumsum(sorted_mass)
  cutoff_idx <- which(cum_mass >= 0.95)[1]
  hdi_indices <- sorted_idx[1:cutoff_idx]

  return(list(
    metric = metric,
    grid = data.frame(x = mid_x, density = pdf_vals),
    stats = c(Mean = post_mean, Median = median_val, SD = post_sd,
              Q2.5 = q2.5, Q97.5 = q97.5,
              HDI_Lo = min(mid_x[hdi_indices]),
              HDI_Hi = max(mid_x[hdi_indices]))
  ))
}

# ==============================================================================
# DEMO WITH AUTO-RANGING
# ==============================================================================

set.seed(2025)
data_vec <- rnorm(30, mean = 50.0, sd = 0.5) # Slightly off-center
LSL <- 44; USL <- 56
target  <- (LSL + USL) / 2
bounds <- c(1.0, 2.0)

# Jeffreys Prior
prior_conj <- create_prior_conjugate()
prior_gen  <- create_prior_generic(function(m, s) -log(s)) # Jeffreys manually

# something random
prior_gen2 <- create_prior_generic(function(m, s) dnorm(m, mean = 0, sd = 10, log = TRUE) + dexp(s, rate = 1, log = TRUE))

fit_mcmc <- qc::bpc(x = data_vec, LSL = LSL, USL = USL, target = target, iter = 100000, chains = 4)
result_mcmc <- sapply(fit_mcmc$metrics, function(metric_values) mean(bounds[1] < metric_values & metric_values < bounds[2]))

# Compare Indices
metrics <- c("Cp", "Cpk", "Cpm", "Cpu", "Cpl", "Cpc")

# FAST: Precompute expensive state once, reuse for all metrics
cached_conj <- NULL  # Conjugate doesn't need caching (already fast)
cached_gen <- precompute_generic_state(data_vec, prior_gen)
cached_gen2 <- precompute_generic_state(data_vec, prior_gen2)

system.time({
  qc::bpc(x = data_vec, LSL = LSL, USL = USL, target = target, iter = 100000, chains = 4)
})
system.time({
  res_prior_conj <- sapply(metrics, function(m) {
    analyze_capability(data_vec, LSL, USL, prior_conj, metric=m, target=target,
                      cached_state=cached_conj)
  })
})

res_prior_conj <- sapply(metrics, function(m) {
  analyze_capability(data_vec, LSL, USL, prior_conj, metric=m, target=target,
                    cached_state=cached_conj)
})
res_prior_gen <- sapply(metrics, function(m) {
  analyze_capability(data_vec, LSL, USL, prior_gen, metric=m, target=target,
                     cached_state=cached_gen)
})
res_prior_gen2 <- sapply(metrics, function(m) {
  analyze_capability(data_vec, LSL, USL, prior_gen2, metric=m, target=target,
                     cached_state=cached_gen2)
})

dd <- density(fit_mcmc$metrics$Cpk)
plot(res_prior_gen[["grid", "Cpk"]]$x, res_prior_gen[["grid", "Cpk"]]$density)
lines(dd$x, dd$y, col="red", lwd=2)
lines(res_prior_conj[["grid", "Cpk"]]$x, res_prior_conj[["grid", "Cpk"]]$density, col="blue", lwd=2)

sss <- summary(fit_mcmc)
res_prior_conj[["stats", "Cpk"]]
res_prior_gen[["stats", "Cpk"]]
sss$summary[sss$summary$metric == "Cpk", ]

# ==============================================================================
# INTERVAL PROBABILITY COMPARISON: Analytical vs MCMC
# ==============================================================================

cat("\n==============================================================================\n")
cat("INTERVAL PROBABILITY COMPARISON: P(", bounds[1], "< Index <", bounds[2], ")\n")
cat("==============================================================================\n\n")

# Compute analytical probabilities for all metrics
results <- data.frame(
  Metric = metrics,
  Conj_Prob = sapply(metrics, function(m) {
    compute_cpk_prob(data_vec, LSL, USL, bounds, prior_conj, metric=m, target=target)
  }),
  Gen_Prob = sapply(metrics, function(m) {
    compute_cpk_prob(data_vec, LSL, USL, bounds, prior_gen, metric=m, target=target,
                     cached_state=cached_gen)
  }),
  MCMC_Prob = result_mcmc[metrics]
)

print(results, digits=6)

cat("\nMax absolute difference (Conjugate vs MCMC):",
    max(abs(results$Conj_Prob - results$MCMC_Prob)), "\n")
cat("Max absolute difference (Generic vs MCMC):",
    max(abs(results$Gen_Prob - results$MCMC_Prob)), "\n")


sss <- summary(fit_mcmc)
sss$interval_summary
intervalSummary <- sss$interval_summary

