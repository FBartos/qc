# ==============================================================================
# Prior Conversion and Divergence Helpers
# ==============================================================================

#' Extract the effective shape parameter alpha controlling tail behavior at sigma -> 0
#'
#' For a prior pi(sigma), alpha characterizes the density near zero:
#' pi(sigma) ~ sigma^(alpha - 1) as sigma -> 0.
#' Finite alpha indicates a potential divergence of \eqn{E[\sigma^{-k}]} for \eqn{k \ge \alpha}.
#' alpha = Inf means the prior is bounded away from zero or decays superexponentially.
#'
#' @param sigma_prior BayesTools prior object or string ("Jeffreys_sigma")
#' @return Numeric scalar: the effective alpha (possibly Inf)
#' @keywords internal
.extract_alpha_parameter <- function(sigma_prior) {
  if (identical(sigma_prior, "Jeffreys_sigma"))
    return(0)

  if (inherits(sigma_prior, "PriorConjugate"))
    return(sigma_prior$alpha0)

  if (!inherits(sigma_prior, "prior"))
    stop("sigma_prior must be a string or BayesTools::prior object")

  dist   <- sigma_prior[["distribution"]]
  params <- sigma_prior[["parameters"]]
  trunc  <- sigma_prior[["truncation"]]
  lower  <- trunc[["lower"]] %||% -Inf
  upper  <- trunc[["upper"]] %||%  Inf

  # Truncation bounded away from zero overrides everything
  if (is.finite(lower) && lower > 0)
    return(Inf)

  alpha <- switch(dist,
    "gamma"    = params[["shape"]],
    "exp"      = 1,
    "invgamma" = Inf,
    "lognormal"= Inf,
    # A normal prior on sigma that reaches 0 has finite positive density there,
    # so its local order matches alpha = 1 unless truncation bounds exclude 0.
    "normal"   = 1,
    "t"        = 1,
    "cauchy"   = 1,
    "point"    = Inf,
    "uniform"  = if (params[["a"]] <= 0) 1 else Inf,
    "beta"     = if (params[["alpha"]] <= 0 || (is.finite(lower) && lower <= 0)) params[["alpha"]] else Inf,
    Inf  # conservative default for unknown families
  )

  alpha
}

.is_point_prior <- function(prior) {
  inherits(prior, "prior") && identical(prior[["distribution"]], "point")
}

.abort_point_prior_integration <- function(parameter) {
  stop(
    "Point priors on `", parameter,
    "` are discrete and are not supported by the integration backend. ",
    "Use `method = 'mcmc'` for fixed-value priors.",
    call. = FALSE
  )
}

#' Check whether posterior moments of a capability metric diverge
#'
#' Uses the analytic decision rules from the divergence analysis: for metrics
#' scaling as \eqn{\sigma^{-1}} (Cp, Cpu, Cpl, Cpk), \eqn{E[C^k] < \infty} iff \eqn{\alpha > k}.
#' For metrics involving \eqn{\sqrt{\sigma^2 + (\mu - T)^2}} (Cpm, Cpc), the singularity
#' is regularised and \eqn{E[C^k] < \infty} iff \eqn{\alpha > k - 1}.
#'
#' @param metric One of "Cp", "Cpu", "Cpl", "Cpk", "Cpm", "Cpc"
#' @param alpha_sigma Effective alpha from .extract_alpha_parameter()
#' @param sigma_prior_label Human-readable label for the sigma prior (for messages)
#' @return List with mean_divergent, sd_divergent, alpha, reason
#' @keywords internal
.check_moment_divergence <- function(metric, alpha_sigma,
                                     sigma_prior_label = "sigma prior") {
  if (metric %in% c("Cp", "Cpu", "Cpl", "Cpk")) {
    mean_threshold <- 1
    var_threshold  <- 2
  } else if (metric %in% c("Cpm", "Cpc")) {
    mean_threshold <- 0
    var_threshold  <- 1
  } else {
    return(list(mean_divergent = FALSE, sd_divergent = FALSE,
                alpha = alpha_sigma, reason = NULL))
  }

  mean_div <- is.finite(alpha_sigma) && alpha_sigma <= mean_threshold
  sd_div   <- is.finite(alpha_sigma) && alpha_sigma <= var_threshold

  reason <- NULL
  if (mean_div) {
    reason <- sprintf(
      "%s has alpha=%.3g; %s requires alpha>%g for a finite mean (and alpha>%g for finite variance)",
      sigma_prior_label, alpha_sigma, metric, mean_threshold, var_threshold
    )
  } else if (sd_div) {
    reason <- sprintf(
      "%s has alpha=%.3g; %s requires alpha>%g for finite variance",
      sigma_prior_label, alpha_sigma, metric, var_threshold
    )
  }

  list(
    mean_divergent = mean_div,
    sd_divergent   = sd_div,
    alpha          = alpha_sigma,
    reason         = reason
  )
}

#' Build a human-readable label for a BayesTools prior
#' @param sigma_prior BayesTools prior object or string
#' @return Character string
#' @keywords internal
.sigma_prior_label <- function(sigma_prior) {
  if (identical(sigma_prior, "Jeffreys_sigma"))
    return("Jeffreys(sigma)")
  if (!inherits(sigma_prior, "prior"))
    return("unknown prior")

  dist   <- sigma_prior[["distribution"]]
  params <- sigma_prior[["parameters"]]
  param_str <- paste(vapply(params, function(p) format(p, digits = 3), character(1)), collapse = ", ")
  sprintf("%s(%s)", dist, param_str)
}

#' Convert parameter-wise BayesTools priors to integration prior format
#' @param mu_prior Prior for mu (string or BayesTools prior)
#' @param sigma_prior Prior for sigma (string or BayesTools prior)
#' @return List with $prior (integration prior object), $case (1-4), and $is_conjugate (logical)
#' @keywords internal

#' @noRd
.classify_mu_prior <- function(prior) {
  # Jeffreys (flat) is conjugate with k0 = 0
  if (identical(prior, "Jeffreys_mu")) {
    return(list(is_conjugate = TRUE, mu0 = 0, k0 = 0))
  }

  # PriorConjugate carries the full NIG hyperparameters
  if (inherits(prior, "PriorConjugate")) {
    return(list(is_conjugate = TRUE, mu0 = prior$mu0, k0 = prior$k0))
  }

  # BayesTools Normal(mean, sd) prior on mu is sigma-INDEPENDENT: N(mean, sd^2).
  # The conjugate NIG model requires N(mu0, sigma^2/k0), which is sigma-DEPENDENT.
  # These are different models, so Normal mu priors are NOT conjugate here.

  # Not conjugate - build log-density function
  list(is_conjugate = FALSE, log_dens_fn = .make_prior_log_dens_fn(prior))
}

#' @noRd
.classify_sigma_prior <- function(prior) {
  # Jeffreys (1/sigma) is conjugate with alpha0 = -0.5, beta0 = 0
  if (identical(prior, "Jeffreys_sigma")) {
    return(list(is_conjugate = TRUE, alpha0 = -0.5, beta0 = 0))
  }

  # PriorConjugate carries the full NIG hyperparameters
  if (inherits(prior, "PriorConjugate")) {
    return(list(is_conjugate = TRUE, alpha0 = prior$alpha0, beta0 = prior$beta0))
  }

  # InvGamma on sigma is NOT conjugate with the Normal likelihood because
  # BayesTools places InvGamma on sigma (not sigma^2). The conjugate model
  # requires InvGamma on sigma^2, which has a different functional form.
  # Gamma on sigma is also NOT conjugate.
  # All other distributions are non-conjugate.
  list(is_conjugate = FALSE, log_dens_fn = .make_prior_log_dens_fn(prior))
}

.bayestools_to_integration_prior <- function(mu_prior, sigma_prior,
                                             mu_label = "prior$mu",
                                             sigma_label = "prior$sigma") {
  if (.is_point_prior(mu_prior)) {
    .abort_point_prior_integration(mu_label)
  }
  if (.is_point_prior(sigma_prior)) {
    .abort_point_prior_integration(sigma_label)
  }

  # Classify each prior
  mu_info <- .classify_mu_prior(mu_prior)
  sigma_info <- .classify_sigma_prior(sigma_prior)
  bayestools_priors <- list(mu = mu_prior, sigma = sigma_prior)

  # Case 1: Full Conjugate (Normal/Jeffreys on mu AND InvGamma/Jeffreys on sigma)
  if (mu_info$is_conjugate && sigma_info$is_conjugate) {
    return(list(
      prior = create_prior_conjugate(mu_info$mu0, mu_info$k0,
                                      sigma_info$alpha0, sigma_info$beta0),
      case = 1L,
      is_conjugate = TRUE
    ))
  }

  # Case 2: Semi-Conjugate Mu (conjugate mu, non-conjugate sigma)
  if (mu_info$is_conjugate && !sigma_info$is_conjugate) {
    return(list(
      prior = create_prior_semi_mu(mu_info$mu0, mu_info$k0,
                                    sigma_info$log_dens_fn, bayestools_priors),
      case = 2L,
      is_conjugate = FALSE
    ))
  }

  # Case 3: Semi-Conjugate Sigma (non-conjugate mu, conjugate sigma)
  if (!mu_info$is_conjugate && sigma_info$is_conjugate) {
    return(list(
      prior = create_prior_semi_sigma(sigma_info$alpha0, sigma_info$beta0,
                                       mu_info$log_dens_fn, bayestools_priors),
      case = 3L,
      is_conjugate = FALSE
    ))
  }

  # Case 4: Generic (both non-conjugate)
  log_dens_fn <- function(mu, sigma) {
    mu_info$log_dens_fn(mu) + sigma_info$log_dens_fn(sigma)
  }
  list(
    prior = create_prior_generic(log_dens_fn, bayestools_priors = bayestools_priors),
    case = 4L,
    is_conjugate = FALSE
  )
}

#' Extract a reasonable starting value and scale from a BayesTools prior
#' @keywords internal
.extract_prior_init <- function(bt_prior) {
  if (is.character(bt_prior)) {
    return(list(value = 1, scale = 10))
  }
  if (!inherits(bt_prior, "prior")) return(list(value = 1, scale = 1))

  dist   <- bt_prior[["distribution"]]
  params <- bt_prior[["parameters"]]
  trunc  <- bt_prior[["truncation"]]
  lower  <- trunc[["lower"]] %||% -Inf
  upper  <- trunc[["upper"]] %||%  Inf

  info <- switch(dist,
    "normal"    = list(value = params[["mean"]], scale = params[["sd"]]),
    "t"         = list(value = params[["location"]], scale = params[["scale"]]),
    "uniform"   = list(value = (params[["a"]] + params[["b"]]) / 2,
                       scale = (params[["b"]] - params[["a"]]) / 4),
    "gamma"     = { sh <- params[["shape"]]; rt <- params[["rate"]]
                    list(value = max((sh - 1) / rt, 0.5 / rt),
                         scale = sqrt(sh) / rt) },
    "invgamma"  = { a <- params[["shape"]]; b <- params[["scale"]]
                    list(value = b / (a + 1), scale = b / a) },
    "lognormal" = { ml <- params[["meanlog"]]; sl <- params[["sdlog"]]
                    list(value = exp(ml - sl^2),
                         scale = exp(ml) * sqrt(exp(sl^2) - 1)) },
    "exp"       = { rt <- params[["rate"]]
                    list(value = 1 / rt, scale = 1 / rt) },
    list(value = 1, scale = 1)
  )

  # Clamp value to truncation bounds
  if (is.finite(lower) && info$value < lower) {
    span <- if (is.finite(upper)) upper - lower else info$scale
    info$value <- lower + 0.1 * min(info$scale, span)
  }
  if (is.finite(upper) && info$value > upper) {
    span <- if (is.finite(lower)) upper - lower else info$scale
    info$value <- upper - 0.1 * min(info$scale, span)
  }
  info
}

#' Extract effective support bounds for a BayesTools prior
#' @keywords internal
.extract_prior_bounds <- function(bt_prior) {
  if (!inherits(bt_prior, "prior")) {
    return(list(lower = -Inf, upper = Inf))
  }

  dist <- bt_prior[["distribution"]]
  params <- bt_prior[["parameters"]]
  trunc <- bt_prior[["truncation"]]

  lower <- trunc[["lower"]] %||% -Inf
  upper <- trunc[["upper"]] %||% Inf

  natural_bounds <- switch(dist,
    "uniform" = list(lower = params[["a"]], upper = params[["b"]]),
    "gamma" = list(lower = 0, upper = Inf),
    "invgamma" = list(lower = 0, upper = Inf),
    "lognormal" = list(lower = 0, upper = Inf),
    "exp" = list(lower = 0, upper = Inf),
    "beta" = list(lower = 0, upper = 1),
    "point" = list(lower = params[["location"]], upper = params[["location"]]),
    list(lower = -Inf, upper = Inf)
  )

  list(
    lower = max(lower, natural_bounds$lower),
    upper = min(upper, natural_bounds$upper)
  )
}

#' Find finite optimizer initials for generic integration problems
#' @keywords internal
.find_feasible_generic_init <- function(init_mu, init_sd, prior, log_post) {
  candidate_pairs <- list(c(init_mu, max(init_sd, 1e-6)))

  bt <- prior$bayestools_priors
  if (!is.null(bt)) {
    mu_init  <- .extract_prior_init(bt$mu)
    sig_init <- .extract_prior_init(bt$sigma)

    mu_candidates <- unique(c(init_mu, mu_init$value, 0))
    sig_candidates <- unique(c(init_sd, sig_init$value, 0.1, 1))

    if (inherits(bt$sigma, "prior")) {
      sig_q <- tryCatch(
        as.numeric(BayesTools::quant(bt$sigma, c(0.25, 0.5, 0.75))),
        error = function(e) numeric(0)
      )
      sig_candidates <- c(sig_candidates, sig_q)
    }

    mu_candidates <- mu_candidates[is.finite(mu_candidates)]
    sig_candidates <- sig_candidates[is.finite(sig_candidates) & sig_candidates > 0]

    for (sigma in sig_candidates) {
      for (mu in mu_candidates) {
        candidate_pairs[[length(candidate_pairs) + 1L]] <- c(mu, sigma)
      }
    }
  }

  candidate_pairs <- c(candidate_pairs, list(c(0, max(init_sd, 1e-6)), c(0, 1)))

  for (candidate in candidate_pairs) {
    value <- tryCatch(log_post(candidate[1], candidate[2]), error = function(e) -Inf)
    if (is.finite(value)) {
      return(list(mu = candidate[1], sigma = candidate[2]))
    }
  }

  stop(
    "Could not find finite initial values for integration with the supplied priors. ",
    "Consider widening the prior support or using method = 'mcmc'."
  )
}

#' Create a log-density function for a BayesTools prior
#' @param prior Prior specification (string or BayesTools prior object)
#' @return Function(x) returning vectorized log-density
#' @keywords internal
.make_prior_log_dens_fn <- function(prior) {

  # Handle string priors (Jeffreys)
  if (is.character(prior)) {
    if (prior == "Jeffreys_mu") {
      return(function(x) rep(0, length(x)))  # Improper flat prior
    } else if (prior == "Jeffreys_sigma") {
      return(function(x) ifelse(x <= 0, -Inf, -log(x)))
    }
    stop("Unknown string prior: ", prior)
  }

  # Handle BayesTools prior objects
  if (!inherits(prior, "prior")) {
    stop("prior must be a string or BayesTools::prior object")
  }

  dist <- prior[["distribution"]]
  params <- prior[["parameters"]]
  trunc <- prior[["truncation"]]

  # Check truncation bounds
  lower <- if (!is.null(trunc[["lower"]])) trunc[["lower"]] else -Inf
  upper <- if (!is.null(trunc[["upper"]])) trunc[["upper"]] else Inf

  # Pre-compute normalization constant for truncation
  log_norm <- 0
  if (is.finite(lower) || is.finite(upper)) {
    log_norm <- .compute_truncation_norm(dist, params, lower, upper)
    if (!is.finite(log_norm)) {
      return(function(x) rep(-Inf, length(x)))
    }
  }


  # Construct efficient closure based on distribution
  base_dens_fn <- switch(dist,
    "point" = {
      loc <- params[["location"]]
      function(x) ifelse(x == loc, 0, -Inf)
    },
    "normal" = {
      mean <- params[["mean"]]
      sd <- params[["sd"]]
      function(x) stats::dnorm(x, mean = mean, sd = sd, log = TRUE)
    },
    "lognormal" = {
      meanlog <- params[["meanlog"]]
      sdlog <- params[["sdlog"]]
      function(x) ifelse(x <= 0, -Inf, stats::dlnorm(x, meanlog = meanlog, sdlog = sdlog, log = TRUE))
    },
    "t" = {
      loc <- params[["location"]]
      scale <- params[["scale"]]
      df <- params[["df"]]
      log_scale <- log(scale)
      function(x) {
        z <- (x - loc) / scale
        stats::dt(z, df = df, log = TRUE) - log_scale
      }
    },
    "gamma" = {
      shape <- params[["shape"]]
      rate <- params[["rate"]]
      # function(x) ifelse(x <= 0, -Inf, stats::dgamma(x, shape = shape, rate = rate, log = TRUE))
      function(x) stats::dgamma(x, shape = shape, rate = rate, log = TRUE)
    },
    "invgamma" = {
      alpha <- params[["shape"]]
      beta <- params[["scale"]]
      log_beta <- log(beta)
      lgamma_alpha <- lgamma(alpha)
      function(x) {
        ifelse(x <= 0, -Inf,
               alpha * log_beta - lgamma_alpha - (alpha + 1) * log(x) - beta / x)
      }
    },
    "uniform" = {
      a <- params[["a"]]
      b <- params[["b"]]
      log_width <- log(b - a)
      function(x) ifelse(x >= a & x <= b, -log_width, -Inf)
    },
    "beta" = {
      alpha <- params[["alpha"]]
      beta <- params[["beta"]]
      function(x) ifelse(x < 0 | x > 1, -Inf, stats::dbeta(x, shape1 = alpha, shape2 = beta, log = TRUE))
    },
    "exp" = {
      rate <- params[["rate"]]
      function(x) ifelse(x < 0, -Inf, stats::dexp(x, rate = rate, log = TRUE))
    },
    stop("Unsupported prior distribution: ", dist)
  )

  # Return function that handles truncation and normalization
  function(x) {
    # Quick bounds check
    in_bounds <- x >= lower & x <= upper

    # If all valid, fast path
    if (all(in_bounds)) {
      return(base_dens_fn(x) - log_norm)
    }

    # If none valid, fast return
    if (!any(in_bounds)) {
      return(rep(-Inf, length(x)))
    }

    # Mixed case
    res <- rep(-Inf, length(x))
    res[in_bounds] <- base_dens_fn(x[in_bounds]) - log_norm
    res
  }
}


#' Compute log normalizing constant for truncated distribution
#' @keywords internal
.truncation_log_cdf <- function(dist, params, x) {
  switch(dist,
    "point" = {
      loc <- params[["location"]]
      ifelse(x < loc, -Inf, 0)
    },
    "normal" = stats::pnorm(x, mean = params[["mean"]], sd = params[["sd"]], log.p = TRUE),
    "lognormal" = stats::plnorm(
      x, meanlog = params[["meanlog"]], sdlog = params[["sdlog"]], log.p = TRUE
    ),
    "t" = stats::pt(
      (x - params[["location"]]) / params[["scale"]],
      df = params[["df"]],
      log.p = TRUE
    ),
    "gamma" = stats::pgamma(
      x, shape = params[["shape"]], rate = params[["rate"]], log.p = TRUE
    ),
    "invgamma" = {
      alpha <- params[["shape"]]
      beta <- params[["scale"]]
      ifelse(
        x <= 0,
        -Inf,
        stats::pgamma(beta / x, shape = alpha, lower.tail = FALSE, log.p = TRUE)
      )
    },
    "uniform" = {
      a <- params[["a"]]
      b <- params[["b"]]
      log_width <- log(b - a)
      res <- rep(-Inf, length(x))
      res[x >= b] <- 0
      inside <- x > a & x < b
      res[inside] <- log(x[inside] - a) - log_width
      res
    },
    "beta" = stats::pbeta(
      x, shape1 = params[["alpha"]], shape2 = params[["beta"]], log.p = TRUE
    ),
    "exp" = stats::pexp(x, rate = params[["rate"]], log.p = TRUE),
    stop("Unsupported prior distribution: ", dist)
  )
}

.truncation_log_ccdf <- function(dist, params, x) {
  switch(dist,
    "point" = {
      loc <- params[["location"]]
      ifelse(x < loc, 0, -Inf)
    },
    "normal" = stats::pnorm(
      x, mean = params[["mean"]], sd = params[["sd"]],
      lower.tail = FALSE, log.p = TRUE
    ),
    "lognormal" = stats::plnorm(
      x, meanlog = params[["meanlog"]], sdlog = params[["sdlog"]],
      lower.tail = FALSE, log.p = TRUE
    ),
    "t" = stats::pt(
      (x - params[["location"]]) / params[["scale"]],
      df = params[["df"]],
      lower.tail = FALSE,
      log.p = TRUE
    ),
    "gamma" = stats::pgamma(
      x, shape = params[["shape"]], rate = params[["rate"]],
      lower.tail = FALSE, log.p = TRUE
    ),
    "invgamma" = {
      alpha <- params[["shape"]]
      beta <- params[["scale"]]
      ifelse(
        x <= 0,
        0,
        stats::pgamma(beta / x, shape = alpha, log.p = TRUE)
      )
    },
    "uniform" = {
      a <- params[["a"]]
      b <- params[["b"]]
      log_width <- log(b - a)
      res <- rep(-Inf, length(x))
      res[x <= a] <- 0
      inside <- x > a & x < b
      res[inside] <- log(b - x[inside]) - log_width
      res
    },
    "beta" = stats::pbeta(
      x, shape1 = params[["alpha"]], shape2 = params[["beta"]],
      lower.tail = FALSE, log.p = TRUE
    ),
    "exp" = stats::pexp(
      x, rate = params[["rate"]],
      lower.tail = FALSE, log.p = TRUE
    ),
    stop("Unsupported prior distribution: ", dist)
  )
}

.compute_truncation_norm <- function(dist, params, lower, upper) {
  if (!is.finite(lower) && !is.finite(upper)) {
    return(0)
  }
  if (dist == "point") {
    loc <- params[["location"]]
    return(if (loc >= lower && loc <= upper) 0 else -Inf)
  }
  if (lower >= upper) {
    return(-Inf)
  }

  log_upper <- if (!is.finite(upper)) {
    if (upper > 0) 0 else -Inf
  } else {
    .truncation_log_cdf(dist, params, upper)
  }
  log_lower <- if (!is.finite(lower)) {
    if (lower < 0) -Inf else 0
  } else {
    .truncation_log_cdf(dist, params, lower)
  }
  log_mass_lower <- log_diff_exp(log_upper, log_lower)

  log_survival_lower <- if (!is.finite(lower)) {
    if (lower < 0) 0 else -Inf
  } else {
    .truncation_log_ccdf(dist, params, lower)
  }
  log_survival_upper <- if (!is.finite(upper)) {
    if (upper > 0) -Inf else 0
  } else {
    .truncation_log_ccdf(dist, params, upper)
  }
  log_mass_upper <- log_diff_exp(log_survival_lower, log_survival_upper)

  if (is.finite(log_mass_lower) && is.finite(log_mass_upper)) {
    return(max(log_mass_lower, log_mass_upper))
  }
  if (is.finite(log_mass_lower)) {
    return(log_mass_lower)
  }
  if (is.finite(log_mass_upper)) {
    return(log_mass_upper)
  }
  -Inf
}
