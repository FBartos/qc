#' Unified prior specifications for \code{bpc()}
#'
#' The \code{prior} argument in \code{\link{bpc}} accepts either a named prior
#' specification or a prior object created by one of these constructors.
#'
#' Named specifications can be supplied as strings or constructors:
#' \code{"DCSI"} / \code{prior_DCSI()},
#' \code{"Jeffreys"} / \code{prior_jeffreys()}, and
#' \code{"unit_information"} / \code{prior_unit_information()}.
#' Use \code{prior_unit_information(n0 = ...)} to adjust its information weight;
#' the string specification uses \code{n0 = 1}.
#'
#' Parameter-wise priors are supplied with \code{prior_independent()}, for
#' example \code{prior_independent(mu = ..., sigma = ...)} for the normal
#' likelihood or \code{prior_independent(mu = ..., sigma = ..., nu = ...)} for
#' the Student-t likelihood. Unspecified parameters use the distribution's
#' default parameter-wise prior.
#'
#' Joint priors on the mean and standard deviation are supplied with
#' \code{prior_conjugate()}, \code{prior_joint()}, or
#' \code{prior_semi_conjugate()}. Other parameters, such as the Student-t
#' degrees of freedom, use their default parameter-wise prior. With
#' \code{method = "mcmc"}, joint priors defined by R log-density functions
#' need their \code{bayestools_priors}.
#'
#' Every prior is available with both methods. When \code{prior} is omitted,
#' \code{bpc()} uses DCSI.
#'
#' DCSI is a Normal-Inverse-Gamma prior centered at the midpoint of the
#' specification limits. For lower limit \eqn{L}, upper limit \eqn{U}, and
#' width \eqn{W = U - L}, DCSI uses
#' \eqn{\mu_0 = (L + U) / 2}, \eqn{k_0 = 2}, \eqn{\nu = 6},
#' \eqn{c_0 = 1.25}, \eqn{\sigma_0 = W / (6 c_0)},
#' \eqn{\alpha_0 = \nu / 2}, and
#' \eqn{\beta_0 = \nu \sigma_0^2 / 2}.
#'
#' @examples
#' prior_DCSI()
#' prior_jeffreys()
#' prior_unit_information()
#' prior_unit_information(n0 = 5)
#' prior_independent(
#'   mu = BayesTools::prior("normal", list(0, 1)),
#'   sigma = "Jeffreys_sigma"
#' )
#' prior_conjugate(mu0 = 0, k0 = 2, alpha0 = 3, beta0 = 1)
#'
#' @name bpc_prior
NULL

#' Default conjugate specification-information prior
#'
#' Requests the default conjugate specification-information prior for the normal
#' likelihood. The prior is resolved by \code{\link{bpc}} from the specification
#' limits: it is centered at the midpoint of \code{LSL} and \code{USL}.
#'
#' @return A prior specification object.
#' @export
prior_DCSI <- function() {
  structure(list(name = "DCSI"), class = c("PriorDCSI", "PriorSpecification"))
}

#' Jeffreys prior specification
#'
#' Requests the default parameter-wise Jeffreys prior specification for the
#' selected likelihood.
#'
#' @return A prior specification object.
#' @export
prior_jeffreys <- function() {
  structure(list(name = "Jeffreys"), class = c("PriorJeffreys", "PriorSpecification"))
}

#' Unit-information prior specification
#'
#' Requests a data-derived unit-information prior for the normal likelihood.
#' The prior is resolved by \code{\link{bpc}} from the supplied observations or
#' normal summary statistics.
#'
#' @param n0 Positive, finite information weight. The default \code{1}
#'   gives the unit-information prior; larger values give a stronger prior,
#'   and values between zero and one give a weaker prior.
#' @details The Normal-Inverse-Gamma prior uses
#'   \eqn{\mu \mid \sigma^2 \sim N(\bar{x}, \sigma^2 / n_0)} and
#'   \eqn{\sigma^2 \sim \mathrm{InvGamma}(n_0 / 2, n_0 s^2 / 2)},
#'   where \eqn{s^2} is the unbiased sample variance and the inverse-gamma
#'   density is proportional to \eqn{v^{-\alpha-1}\exp(-\beta/v)}.
#' @return A prior specification object.
#' @export
prior_unit_information <- function(n0 = 1) {
  if (!is.numeric(n0) || length(n0) != 1L || !is.finite(n0) || n0 <= 0) {
    stop(gettext("`n0` must be a positive, finite number."), call. = FALSE)
  }
  structure(list(name = "unit_information", n0 = n0), class = c("PriorUnitInformation", "PriorSpecification"))
}

#' Independent parameter prior specification
#'
#' @param mu Prior for the location parameter.
#' @param sigma Prior for the process standard deviation parameter.
#' @param nu Prior for the Student-t degrees-of-freedom parameter.
#' @param ... Additional named distribution-specific parameter priors.
#' @return A prior specification object.
#' @export
prior_independent <- function(mu, sigma, nu, ...) {
  params <- list(...)
  if (is.null(names(params))) {
    names(params) <- rep("", length(params))
  }
  if (any(!nzchar(names(params)))) {
    stop("Additional priors in `...` must be named.", call. = FALSE)
  }

  formals <- list()
  if (!missing(mu)) {
    formals$mu <- mu
  }
  if (!missing(sigma)) {
    formals$sigma <- sigma
  }
  if (!missing(nu)) {
    formals$nu <- nu
  }

  duplicated_parameters <- intersect(names(formals), names(params))
  if (length(duplicated_parameters) > 0L) {
    stop(
      sprintf(
        "Prior parameter(s) supplied more than once: %s.",
        paste(duplicated_parameters, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  structure(
    list(parameters = c(formals, params)),
    class = c("PriorIndependent", "PriorSpecification")
  )
}

#' Conjugate normal likelihood prior
#'
#' @param mu0 Prior mean location.
#' @param k0 Prior precision multiplier. Use \code{0} for the default
#'   noninformative location component.
#' @param alpha0 Prior shape for \eqn{\sigma^2}.
#' @param beta0 Prior rate for \eqn{\sigma^2}.
#' @return A \code{PriorConjugate} object.
#' @export
prior_conjugate <- function(mu0 = 0, k0 = 0, alpha0 = -0.5, beta0 = 0) {
  create_prior_conjugate(mu0 = mu0, k0 = k0, alpha0 = alpha0, beta0 = beta0)
}

#' Joint normal likelihood prior
#'
#' @param log_dens_fn Function \code{function(mu, sigma)} returning the joint
#'   log prior density.
#' @param bayestools_priors Optional list preserving original BayesTools prior
#'   objects for diagnostics.
#' @return A \code{PriorGeneric} object.
#' @export
prior_joint <- function(log_dens_fn, bayestools_priors = NULL) {
  create_prior_generic(log_dens_fn = log_dens_fn, bayestools_priors = bayestools_priors)
}

#' Semi-conjugate normal likelihood prior
#'
#' @param conjugate Which component uses the conjugate form: \code{"mu"} or
#'   \code{"sigma"}.
#' @param ... Arguments for the selected semi-conjugate prior. Use
#'   \code{mu0}, \code{k0}, and \code{log_dens_sigma} when
#'   \code{conjugate = "mu"}; use \code{alpha0}, \code{beta0}, and
#'   \code{log_dens_mu} when \code{conjugate = "sigma"}.
#' @return A semi-conjugate prior object.
#' @export
prior_semi_conjugate <- function(conjugate = c("mu", "sigma"), ...) {
  conjugate <- match.arg(conjugate)
  args <- list(...)

  if (identical(conjugate, "mu")) {
    return(do.call(create_prior_semi_mu, args))
  }

  do.call(create_prior_semi_sigma, args)
}

.prior_DCSI <- function(LSL, USL) {
  W <- USL - LSL
  midpoint <- LSL + W / 2
  k0 <- 2
  nu <- 6
  c0 <- 1.25
  sigma0 <- W / (6 * c0)

  create_prior_conjugate(
    mu0 = midpoint,
    k0 = k0,
    alpha0 = nu / 2,
    beta0 = nu * sigma0^2 / 2
  )
}

.is_prior_specification <- function(prior) {
  inherits(prior, "PriorSpecification") || .is_integration_prior(prior)
}

.normalize_prior_specification <- function(prior) {
  if (.is_prior_specification(prior)) {
    return(prior)
  }

  if (is.character(prior)) {
    BayesTools::check_char(prior, name = "prior", check_length = 1)
    key <- tolower(gsub("[-_[:space:]]", "", prior))
    return(switch(
      key,
      "dcsi" = prior_DCSI(),
      "jeffreys" = prior_jeffreys(),
      "unitinformation" = prior_unit_information(),
      stop(
        "`prior` must be one of \"DCSI\", \"Jeffreys\", \"unit_information\", ",
        "or a prior specification object.",
        call. = FALSE
      )
    ))
  }

  if (inherits(prior, "prior")) {
    stop(
      "BayesTools prior objects must be wrapped in `prior_independent()`, ",
      "for example `prior = prior_independent(mu = ..., sigma = ...)`.",
      call. = FALSE
    )
  }

  stop(
    "`prior` must be a character prior name or a prior specification object.",
    call. = FALSE
  )
}

.prior_unit_information_from_state <- function(cached_state, n0 = 1) {
  n <- cached_state$n %||% 0L
  if (n < 2L) {
    stop(
      "`prior = \"unit_information\"` requires at least two observations or ",
      "normal summary statistics with `N >= 2`.",
      call. = FALSE
    )
  }

  s2 <- cached_state$sse / (n - 1L)
  create_prior_conjugate(
    mu0 = cached_state$x_bar,
    k0 = n0,
    alpha0 = n0 / 2,
    beta0 = n0 * s2 / 2
  )
}

.resolve_independent_prior_map <- function(distribution, parameters) {
  prior_parameters <- .qc_distribution_parameter_names(distribution, type = "prior")
  unknown_parameters <- setdiff(names(parameters), prior_parameters)

  if (length(unknown_parameters) > 0L) {
    stop(
      sprintf(
        "%s %s not a prior parameter for `distribution = \"%s\"`.",
        paste(sprintf("`%s`", unknown_parameters), collapse = ", "),
        if (length(unknown_parameters) == 1L) "is" else "are",
        distribution
      ),
      call. = FALSE
    )
  }

  prior_map <- .qc_distribution_prior_defaults(distribution)
  for (parameter in names(parameters)) {
    prior_map[parameter] <- list(parameters[[parameter]])
  }

  prior_map[prior_parameters]
}

.resolve_bpc_prior <- function(prior, distribution, method,
                               LSL, USL, cached_state) {
  prior_info <- .resolve_bpc_prior_info(
    prior = prior,
    distribution = distribution,
    method = method,
    LSL = LSL,
    USL = USL,
    cached_state = cached_state
  )

  if (identical(method, "mcmc")) {
    prior_info$stan_prior_map <- .stan_prior_map(prior_info, distribution)
  }

  prior_info
}

.has_location_scale_priors <- function(distribution) {
  all(c("mu", "sigma") %in% .qc_distribution_parameter_names(distribution, type = "prior"))
}

.resolve_bpc_prior_info <- function(prior, distribution, method,
                                    LSL, USL, cached_state) {
  specification <- .normalize_prior_specification(prior)

  if (inherits(specification, "PriorDCSI")) {
    if (!.has_location_scale_priors(distribution)) {
      stop(sprintf("`prior = \"DCSI\"` is not available for `distribution = \"%s\"`.", distribution), call. = FALSE)
    }

    prior_obj <- .prior_DCSI(LSL = LSL, USL = USL)
    return(list(
      specification = specification,
      prior = prior_obj,
      prior_map = NULL,
      kind = "DCSI"
    ))
  }

  if (inherits(specification, "PriorUnitInformation")) {
    if (!.has_location_scale_priors(distribution)) {
      stop(sprintf("`prior = \"unit_information\"` is not available for `distribution = \"%s\"`.", distribution), call. = FALSE)
    }

    prior_obj <- .prior_unit_information_from_state(cached_state, n0 = specification$n0)
    return(list(
      specification = specification,
      prior = prior_obj,
      prior_map = NULL,
      kind = "unit_information"
    ))
  }

  if (inherits(specification, "PriorJeffreys")) {
    specification <- prior_independent()
  }

  if (inherits(specification, "PriorIndependent")) {
    prior_map <- .resolve_independent_prior_map(
      distribution = distribution,
      parameters = specification$parameters
    )

    if (identical(method, "integration")) {
      if (!identical(distribution, "normal")) {
        stop("The integration method currently only supports distribution = 'normal'.", call. = FALSE)
      }

      prior_info <- .bayestools_to_integration_prior(prior_map$mu, prior_map$sigma)
      return(list(
        specification = specification,
        prior = prior_info$prior,
        prior_map = prior_map,
        kind = "independent",
        integration_case = prior_info$case,
        is_conjugate = prior_info$is_conjugate
      ))
    }

    return(list(
      specification = specification,
      prior = specification,
      prior_map = prior_map,
      kind = "independent"
    ))
  }

  if (.is_integration_prior(specification)) {
    if (!.has_location_scale_priors(distribution)) {
      stop(sprintf("Joint location-scale priors are not available for `distribution = \"%s\"`.", distribution), call. = FALSE)
    }

    return(list(
      specification = specification,
      prior = specification,
      prior_map = NULL,
      kind = "joint"
    ))
  }

  stop("Unsupported prior specification.", call. = FALSE)
}

# Translates a resolved prior into the parameter-wise map used by the Stan
# backend. Joint priors use the same (mu, sigma) density as the integration
# backend; the remaining parameters (e.g., nu) use the distribution defaults.
.stan_prior_map <- function(prior_info, distribution) {
  prior_map <- prior_info$prior_map

  if (!is.null(prior_map)) {
    uses_conjugate <- vapply(prior_map, inherits, logical(1), what = "PriorConjugate")
    if (!any(uses_conjugate)) {
      return(prior_map)
    }
    location_scale <- .stan_location_scale_priors(
      .bayestools_to_integration_prior(prior_map$mu, prior_map$sigma)$prior
    )
  } else {
    prior_map <- .qc_distribution_prior_defaults(distribution)
    location_scale <- .stan_location_scale_priors(prior_info$prior)
  }

  prior_map$mu <- location_scale$mu
  prior_map$sigma <- location_scale$sigma
  prior_map
}

.stan_location_scale_priors <- function(prior) {
  # sigma^(-(2 * alpha0 + 1)) * exp(-beta0 / sigma^2) is InvGamma(alpha0, beta0)
  # on sigma^2 expressed as a density on sigma
  conjugate_power <- function(alpha0) 2 * alpha0 + 1

  if (inherits(prior, "PriorConjugate")) {
    # the flat k0 = 0 limit keeps the sigma^-1 normalizing term of N(mu0, sigma^2 / k0)
    return(list(
      mu = .stan_conditional_mu_prior(prior$mu0, prior$k0),
      sigma = .stan_conjugate_sigma_prior(conjugate_power(prior$alpha0) + (prior$k0 == 0), prior$beta0)
    ))
  }

  bayestools_priors <- prior$bayestools_priors
  required <- switch(
    class(prior)[1],
    "PriorSemiConjugateMu" = "sigma",
    "PriorSemiConjugateSigma" = "mu",
    "PriorGeneric" = c("mu", "sigma")
  )
  if (!all(required %in% names(bayestools_priors))) {
    stop(
      "`method = \"mcmc\"` requires priors that can be evaluated by Stan. ",
      "Joint priors defined only through R log-density functions need ",
      "`bayestools_priors` with BayesTools priors for ",
      paste(sprintf("`%s`", required), collapse = " and "),
      ".",
      call. = FALSE
    )
  }

  if (inherits(prior, "PriorSemiConjugateMu")) {
    return(list(
      mu = .stan_conditional_mu_prior(prior$mu0, prior$k0),
      sigma = bayestools_priors$sigma
    ))
  }

  if (inherits(prior, "PriorSemiConjugateSigma")) {
    # the semi-conjugate backend treats (alpha0, beta0) = (-0.5, 0) as Jeffreys 1 / sigma
    jeffreys <- prior$alpha0 == -0.5 && prior$beta0 == 0
    return(list(
      mu = bayestools_priors$mu,
      sigma = .stan_conjugate_sigma_prior(conjugate_power(prior$alpha0) + jeffreys, prior$beta0)
    ))
  }

  list(mu = bayestools_priors$mu, sigma = bayestools_priors$sigma)
}
