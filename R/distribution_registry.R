.qc_distribution_registry <- local({
  specs <- list(
    "normal" = list(
      name = "normal",
      display_name = "Normal",
      supports = list(
        pc = TRUE,
        mcmc = TRUE,
        integration = TRUE
      ),
      supports_summary_statistics = TRUE,
      sample_class = "normal",
      sample_parameter_names = c("mu", "sigma"),
      prior_parameter_names = c("mu", "sigma"),
      prior_defaults = list(
        mu = "Jeffreys_mu",
        sigma = "Jeffreys_sigma"
      ),
      s3_methods = list(
        percentiles = "samples_to_percentiles",
        posterior_predictives = "samples_to_posterior_predictives",
        E_abs_dev = "samples_to_E_abs_dev",
        pc_fit = "pc_fit_distribution"
      ),
      stan_model = "normal",
      integration_backend = ".bpc_fit_integration_normal"
    ),
    "t" = list(
      name = "t",
      display_name = "Student-t",
      supports = list(
        pc = TRUE,
        mcmc = TRUE,
        integration = FALSE
      ),
      supports_summary_statistics = FALSE,
      sample_class = "t",
      sample_parameter_names = c("mu", "scale", "nu"),
      prior_parameter_names = c("mu", "sigma", "nu"),
      prior_defaults = list(
        mu = "Jeffreys_mu",
        sigma = "Jeffreys_sigma",
        nu = "uniform_nu"
      ),
      s3_methods = list(
        percentiles = "samples_to_percentiles",
        posterior_predictives = "samples_to_posterior_predictives",
        E_abs_dev = "samples_to_E_abs_dev",
        pc_fit = "pc_fit_distribution"
      ),
      stan_model = "t",
      integration_backend = NULL
    )
  )

  function() specs
})

.qc_distribution_names <- function(method = NULL) {
  specs <- .qc_distribution_registry()

  if (is.null(method)) {
    return(names(specs))
  }

  method <- match.arg(method, choices = c("pc", "mcmc", "integration"))
  names(Filter(function(spec) isTRUE(spec$supports[[method]]), specs))
}

.qc_distribution_summary_stat_names <- function() {
  specs <- .qc_distribution_registry()
  names(Filter(function(spec) isTRUE(spec$supports_summary_statistics), specs))
}

.qc_distribution_spec <- function(distribution, method = NULL, arg = "distribution") {
  specs <- .qc_distribution_registry()
  BayesTools::check_char(
    distribution,
    name = arg,
    check_length = 1,
    allow_values = names(specs)
  )

  spec <- specs[[distribution]]

  if (!is.null(method)) {
    method <- match.arg(method, choices = c("pc", "mcmc", "integration"))
    if (!isTRUE(spec$supports[[method]])) {
      if (identical(method, "integration") && identical(distribution, "t")) {
        stop(
          "The integration method currently only supports distribution = 'normal'. ",
          "Use method = 'mcmc' for t-distribution.",
          call. = FALSE
        )
      }

      stop(
        sprintf(
          "`%s = \"%s\"` does not support `method = \"%s\"`.",
          arg,
          distribution,
          method
        ),
        call. = FALSE
      )
    }
  }

  spec
}

.qc_distribution_display_name <- function(distribution) {
  .qc_distribution_spec(distribution)$display_name
}

.qc_distribution_supports_summary_statistics <- function(distribution) {
  isTRUE(.qc_distribution_spec(distribution)$supports_summary_statistics)
}

.qc_distribution_parameter_names <- function(distribution, type = c("sample", "prior")) {
  spec <- .qc_distribution_spec(distribution)
  type <- match.arg(type)

  switch(
    type,
    "sample" = spec$sample_parameter_names,
    "prior" = spec$prior_parameter_names
  )
}

.qc_distribution_supports_parameter <- function(distribution, parameter,
                                                type = c("sample", "prior")) {
  parameter %in% .qc_distribution_parameter_names(
    distribution = distribution,
    type = match.arg(type)
  )
}

.qc_distribution_prior_defaults <- function(distribution) {
  prior_parameters <- .qc_distribution_parameter_names(
    distribution = distribution,
    type = "prior"
  )
  defaults <- .qc_distribution_spec(distribution)$prior_defaults %||% list()
  unknown_defaults <- setdiff(names(defaults), prior_parameters)

  if (length(unknown_defaults) > 0L) {
    stop(
      sprintf(
        "`distribution = \"%s\"` declares prior defaults for unknown parameter(s): %s.",
        distribution,
        paste(unknown_defaults, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  out <- stats::setNames(vector("list", length(prior_parameters)), prior_parameters)
  out[names(defaults)] <- defaults
  out
}

.qc_distribution_prior_map <- function(distribution, parameters = list()) {
  prior_parameters <- .qc_distribution_parameter_names(
    distribution = distribution,
    type = "prior"
  )
  prior_map <- .qc_distribution_prior_defaults(distribution)

  if (inherits(parameters, "PriorIndependent")) {
    parameters <- parameters$parameters
  }
  if (is.null(parameters)) {
    parameters <- list()
  }
  if (!is.list(parameters)) {
    stop("`parameters` must be a named list or a `PriorIndependent` object.", call. = FALSE)
  }

  parameter_names <- names(parameters) %||% rep("", length(parameters))
  if (length(parameters) == 0L) {
    return(prior_map[prior_parameters])
  }
  if (any(!nzchar(parameter_names))) {
    stop("`parameters` must be named.", call. = FALSE)
  }

  unknown_parameters <- setdiff(parameter_names, prior_parameters)
  if (length(unknown_parameters) > 0L) {
    supported <- paste(sprintf("`%s`", prior_parameters), collapse = ", ")
    stop(
      sprintf(
        "%s %s not a prior parameter for `distribution = \"%s\"`; supported prior parameters are %s.",
        paste(sprintf("`%s`", unknown_parameters), collapse = ", "),
        if (length(unknown_parameters) == 1L) "is" else "are",
        distribution,
        supported
      ),
      call. = FALSE
    )
  }

  for (parameter in parameter_names) {
    prior_map[parameter] <- list(parameters[[parameter]])
  }

  prior_map[prior_parameters]
}

.qc_distribution_stan_model <- function(distribution) {
  .qc_distribution_spec(distribution)$stan_model
}

.qc_distribution_integration_backend <- function(distribution) {
  spec <- .qc_distribution_spec(distribution)
  backend <- spec$integration_backend %||% NULL

  if (isTRUE(spec$supports[["integration"]]) &&
      (is.null(backend) || !nzchar(backend))) {
    stop(
      sprintf(
        "No integration backend is registered for `distribution = \"%s\"`.",
        distribution
      ),
      call. = FALSE
    )
  }

  backend
}

.qc_distribution_s3_method_name <- function(generic, distribution) {
  sprintf("%s.%s", generic, distribution)
}

.qc_distribution_s3_method_exists <- function(generic, distribution) {
  exists(
    .qc_distribution_s3_method_name(generic, distribution),
    envir = asNamespace("qc"),
    mode = "function",
    inherits = FALSE
  )
}

.qc_distribution_resolve_s3_method <- function(generic, distribution,
                                               allow_default = FALSE) {
  method_name <- .qc_distribution_s3_method_name(generic, distribution)
  fn <- get0(
    method_name,
    envir = asNamespace("qc"),
    mode = "function",
    inherits = FALSE
  )

  used_default <- FALSE
  if (is.null(fn) && allow_default) {
    default_name <- .qc_distribution_s3_method_name(generic, "default")
    fn <- get0(
      default_name,
      envir = asNamespace("qc"),
      mode = "function",
      inherits = FALSE
    )
    if (!is.null(fn)) {
      method_name <- default_name
      used_default <- TRUE
    }
  }

  list(
    generic = generic,
    method = method_name,
    fn = fn,
    used_default = used_default,
    available = !is.null(fn)
  )
}

.qc_distribution_has_registered_stan_model <- function(distribution) {
  stan_model_name <- .qc_distribution_stan_model(distribution)
  stan_registry <- get0("stanmodels", envir = asNamespace("qc"), inherits = FALSE)

  !is.null(stan_registry) &&
    stan_model_name %in% names(stan_registry) &&
    !is.null(stan_registry[[stan_model_name]])
}

.qc_distribution_resolve_integration_backend <- function(distribution) {
  backend_name <- .qc_distribution_integration_backend(distribution)
  backend <- get0(
    backend_name,
    envir = asNamespace("qc"),
    mode = "function",
    inherits = FALSE
  )

  list(
    backend_name = backend_name,
    backend = backend,
    available = !is.null(backend)
  )
}

.qc_distribution_has_registered_integration_backend <- function(distribution) {
  .qc_distribution_resolve_integration_backend(distribution)$available
}

.qc_distribution_adapter <- function(distribution, method = NULL,
                                     validate = TRUE) {
  spec <- .qc_distribution_spec(distribution, method = method)

  s3_generics <- spec$s3_methods %||% list(
    percentiles = "samples_to_percentiles",
    posterior_predictives = "samples_to_posterior_predictives",
    E_abs_dev = "samples_to_E_abs_dev",
    pc_fit = "pc_fit_distribution"
  )
  s3 <- list(
    percentiles = .qc_distribution_resolve_s3_method(
      s3_generics$percentiles,
      spec$name
    ),
    posterior_predictives = .qc_distribution_resolve_s3_method(
      s3_generics$posterior_predictives,
      spec$name
    ),
    E_abs_dev = .qc_distribution_resolve_s3_method(
      s3_generics$E_abs_dev,
      spec$name,
      allow_default = TRUE
    ),
    pc_fit = .qc_distribution_resolve_s3_method(
      s3_generics$pc_fit,
      spec$name
    )
  )

  stan_model_name <- spec$stan_model
  stan_registry <- get0("stanmodels", envir = asNamespace("qc"), inherits = FALSE)
  stan_model <- if (!is.null(stan_registry) &&
                    stan_model_name %in% names(stan_registry)) {
    stan_registry[[stan_model_name]]
  } else {
    NULL
  }

  integration <- if (isTRUE(spec$supports[["integration"]])) {
    .qc_distribution_resolve_integration_backend(spec$name)
  } else {
    list(backend_name = NULL, backend = NULL, available = FALSE)
  }

  adapter <- structure(
    list(
      name = spec$name,
      display_name = spec$display_name,
      supports = spec$supports,
      supports_summary_statistics = spec$supports_summary_statistics,
      parameters = list(
        sample = spec$sample_parameter_names,
        prior = spec$prior_parameter_names
      ),
      prior_defaults = .qc_distribution_prior_defaults(spec$name),
      sample_class = spec$sample_class %||% spec$name,
      s3 = s3,
      mcmc = list(
        stan_model_name = stan_model_name,
        stan_model = stan_model,
        available = !is.null(stan_model)
      ),
      integration = integration
    ),
    class = "qc_distribution_adapter"
  )

  if (validate) {
    .qc_distribution_validate_adapter(adapter, method = method, error = TRUE)
  }

  adapter
}

.qc_distribution_contract <- function(distribution, method = NULL) {
  adapter <- .qc_distribution_adapter(
    distribution = distribution,
    method = method,
    validate = FALSE
  )

  include_pc <- (is.null(method) || identical(method, "pc")) &&
    isTRUE(adapter$supports[["pc"]])
  include_mcmc <- (is.null(method) || identical(method, "mcmc")) &&
    isTRUE(adapter$supports[["mcmc"]])
  include_integration <- (is.null(method) || identical(method, "integration")) &&
    isTRUE(adapter$supports[["integration"]])

  rows <- list(
    list(
      capability = "sample percentiles",
      kind = "s3",
      symbol = adapter$s3$percentiles$method,
      required = include_pc || include_mcmc,
      available = adapter$s3$percentiles$available
    ),
    list(
      capability = "posterior predictive",
      kind = "s3",
      symbol = adapter$s3$posterior_predictives$method,
      required = include_mcmc,
      available = adapter$s3$posterior_predictives$available
    ),
    list(
      capability = "E_abs_dev",
      kind = "s3",
      symbol = if (adapter$s3$E_abs_dev$used_default) {
        paste0(
          .qc_distribution_s3_method_name("samples_to_E_abs_dev", distribution),
          " (or samples_to_E_abs_dev.default)"
        )
      } else {
        adapter$s3$E_abs_dev$method
      },
      required = include_pc || include_mcmc,
      available = adapter$s3$E_abs_dev$available
    )
  )

  if (isTRUE(adapter$supports[["pc"]]) && include_pc) {
    rows[[length(rows) + 1L]] <- list(
      capability = "pc fit",
      kind = "s3",
      symbol = adapter$s3$pc_fit$method,
      required = TRUE,
      available = adapter$s3$pc_fit$available
    )
  }

  if (isTRUE(adapter$supports[["mcmc"]]) && include_mcmc) {
    rows[[length(rows) + 1L]] <- list(
      capability = "stan model",
      kind = "object",
      symbol = sprintf(
        "stanmodels[[\"%s\"]]",
        adapter$mcmc$stan_model_name
      ),
      required = TRUE,
      available = adapter$mcmc$available
    )
  }

  if (isTRUE(adapter$supports[["integration"]]) && include_integration) {
    rows[[length(rows) + 1L]] <- list(
      capability = "integration fit",
      kind = "function",
      symbol = adapter$integration$backend_name,
      required = TRUE,
      available = adapter$integration$available
    )
  }

  data.frame(
    capability = vapply(rows, `[[`, character(1), "capability"),
    kind = vapply(rows, `[[`, character(1), "kind"),
    symbol = vapply(rows, `[[`, character(1), "symbol"),
    required = vapply(rows, `[[`, logical(1), "required"),
    available = vapply(rows, `[[`, logical(1), "available"),
    stringsAsFactors = FALSE
  )
}

.qc_distribution_validate_adapter <- function(adapter, method = NULL,
                                              error = TRUE) {
  contract <- .qc_distribution_contract(adapter$name, method = method)
  missing <- contract[contract$required & !contract$available, , drop = FALSE]

  if (nrow(missing) > 0L && error) {
    stop(
      sprintf(
        "`distribution = \"%s\"` is missing required adapter capability: %s.",
        adapter$name,
        paste(missing$symbol, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  missing
}

.qc_distribution_fit_integration <- function(distribution, ...) {
  adapter <- .qc_distribution_adapter(distribution, method = "integration")
  fit_fun <- adapter$integration$backend

  if (is.null(fit_fun)) {
    stop(
      sprintf(
        "No integration backend is registered for `distribution = \"%s\"`.",
        distribution
      ),
      call. = FALSE
    )
  }

  fit_fun(...)
}

.qc_distribution_missing_requirements <- function(distribution, method = NULL) {
  contract <- .qc_distribution_contract(distribution, method = method)
  contract$symbol[contract$required & !contract$available]
}
