.new_qc_suff_stats_state <- function(n, x_bar, sse) {
  structure(
    list(
      kind = "suff_stats",
      state_version = 1L,
      n = as.integer(n),
      x_bar = x_bar,
      sse = sse
    ),
    class = c("qc_suff_stats_state", "list")
  )
}

.integration_missing_state_fields <- function(cached_state, required_fields,
                                             allow_null = character()) {
  required_fields[vapply(
    required_fields,
    function(name) {
      !(name %in% names(cached_state)) ||
        (is.null(cached_state[[name]]) && !name %in% allow_null)
    },
    logical(1)
  )]
}

.as_qc_suff_stats_state <- function(data = NULL, cached_state = NULL,
                                    arg = "cached_state") {
  if (!is.null(cached_state)) {
    required_fields <- c("n", "x_bar", "sse")
    missing_fields <- .integration_missing_state_fields(cached_state, required_fields)
    if (length(missing_fields) > 0L) {
      stop(
        arg, " must contain `n`, `x_bar`, and `sse` when provided. Missing: ",
        paste(missing_fields, collapse = ", "),
        call. = FALSE
      )
    }

    return(.new_qc_suff_stats_state(
      n = cached_state$n,
      x_bar = cached_state$x_bar,
      sse = cached_state$sse
    ))
  }

  data <- data %||% numeric(0)
  n <- length(data)
  if (n > 0L) {
    x_bar <- mean(data)
    return(.new_qc_suff_stats_state(
      n = n,
      x_bar = x_bar,
      sse = sum((data - x_bar)^2)
    ))
  }

  .new_qc_suff_stats_state(n = 0L, x_bar = 0, sse = 0)
}

.integration_is_prior_only <- function(data = NULL, cached_state = NULL) {
  .as_qc_suff_stats_state(data = data, cached_state = cached_state)$n == 0L
}

.new_qc_integration_request <- function(data,
                                        LSL,
                                        USL,
                                        prior,
                                        metric = "Cpk",
                                        target = NULL,
                                        cached_state = NULL,
                                        sigma_level = 3) {
  metric <- .validate_metric_name(metric)
  resolved_target <- .integration_resolve_target(metric, target, LSL, USL)
  .validate_capability_request(
    LSL = LSL,
    USL = USL,
    target = resolved_target,
    sigma_level = sigma_level,
    metric = metric,
    target_required = !is.null(resolved_target),
    sigma_name = "sigma_level"
  )

  structure(
    list(
      kind = "integration_request",
      request_version = 1L,
      data = data %||% numeric(0),
      LSL = LSL,
      USL = USL,
      prior = prior,
      metric = metric,
      requested_target = target,
      target = resolved_target,
      cached_state = cached_state,
      sigma_level = sigma_level
    ),
    class = c("qc_integration_request", "list")
  )
}

.validate_qc_integration_request <- function(request, arg = "request") {
  required_fields <- c(
    "data", "LSL", "USL", "prior", "metric", "target",
    "cached_state", "sigma_level"
  )
  missing_fields <- .integration_missing_state_fields(
    request,
    required_fields,
    allow_null = c("target", "cached_state")
  )
  if (length(missing_fields) > 0L) {
    stop(
      arg, " is missing required integration request fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  .validate_capability_request(
    LSL = request$LSL,
    USL = request$USL,
    target = request$target,
    sigma_level = request$sigma_level,
    metric = request$metric,
    target_required = !is.null(request$target),
    sigma_name = "sigma_level"
  )

  invisible(request)
}

.as_qc_integration_request <- function(request = NULL,
                                       data = NULL,
                                       LSL = NULL,
                                       USL = NULL,
                                       prior = NULL,
                                       metric = "Cpk",
                                       target = NULL,
                                       cached_state = NULL,
                                       sigma_level = 3,
                                       arg = "request") {
  if (!is.null(request)) {
    .validate_qc_integration_request(request, arg = arg)
    return(request)
  }

  .new_qc_integration_request(
    data = data,
    LSL = LSL,
    USL = USL,
    prior = prior,
    metric = metric,
    target = target,
    cached_state = cached_state,
    sigma_level = sigma_level
  )
}

.integration_request_update <- function(request,
                                        data = request$data,
                                        LSL = request$LSL,
                                        USL = request$USL,
                                        prior = request$prior,
                                        metric = request$metric,
                                        target = request$requested_target,
                                        cached_state = request$cached_state,
                                        sigma_level = request$sigma_level) {
  request <- .as_qc_integration_request(request = request)

  .new_qc_integration_request(
    data = data,
    LSL = LSL,
    USL = USL,
    prior = prior,
    metric = metric,
    target = target,
    cached_state = cached_state,
    sigma_level = sigma_level
  )
}

.new_generic_backend_state <- function(...) {
  structure(
    c(
      list(
        kind = "generic_cached_state",
        state_version = 1L
      ),
      list(...)
    ),
    class = c("qc_generic_cached_state", "qc_suff_stats_state", "list")
  )
}

.validate_qc_generic_cached_state <- function(state, arg = "cached_state") {
  required_fields <- c(
    "n", "x_bar", "sse",
    "log_post", "log_post_vec", "h_max", "int_2d", "Z",
    "map_mu", "map_sig", "uni_s",
    "mu_lower", "mu_upper", "sigma_lower", "sigma_upper"
  )
  missing_fields <- .integration_missing_state_fields(state, required_fields)
  if (length(missing_fields) > 0L) {
    stop(
      arg, " is missing required generic integration fields: ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  invisible(state)
}

.new_integration_suff_state <- .new_qc_suff_stats_state
.integration_require_suff_state <- .as_qc_suff_stats_state
