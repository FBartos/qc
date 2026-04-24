
#' Process Capability
#'
#' @param bootstrap whether to run the bootstrap to compute CI
#' @param samples number of bootstrap samples for calculating CI
#' @param parallel whether to run the bootstrap in parallel
#' @param control a list of control settings
#' @param seed a random seed for reproducibility, defaults to NULL
#' @param silent whether to suppress output during fitting, defaults to TRUE
#' @param ... additional arguments
#' @inheritParams bpc
#'
#' @examples
#' data("pistonrings", package = "qc")
#' diameter <- subset(pistonrings, trial)[["diameter"]]
#'
#' # Point estimates only.
#' fit <- pc(
#'   diameter,
#'   LSL = 73.98, target = 74.00, USL = 74.02,
#'   bootstrap = FALSE
#' )
#' coef(fit)
#' summary(fit)
#'
#' # Bootstrap intervals and interval probabilities.
#' fit_boot <- pc(
#'   diameter,
#'   LSL = 73.98, target = 74.00, USL = 74.02,
#'   bootstrap = TRUE, samples = 200, seed = 1
#' )
#' summary(fit_boot)
#'
#' # Re-summarize the same fitted process with wider specification limits.
#' summary(fit_boot, LSL = 73.975, target = 74.00, USL = 74.025)
#'
#' @export
pc <- function(
    x,
    LSL, target, USL,

    # prior settings
    distribution = "normal",

    # bootstrap control settings
    bootstrap = TRUE, samples = 1000, cores = NULL, parallel = FALSE,
    control = NULL,
    seed = NULL, silent = TRUE,

    # capability metrics settings
    sigma = 3,
    ...) {

  # check input for capability metrics calculation to fail fast before estimation
  .validate_capability_request(
    LSL = LSL,
    USL = USL,
    target = target,
    sigma_level = sigma,
    sigma_name = "sigma"
  )
  .qc_distribution_spec(distribution, method = "pc")

  # prepare data
  data <- .bpc_data(distribution = distribution, x = x)

  # collect control settings
  control <- .optim_check_and_list_fit_settings(
    bootstrap = bootstrap, samples = samples, parallel = parallel, cores = cores, seed = seed, control = control
  )

  fit <- .pc_fit(distribution = distribution, data = data, control = control)

  object <- .new_pc_fit(
    call = match.call(),
    distribution = distribution,
    data = data,
    control = control,
    fit = fit,
    sigma = sigma
  )

  # compute capability metrics
  object$metrics      <- .compute_capability_metrics(object, LSL = LSL, USL = USL, target = target, sigma = sigma, bootstrap = FALSE)
  object$metrics_boot <- if (isTRUE(object$control[["bootstrap"]])) {
    .compute_capability_metrics(object, LSL = LSL, USL = USL, target = target, sigma = sigma, bootstrap = TRUE)
  } else {
    NULL
  }

  # add coefficients
  object$coefficients <- unlist(object$metrics)

  return(object)
}


.pc_fit           <- function(distribution, data, control) {

  if(!is.null(control[["seed"]]))
    set.seed(control[["seed"]])

  fit <- .pc_single_fit(distribution = distribution, data = data, control = control)
  attr(fit, "distribution") <- distribution

  if (control[["bootstrap"]]) {
    boot_fit <- .pc_bootstrap_fit(distribution = distribution, data = data, control = control)
    attr(boot_fit, "distribution") <- distribution
  }

  return(list(
    fit      = fit,
    boot_fit = if (control[["bootstrap"]]) boot_fit else NULL
  ))
}
.pc_single_fit    <- function(distribution, data, control) {
  distribution <- structure(distribution, class = distribution)
  pc_fit_distribution(distribution, data, control)
}
.pc_bootstrap_fit <- function(distribution, data, control) {

  if (!is.null(control[["seed"]]))
    set.seed(control[["seed"]])


  if (control[["parallel"]]) {

    # TODO: this should not create all bootstrap datasets at once, which is memory inefficient, but I guess this is easier to implement for now.
    data <- lapply(seq_len(control[["samples"]]), function(i) {
      list(x = sample(data$x, size = length(data$x), replace = TRUE))
    })
    cl <- parallel::makeCluster(control[["cores"]])
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterEvalQ(cl, {library("qc")})
    parallel::clusterExport(cl, c("distribution", "data", "control"), envir = environment())
    out <- parallel::parLapplyLB(cl, seq_len(control[["samples"]]), function(i) {
      .pc_single_fit(distribution = distribution, data = data[[i]], control = control)
    })

  } else {

    out  <- vector("list", control[["samples"]])
    for (i in seq_len(control[["samples"]])) {
      data_i <- list(x = sample(data$x, size = length(data$x), replace = TRUE))
      out[[i]]  <- .pc_single_fit(distribution = distribution, data = data_i, control = control)
    }

  }

  out <- do.call(rbind.data.frame, out)

  return(out)
}

pc_fit_distribution        <- function(distribution, data, control) {
  UseMethod("pc_fit_distribution", distribution)
}
#' @export
pc_fit_distribution.normal <- function(distribution, data, control) {
  return(list(
    mu    = mean(data[["x"]]),
    sigma = stats::sd(data[["x"]])
  ))
}

lpdf_scaled_t <- function(x, df, mu, sigma) {
  # base R version that underflows more slowly than extraDistr::dlst
  stats::dt((x - mu) / sigma, df = df, log = TRUE) - log(sigma)
}
# reference version for correctness checks
# lpdf_scaled_t2 <- function(x, df, mu, sigma) {
#   extraDistr::dlst(x, df = df, mu = mu, sigma = sigma, log = TRUE)
# }

#' @export
pc_fit_distribution.t      <- function(distribution, data, control) {

  fit <- try(stats::optim(
    par     = c(df = 30, mu = mean(data[["x"]]), sigma = stats::sd(data[["x"]])),
    fn      = function(par, x) {
      ret_val <- -sum(lpdf_scaled_t(x, par[["df"]], par[["mu"]], par[["sigma"]]))
      # if (is.infinite(ret_val) || is.na(ret_val)) browser()
      return(ret_val)
    },
    x       = data[["x"]],
    lower   = c(2, -Inf, stats::sd(data[["x"]]) / 1e3),
    method  = "L-BFGS-B"
  ))

  if (inherits(fit, "try-error"))
    return(list(
      mu    = NA,
      scale = NA,
      nu    = NA
    ))

  return(list(
    mu    = fit$par[["mu"]],
    scale = fit$par[["sigma"]],
    nu    = fit$par[["df"]]
  ))
}

.optim_check_and_list_fit_settings  <- function(bootstrap, samples, parallel, cores, seed, control, call = ""){

  BayesTools::check_int(samples,      "samples",   lower = 10,        call = call)
  BayesTools::check_bool(parallel,    "parallel",                     call = call)
  BayesTools::check_int(cores,        "cores",     allow_NULL = TRUE, lower = 1,  call = call)
  BayesTools::check_bool(bootstrap,   "bootstrap",                    call = call)
  BayesTools::check_int(seed,         "seed",      allow_NULL = TRUE, call = call)

  if(!parallel){
    cores <- 1
  }else if(is.null(cores)){
    cores <- qc.get_option("max_cores")
  }

  return(invisible(list(
    bootstrap  = bootstrap,
    samples    = samples,
    parallel   = parallel,
    cores      = cores,
    seed       = seed,
    control    = control
  )))
}

.pc_query_metrics <- function(object,
                              limits = NULL,
                              sigma = object$sigma %||% 3) {
  has_bootstrap <- isTRUE(object$control[["bootstrap"]]) &&
    !is.null(object$fit[["boot_fit"]])

  if (is.null(limits)) {
    return(list(
      metrics = object$metrics,
      metrics_boot = object$metrics_boot
    ))
  }

  list(
    metrics = .compute_capability_metrics(
      object,
      LSL = limits$LSL,
      USL = limits$USL,
      target = limits$target,
      sigma = sigma,
      bootstrap = FALSE
    ),
    metrics_boot = if (has_bootstrap) {
      .compute_capability_metrics(
        object,
        LSL = limits$LSL,
        USL = limits$USL,
        target = limits$target,
        sigma = sigma,
        bootstrap = TRUE
      )
    } else {
      NULL
    }
  )
}

.pc_summary_tables <- function(metrics, metrics_boot, ci.level, interval_probability) {
  if (is.null(metrics_boot)) {
    summary <- tibble::tibble(
      metric = names(metrics),
      mean = as.numeric(metrics),
      median = NA_real_,
      sd = NA_real_,
      lower = NA_real_,
      upper = NA_real_
    )

    interval_columns <- levels(cut(
      0,
      breaks = c(-Inf, interval_probability, Inf),
      include.lowest = TRUE
    ))
    interval_summary <- tibble::as_tibble(matrix(
      NA_real_,
      nrow = length(metrics),
      ncol = length(interval_columns),
      dimnames = list(names(metrics), interval_columns)
    ), rownames = "metric")

    return(list(
      summary = summary,
      interval_summary = interval_summary
    ))
  }

  distributions <- .as_qc_metric_distributions(metrics_boot, what = names(metrics))
  tables <- .qc_metric_distributions_summary_bundle(
    distributions,
    ci_level = ci.level,
    interval_probability = interval_probability,
    mean_override = metrics,
    include_divergence = FALSE
  )

  list(
    summary = tables$summary,
    interval_summary = tables$interval_summary
  )
}

### print and summary functions ----
#' @export
print.pc <- function(x, ...) {

  cat("\nCall:\n")
  print(x$call)
  cat("\nProcess Capability:\n")
  print(round(coef(x), 4))

  invisible(x)
}

#' @export
summary.pc <- function(object, LSL, target, USL, sigma = 3, ci.level = 0.95, interval_probability = c(1.00, 1.33, 1.50, 2.00), ...) {

  if (missing(sigma)) {
    sigma <- object$sigma %||% 3
  }
  .validate_sigma_level(sigma, name = "sigma")
  BayesTools::check_real(ci.level, name = "ci.level", check_length = 1, allow_NA = FALSE, lower = 0, upper = 1)
  BayesTools::check_real(interval_probability, name = "interval_probability", check_length = 0, allow_NA = FALSE)

  limits <- .qc_requested_limits(
    LSL = LSL,
    target = target,
    USL = USL,
    LSL_missing = missing(LSL),
    target_missing = missing(target),
    USL_missing = missing(USL)
  )
  query <- .pc_query_metrics(
    object,
    limits = limits,
    sigma = sigma
  )
  metrics <- query$metrics
  metrics_boot <- query$metrics_boot

  tables <- .pc_summary_tables(
    metrics = metrics,
    metrics_boot = metrics_boot,
    ci.level = ci.level,
    interval_probability = interval_probability
  )

  .new_pc_summary(
    call             = object$call,
    metrics          = metrics,
    metrics_boot     = metrics_boot,
    summary          = tables$summary,
    interval_summary = tables$interval_summary
  )
}

#' @export
print.pc_summary <- function(x, ...) {
  cat("\nCall:\n")
  print(x$call)

  cat("\nProcess Capability:\n")
  print(as.data.frame(round(x$summary[,-1], 4)), quote = FALSE, right = TRUE, row.names = unlist(x$summary[,1]))

  cat("\nInterval Probability:\n")
  print(as.data.frame(round(x$interval_summary[,-1], 4)), quote = FALSE, right = TRUE, row.names = unlist(x$interval_summary[,1]))

  invisible(x$summary)
}

