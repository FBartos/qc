
#' Process Capability
#'
#' @param bootstrap whether to run the bootstrap to compute CI
#' @param samples number of bootstrap samples for calculating CI
#' @param parallel whether to run the bootstrap in parallel
#' @param control a list of control settings
#' @param seed a random seed for reproducibility, defaults to NULL
#' @param silent whether to suppress output during fitting, defaults to TRUE
#' @param ... additional arguments
#' @inherit bpc
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
    sigma = 3, force_normal = FALSE,
    ...) {

  ### Instructions for adding a new distribution
  # the following functions need to be created:
  # samples_to_mu_and_sigma.<distribution> (S3 Class)
  # samples_to_percentiles.<distribution> (S3 Class)
  # samples_to_posterior_predictives.<distribution> (S3 Class)
  # samples_to_E_abs_dev.<distribution> (S3 Class)
  # pc_fit_distribution.<distribution> (S3 Class)
  # the following functions need to be extended:


  dots         <- list(...)
  object       <- list()
  object$call  <- match.call()

  # check input for capability metrics calculation to fail fast before estimation
  .validate_LSL_USL_target(LSL = LSL, USL = USL, target = target)
  BayesTools::check_char(distribution, name = "distribution", check_length = 1, allow_values = c("normal", "t"))
  BayesTools::check_real(sigma, name = "sigma", check_length = 1, lower = 0, allow_NA = FALSE)
  BayesTools::check_bool(force_normal, name = "force_normal", check_length = 1, allow_NA = FALSE)

  # prepare data
  object$data <- .bpc_data(x = x)

  # collect control settings
  object$control <- .optim_check_and_list_fit_settings(
    bootstrap = bootstrap, samples = samples, parallel = parallel, cores = cores, seed = seed, control = control
  )

  # fit stan model
  object$distribution <- distribution
  object$fit          <- .pc_fit(distribution = distribution, data = object$data, control = object$control)

  # add class to the object because it's required for dispatching in compute_capability_metrics
  class(object) <- "pc"

  # compute capability metrics
  object$metrics      <- .compute_capability_metrics(object, LSL = LSL, USL = USL, target = target, sigma = sigma, force_normal = force_normal, bootstrap = FALSE)
  object$metrics_boot <- .compute_capability_metrics(object, LSL = LSL, USL = USL, target = target, sigma = sigma, force_normal = force_normal, bootstrap = TRUE)

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

  # bootstrap the data
  data <- lapply(seq_len(control[["samples"]]), function(i) {
    list(x = sample(data$x, size = length(data$x), replace = TRUE))
  })

  if (control[["parallel"]]) {

    cl <- parallel::makeCluster(control[["cores"]])
    parallel::clusterEvalQ(cl, {library("qc")})
    parallel::clusterExport(cl, c("distribution", "data", "control"), envir = environment())
    out <- parallel::parLapplyLB(cl, seq_len(control[["samples"]]), function(i) {
      .pc_single_fit(distribution = distribution, data = data[[i]], control = control)
    })

  } else {

    out  <- vector("list", control[["samples"]])
    for (i in seq_len(control[["samples"]])) {
      out[[i]]  <- .pc_single_fit(distribution = distribution, data = data[[i]], control = control)
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
#' @export
pc_fit_distribution.t      <- function(distribution, data, control) {

  fit <- try(stats::optim(
    par     = c(df = 30, mu = mean(data[["x"]]), sigma = stats::sd(data[["x"]])),
    fn      = function(par, x) -sum(extraDistr::dlst(x, df = par["df"], mu = par["mu"], sigma = par["sigma"], log = TRUE)),
    x       = data[["x"]],
    lower   = c(2, -Inf, stats::sd(x) / 1e3),
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
summary.pc <- function(object, LSL, target, USL, sigma = 3, force_normal = FALSE, ci.level = 0.95, interval_probability = c(1.00, 1.33, 1.50, 2.00), ...) {

  BayesTools::check_real(sigma, name = "sigma", check_length = 1, lower = 0, allow_NA = FALSE)
  BayesTools::check_real(ci.level, name = "ci.level", check_length = 1, allow_NA = FALSE, lower = 0, upper = 1)
  BayesTools::check_real(interval_probability, name = "interval_probability", check_length = 0, allow_NA = FALSE)
  BayesTools::check_bool(force_normal, name = "force_normal", check_length = 1, allow_NA = FALSE)

  # recompute capability metrics if specification limits are set
  if (!missing(LSL) && !missing(target) && !missing(USL)) {
    metrics      <- .compute_capability_metrics(object, LSL = LSL, USL = USL, target = target, sigma = sigma, force_normal = force_normal, bootstrap = FALSE)
    metrics_boot <- .compute_capability_metrics(object, LSL = LSL, USL = USL, target = target, sigma = sigma, force_normal = force_normal, bootstrap = TRUE)
  } else {
    metrics      <- object$metrics
    metrics_boot <- object$metrics_boot
  }

  ### compute mean, median, and credible intervals for the metrics
  h     <- (1 - ci.level) / 2
  probs <- c(h, .5, 1 - h)

  summary <- cbind(mean = as.numeric(metrics) ,t(vapply(metrics_boot, function(x) {
    quantiles <- unname(stats::quantile(x, probs = probs, na.rm = TRUE))
    c(median = quantiles[2], sd = stats::sd(x), lower = quantiles[1], upper = quantiles[3])
  }, numeric(4L))))
  summary <- tibble::as_tibble(summary, rownames = "metric")

  ### compute interval summaries
  interval_summary <- t(vapply(metrics_boot, FUN = function(x) {
    table(cut(x, breaks = c(-Inf, interval_probability, Inf), include.lowest = TRUE)) / length(x)
  }, FUN.VALUE = numeric(length(interval_probability) + 1L)))
  interval_summary <- tibble::as_tibble(interval_summary, rownames = "metric")

  out <- list(
    call             = object$call,
    summary          = summary,
    interval_summary = interval_summary,
    metrics          = metrics
  )
  class(out) <- "pc_summary"

  return(out)
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

