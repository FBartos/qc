#' @title Options for the 'qc' package
#'
#' @description Get or set package-level options.
#'
#' @param name the name of the option to get the current value of - for a list of
#' available options, see details below.
#' @param ... named option(s) to change - for a list of available options, see
#' details below. Available options are \code{max_cores}, the maximum number
#' of cores used when a parallel backend is requested, and
#' \code{prior_mc_samples}, the number of Monte Carlo draws used for prior-only
#' integration summaries that require simulation.
#'
#' @return The current value of all available 'qc' options (after applying any
#' changes specified) is returned invisibly as a named list.
#'
#' @export qc.options
#' @export qc.get_option
#' @name qc_options
#' @aliases qc_options qc.options qc.get_option
NULL


#' @rdname qc_options
qc.options    <- function(...){

  opts <- list(...)
  option_names <- qc.private$.option_names

  for(i in seq_along(opts)){

    if(!names(opts)[i] %in% option_names)
      stop(paste("Unmatched or ambiguous option '", names(opts)[i], "'", sep=""))

    if (names(opts)[i] == "max_cores") {
      BayesTools::check_int(opts[[i]], "max_cores", lower = 1)
    }
    if (names(opts)[i] == "prior_mc_samples") {
      BayesTools::check_int(opts[[i]], "prior_mc_samples", lower = 1000)
    }

    assign(names(opts)[i], opts[[i]] , envir = qc.private)
  }

  values <- setNames(lapply(option_names, qc.get_option), option_names)
  return(invisible(values))
}

#' @rdname qc_options
qc.get_option <- function(name){
  option_names <- qc.private$.option_names

  if(length(name)!=1)
    stop("Only 1 option can be retrieved at a time")

  if(!name %in% option_names)
    stop(paste("Unmatched or ambiguous option '", name, "'", sep=""))

  # Use eval as some defaults are put in using 'expression' to avoid evaluating at load time:
  value <- qc.private[[name]]
  if (inherits(value, "expression")) {
    return(eval(value))
  }

  return(value)
}


qc.private <- new.env()
assign(".option_names", c("max_cores", "prior_mc_samples"), envir = qc.private)
detected_cores <- parallel::detectCores(logical = TRUE)
default_max_cores <- if (is.na(detected_cores)) 1L else max(1L, detected_cores - 1L)
assign("max_cores", default_max_cores, envir = qc.private)
assign("prior_mc_samples", 200000L, envir = qc.private)
