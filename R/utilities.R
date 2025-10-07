#' @title Options for the 'qc' package
#'
#' @description A placeholder object and functions for the 'qc' package.
#'
#' @param name the name of the option to get the current value of - for a list of
#' available options, see details below.
#' @param ... named option(s) to change - for a list of available options, see
#' details below.
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

  for(i in seq_along(opts)){

    if(!names(opts)[i] %in% names(qc.private))
      stop(paste("Unmatched or ambiguous option '", names(opts)[i], "'", sep=""))

    assign(names(opts)[i], opts[[i]] , envir = qc.private)
  }

  return(invisible(qc.private$options))
}

#' @rdname qc_options
qc.get_option <- function(name){

  if(length(name)!=1)
    stop("Only 1 option can be retrieved at a time")

  if(!name %in% names(qc.private))
    stop(paste("Unmatched or ambiguous option '", name, "'", sep=""))

  # Use eval as some defaults are put in using 'expression' to avoid evaluating at load time:
  return(eval(qc.private[[name]]))
}


qc.private <- new.env()
# Use 'expression' for functions to avoid having to evaluate before the package is fully loaded:
assign("default_toptions",
       list(envir = qc.private)
)

assign("options",   qc.private$default_options,                 envir = qc.private)
assign("max_cores", parallel::detectCores(logical = TRUE) - 1,  envir = qc.private)
