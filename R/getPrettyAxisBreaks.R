# file copied fom https://github.com/jasp-stats/jaspGraphs/blob/8be6081bc7b1719a70538e563be853e11c4aedbc/R/getPrettyAxisBreaks.R#L1-L33
# to avoid a dependency on a GitHub package

#' @title Compute axis breaks
#' @param x the object to compute axis breaks for
#'
#' @param ... if x is numeric, this is passed to pretty
#' @details this is just a wrapper for pretty.
#'
#' @export
getPrettyAxisBreaks <- function(x, ...) {
  force(x)
  UseMethod("getPrettyAxisBreaks", x)
}

#' @export
getPrettyAxisBreaks.numeric <- function(x, ...) {
  return(base::pretty(x, ...))
}

#' @export
getPrettyAxisBreaks.factor <- function(x, ...) {
  return(unique(x))
}

#' @export
getPrettyAxisBreaks.character <- function(x, ...) {
  return(unique(x))
}

#' @export
getPrettyAxisBreaks.default <- function(x, ...) {
  if (is.numeric(x))
    return(pretty(x, ...))
  else return(unique(x))
}