#' @noRd
getPrettyAxisBreaks <- function(x, ...) {
  force(x)
  UseMethod("getPrettyAxisBreaks", x)
}

#' @noRd
getPrettyAxisBreaks.numeric <- function(x, ...) {
  return(base::pretty(x, ...))
}

#' @noRd
getPrettyAxisBreaks.factor <- function(x, ...) {
  return(unique(x))
}

#' @noRd
getPrettyAxisBreaks.character <- function(x, ...) {
  return(unique(x))
}

#' @noRd
getPrettyAxisBreaks.default <- function(x, ...) {
  if (is.numeric(x))
    return(pretty(x, ...))
  else return(unique(x))
}
