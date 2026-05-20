#' The 'qc' package.
#'
#' @description Frequentist and Bayesian process-capability tools for quality
#'   control workflows.
#'
#' @name qc-package
#' @aliases qc
#' @useDynLib qc, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#' @importFrom rstantools rstan_config
#' @importFrom RcppParallel RcppParallelLibs
#' @importFrom BayesTools is.prior is.prior.none is.prior.point is.prior.simple
#' @importFrom stats coef setNames
#' @importFrom utils tail
#'
#' @references
#' Stan Development Team (NA). RStan: the R interface to Stan. R package version 2.32.7. https://mc-stan.org
#'
"_PACKAGE"

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("metric", "region", "type", "x", "xmax", "xmin", "y", "yend"))
}
