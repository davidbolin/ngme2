#' ngme2
#'
#' Inference and prediction for mixed effects models with flexible non-Gaussian and Gaussian distributions.
#'
#' @docType package
#' @author David Bolin <davidbolin@gmail.com>
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @importFrom stats simulate delete.response dnorm formula model.matrix rnorm sd terms terms.formula
#' @importFrom methods as
#' @importFrom utils str head tail modifyList
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab guides labs
#' @importFrom rlang .data
#' @importFrom fmesher fm_mesh_1d fm_mesh_2d fm_basis
#' @importFrom graphics hist
#' @importFrom stats median quantile as.formula
#' @useDynLib ngme2, .registration = TRUE
#' @name ngme2
NULL

.onLoad <- function(libname, pkgname) {
  version <- utils::packageVersion("ngme2")
   packageStartupMessage(
    "This is ngme2 of version ", version, "\n",
    "- See our homepage: https://davidbolin.github.io/ngme2 for more details."
  )
}
