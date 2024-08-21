#' Inference and prediction for mixed effects models with flexible non-Gaussian and Gaussian distributions.
#'
#' @name ngme2
#' @keywords internal 
"_PACKAGE"
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @importFrom stats simulate delete.response dnorm formula model.matrix rnorm sd terms terms.formula
#' @importFrom methods as
#' @importFrom utils str head tail modifyList
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab guides labs geom_hline
#' @importFrom rlang .data
#' @importFrom fmesher fm_mesh_1d fm_mesh_2d fm_basis fm_fem
#' @importFrom graphics hist
#' @importFrom stats median quantile as.formula dist
#' @useDynLib ngme2, .registration = TRUE
#' @name ngme2
NULL

.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion("ngme2")
   packageStartupMessage(
    "This is ngme2 of version ", version, "\n",
    "- See our homepage: https://davidbolin.github.io/ngme2 for more details."
  )
}
