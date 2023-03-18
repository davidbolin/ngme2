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
#' @importFrom utils str
#' @importFrom ggplot2 ggplot aes geom_line xlab ylab guides labs
#' @importFrom rlang .data
#' @importFrom graphics hist
#' @importFrom stats median quantile
#' @useDynLib ngme2, .registration = TRUE
#' @name ngme2
NULL
