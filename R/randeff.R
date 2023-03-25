# implementation of random effect structure

#' Build random effect structure
#' @param formula formula
#' @param data data
#' @param family distribution of random effect
#' @return a list of random effect structure
#' @export
randeff <- function(formula, data, family) {
  mat <- model.matrix(formula, data)

  # num. of covariance matrix
  parameter <- rep(0, ncol(mat)^2)

  # how to build covariance
  structure(
    list(
      family = family,
      mat = mat,
      parameter = parameter,
      n_params = length(parameter)
    ),
    class = "randeff"
  )
}