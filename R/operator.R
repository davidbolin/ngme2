ngme_operator <- function(
  map,
  mesh,
  n_rep,
  model,
  K,
  h,
  theta_K = NULL,
  zero_trace = FALSE,
  symmetric = FALSE,
  ...
) {

  structure(
    list(
      map = map,
      mesh = mesh,
      n_rep = n_rep,
      model = model,
      K = K,
      h = h,
      theta_K = theta_K,
      n_theta_K = length(theta_K),
      zero_trace = zero_trace,
      symmetric = symmetric,
      ...
    ),
    class = "ngme_operator"
  )
}


#' Print ngme operator
#'
#' @param x ngme operator object
#' @param padding number of white space padding in front
#' @param ... ...
#'
#' @return a list (operator specifications)
#' @export
print.ngme_operator <- function(x, padding = 0, ...) {
  operator <- x
  pad_space <- paste(rep(" ", padding), collapse = "")
  pad_add4_space <- paste(rep(" ", padding + 4), collapse = "")

  model_name <- switch(operator$model,
    ar1 = "AR(1)",
    matern = "Matern",
    tp  = "Tensor product (with 2 sub-models)",
    bv  = "Bivariate model (with 2 sub-models)",
    iid = "IID model",
    rw1 = "Random walk (order 1)",
    rw2 = "Random walk (order 2)",
    ou  = "Ornstein-Uhlenbeck",
    "Unknown"
  )

  model_name <- paste("Model type: ", model_name, "\n", sep = "")
  cat(pad_space, model_name, sep="")

  parameter <- with(operator, switch(model,
    ar1 = cat(pad_add4_space, "alpha = ", format(ar1_th2a(theta_K), digits=3), "\n", sep=""),
    matern = if (length(theta_K) > 1)
      cat(pad_add4_space, "theta_K = ", format(theta_K, digits=3), "\n", sep="") else
      cat(pad_add4_space, "kappa = ", format(exp(theta_K), digits=3), "\n", sep=""),
    tp = {
      print(operator$first,  padding = padding + 4)
      print(operator$second, padding = padding + 4)
    },
    bv = {
      cat(pad_add4_space, "zeta = ", format(theta_K[1], digits=3), "\n", sep="")
      cat(pad_add4_space, "rho = ", format(theta_K[2], digits=3), "\n", sep="")
      print(operator$first,  padding = padding + 4)
      print(operator$second, padding = padding + 4)
    },
    cat(pad_add4_space, "No parameter.", "\n", sep="")
  ))

  invisible(operator)
}

# build_noise <- function(
#   noise,
#   operator
# ) {
#   stopifnot(inherits(operator, "ngme_operator"))
# }