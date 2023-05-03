ngme_operator <- function(
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
    tp  = "Tensor product",
    bv  = "Bivariate model",
    "Unknown"
  )

  model_name <- paste("Model type: ", model_name, "\n", sep = "")
  cat(pad_space, model_name, sep="")

  parameter <- with(operator, switch(model,
    ar1 = cat(pad_add4_space, "alpha = ", format(ar1_th2a(theta_K), digits=3), "\n", sep=""),
    matern = cat(pad_add4_space, "theta_K = ", format(theta_K, digits=3), "\n", sep=""),
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
    "Unknown"
  ))

  # model_string <- model$model
  #   if (model_string == "rw" && model$rw_order==1) model_string <- "rw1"
  #   if (model_string == "rw" && model$rw_order==2) model_string <- "rw2"
  # cat(pad_space); cat("Ngme model: "); cat(model_string); cat("\n")

  # cat(pad_space); cat("Model parameters: \n")
  # params <- with(model, {
  #   switch(model,
  #     "ar1"         = paste0(pad_add4_space, ngme_format("K", theta_K, "ar1")),
  #     "ou"          = paste0(pad_add4_space, ngme_format("K", theta_K, "ou")),
  #     "re"          = {cat(paste0("  Covariance Matrix = \n"));
  #                      ngme_format("K", theta_K, "re", x$W_size); NULL},
  #     "matern"      = paste0(pad_add4_space, ngme_format("K", theta_K, "matern")),
  #     "rw"          = paste0(pad_add4_space, "No parameter."),
  #     "unkown"      = paste0(pad_add4_space, "No parameter."),
  #     "tp" = paste0(pad_add4_space, "left - ", left$model, ": ",
  #       ngme_format("K", left$theta_K, left$model), "\n",
  #       pad_add4_space, "right - ", right$model, ": ",
  #       ngme_format("K", right$theta_K, right$model)
  #     ),
  #     "bv" = paste0(pad_add4_space, "first - ", m1$model, ": ",
  #       ngme_format("K", m1$theta_K, m1$model), "\n",
  #       pad_add4_space, "second - ", m2$model, ": ",
  #       ngme_format("K", m2$theta_K, m2$model), "\n",
  #       "theta : ", theta_K[1], "\n",
  #       "rho : ", theta_K[2]
  #     ),
  #     paste0(pad_add4_space, "Not implemented yet!")
  #   )
  # })
  # cat(params);
  # cat("\n\n")

  # print.ngme_noise(model$noise, padding = padding)
  invisible(operator)
}

# build_noise <- function(
#   noise,
#   operator
# ) {
#   stopifnot(inherits(operator, "ngme_operator"))
# }