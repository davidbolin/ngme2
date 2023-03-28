# implementation of random effect structure

#' Build random effect structure
#' @param formula formula
#' @param data data
#' @param family distribution of random effect
#' @return a list of random effect structure
#' @export
randeff <- function(formula, data, effect_type, name) {
  stopifnot(effect_type %in% ngme_randeff_types())

  # Design Matrix and Covariance matrix
  B_reff <- model.matrix(formula, data)
  n_reff <- ncol(B_reff)

  Sigma  <- diag(n_reff)
  n_cov_params <- sum(1:n_reff)

  mix_var <- if (effect_type == "normal") noise_normal() else noise_nig()

  structure(
    list(
      formula      = formula,
      effect_type  = effect_type,
      B_reff       = B_reff,
      Sigma        = Sigma,
      n_cov_params = n_cov_params,
      n_params     = n_cov_params + n_reff + 1,
      name         = name,
      mix_var      = mix_var
    ),
    class = "randeff"
  )
}

#' Print ngme random effect
#'
#' @param x random effect object
#' @param padding number of white space padding in front
#' @param ... ...
#'
#' @return radom effect model
#' @export
print.randeff <- function(x, padding = 0, ...) {
  reff <- x
  pad_space <- paste(rep(" ", padding), collapse = "")
  pad_add4_space <- paste(rep(" ", padding + 4), collapse = "")

  if (is.null(reff)) {
    cat(pad_space); cat("Effect type - "); cat("NULL"); cat("\n")
  } else {
    cat(pad_space); cat("Effect type - "); cat(reff$effect_type); cat("\n")

    cat(pad_space); cat("Effect formula: ");
    cat(paste(as.character(reff$formula))); cat("\n")

    cat(pad_space); cat("Effect covariance: \n")
    cat(pad_space); print(reff$Sigma)
  }
  cat("\n")
  invisible(reff)
}