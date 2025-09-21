#' Compute Gaussian log-likelihood
#'
#' Creates the underlying C++ block models and evaluates their Gaussian log-likelihood when applicable.
#'
#' @param ngme_model An object created by `ngme()`.
#'
#' @return A numeric scalar. Returns `0` if any replicate is non-Gaussian.
#' @export
compute_log_like <- function(ngme_model) {
  stopifnot(inherits(ngme_model, "ngme"))
  estimation_enabled <- attr(ngme_model, "estimation_enabled")
  if (!is.null(estimation_enabled) && !isTRUE(estimation_enabled)) {
    warning("Model was created with control_opt$estimation = FALSE; fit the model before computing log-likelihood.")
    return(NA_real_)
  }

  tryCatch(
    compute_log_like_cpp(ngme_model),
    error = function(err) {
      warning("Failed to compute log-likelihood: ", conditionMessage(err))
      NA_real_
    }
  )
}
