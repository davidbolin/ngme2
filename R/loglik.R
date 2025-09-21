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
  was_estimated <- any(vapply(ngme_model$replicates, function(rep) !is.null(rep$log_likelihood), logical(1)))
  if (!was_estimated) {
    warning("Model has not been estimated yet. Run ngme() with control_opt$estimation = TRUE first.")
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
