# f(idx, )

# prediction function for ngme block model
# predict.block_model <- function(
#   ngme,
#   locs = NULL,
#   type = "fe"
# ) {

# }


# predict for the latent model
predict.ngme_model <- function(
  object,
  ...
) {
  model <- object
  if (is.null(model$A_pred))
    stop("Make sure model has A_pred field!")
  if (is.null(model$W))
    stop("Make sure model has W field!")

  with(model,{
    A_pred %*% W
  })
}

#' Predict function of ngme2
#' predict using ngme after estimation
#'
#' @param ngme a ngme object
#' @param data a data.frame of covariates
#' @param locs the predicted locations
#' @param type character, c("fe", "lp", "filed1", "cv")
#' @param cv_type cross-validation, c("loo", "k-fold")
#'
#' @return a list of outputs contains estimation of operator paramters, noise parameters
#' @export
#'
predict <- function(
  ngme,
  data = NULL,
  locs = NULL,
  type = NULL,
  cv_type = NULL
) {

  if (type == "fe") {
    # X_pred?

    fe <- drop(X_pred %*% ngme$beta)
  }

  # cv
  # compute MAE, MSE, CRPS, sCRPS

}
