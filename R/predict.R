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
#' @param type character, c("fe", "lp", "field1", "cv")
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
  # names <- c("fe", "field1", "field2")
  names <- if (type == "lp") c("fe", names(ngme$latents)) else type

  # ngme after estimation, has W, beta
  for (i in seq_along(names)) {
    name <- names[[i]]
    preds <- double(length(names))

    if (name == "fe") {
      # X_pred?
      if (is.null(ngme$X_pred)) {
        X_pred <- as.matrix(data) # from formula
      } else {
        X_pred <- ngme$X_pred
      }
      preds[[i]] <- drop(X_pred %*% ngme$beta)
    } else if (type %in% names(ngme$latents)) {
      model <- ngme$latents[[type]]
      # A_pred?
      if (is.null(model$A_pred))
        stop("A_pred not available")
      else
        A_pred <- model$A_pred

      preds[[i]] <- drop(A_pred %*% model$W)
    }
  }

  # lp case is just fe + A1 * W1 + A2 * W2

  # cv
  # compute MAE, MSE, CRPS, sCRPS
  sum(preds)
}

cv <- function(ngme) {

  # cv
  # compute MAE, MSE, CRPS, sCRPS
}