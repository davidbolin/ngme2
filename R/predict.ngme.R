#' prediction of ngme
#'
#' @param ngme ngme object
#' @return a list of outputs contains estimation of operator paramters, noise parameters
#' @export
#'
predict.ngme <- function(
  ngme,
  ...
) {
  linear_predictor <- double(ngme$W_sizes)

  AW_pred <- 0; AW_data <- 0
  for (i in seq_along(latents_in)) {
    W <- ngme$latents[[i]]$W
    AW_pred <- AW_pred + drop(latents_in[[i]]$A_pred %*% W)
    AW_data <- AW_data + drop(latents_in[[i]]$A %*% W)
  }

  # fixed effects. watch out! Xb could be double(0)
  X_pred <- split_data$X_NA;
  Xb_pred <- drop(X_pred %*% ngme$beta)
  Xb_data <- drop(X_data %*% ngme$beta)

  # ngme_response[split_data$index_NA] <- if (length(Xb_pred) == 0) AW_pred else AW_pred + Xb_pred
  #
  linear_predictor[split_data$index_NA]   <- if (length(Xb_pred) == 0) AW_pred else AW_pred + Xb_pred
  linear_predictor[split_data$index_data] <- if (length(Xb_data) == 0) AW_data else AW_data + Xb_data

  linear_predictor
}
