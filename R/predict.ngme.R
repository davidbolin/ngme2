# #' prediction of ngme
# #'
# #' @param ngme ngme object
# #' @return a list of outputs contains estimation of operator paramters, noise parameters
# #' @export
# #'
# predict.ngme <- function(
#   ngme,
#   ...
# ) {
#   linear_predictor <- double(ngme$W_sizes)

#   AW_pred <- 0; AW_data <- 0
#   for (i in seq_along(models_in)) {
#     W <- ngme$models[[i]]$W
#     AW_pred <- AW_pred + drop(models_in[[i]]$A_pred %*% W)
#     AW_data <- AW_data + drop(models_in[[i]]$A %*% W)
#   }

#   # fixed effects. watch out! Xb could be double(0)
#   X_pred <- split_data$X_NA;
#   Xb_pred <- drop(X_pred %*% ngme$feff)
#   Xb_data <- drop(X_data %*% ngme$feff)

#   # ngme_response[split_data$index_NA] <- if (length(Xb_pred) == 0) AW_pred else AW_pred + Xb_pred
#   #
#   linear_predictor[split_data$index_NA]   <- if (length(Xb_pred) == 0) AW_pred else AW_pred + Xb_pred
#   linear_predictor[split_data$index_data] <- if (length(Xb_data) == 0) AW_data else AW_data + Xb_data

#   linear_predictor
# }


# #' prediction of ngme given X_pred and A_pred
# #'
# #' @param ngme ngme object
# #' @param X_pred a matrix for covariate at unknown location
# #' @param A_pred a list of observation matrix for each process
# #'
# #' @return linear prediction at unknown locations
# #' @export
# #'
# predict.ngme <- function(
#   ngme,
#   X_pred,
#   A_pred,
#   sampling_iteration = 300
# ) {
#   # posterior sampling
#   ngme <- post_sampling(ngme, sampling_iteration)

#   linear_predictor <- X_pred %*% ngme$feff
#     for (latent in ngme$models)
#       linear_predictor <- linear_predictor + A_pred[i] %*% latent$W

#   linear_predictor
# }
