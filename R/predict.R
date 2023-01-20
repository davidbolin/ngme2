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

# compute a list of criterions
# loc = all
# k-fold location at unknown
model_validation <- function(ngme, N=100, loc=NULL, test_at_loc=NULL) {
  # sampling Y by, Y = X beta + (block_A %*% block_W) + eps

  # AW_N[[1]] is concat(A1 W1, A2 W2, ..)

  # option 1. AW comes from 1 chain
  # turn into df of dim: n_obs * N
  AW_N <- sampling_cpp(ngme, n = N, posterior = TRUE)[["AW"]]
  AW_N <- as.data.frame(AW_N)
  names(AW_N) <- 1:N

  AW2_N <- sampling_cpp(ngme, n = N, posterior = TRUE)[["AW"]]
  AW2_N <- as.data.frame(AW2_N)
  names(AW_N) <- 1:N

  # option 2. AW comes from N chains
  # to-do

  fe <- with(ngme, as.numeric(X %*% beta))
  fe_N <- matrix(rep(fe, N), ncol=N, byrow=F)

  mn_N <- sapply(1:N, function(x) simulate(ngme$noise))
  mn2_N <- sapply(1:N, function(x) simulate(ngme$noise))

  mu_N <- fe_N + AW_N
  Y_N <- fe_N + AW_N + mn_N
  Y2_N <- fe_N + AW2_N + mn2_N

  # Now Y is of dim n_obs * N
  y_data <- ngme$Y; n_obs <- length(y_data)
  E3 <- E2 <- E1 <- double(length(y_data))
  for (i in 1:n_obs) {
    # turn row of df into numeric vector.
    yi <- as.numeric(Y_N[i, ])
    yi2 <- as.numeric(Y2_N[i, ])

    # estimate E(| Y_i - y_data |). y_data is observation
    E1[[i]] <- mean(abs(yi - y_data[i]))
    # estimate E(| Y_i - y_data |)
    E2[[i]] <- mean(abs(yi - yi2))
    # estimate E(| Y_i - y_data |^2)
    E3[[i]] <- mean((yi - y_data[i])^2)
  }

  # compute MSE, MAE, CRPS, sCRPS
  list(
    MAE = mean(E1),
    MSE = mean(E3),
    CRPS = mean(0.5 * E2 - E1),
    sCRPS = mean(-E2 / E1 - 0.5 * log(E2))
  )
}
