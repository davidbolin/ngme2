# This file contains function related to model predict function


#' Predict function of ngme2
#' predict using ngme after estimation
#'
#' @param object a ngme object
#' @param data a data.frame of covariates (used for fixed effects)
#' @param loc a list of the locations to make the prediction
#'  corresponding to each latent model
#'  vector or matrix (n * 2) for spatial coords
#' @param type what type of prediction, c("fe", "lp", "field1", "cv")
#' @param cv_type cross-validation, c("loo", "k-fold", "")
#' @param ... extra args
#'
#' @return a list of outputs contains estimation of operator paramters, noise parameters
#' @export
#'
predict.ngme <- function(
  object,
  data = NULL,
  loc = NULL,
  type = "lp",
  cv_type = NULL,
  ...
) {
  ngme <- object
  # maybe not necessary...
  # 0. length of output (loc 2d, loc 1d, )
  # if (!is.null(loc) && (is.matrix(loc) || is.data.frame(loc)))
  #   n_pred <- nrow(loc)
  # else if (!is.null(loc))
  #   n_pred <- length(loc)
  # else
  #   stop("not implement for loc=NULL yet")

  # 1. If provide loc, make A_pred for each model
  if (!is.null(loc)) {
    # make sure it's list
    if (!is.list(loc)) loc <- list(loc)

    # 1. Make A_pred at new loc for each latent model!!!
    for (i in seq_along(ngme$latents)) {
      mesh <- ngme$latents[[i]]$mesh
      if (!is.null(mesh) &&
          inherits(mesh, c("inla.mesh", "inla.mesh.1d"))) {
        # model using INLA mesh
        ngme$latents[[i]]$A_pred <- INLA::inla.spde.make.A(
          mesh = mesh,
          loc = loc[[i]]
        )
      } else if (!is.null(mesh)) {
        # time series model using
        # watch out! to-do
        ngme$latents[[i]]$A_pred <- ngme_ts_make_A(loc=loc[[i]])
      } else {
        stop("mesh is null, don't know how to make A")
      }
    }
  }

  # 2. Every latent has A_pred, now let's predict by Xbeta + AW
  # preds <- double(n_pred)
  preds <- 0

  names <- if (type == "lp") c("fe", names(ngme$latents)) else type
  # names <- c("fe", "field1", "field2")

  # ngme after estimation, has W, beta
  for (i in seq_along(names)) {
    name <- names[[i]]

    if (name == "fe") {
      # compute fe
      if (length(ngme$beta) > 0) {
        # make X_pred
        if (!is.null(data)) {
          X_pred <- as.matrix(data) # from formula
        } else {
          X_pred <- ngme$X_pred
        }
        # Add intercept if has
        # If only has intercept, then X is of dim 1*1, which is fine
        # (if we have further AW)
        if ((is.null(X_pred) && length(ngme$beta) == 1) ||
            (ncol(X_pred) == 1 + length(ngme$beta)))
          X_pred <- cbind(1, X_pred)

        preds <- preds + as.numeric(X_pred %*% ngme$beta)
      }
    } else if (name %in% names(ngme$latents)) {
      model <- ngme$latents[[name]]
      # A_pred?
      if (is.null(model$A_pred))
        stop("A_pred not available")
      else
        A_pred <- model$A_pred

    # browser()
      preds <- preds + as.numeric(A_pred %*% model$W)
    }
  }
  # lp case is just fe + A1 * W1 + A2 * W2
  stopifnot(length(preds) > 1)
  preds
}

# helper function to compute MSE, MAE, ..
compute_indices <- function(ngme, test_idx, N = 100) {
  stopifnot("idx wrong" = all(test_idx %in% seq_along(ngme$Y)))

  # 1. Make A1_pred, An_pred
  # A_pred_blcok  %*% block_W[[1 .. N]]
  # Y_N_pred <- AW_N_pred[1..N] + X_pred %*% beta + eps_N_pred
  # compute criterion (Y_N_pred, Y_data)

  # option 1. AW comes from 1 chain
  #   turn into df of dim: n_obs * N
  # option 2. AW comes from N chains
  #   to-do
  y_data <- ngme$Y[test_idx]
  n_obs <- length(y_data)

  new_model <- modify_ngme_with_idx_NA(ngme, idx_NA = test_idx)

  # A_pred_blcok <- [A1_pred .. An_pred]
  # extract A and cbind!
  A_pred_block <- Reduce(cbind, x = sapply(
    seq_along(new_model$latents),
    function(i) new_model$latents[[i]]$A_pred
  ))
  Ws_block <- sampling_cpp(ngme, n=N, posterior=TRUE)[["W"]]
  W2s_block <- sampling_cpp(ngme, n=N, posterior=TRUE)[["W"]]
  AW_N <- Reduce(cbind, sapply(Ws_block, function(W) A_pred_block %*% W))
  # AW_N <- as.data.frame(AW_N)
  # names(AW_N) <- 1:N

  AW2_N <- Reduce(cbind, sapply(W2s_block, function(W) A_pred_block %*% W))

  # sampling Y by, Y = X beta + (block_A %*% block_W) + eps
  # AW_N[[1]] is concat(A1 W1, A2 W2, ..)

  fe <- with(new_model, as.numeric(X_pred %*% beta))
  fe_N <- matrix(rep(fe, N), ncol=N, byrow=F)

  mn_N <- sapply(1:N, function(x)
    simulate(new_model$noise, n_noise = length(y_data)))
  mn2_N <- sapply(1:N, function(x)
    simulate(new_model$noise, n_noise = length(y_data)))
# browser()

  mu_N <- fe_N + AW_N
  Y_N <- fe_N + AW_N + mn_N
  Y2_N <- fe_N + AW2_N + mn2_N

  # Now Y is of dim n_obs * N
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

#' Compute the cross-validation for the ngme model
#'
#' @param ngme a ngme object
#' @param type character, cv type, in c("k-fold", "loo", "lpo")
#'  loo is leave-one-out, and lpo is leave-percent-out
#' @param k integer, if using "k-fold" cv
#' @param percent from 1 to 100, if using leave percent off cv ("lpo")
#' @param N integer, number of samplings (the higher, the better)
#' @param seed random seed
#' @param times run how many times (only for lpo type)
#' @param group group of indices of test set
#'  Can be the result of last CV function, or a list of indices
#' @param print print information along the process
#'
#' @return a list of MSE, MAE, CRPS, sCRPS
#' @export
cross_validation <- function(
  ngme,
  type = "k-fold",
  k = 5,
  N = 100,
  seed = 1,
  percent = 50,
  times = 10,
  group = NULL,
  print = FALSE
) {
  stopifnot(type %in% c("k-fold", "loo", "lpo"))

  # cut the group if not
  if (is.null(group)) {
    if (type == "k-fold") {
      # split idx into k
      idx <- seq_along(ngme$Y)
      folds <- cut(sample(idx), breaks = k, label = FALSE)
      group <- lapply(1:10, function(x) {which(folds == x, arr.ind = TRUE)})
    } else if (type == "loo") {
      return(cross_validation(ngme, "k-fold", k = length(ngme$Y), seed=seed))
    } else if (type == "lpo") {
      n_Y <- length(ngme$Y)
      for (i in 1:times) {
        group[[i]] <- sample(1:n_Y, size = (percent/100) * n_Y)
      }
    } else {
      stop("This cross-validation is not implemented!")
    }
  } else {
    if (!is.null(attr(ngme, "group")))
      group <- attr(ngme, "group")
    else if (!is.list(group))
      stop("Unkownn group")
  }

  # compute for each group
  crs <- NULL
  for (i in seq_along(group)) {
    crs[[i]] <- compute_indices(ngme, group[[i]], N=N)
if (print) {
  cat(paste("In group", i, ": \n"))
  print(as.data.frame(crs[[i]]))
  cat("\n")
}
  }

cat("The average of indices computed: \n")
  ret <- mean_list(crs)
  print(as.data.frame(ret))
  attr(ret, "group") <- group
  invisible(ret)
}