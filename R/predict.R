# This file contains function related to model predict function

#' Predict function of ngme2
#' predict using ngme after estimation
#'
#' @param object a ngme object
#' @param map a named list (or dataframe) of the locations to make the prediction
#' @param data a data.frame or matrix of covariates (used for fixed effects)
#'  names(loc) corresponding to the name each latent model
#'  vector or matrix (n * 2) for spatial coords
#' @param type what type of prediction, c("fe", "lp", <model_name>)
#' "fe" is fixed effect prediction
#' <model_name> is prediction of a specific model
#' "lp" is linear predictor (including fixed effect and all sub-models)
#' @param estimator what type of estimator, c("mean", "median", "mode", "quantile")
#' @param sampling_size size of posterior sampling
#' @param seed random seed
#' @param q quantile if using "quantile"
#' @param ... extra argument from 0 to 1 if using "quantile"
#'
#' @return a list of outputs contains estimation of operator paramters, noise parameters
#' @export
predict.ngme <- function(
  object,
  map,
  data = NULL,
  type = "lp",
  estimator = c("mean", "sd", "5quantile", "95quantile", "median", "mode"),
  sampling_size = 30,
  q = NULL,
  seed = Sys.time(),
  ...
) {
  fm <- attr(object, "fitting")$formula
  ngme <- object$replicate[[1]]
  stopifnot(sampling_size > 0)
  samples_W <- sampling_cpp(ngme, n=sampling_size, posterior=TRUE, seed=seed)[["W"]]

  ret <- NULL
  for (estimator in estimator) {
    post_W <- switch(estimator,
      "mean"      = mean_list(samples_W),
      "median"    = apply(as.data.frame(samples_W), 1, median),
      "sd"        = apply(as.data.frame(samples_W), 1, sd),
      "mode"      = apply(as.data.frame(samples_W), 1, emprical_mode),
      "5quantile" = apply(as.data.frame(samples_W), 1, function(x) {quantile(x, 0.05)}),
      "95quantile" = apply(as.data.frame(samples_W), 1, function(x) {quantile(x, 0.95)}),
      "quantile"  = {
        stopifnot("please provide quantile argument q between 0 to 1"
          = !is.null(q) && length(q) == 1 && q > 0 && q < 1)
        apply(as.data.frame(samples_W), 1, function(x) {quantile(x, q)})
      },
      stop("No such estimator available")
    )

    # update post W (notice here W is concated)
    j <- 1
    for (i in seq_along(ngme$models)) {
      sz <- ngme$models[[i]]$W_size
      ngme$models[[i]][[estimator]] <- post_W[j:(j + sz - 1)]
      j <- j + sz
    }

    if (!is.null(map)) {
      stopifnot(
        "map should be a named list (name for each model)"
          = is.list(map) && !is.null(names(map))
      )
      names <- names(map)
      stopifnot(length(names) == length(ngme$models))

      AW <- list()
      for (i in seq_along(ngme$models)) {
        loc <- map[[ngme$models[[i]]$name]]
    if (is.null(loc)) stop("The loction for model ", ngme$models[[i]]$name, " is not provided")
        loc <- as.matrix(loc)
        AW[[ngme$models[[i]]$name]] <- with(ngme$models[[i]], {
          mesh <- operator$mesh
          W <- ngme$models[[i]][[estimator]]
          A <- fmesher::fm_basis(mesh, loc = loc)
          if (model == "bv") A <- Matrix::bdiag(A, A)
          stopifnot(ncol(A) == length(W))
          as.numeric(A %*% W)
        })
      }
    }

    # lp case is just fe + A1 * W1 + A2 * W2
    # e.g. names <- c("fe", "field1", "field2")
    type_names <- if (type == "lp") c("fe", names(ngme$models)) else type

    preds <- 0
    for (i in seq_along(type_names)) {
      name <- type_names[[i]]
      if (name == "fe" && length(ngme$feff) > 0) {
          X_pred <- if (is.null(data) && attr(terms(fm), "intercept")) {
            matrix(1, nrow = length(AW[[1]]), ncol = 1)
          } else {
            # build plain_fm
            tf <- terms.formula(fm, specials = c("f"))
            terms <- attr(tf, "term.labels")
            intercept <- attr(tf, "intercept")
            spec_order <- attr(tf, "specials")$f - 1
            fixf <- if (length(spec_order) == 0) terms else terms[-spec_order]
            plain_fm_str <- paste("~", intercept, paste(c("", fixf), collapse = " + "))
            plain_fm <- formula(plain_fm_str)
            model.matrix(plain_fm, data = data)
          }
          preds <- preds + as.numeric(X_pred %*% ngme$feff)
      } else if (name %in% names(ngme$models)) {
        preds <- preds + AW[[name]]
      }
    }
    ret[[estimator]] <- preds
  }
  ret
}

# ------- CV -------

# helper function to compute MSE, MAE, ... for each subset of target / data
# assume test_idx and train_idx belongs to same replicate
compute_err_1rep <- function(
  ngme_1rep,
  bool_test_idx,
  bool_train_idx,
  N = 100,
  seed=Sys.time()
) {
  stopifnot(
    "bool_<..>_idx should be a logical vector" =
      is.logical(bool_test_idx) && is.logical(bool_train_idx)
  )
  # if test_idx and train_idx overlap, warning message
  if (sum(bool_test_idx & bool_train_idx) > 0) {
    warning("Notice that test_idx and train_idx overlap!")
  }

  # 1. Make A1_pred, An_pred
  # A_pred_blcok  %*% block_W[[1 .. N]]
  # Y_N_pred <- AW_N_pred[1..N] + X_pred %*% feff + eps_N_pred
  # compute criterion (Y_N_pred, Y_data)

  # option 1. AW comes from 1 chain
  #   turn into df of dim: n_obs * N
  # option 2. AW comes from N chains
  #   to-do

  # Subset noise[test_idx, ] for test location
  y_data <- ngme_1rep$Y[bool_test_idx]
  n_obs <- length(y_data)
  X_pred <- ngme_1rep$X[bool_test_idx,, drop=FALSE]
  noise_test_idx <- subset_noise(ngme_1rep$noise, sub_idx = bool_test_idx)

  # Subset noise, X, Y in train location
  ngme_1rep$X <- ngme_1rep$X[bool_train_idx,, drop=FALSE]
  ngme_1rep$Y <- ngme_1rep$Y[bool_train_idx]
  ngme_1rep$noise <- subset_noise(ngme_1rep$noise, sub_idx = bool_train_idx)

  # Subset A for test and train location
  A_preds <- list();
  for (i in seq_along(ngme_1rep$models)) {
    A_preds[[i]] <- ngme_1rep$models[[i]]$A[bool_test_idx, ,drop=FALSE]
    ngme_1rep$models[[i]]$A <- ngme_1rep$models[[i]]$A[bool_train_idx,,drop=FALSE]
  }

  # A_pred_blcok <- [A1_pred .. An_pred]
  # extract A and cbind!
  A_pred_block <- Reduce(cbind, x = A_preds)
# print(A_pred_block)

  Ws <- sampling_cpp(ngme_1rep, n=2*N, posterior=TRUE, seed=seed)[["W"]]
  Ws_block <- head(Ws, N); W2s_block <- tail(Ws, N)

  AW_N <- Reduce(cbind, sapply(Ws_block, function(W) A_pred_block %*% W))
  AW2_N <- Reduce(cbind, sapply(W2s_block, function(W) A_pred_block %*% W))

  # sampling Y by, Y = X feff + (block_A %*% block_W) + eps
  # AW_N[[1]] is concat(A1 W1, A2 W2, ..)

  # generate fixed effect
  fe <- with(ngme_1rep, as.numeric(X_pred %*% feff))
  fe_N <- matrix(rep(fe, N), ncol=N, byrow=F)

  # simulate measurement noise
  mn_N  <- sapply(1:N, function(x) simulate(noise_test_idx))
  mn2_N <- sapply(1:N, function(x) simulate(noise_test_idx))

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
#' Perform cross-validation for ngme model
#' first into sub_groups (a list of target, and train data)
#'
#' @param ngme a ngme object
#' @param type character, in c("k-fold", "loo", "lpo", "custom")
#' k-fold is k-fold cross-validation, provide \code{k}
#' loo is leave-one-out,
#' lpo is leave-percent-out, provide \code{percent} from 1 to 100
#' custom is user-defined group, provide \code{target} and \code{data}
#' @param seed random seed
#' @param N integer, number of MC sampling (i.e. sampling N process and average)
#' @param k integer (only for k-fold type)
#' @param print print information during computation
#' @param percent from 1 to 100 (only for lpo type)
#' @param times run how many times (only for lpo type)
#' @param test_idx a list of indices of the data (which data points to be predicted) (only for custom type)
#' @param train_idx  a list of indices of the data (which data points to be used for re-sampling (not re-estimation)) (only for custom type)
#'
#' @return a list of criterions: MSE, MAE, CRPS, sCRPS
#' @export
cross_validation <- function(
  ngme,
  type = "k-fold",
  seed = Sys.time(),
  print = FALSE,
  N = 100,
  k = 5,
  percent = 50,
  times = 10,
  test_idx = NULL,
  train_idx = NULL
) {
  stopifnot(
    "type should be in c('k-fold', 'loo', 'lpo', 'custom')"
      = type %in% c("k-fold", "loo", "lpo", "custom"),
    "ngme is a ngme object"
      = inherits(ngme, "ngme")
  )
  n_data <- attr(ngme, "fitting")$n_data

  # 1. compute indices of tartget and train if not custom type
  if (type == "k-fold") {
    # split idx into k
    idx <- seq_len(n_data)
    folds <- cut(sample(idx), breaks = k, label = FALSE)
    test_idx <- lapply(1:k, function(x) {which(folds == x, arr.ind = TRUE)})
    train_idx <- lapply(1:k, function(x) {which(folds != x, arr.ind = TRUE)})
  } else if (type == "loo") {
    return(cross_validation(ngme, "k-fold", k = n_data, seed=seed))
  } else if (type == "lpo") {
    for (i in 1:times) {
      test_idx[[i]] <- sample(1:n_data, size = (percent/100) * n_data)
      train_idx[[i]] <- setdiff(1:n_data, test_idx[[i]])
    }
  } else {
    # check if test_idx and train_idx is provided and of same length
    stopifnot(
      "test_idx and train_idx should be provided"
        = !is.null(test_idx) && !is.null(train_idx),
      "test_idx and train_idx should be a list"
        = is.list(test_idx) && is.list(train_idx),
      "test_idx and train_idx should be of same length"
        = length(test_idx) == length(train_idx)
    )
  }

  # 2. loop over each test_idx and train_idx, and compute the criterion
  crs <- NULL
  for (i in seq_along(test_idx)) {
    crs[[i]] <- compute_err_reps(
      ngme,
      test_idx[[i]],
      train_idx[[i]],
      N=N,
      seed=seed
    )
    if (print) {
      cat(paste("In test_idx", i, ": \n"))
      print(as.data.frame(crs[[i]]))
      cat("\n")
    }
  }

  # 3. take (weighted by train_data?) average over all groups
  ret <- mean_list(crs)

  # weights <- sapply(ngme$replicates, function(x) length(x$Y))
  # weights <- weights / sum(weights)
  # ret <- list()
  # for (i in seq_along(ngme$replicates)) {
  #   ret[[i]] <- cross_validation(ngme$replicates[[i]], type=type, k=k, N=N, percent=percent,
  #   times=times, test_idx=test_idx, print=print, seed=seed)
  # }
  # ret <- mean_list(ret, weights=weights)
  # cat("\n")
  # cat("The final result averaged over replicates: \n")

  print(as.data.frame(ret))
  return(invisible(ret))
}

# helper function to dispatch over reps
compute_err_reps <- function(
  ngme, test_idx, train_idx, N=100, seed=Sys.time()
) {
  stopifnot("Not a ngme object." = inherits(ngme, "ngme"))
  repls <- attr(ngme, "fitting")$replicate
  uni_repl <- unique(repls)

  crs <- NULL; weight <- NULL; n_crs <- 0
  for (i in seq_along(uni_repl)) {
    data_idx <- ngme$replicates[[i]]$data_idx
    bool_train_idx <- data_idx %in% train_idx
    bool_test_idx  <- data_idx %in% test_idx

    # skip this replicate if no target or train data
    if (sum(bool_train_idx) == 0 || sum(bool_test_idx) == 0) next
    n_crs <- n_crs + 1
    crs[[n_crs]] <- compute_err_1rep(
      ngme$replicates[[i]],
      bool_train_idx = bool_train_idx,
      bool_test_idx = bool_test_idx,
      N=N,
      seed=seed
    )

    which_repl <- which(repls == uni_repl[i])
    weight <- c(weight, length(which_repl))
  }

  # take weighted average over replicates
  mean_list(crs, weight)
}


# questions:
# 1. loop over replicates, and do CV for each replciate, and average
# 2. partition first, then do CV for each partition (with many replicates), and average