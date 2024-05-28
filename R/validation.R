# ------- CV -------
#' Compute the cross-validation for the ngme model
#' Perform cross-validation for ngme model
#' first into sub_groups (a list of target, and train data)
#'
#' @param ngme a ngme object, or a list of ngme object (if comparing multiple models)
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
#' @param times how many test cases (only for lpo type)
#' @param n_gibbs_samples number of gibbs samples
#' @param test_idx a list of indices of the data (which data points to be predicted) (only for custom type)
#' @param train_idx  a list of indices of the data (which data points to be used for re-sampling (not re-estimation)) (only for custom type)
#' @param keep_test_information logical, keep test information (pred_1, pred_2) in the return (as attributes), pred_1 and pred_2 are the prediction of the two chains
#'
#' @return a list of criterions: MSE, MAE, CRPS, sCRPS
#' @export
cross_validation <- function(
  ngme,
  type = "k-fold",
  seed = NULL,
  print = FALSE,
  N = 100,
  k = 5,
  percent = 50,
  times = 10,
  n_gibbs_samples = 50,
  test_idx = NULL,
  train_idx = NULL,
  keep_test_information = FALSE
) {
  if (inherits(ngme, "ngme")) ngme <- list(ngme)
  if (is.null(names(ngme))) names(ngme) <- paste("model", seq_along(ngme), sep = "_" )

  n_data <- attr(ngme[[1]], "fit")$n_data
  if (is.null(n_data)) stop("Please provide ngme object or a list of ngme object")

  if (!is.null(seed)) set.seed(seed)
  stopifnot(
    "type should be in c('k-fold', 'loo', 'lpo', 'custom')"
      = type %in% c("k-fold", "loo", "lpo", "custom")
  )

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
  # Alternative. do not distinguish between replicates?
  # But the internal mesh may not be the same for each replicate....

  # 2. loop over each test_idx and train_idx, and compute the criterion
  final_crs <- list(); pred_1 <- list(); pred_2 <- list();
  for (idx in seq_along(ngme)) {
    crs <- NULL
    for (i in seq_along(test_idx)) {
      test_idx[[i]] <- sort(test_idx[[i]])
      result <- compute_err_reps(
        ngme[[idx]],
        test_idx[[i]],
        train_idx[[i]],
        N=N,
        n_gibbs_samples = n_gibbs_samples,
        seed=seed
      )
      crs[[i]] <- result$score
      pred_1[[i]] <- result$pred_1
      pred_2[[i]] <- result$pred_2
      # Y_1[[i]] <- result$Y_1
      # Y_2[[i]] <- result$Y_2

      if (print) {
        cat(paste("In test_idx", i, ": \n"))
        print(as.data.frame(crs[[i]]))
        cat("\n")
      }
    }
    final_crs[[idx]] <- as.data.frame(mean_list(crs))
  }
  ret <- do.call(rbind, final_crs)

  rownames(ret) <- if (length(rownames(final_crs[[1]])) == 1)
      names(ngme)
    else {
      names_list = lapply(names(ngme), function(x) paste(x, rownames(final_crs[[1]]), sep = "_"))
      do.call(c, names_list)
    }

  # 3. take (weighted by train_data?) average over all groups

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

  print(ret)

  if (!keep_test_information) {
    invisible(ret)
  } else {
    structure(
      invisible(ret),
      train_idx = train_idx,
      test_idx = test_idx,
      pred_1 = pred_1,
      pred_2 = pred_2
      # Y_1 = Y_1,
      # Y_2 = Y_2
    )
  }
}

# helper function to dispatch over reps
compute_err_reps <- function(
  ngme,
  test_idx,
  train_idx,
  N = 100,
  n_gibbs_samples = 50,
  seed = NULL
) {
  stopifnot("Not a ngme object." = inherits(ngme, "ngme"))
  repls <- attr(ngme, "fit")$replicate
  uni_repl <- unique(repls)

  crs <- NULL; weight <- NULL; n_crs <- 0;
  pred_1 <- double(length = length(test_idx)); pred_2 <- double(length = length(test_idx))

  for (i in seq_along(uni_repl)) {
    data_idx_rep <- ngme$replicates[[i]]$data_idx
    bool_train_idx <- data_idx_rep %in% train_idx # current rep has train
    bool_test_idx  <- data_idx_rep %in% test_idx  # current rep has test

    # skip this replicate if no target or train data
    if (sum(bool_train_idx) == 0 || sum(bool_test_idx) == 0) next
    n_crs <- n_crs + 1
    result_1rep <- compute_err_1rep(
      ngme$replicates[[i]],
      bool_train_idx = bool_train_idx,
      bool_test_idx = bool_test_idx,
      N=N,
      n_gibbs_samples = n_gibbs_samples,
      seed=seed
    )
    crs[[n_crs]] <- result_1rep$scores

    # Assume which_idx_pred ordered (order test idx)
    which_idx_pred <- data_idx_rep[bool_test_idx]
    pred_1[test_idx %in% which_idx_pred] <- result_1rep$pred_1
    pred_2[test_idx %in% which_idx_pred] <- result_1rep$pred_2
    # Y_1[test_idx %in% which_idx_pred] <- result_1rep$Y_1
    # Y_2[test_idx %in% which_idx_pred] <- result_1rep$Y_2

    which_repl <- which(repls == uni_repl[i])
    weight <- c(weight, length(which_repl))
  }

  # take weighted average over replicates
  list(
    score = mean_list(crs, weight),
    pred_1 = pred_1,
    pred_2 = pred_2
    # Y_1 = Y_1,
    # Y_2 = Y_2
  )
}


# helper function to compute MSE, MAE, ... for each subset of target / data
# assume test_idx and train_idx belongs to same replicate
compute_err_1rep <- function(
  ngme_1rep,
  bool_test_idx,
  bool_train_idx,
  N = 100,
  n_gibbs_samples = 50,
  seed = NULL
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
  group_data <- ngme_1rep$group[bool_test_idx]
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

  if (is.null(seed)) seed <- Sys.time()

  # increase n_gibbs_samples
  ngme_1rep$control_ngme$n_gibbs_samples=n_gibbs_samples

  Ws <- sampling_cpp(
    ngme_1rep, 
    n=2*N, 
    n_burnin = N,
    posterior=TRUE, 
    seed=seed
  )[["W"]]
  
  Ws_block <- head(Ws, N); W2s_block <- tail(Ws, N)

# Note: Ws_block is a list of N realizations of W of current replicate
# Note: AW_N_1 is a matrix of n_test * N
  AW_N_1 <- Reduce(cbind, sapply(Ws_block, function(W) A_pred_block %*% W))
  AW_N_2 <- Reduce(cbind, sapply(W2s_block, function(W) A_pred_block %*% W))

  # sampling Y by, Y = X feff + (block_A %*% block_W) + eps
  # AW_N_1[[1]] is concat(A1 W1, A2 W2, ..)

  # generate fixed effect
  fe <- with(ngme_1rep, as.numeric(X_pred %*% feff))
  fe_N <- matrix(rep(fe, N), ncol=N, byrow=F)

  pred_N_1 <- fe_N + AW_N_1
  pred_N_2 <- fe_N + AW_N_2

  # simulate measurement noise
  # mn_N_1 <- sapply(1:N, function(x) simulate(noise_test_idx)[[1]])
  # mn_N_2 <- sapply(1:N, function(x) simulate(noise_test_idx)[[1]])
  # Y_N_1 <- pred_N_1 + mn_N_1
  # Y_N_2 <- pred_N_2 + mn_N_2

  pred <- 0.5*(rowMeans(as.matrix(pred_N_1)) + rowMeans(as.matrix(pred_N_2)))

  # Now Y is of dim n_obs * N
  E_pred_data <- E_pred_pred <- double(length(y_data))
  for (i in 1:n_obs) {
    # turn row of df into numeric vector.
    pred_1 <- as.numeric(pred_N_1[i, ])
    pred_2 <- as.numeric(pred_N_2[i, ])

    # estimate E(| X_i - y_data |). y_data is observation, X_i ~ predictive distribution at i
    E_pred_data[[i]] <- mean(abs(pred_1 - y_data[i]))

    # estimate E(| X_i - Y_i |) , X_i, Y_i ~ predictive distribution at i
    E_pred_pred[[i]] <- mean(abs(pred_1 - pred_2))
  }

  # compute MSE, MAE, CRPS, sCRPS within each group
  pred_each_group   <- split(pred, group_data)
  y_data_each_group <- split(y_data, group_data)

  CRPS  <- split(0.5 * E_pred_pred - E_pred_data, group_data)
  sCRPS <- split(
    -E_pred_data / E_pred_pred - 0.5 * log(E_pred_pred),
    group_data
  )

  # Compute MAE and MSE within each group
  MAE = MSE = double(length(pred_each_group))
  for (j in seq_along(pred_each_group)) {
    MAE[[j]] <- mean(abs(pred_each_group[[j]] - y_data_each_group[[j]]))
    MSE[[j]] <- mean((pred_each_group[[j]] - y_data_each_group[[j]])^2)
  }

  scores <- data.frame(
    MAE   = MAE,
    MSE   = MSE,
    negCRPS  = -sapply(CRPS, mean), # mean over 1:n_obs_test within each group
    negsCRPS = -sapply(sCRPS, mean) # same
  )

  # scores results and 2 predictions
  list(
    scores = scores,
    # Y_1 = rowMeans(as.matrix(Y_N_1)),
    # Y_2 = rowMeans(as.matrix(Y_N_2)),
    pred_1 = rowMeans(as.matrix(pred_N_1)),
    pred_2 = rowMeans(as.matrix(pred_N_2))
  )
}



# questions:
# 1. loop over replicates, and do CV for each replicate, and average
# 2. partition first, then do CV for each partition (with many replicates), and average