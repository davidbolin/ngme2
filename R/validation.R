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
#' @param N_sim integer, number of simulations (e.g., estimate MAE, MSE, .. N times)
#' @param k integer (only for k-fold type)
#' @param print print information during computation
#' @param percent how many percent for testing? from 0 to 1 (for lpo type)
#' @param times how many test cases (only for lpo type)
#' @param transform a function to transform the data (e.g., log, exp, ...)
#' e.g., the MAE will be computed as |transform(Y) - transform(Y_pred)|
#' @param n_gibbs_samples number of gibbs samples of latent process, used for computing CRPS, sCRPS
#' @param n_burnin number of burnin
#' @param test_idx a list of indices of the data (which data points to be predicted) (only for custom type)
#' @param train_idx  a list of indices of the data (which data points to be used for re-sampling (not re-estimation)) (only for custom type)
#' @param keep_pred logical, keep test information (pred_1, pred_2) in the return (as attributes), pred_1 and pred_2 are the prediction of the two chains
#' @param thining_gap integer, the gap between samples for thinning, if 0, then no thinning, if 1, then keep 50% of the samples for CRPS, sCRPS, etc.
#' @param parallel logical, run in parallel mode
#' @param cores_layer1 integer, number of cores for the first layer (over testing samples)
#' @param cores_layer2 integer, number of cores for the second layer (over computing scores for N_sim simulations)
#' @return 
#'  1. mean of N_sim estimations of 4 criterions: MSE, MAE, CRPS, sCRPS
#'  2. standard deviation of N_sim estimations of 4 criterions: MSE, MAE, CRPS, sCRPS
#' @export
cross_validation <- function(
  ngme,
  type = "k-fold",
  seed = NULL,
  print = FALSE,
  N_sim = 5,
  n_gibbs_samples = 500,
  n_burnin = 100,
  k = 5,
  percent = 0.2,
  times = 10,
  transform = identity,
  test_idx = NULL,
  train_idx = NULL,
  keep_pred = FALSE,
  parallel = FALSE,
  thining_gap = 1, # Used for computing CRPS, sCRPS, the gap between samples for thinning, if 0, then no thinning, if 1, then keep 50% of the samples for CRPS, sCRPS, etc.
  # merge_replicates = FALSE, # remove this option
  cores_layer1 = if (parallel) min(parallel::detectCores(), 2) else 1,  # Limit to 2 cores for safety
  cores_layer2 = if (parallel) min(parallel::detectCores(), 2) else 1   # Limit to 2 cores for safety
) {
  merge_replicates = FALSE

  if (!requireNamespace("parallel", quietly = TRUE)) {
    message("Parallel package not available. Running in serial mode. You can install `parallel` package to speed up the computation.")
    parallel <- FALSE
  }
  
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
    stopifnot(
      "percent should be between 0 and 1" = percent > 0 && percent < 1,
      "times should be a positive integer" = times > 0
    )
    for (i in 1:times) {
      test_idx[[i]] <- sample(1:n_data, size = percent * n_data)
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
  final_mean <- final_sd <- list(); 
  pred_1 <- list(); pred_2 <- list();
  Y_1 <- list(); Y_2 <- list();

  compute_err <- if (merge_replicates) compute_err_merged_reps else compute_err_reps

  for (idx in seq_along(ngme)) {
    if (print) cat(paste0("Model ", names(ngme)[[idx]], ": \n\n"))
    scores <- sd_scores <- NULL

    # loop over each test_idx and train_idx
    if (parallel && requireNamespace("parallel", quietly = TRUE)) {
      if (print) cat("Running in parallel mode. \n")
      
      ngme_list <- list()
      for (i in 1:length(test_idx)) {
        # ngme_list[[i]] <- rlang::duplicate(ngme[[idx]])
        ngme_list[[i]] <- ngme[[idx]]
      }

      # Add error handling for parallel execution
      tryCatch({
        ret = parallel::mclapply(
          seq_along(test_idx),
          function(i) {
            tryCatch({
              result <- compute_err(
                ngme_list[[i]],
                test_idx[[i]],
                train_idx[[i]],
                N_sim=N_sim,
                n_gibbs_samples = n_gibbs_samples,
                seed=seed,
                keep_pred=keep_pred,
                parallel = TRUE,
                transform = transform,
                num_cores = cores_layer2,
                thining_gap = thining_gap
              )
              if (print) {
                cat(paste("In test batch", i, ": \n"))
                print(as.data.frame(result$scores))
                cat("\n")
              }
              result
            }, error = function(e) {
              warning(paste("Error in test batch", i, ":", e$message))
              NULL
            })
          },
          mc.cores = cores_layer1
        )
        
        # Filter out NULL results
        ret <- ret[!sapply(ret, is.null)]
        
        if (length(ret) == 0) {
          warning("All parallel attempts failed, falling back to sequential computation")
          parallel <- FALSE
        } else {
          scores <- lapply(ret, function(x) x$scores)
          sd_scores <- lapply(ret, function(x) x$sd_scores)
        }
      }, error = function(e) {
        warning(paste("Error in parallel execution:", e$message))
        parallel <- FALSE
      })
    }
    
    # Fall back to sequential if parallel failed or was not requested
    if (!parallel || length(scores) == 0) {
      for (i in seq_along(test_idx)) {
        tryCatch({
          result <- compute_err(
            ngme[[idx]],
            test_idx[[i]],
            train_idx[[i]],
            N_sim=N_sim,
            n_gibbs_samples = n_gibbs_samples,
            seed=seed,
            keep_pred=keep_pred,
            parallel = FALSE,
            transform = transform,
            thining_gap = thining_gap
          )
          scores[[i]] <- result$scores
          sd_scores[[i]] <- result$sd_scores

          if (print) {
            cat(paste("In test batch", i, ": \n"))
            print(as.data.frame(scores[[i]]))
            cat("\n")
          }
        }, error = function(e) {
          warning(paste("Error in sequential computation for test batch", i, ":", e$message))
        })
      }
    }
    
    # Check if we have any valid scores
    if (length(scores) == 0) {
      warning(paste("Failed to compute any valid scores for model", names(ngme)[[idx]]))
      next
    }

    final_mean[[idx]] <- as.data.frame(mean_list(scores))
    final_sd[[idx]] <- as.data.frame(mean_list(sd_scores))
  }
  
  # Check if we have any results
  if (length(final_mean) == 0) {
    stop("Failed to compute any valid scores for any model")
  }
  
  # convert to data.frame
  final_mean <- do.call(rbind, final_mean)
  final_sd <- do.call(rbind, final_sd)

  if (nrow(scores[[1]]) == 1) {
    # univariate model
    rownames(final_mean) <- names(ngme) 
  } else {
    # bivariate model
    names_list = lapply(names(ngme), function(x) paste(x, rownames(final_mean)[1:2], sep = "_"))
    rownames(final_mean) <- do.call(c, names_list)
  }

  rownames(final_sd) <- rownames(final_mean)

  ret = list(
    mean.scores = final_mean,
    sd.scores = final_sd
  )
  
  if (!keep_pred) {
    invisible(ret)
  } else {
    structure(
      invisible(ret),
      train_idx = train_idx,
      test_idx = test_idx,
      pred_1 = pred_1,
      pred_2 = pred_2,
      Y_1 = Y_1,
      Y_2 = Y_2
    )
  }
}

# Compute error with merged 1 replicate using helper function merge_replicates
compute_err_merged_reps <- function(
  ngme,
  test_idx,
  train_idx,
  N_sim,
  n_gibbs_samples = 100,
  n_burnin = 100,
  seed = NULL,
  keep_pred = FALSE,
  parallel = TRUE,
  transform = identity,
  num_cores = 1,
  thining_gap = 0
) {
  if (is.null(seed)) seed <- Sys.time()
  stopifnot("Not a ngme object." = inherits(ngme, "ngme"))

  test_idx <- sort(test_idx)
  repls <- attr(ngme, "fit")$replicate
  uni_repl <- unique(repls)

  scores <- sd_scores <- weight <- NULL; 
  n_scores <- 0;
  pred_1 <- double(length = length(test_idx)); pred_2 <- double(length = length(test_idx))
  Y_1 <- double(length = length(test_idx)); Y_2 <- double(length = length(test_idx))

  ret <- merge_replicates(ngme, train_idx, test_idx)
  merged_rep <- ret$merged_rep
  
  scores <- list()
  for (sim_n in 1:N_sim) {
    scores[[sim_n]] <- compute_scores(
      merged_rep, 
      n_gibbs_samples, 
      n_burnin, 
      seed + sim_n, 
      ret$test_A_block, 
      ret$test_noise, 
      ret$test_Y, 
      ret$test_group, 
      ret$test_X, 
      transform,
      thining_gap = thining_gap
    )
  }

  # compute mean and sd of scores
  n_group = length(levels(ret$test_group))
  array_3d <- array(unlist(scores), 
    dim = c(n_group, 4, length(scores)))

  mean_array <- apply(array_3d, c(1, 2), mean)
  sd_array <- apply(array_3d, c(1, 2), sd)

  colnames(mean_array) <- colnames(sd_array) <- names(scores[[1]])
  rownames(mean_array) <- rownames(sd_array) <- rownames(scores[[1]])

  list(
    scores = mean_array,
    sd_scores = sd_array
  )
}



# Compute error with multiple replicates using helper function compute_err_1rep
compute_err_reps <- function(
  ngme,
  test_idx,
  train_idx,
  N_sim,
  n_gibbs_samples = 100,
  n_burnin = 100,
  seed = NULL,
  keep_pred = FALSE,
  parallel = TRUE,
  transform = identity,
  num_cores = 1,
  thining_gap = 1
) {
  test_idx <- sort(test_idx)
  stopifnot("Not a ngme object." = inherits(ngme, "ngme"))
  repls <- attr(ngme, "fit")$replicate
  uni_repl <- unique(repls)

  scores <- sd_scores <- weight <- NULL; 
  n_scores <- 0;
  pred_1 <- double(length = length(test_idx)); pred_2 <- double(length = length(test_idx))
  Y_1 <- double(length = length(test_idx)); Y_2 <- double(length = length(test_idx))

  for (i in seq_along(uni_repl)) {
    data_idx_rep <- ngme$replicates[[i]]$data_idx
    bool_train_idx <- data_idx_rep %in% train_idx # current rep has train
    bool_test_idx  <- data_idx_rep %in% test_idx  # current rep has test

    # skip this replicate if no target or train data
    if (sum(bool_train_idx) == 0 || sum(bool_test_idx) == 0) next
    n_scores <- n_scores + 1

    ngme_1rep = ngme$replicates[[i]]
    result_1rep <- compute_err_1rep(
      ngme_1rep,
      bool_train_idx = bool_train_idx,
      bool_test_idx = bool_test_idx,
      N_sim=N_sim,
      n_gibbs_samples = n_gibbs_samples,
      seed=seed,
      keep_pred=keep_pred,
      parallel = parallel,
      transform = transform,
      num_cores = num_cores,
      thining_gap = thining_gap
    )
    scores[[n_scores]] <- result_1rep$mean_scores
    sd_scores[[n_scores]] <- result_1rep$sd_scores

    # Assume which_idx_pred ordered (order test idx)
    which_idx_pred <- data_idx_rep[bool_test_idx]
    
    if (keep_pred) {
      pred_1[test_idx %in% which_idx_pred] <- result_1rep$pred_1
      pred_2[test_idx %in% which_idx_pred] <- result_1rep$pred_2
      Y_1[test_idx %in% which_idx_pred] <- result_1rep$Y_1
      Y_2[test_idx %in% which_idx_pred] <- result_1rep$Y_2
    }

    which_repl <- which(repls == uni_repl[i])
    weight <- c(weight, length(which_repl))
  }

  # take weighted average over replicates
  list(
    scores = mean_list(scores, weight),
    sd_scores = mean_list(sd_scores, weight),
    pred_1 = pred_1,
    pred_2 = pred_2,
    Y_1 = Y_1,
    Y_2 = Y_2
  )
}


# helper function to compute MSE, MAE, ... for each subset of target / data
# assume test_idx and train_idx belongs to same replicate
compute_err_1rep <- function(
  ngme_1rep,
  bool_test_idx,
  bool_train_idx,
  N_sim,
  n_gibbs_samples = 50,
  n_burnin = 10,
  seed = NULL,
  keep_pred = FALSE,
  parallel = TRUE,
  transform = identity,
  num_cores = 1,  # Default to 1 core to avoid potential issues
  thining_gap = 1
) {
  stopifnot(
    "bool_<..>_idx should be a logical vector" =
      is.logical(bool_test_idx) && is.logical(bool_train_idx)
  )
  # if test_idx and train_idx overlap, warning message
  if (sum(bool_test_idx & bool_train_idx) > 0) {
    warning("Notice that test_idx and train_idx overlap!")
  }

# Since we revert the order of Y, now we need to 
# revert the train and test idx to match
# NOT REALLY, I DID IT in the outside function!!!

  # Subset noise[test_idx, ] for test location
  y_data <- ngme_1rep$Y[bool_test_idx]
  group_data <- ngme_1rep$group[bool_test_idx]
  X_pred <- ngme_1rep$X[bool_test_idx,, drop=FALSE]
  noise_test_idx <- subset_noise(
    ngme_1rep$noise, sub_idx = bool_test_idx, compute_corr = TRUE
  )

  # Subset noise, X, Y in train location
  ngme_1rep$X <- ngme_1rep$X[bool_train_idx,, drop=FALSE]
  ngme_1rep$Y <- ngme_1rep$Y[bool_train_idx]
  
  ngme_1rep$noise <- subset_noise(
    ngme_1rep$noise, sub_idx = bool_train_idx, compute_corr = TRUE
  )

  # Subset A for test and train location
  A_preds <- list();
  for (i in seq_along(ngme_1rep$models)) {
    A_preds[[i]] <- ngme_1rep$models[[i]]$A[bool_test_idx, ,drop=FALSE]
    ngme_1rep$models[[i]]$A <- ngme_1rep$models[[i]]$A[bool_train_idx,,drop=FALSE]
  }

  # A_pred_blcok <- [A1_pred .. An_pred]
  # extract A and cbind!
  A_pred_block <- Reduce(cbind, x = A_preds)

  if (is.null(seed)) seed <- as.integer(Sys.time())

  scores <- pred_1 <- pred_2 <- Y_1 <- Y_2 <- list()

  # Fall back to sequential if parallel failed or was not requested
  # if (!parallel || length(scores) == 0) {
    for (nn in 1:N_sim) {
      tryCatch({
        scores[[nn]] <- compute_scores(
          ngme_1rep, n_gibbs_samples, n_burnin, seed+nn, A_pred_block, 
          noise_test_idx, y_data, group_data, X_pred, transform, thining_gap
        )
      }, error = function(e) {
        warning(paste("Error in sequential computation:", e$message))
      })
    }
  # }
  
  # Check if we have any valid scores
  if (length(scores) == 0) {
    stop("Failed to compute any valid scores")
  }

  # compute extra pred_N?
  if (keep_pred) {
    tmp = compute_pred_N(
      ngme_1rep, n_gibbs_samples, n_burnin, seed, A_pred_block, noise_test_idx, y_data, group_data, X_pred
    )
    pred_1 <- tmp$pred_N_1; pred_2 <- tmp$pred_N_2
    Y_1 <- tmp$Y_N_1; Y_2 <- tmp$Y_N_2
  }

  n_group = length(levels(group_data))
  # compute mean and sd
  array_3d <- array(unlist(scores), 
    dim = c(n_group, 4, length(scores)))

  mean_array <- apply(array_3d, c(1, 2), mean)
  sd_array <- apply(array_3d, c(1, 2), sd)

  colnames(mean_array) <- colnames(sd_array) <- names(scores[[1]])
  rownames(mean_array) <- rownames(sd_array) <- rownames(scores[[1]])

  # scores results and 2 predictions
  list(
    mean_scores = mean_array,
    sd_scores = sd_array,
    Y_1 = Y_1,
    Y_2 = Y_2,
    pred_1 = pred_1,
    pred_2 = pred_2
  )
}



# questions:
# 1. loop over replicates, and do CV for each replicate, and average
# 2. partition first, then do CV for each partition (with many replicates), and average


# helper function to compute the scores
compute_scores <- function(
  ngme_1rep, 
  n_gibbs_samples, 
  n_burnin, 
  seed, 
  A_pred_block, 
  test_noise, 
  y_data, 
  group_data, 
  X_pred, 
  transform,
  thining_gap # keep 50% of the thinning samples for CRPS, sCRPS
) {
  # Add error handling for the C++ sampling_cpp call
  tryCatch({
    # Ensure seed is an integer and within valid range
    seed_int <- as.integer(abs(seed) %% 2147483647)
    
    # Add diagnostic message
    # message(sprintf("Calling sampling_cpp with n=%d, n_burnin=%d, seed=%d", 
    #                n_gibbs_samples * 2, n_burnin, seed_int))
    
    # Call sampling_cpp with error handling
    Ws_result <- tryCatch({
      sampling_cpp(
        ngme_1rep, 
        n = n_gibbs_samples * 2, 
        n_burnin = n_burnin,
        posterior = TRUE, 
        seed = seed_int
      )
    }, error = function(e) {
      message("Error in sampling_cpp: ", e$message)
      # Try again with different seed
      message("Retrying with different seed...")
      sampling_cpp(
        ngme_1rep, 
        n = n_gibbs_samples * 2, 
        n_burnin = n_burnin,
        posterior = TRUE, 
        seed = seed_int
      )
    })
    
    Ws <- Ws_result[["W"]]
    # Ws <- Ws_result[["cond_W"]]
    
    if (is.null(Ws) || length(Ws) == 0) {
      stop("sampling_cpp returned NULL or empty result")
    }
    
    # Rest of the function remains the same
    Ws_block  <- head(Ws, n_gibbs_samples)
    W2s_block <- tail(Ws, n_gibbs_samples)
    
    # keep 50% of the thinning samples for CRPS, sCRPS
    Ws_block_thin <- Ws_block[seq(1, n_gibbs_samples, by=thining_gap+1)]
    W2s_block_thin <- W2s_block[seq(1, n_gibbs_samples, by=thining_gap+1)]
    
    # Note: Ws_block is a list of N realizations of W of current replicate
    # Note: AW_N_1 is a matrix of n_test * N
    AW_N_1 <- if (length(ngme_1rep$models) == 0) 0 else
      Reduce(cbind, sapply(Ws_block, function(W) A_pred_block %*% W))
    AW_N_2 <- if (length(ngme_1rep$models) == 0) 0 else
      Reduce(cbind, sapply(W2s_block, function(W) A_pred_block %*% W))
    
    AW_N_1_thin <- if (length(ngme_1rep$models) == 0) 0 else
      Reduce(cbind, sapply(Ws_block_thin, function(W) A_pred_block %*% W))
    AW_N_2_thin <- if (length(ngme_1rep$models) == 0) 0 else
      Reduce(cbind, sapply(W2s_block_thin, function(W) A_pred_block %*% W))
    
    # sampling Y by, Y = X feff + (block_A %*% block_W) + eps
    # AW_N_1[[1]] is concat(A1 W1, A2 W2, ..)
    
    # generate fixed effect
    fe <- with(ngme_1rep, as.numeric(X_pred %*% feff))
    fe_N <- matrix(
      rep(fe, n_gibbs_samples), 
      ncol = n_gibbs_samples, 
      byrow=F
    )
    
    n_thin <- length(Ws_block_thin)
    fe_N_thin <- matrix(
      rep(fe, n_thin), 
      ncol = n_thin, 
      byrow=F
    )
    
    # prediction using all samples
    pred_N_1 <- fe_N + AW_N_1
    pred_N_2 <- fe_N + AW_N_2
    
    # prediction using thinning samples
    pred_N_1_thin <- fe_N_thin + AW_N_1_thin
    pred_N_2_thin <- fe_N_thin + AW_N_2_thin
    
    # simulate measurement noise
    seed_group_1 <- seed_int + 1:n_thin
    seed_group_2 <- seed_int + 1:n_thin + n_thin
    mn_N_1 <- sapply(1:n_thin, function(x) simulate(test_noise, seed = seed_group_1[x])[[1]])
    mn_N_2 <- sapply(1:n_thin, function(x) simulate(test_noise, seed = seed_group_2[x])[[1]])
    
    # simulate y using thinning samples
    Y_N_1_thin <- pred_N_1_thin + mn_N_1
    Y_N_2_thin <- pred_N_2_thin + mn_N_2
    
    compute_score_given_pred(
      transform(pred_N_1), transform(pred_N_2),  # for MAE, MSE
      transform(Y_N_1_thin), transform(Y_N_2_thin), # for CRPS, sCRPS
      transform(y_data), 
      group_data
    )
  }, error = function(e) {
    warning(paste("Error in compute_scores:", e$message))
    # Return a default score structure with NAs
    n_group <- length(levels(group_data))
    if (n_group == 0) n_group <- 1  # Default to 1 if no groups
    
    scores <- matrix(NA, nrow=n_group, ncol=4)
    colnames(scores) <- c("MAE", "MSE", "neg.CRPS", "neg.sCRPS")
    if (n_group == 1) {
      rownames(scores) <- "all"
    } else {
      rownames(scores) <- levels(group_data)
    }
    scores
  })
}



# helper function to generate predictions and simulate Y
compute_pred_N <- function(
  ngme_1rep, n_gibbs_samples, n_burnin, seed, A_pred_block, noise_test_idx, y_data, group_data, X_pred
) {
  seed_int <- as.integer(abs(seed) %% 2147483647)

  Ws <- sampling_cpp(
    ngme_1rep, 
    n = n_gibbs_samples * 2, 
    n_burnin = n_burnin,
    posterior = TRUE, 
    seed=seed_int
  )[["W"]]
  # )[["cond_W"]]
  
  Ws_block  <- head(Ws, n_gibbs_samples)
  W2s_block <- tail(Ws, n_gibbs_samples)

# Note: Ws_block is a list of N realizations of W of current replicate
# Note: AW_N_1 is a matrix of n_test * N
  AW_N_1 <- if (length(ngme_1rep$models) == 0) 0 else
    Reduce(cbind, sapply(Ws_block, function(W) A_pred_block %*% W))
  AW_N_2 <- if (length(ngme_1rep$models) == 0) 0 else
    Reduce(cbind, sapply(W2s_block, function(W) A_pred_block %*% W))

  # sampling Y by, Y = X feff + (block_A %*% block_W) + eps
  # AW_N_1[[1]] is concat(A1 W1, A2 W2, ..)

  # generate fixed effect
  fe <- with(ngme_1rep, as.numeric(X_pred %*% feff))
  fe_N <- matrix(
    rep(fe, n_gibbs_samples), 
    ncol = n_gibbs_samples, 
    byrow=F
  )

  pred_N_1 <- fe_N + AW_N_1
  pred_N_2 <- fe_N + AW_N_2
  
  # simulate measurement noise
  seed_group_1 <- seed_int + 1:n_thin
  seed_group_2 <- seed_int + 1:n_thin + n_thin
  mn_N_1 <- sapply(1:n_thin, function(x) simulate(test_noise, seed = seed_group_1[x])[[1]])
  mn_N_2 <- sapply(1:n_thin, function(x) simulate(test_noise, seed = seed_group_2[x])[[1]])

  # simulate y
  Y_N_1 <- pred_N_1 + mn_N_1
  Y_N_2 <- pred_N_2 + mn_N_2

  list(
    pred_N_1 = rowMeans(as.matrix(pred_N_1)),
    pred_N_2 = rowMeans(as.matrix(pred_N_2)),
    Y_N_1 = rowMeans(as.matrix(Y_N_1)),
    Y_N_2 = rowMeans(as.matrix(Y_N_2))
  )
}



#' Compute the scores given the prediction
#'
#' @param pred_N_1 a matrix of n_obs * N
#' @param pred_N_2 a matrix of n_obs * N
#' @param Y_N_1 a matrix of n_obs * N
#' @param Y_N_2 a matrix of n_obs * N
#' @param y_data a vector of length n_obs
#' @param group_data a vector of length n_obs
compute_score_given_pred <- function(
  pred_N_1, pred_N_2, 
  Y_N_1_thin, Y_N_2_thin,
  y_data, group_data
) {
  # Use ensemble prediction (buggy)
  pred <- 0.5*(rowMeans(as.matrix(pred_N_1)) + rowMeans(as.matrix(pred_N_2)))

  # Use only 1 realization
  # pred <- pred_N_1[, 1] 

  # Now Y is of dim n_obs * N
  n_obs <- length(y_data)
  E_sim_data_thin <- E_sim_sim_thin <- double(length(y_data))
  
  for (i in 1:n_obs) {
    # turn row of df into numeric vector.
    y_sim1_thin <- as.numeric(Y_N_1_thin[i, ])
    y_sim2_thin <- as.numeric(Y_N_2_thin[i, ])

    # estimate E(| X_i - y_data |). y_data is observation, X_i ~ predictive distribution at i
    E_sim_data_thin[[i]] <- mean(abs(y_sim1_thin - y_data[i]))

    # estimate E(| X_i - Y_i |) , X_i, Y_i ~ predictive distribution at i
    E_sim_sim_thin[[i]] <- mean(abs(y_sim1_thin - y_sim2_thin))
  }

  # compute MSE, MAE, CRPS, sCRPS within each group
  pred_each_group   <- split(pred, group_data)
  y_data_each_group <- split(y_data, group_data)

  CRPS_thin <- split(0.5 * E_sim_sim_thin - E_sim_data_thin, group_data)
  sCRPS_thin <- split(
    -E_sim_data_thin / E_sim_sim_thin - 0.5 * log(E_sim_sim_thin),
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
    neg.CRPS  = -sapply(CRPS_thin, mean), # mean over 1:n_obs_test within each group
    neg.sCRPS = -sapply(sCRPS_thin, mean) # same
  )

  scores
}


#' Merge model of replicates into model of 1 replicate given train_idx and 
#' test_idx, the merged model contains all the information of train_idx from 
#' different replicates.
#'
#' @param ngme a ngme object
#' @param train_idx a vector of indices of train data
#' @param test_idx a vector of indices of test data
merge_replicates <- function(
  ngme, 
  train_idx, 
  test_idx
) {
  repls <- attr(ngme, "fit")$replicate
  uni_repl <- unique(repls)
  merged_rep <- ngme$replicates[[1]]
  n_latent <- length(merged_rep$models)
  

  train_X <- train_Y <- train_noise <- list()
  test_X <- test_Y <- test_noise <- test_group <- list()
  A_preds_block <- list()
  A_train <- vector("list", n_latent)
  # Loop over each replicate
  for (i in seq_along(uni_repl)) {
    data_idx_rep <- ngme$replicates[[i]]$data_idx
    bool_train_idx <- data_idx_rep %in% train_idx # current rep has train
    bool_test_idx  <- data_idx_rep %in% test_idx  # current rep has test

    ngme_1rep <- ngme$replicates[[i]]
    train_X[[i]]     <- ngme_1rep$X[bool_train_idx,, drop=FALSE]
    train_Y[[i]]     <- ngme_1rep$Y[bool_train_idx]
    train_noise[[i]] <- subset_noise(
      ngme_1rep$noise, sub_idx = bool_train_idx, compute_corr = TRUE
    )

    test_X[[i]]     <- ngme_1rep$X[bool_test_idx,, drop=FALSE]
    test_Y[[i]]     <- ngme_1rep$Y[bool_test_idx]
    test_group[[i]] <- ngme_1rep$group[bool_test_idx]
    test_noise[[i]] <- subset_noise(
      ngme_1rep$noise, sub_idx = bool_test_idx, compute_corr = TRUE
    )

    # Subset A for test and train location
    A_preds <- list();
    for (j in seq_len(n_latent)) {
      A_train[[j]][[i]] <- ngme_1rep$models[[j]]$A[bool_train_idx,,drop=FALSE]
      A_preds[[j]] <- ngme_1rep$models[[j]]$A[bool_test_idx, ,drop=FALSE]
    }
    A_preds_block[[i]] <- Reduce(cbind, x = A_preds)
  }

  # Merge train data
  merged_rep$X     <- Reduce(rbind, train_X)
  merged_rep$Y     <- Reduce(c, train_Y)
  merged_rep$noise <- Reduce(merge_noise, train_noise)
  for (j in seq_len(n_latent)) {
    merged_rep$models[[j]]$A <- Reduce(rbind, A_train[[j]])
  }

  # organize test data
  test_Y       <- Reduce(c, test_Y)
  test_group   <- Reduce(c, test_group)
  test_X       <- Reduce(rbind, test_X)
  test_noise   <- Reduce(merge_noise, test_noise)
  test_A_block <- Reduce(rbind, x = A_preds_block)

  return(list(
    merged_rep = merged_rep,
    test_Y = test_Y,
    test_X = test_X,
    test_noise = test_noise,
    test_A_block = test_A_block,
    test_group = test_group
  ))
}
