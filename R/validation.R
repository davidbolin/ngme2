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
#' @param metric Optional function or list of functions (one per model) that maps the group-wise observations/predictions for a single location to the quantity that should be scored.
#'   The function receives a list containing at least `y` (named numeric vector of group values) and may optionally use `samples1`/`samples2` (matrices with rows named by group and columns indexing posterior draws).
#'   The function must return either a numeric scalar (optionally named) or a list with components `y` (scalar), and optionally `samples1`/`samples2` (numeric vectors matching the posterior draw count) and `label` (character).
#'   When `NULL`, the original per-group scores are computed.
#'   For example, to compare a linear combination of two fields you can use
#'   `metric = function(data) { res <- 2 * data$y["A"] + data$y["B"]; names(res) <- "combo"; res }`. To simply sum all group values, return `sum(data$y)`.
#' @param n_gibbs_samples number of gibbs samples of latent process, used for computing CRPS, sCRPS
#' @param n_burnin number of burnin
#' @param test_idx a list of indices of the data (which data points to be predicted) (only for custom type)
#' @param train_idx  a list of indices of the data (which data points to be used for re-sampling (not re-estimation)) (only for custom type)
#' @param keep_pred logical, keep test information (pred_1, pred_2) in the return (as attributes), pred_1 and pred_2 are the prediction of the two chains
#' @param thining_gap integer, the gap between samples for thinning, if 0, then no thinning, if 1, then keep 50% of the samples for CRPS, sCRPS, etc.
#' @param parallel logical, run in parallel mode
#' @param cores_layer1 integer, number of cores for the first layer (over testing samples)
#' @param cores_layer2 integer, number of cores for the second layer (over computing scores for N_sim simulations)
#' @param merge_groups logical, if TRUE, merge groups as vector components (e.g., for vector-valued wind data with north_wind, east_wind).
#' MAE becomes Euclidean distance, MSE becomes squared Euclidean distance, etc.
#' @param merged_group_name character, name for the merged group when merge_groups=TRUE.
#' If NULL, uses "group1_group2" format (default: NULL)
#' @param data optional data.frame used to replace the original fitting data before running CV.
#' If `NULL`, the data stored in `ngme` is used. If provided, the model is rebuilt on `data`
#' while reusing fitted parameters from `ngme`.
#' @param chain_combine how to combine multiple optimization chains:
#'   `"param_mean"` uses the fitted object directly (default), while
#'   `"predictive_average"` computes predictions from each optimization chain
#'   and averages at the predictive level.
#' @return A list containing:
#'   \itemize{
#'     \item mean.scores - mean of N_sim estimations of 4 criterions: MSE, MAE, CRPS, sCRPS
#'     \item sd.scores - standard deviation of N_sim estimations of 4 criterions: MSE, MAE, CRPS, sCRPS
#'   }
#' @export
cross_validation <- function(
    ngme,
    type = "k-fold",
    seed = NULL,
    print = TRUE,
    N_sim = 5,
    n_gibbs_samples = 500,
    n_burnin = 100,
    k = 5,
    percent = 0.2,
    times = 10,
    metric = NULL,
    test_idx = NULL,
    train_idx = NULL,
    keep_pred = FALSE,
    parallel = FALSE,
    thining_gap = 1, # Used for computing CRPS, sCRPS, the gap between samples for thinning, if 0, then no thinning, if 1, then keep 50% of the samples for CRPS, sCRPS, etc.
    # merge_replicates = FALSE, # remove this option
    cores_layer1 = if (parallel) min(parallel::detectCores(), 2) else 1, # Limit to 2 cores for safety
    cores_layer2 = if (parallel) min(parallel::detectCores(), 2) else 1, # Limit to 2 cores for safety
    merge_groups = FALSE, # New parameter for merging bivariate groups
    merged_group_name = NULL, # New parameter for custom merged group name
    data = NULL,
    chain_combine = c("param_mean", "predictive_average")
    ) {
  chain_combine <- match.arg(chain_combine)
  merge_replicates <- FALSE

  if (!requireNamespace("parallel", quietly = TRUE)) {
    message("Parallel package not available. Running in serial mode. You can install `parallel` package to speed up the computation.")
    parallel <- FALSE
  }

  if (inherits(ngme, "ngme")) ngme <- list(ngme)
  if (is.null(names(ngme))) names(ngme) <- paste("model", seq_along(ngme), sep = "_")

  if (!is.null(data)) {
    stopifnot("`data` must be a data.frame." = is.data.frame(data))
  }

  ngme_chain_sets <- NULL
  if (chain_combine == "predictive_average") {
    ngme_chain_sets <- lapply(ngme, function(model_fit) {
      chain_models <- resolve_ngme_chain_fits(model_fit)
      if (!is.null(data)) {
        chain_models <- lapply(chain_models, rebuild_cv_model_with_data, data = data)
      }
      chain_models
    })
    ngme <- lapply(ngme_chain_sets, function(chain_models) chain_models[[1]])
  } else if (!is.null(data)) {
    ngme <- lapply(ngme, rebuild_cv_model_with_data, data = data)
  }

  # Handle metric argument: allow NULL, function, or list of functions (one per model)
  if (is.null(metric)) {
    metric <- vector("list", length(ngme))
  } else if (is.list(metric)) {
    if (length(metric) != length(ngme)) {
      stop("If metric is a list, its length must equal number of models (length(ngme)).")
    }
  } else if (is.function(metric)) {
    metric <- rep(list(metric), length(ngme))
  } else {
    stop("metric must be NULL, a function, or a list of functions")
  }

  metric_supplied <- any(vapply(metric, function(fun) !is.null(fun), logical(1)))
  if (metric_supplied && merge_groups) {
    stop("merge_groups cannot be used together with a custom metric; please encode the aggregation inside the metric function.")
  }

  n_data <- attr(ngme[[1]], "fit")$n_data
  if (is.null(n_data)) stop("Please provide ngme object or a list of ngme object")

  if (!is.null(seed)) set.seed(seed)
  stopifnot(
    "type should be in c('k-fold', 'loo', 'lpo', 'custom')" = type %in% c("k-fold", "loo", "lpo", "custom")
  )

  # 1. compute indices of tartget and train if not custom type
  if (type == "k-fold") {
    # Check if we should use paired CV splits for grouped data
    if (merge_groups && !is.null(ngme[[1]]$group)) {
      stop("Paired CV splits are not supported for grouped data. Please use custom type. Check ?create_paired_cv_splits for more details.")
    } else {
      # Regular k-fold split
      idx <- seq_len(n_data)
      folds <- cut(sample(idx), breaks = k, label = FALSE)
      test_idx <- lapply(1:k, function(x) {
        which(folds == x, arr.ind = TRUE)
      })
      train_idx <- lapply(1:k, function(x) {
        which(folds != x, arr.ind = TRUE)
      })
    }
  } else if (type == "loo") {
    idx <- seq_len(n_data)
    test_idx <- lapply(idx, function(i) i)
    train_idx <- lapply(idx, function(i) idx[idx != i])
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
      "test_idx and train_idx should be provided" = !is.null(test_idx) && !is.null(train_idx),
      "test_idx and train_idx should be a list" = is.list(test_idx) && is.list(train_idx),
      "test_idx and train_idx should be of same length" = length(test_idx) == length(train_idx)
    )
  }
  # Alternative. do not distinguish between replicates?
  # But the internal mesh may not be the same for each replicate....

  # 2. loop over each test_idx and train_idx, and compute the criterion
  final_mean <- final_sd <- list()
  pred_1 <- list()
  pred_2 <- list()
  Y_1 <- list()
  Y_2 <- list()

  compute_err <- if (merge_replicates) compute_err_merged_reps else compute_err_reps

  for (idx in seq_along(ngme)) {
    if (print) cat(paste0("Model ", names(ngme)[[idx]], ": \n\n"))
    scores <- sd_scores <- NULL
    chain_models <- if (is.null(ngme_chain_sets)) NULL else ngme_chain_sets[[idx]]

    # loop over each test_idx and train_idx
    if (parallel && requireNamespace("parallel", quietly = TRUE)) {
      if (print) cat("Running in parallel mode. \n")

      ngme_list <- list()
      for (i in 1:length(test_idx)) {
        ngme_list[[i]] <- ngme[[idx]]
      }

      # Add error handling for parallel execution
      tryCatch(
        {
          ret <- parallel::mclapply(
            seq_along(test_idx),
            function(i) {
              tryCatch(
                {
                  result <- compute_err(
                    ngme_list[[i]],
                    test_idx[[i]],
                    train_idx[[i]],
                    N_sim = N_sim,
                    n_gibbs_samples = n_gibbs_samples,
                    seed = seed,
                    keep_pred = keep_pred,
                    parallel = TRUE,
                    metric = metric[[idx]],
                    num_cores = cores_layer2,
                    thining_gap = thining_gap,
                    merge_groups = merge_groups,
                    merged_group_name = merged_group_name,
                    chain_models = chain_models,
                    chain_combine = chain_combine
                  )
                  if (print) {
                    cat(paste("In test batch", i, ": \n"))
                    print(as.data.frame(result$scores))
                    cat("\n")
                  }
                  result
                },
                error = function(e) {
                  warning(paste("Error in test batch", i, ":", e$message))
                  NULL
                }
              )
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
        },
        error = function(e) {
          warning(paste("Error in parallel execution:", e$message))
          parallel <- FALSE
        }
      )
    }

    # Fall back to sequential if parallel failed or was not requested
    if (!parallel || length(scores) == 0) {
      for (i in seq_along(test_idx)) {
        tryCatch(
          {
            result <- compute_err(
              ngme[[idx]],
              test_idx[[i]],
              train_idx[[i]],
              N_sim = N_sim,
              n_gibbs_samples = n_gibbs_samples,
              seed = seed,
              keep_pred = keep_pred,
              parallel = FALSE,
              metric = metric[[idx]],
              thining_gap = thining_gap,
              merge_groups = merge_groups,
              merged_group_name = merged_group_name,
              chain_models = chain_models,
              chain_combine = chain_combine
            )
            scores[[i]] <- result$scores
            sd_scores[[i]] <- result$sd_scores

            if (print) {
              cat(paste("In test batch", i, ": \n"))
              print(as.data.frame(scores[[i]]))
              cat("\n")
            }
          },
          error = function(e) {
            warning(paste("Error in sequential computation for test batch", i, ":", e$message))
          }
        )
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

  # Check if merging actually succeeded by looking at the number of rows
  # If merge_groups=TRUE but merging failed, we'll have 2 rows per model instead of 1
  n_models <- length(ngme)
  merging_succeeded <- merge_groups && (nrow(final_mean) == n_models)

  if (merging_succeeded) {
    # For successfully merged groups, use simplified naming
    rownames(final_mean) <- names(ngme)
    rownames(final_sd) <- names(ngme)
  } else if (nrow(final_mean) == n_models) {
    # univariate model (1 row per model)
    rownames(final_mean) <- names(ngme)
    rownames(final_sd) <- names(ngme)
  } else {
    # bivariate model (2 rows per model)
    if (nrow(final_mean) > 0) {
      n_groups_per_model <- nrow(final_mean) / n_models
      if (n_groups_per_model == 2) {
        # Extract group names from the first set of results
        group_names <- rownames(final_mean)[1:2]
        names_list <- lapply(names(ngme), function(x) paste(x, group_names, sep = "_"))
        rownames(final_mean) <- do.call(c, names_list)
        rownames(final_sd) <- rownames(final_mean)
      }
    }
  }

  ret <- list(
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


resolve_ngme_chain_fits <- function(object) {
  stopifnot(inherits(object, "ngme"))

  if (exists("get_ngme_chain_fits", mode = "function", inherits = TRUE)) {
    return(get_ngme_chain_fits(object))
  }

  chain_outputs <- attr(object, "chain_outputs")
  if (!is.list(chain_outputs) || length(chain_outputs) == 0) {
    return(list(object))
  }

  lapply(chain_outputs, function(chain_output) {
    chain_fit <- object
    if (!is.list(chain_output) || length(chain_fit$replicates) != length(chain_output)) {
      return(object)
    }

    for (i in seq_along(chain_fit$replicates)) {
      chain_fit$replicates[[i]] <- update_ngme_est(
        chain_fit$replicates[[i]],
        chain_output[[i]]
      )
    }

    attr(chain_fit, "chain_outputs") <- NULL
    attr(chain_fit, "chain_fits") <- NULL
    chain_fit
  })
}


rebuild_cv_model_with_data <- function(ngme_obj, data) {
  stopifnot("`ngme_obj` must be of class 'ngme'." = inherits(ngme_obj, "ngme"))
  fit <- attr(ngme_obj, "fit")
  stopifnot(
    "Missing fit metadata in `ngme_obj`." = !is.null(fit),
    "Missing formula in fit metadata." = !is.null(fit$formula)
  )

  replicate_new <- resolve_cv_partition_for_new_data(
    old_data = fit$data,
    old_values = fit$replicate,
    new_data = data,
    arg_name = "replicate"
  )

  group_new <- resolve_cv_partition_for_new_data(
    old_data = fit$data,
    old_values = fit$group,
    new_data = data,
    arg_name = "group"
  )

  standardize_fixed <- if (is.null(ngme_obj$replicates[[1]]$standardize)) {
    TRUE
  } else {
    isTRUE(ngme_obj$replicates[[1]]$standardize)
  }

  build_args <- list(
    formula = fit$formula,
    data = data,
    family = fit$family,
    control_opt = control_opt(
      estimation = FALSE,
      standardize_fixed = standardize_fixed
    ),
    control_ngme = ngme_obj$replicates[[1]]$control_ngme,
    group = group_new,
    replicate = replicate_new
  )
  rebuild_env <- build_cv_rebuild_env(
    formula = fit$formula,
    ngme_obj = ngme_obj,
    data = data
  )

  tryCatch(
    {
      do.call("ngme", c(build_args, list(start = ngme_obj)), envir = rebuild_env)
    },
    error = function(e) {
      if (grepl("length of [WV] should be the same", e$message)) {
        # Fallback for new data with different latent/noise state dimension:
        # rebuild without start and transplant only compatible hyperparameters.
        rebuilt <- do.call("ngme", build_args, envir = rebuild_env)
        rebuilt <- transplant_cv_hyperparameters(rebuilt, ngme_obj)
        return(rebuilt)
      }
      stop(
        "Failed to rebuild ngme object with supplied `data`: ",
        e$message,
        call. = FALSE
      )
    }
  )
}


build_cv_rebuild_env <- function(formula, ngme_obj, data) {
  env <- new.env(parent = parent.frame())
  env$ngme <- ngme

  formula_env <- environment(formula)
  symbols <- setdiff(unique(all.vars(formula)), names(data))
  if (length(symbols) == 0) return(env)
  if (length(ngme_obj$replicates) == 0) return(env)
  if (length(ngme_obj$replicates[[1]]$models) == 0) return(env)

  rep1 <- ngme_obj$replicates[[1]]
  model1 <- rep1$models[[1]]
  operator <- model1$operator
  model_noise <- model1$noise

  for (sym in symbols) {
    if (exists(sym, envir = env, inherits = TRUE)) next
    val <- resolve_cv_symbol_binding(
      sym = sym,
      formula_env = formula_env,
      operator = operator,
      model_noise = model_noise
    )
    if (!is.null(val)) assign(sym, val, envir = env)
  }

  env
}


resolve_cv_symbol_binding <- function(sym, formula_env, operator, model_noise) {
  if (is.environment(formula_env) && exists(sym, envir = formula_env, inherits = TRUE)) {
    return(get(sym, envir = formula_env, inherits = TRUE))
  }

  if (sym == "mesh" && !is.null(operator$mesh)) return(operator$mesh)
  if (sym %in% c("B", "B_kappa", "B_K") && !is.null(operator$B_K)) return(operator$B_K)
  if (sym %in% c("B_sigma", "B_mu", "B_nu") && !is.null(model_noise[[sym]])) return(model_noise[[sym]])
  if (grepl("^n_?basis$", sym) && !is.null(operator$B_K)) return(ncol(operator$B_K))

  NULL
}


transplant_cv_hyperparameters <- function(new_obj, old_obj) {
  stopifnot(
    inherits(new_obj, "ngme"),
    inherits(old_obj, "ngme")
  )

  n_repl <- min(length(new_obj$replicates), length(old_obj$replicates))
  if (n_repl == 0) return(new_obj)

  for (i in seq_len(n_repl)) {
    rep_new <- new_obj$replicates[[i]]
    rep_old <- old_obj$replicates[[i]]

    if (length(rep_new$feff) == length(rep_old$feff)) {
      rep_new$feff <- as.numeric(rep_old$feff)
      names(rep_new$feff) <- names(new_obj$replicates[[i]]$feff)
    }

    rep_new$noise <- transplant_noise_hyperparameters(rep_new$noise, rep_old$noise)

    n_models <- min(length(rep_new$models), length(rep_old$models))
    if (n_models > 0) {
      for (j in seq_len(n_models)) {
        model_new <- rep_new$models[[j]]
        model_old <- rep_old$models[[j]]

        if (!identical(model_new$model, model_old$model)) next

        theta_new <- model_new$operator$theta_K
        theta_old <- model_old$operator$theta_K

        if (!is.null(names(theta_new)) && !is.null(names(theta_old))) {
          common <- intersect(names(theta_new), names(theta_old))
          if (length(common) > 0) theta_new[common] <- theta_old[common]
        } else if (length(theta_new) == length(theta_old)) {
          theta_new <- theta_old
        }

        model_new$operator$theta_K <- theta_new
        model_new$theta_K <- theta_new
        model_new$operator$K <- model_new$operator$update_K(theta_new)

        model_new$noise <- transplant_noise_hyperparameters(model_new$noise, model_old$noise)
        rep_new$models[[j]] <- model_new
      }
    }

    new_obj$replicates[[i]] <- rep_new
  }

  new_obj
}


transplant_noise_hyperparameters <- function(noise_new, noise_old) {
  old_no_state <- noise_old
  old_no_state$V <- NULL
  update_noise(noise_new, new_noise = old_no_state)
}


resolve_cv_partition_for_new_data <- function(old_data, old_values, new_data, arg_name) {
  if (is.null(old_values)) return(NULL)

  if (length(old_values) == nrow(new_data)) {
    return(old_values)
  }

  source_col <- infer_cv_partition_column(old_data, old_values, arg_name)
  if (!is.null(source_col) && source_col %in% names(new_data)) {
    values <- new_data[[source_col]]
    if (is.factor(old_values)) {
      old_levels <- levels(old_values)
      new_levels <- unique(as.character(values[!is.na(values)]))
      ordered_levels <- c(intersect(old_levels, new_levels), setdiff(new_levels, old_levels))
      values <- factor(values, levels = ordered_levels)
    }
    if (anyNA(values)) {
      stop(
        "Detected NA while preparing `", arg_name, "` from `data$",
        source_col, "`. Please ensure there are no missing values."
      )
    }
    return(values)
  }

  if (length(unique(as.character(old_values))) == 1) {
    return(rep(old_values[1], nrow(new_data)))
  }

  stop(
    "Unable to infer `", arg_name, "` for supplied `data`. ",
    "Please include the same `", arg_name, "` column used during fitting."
  )
}


infer_cv_partition_column <- function(old_data, old_values, arg_name) {
  if (!is.data.frame(old_data) || nrow(old_data) != length(old_values)) {
    return(NULL)
  }

  target <- as.character(old_values)

  if (arg_name %in% names(old_data) &&
      identical(as.character(old_data[[arg_name]]), target)) {
    return(arg_name)
  }

  matched_cols <- names(old_data)[vapply(old_data, function(col) {
    length(col) == length(target) && identical(as.character(col), target)
  }, logical(1))]

  if (length(matched_cols) == 0) return(NULL)
  matched_cols[[1]]
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
    metric = NULL,
    num_cores = 1,
    thining_gap = 0,
    merge_groups = FALSE,
    merged_group_name = NULL,
    chain_models = NULL,
    chain_combine = "param_mean") {
  if (is.null(seed)) seed <- Sys.time()
  stopifnot("Not a ngme object." = inherits(ngme, "ngme"))

  test_idx <- sort(test_idx)
  repls <- attr(ngme, "fit")$replicate
  uni_repl <- unique(repls)

  scores <- sd_scores <- weight <- NULL
  n_scores <- 0
  pred_1 <- double(length = length(test_idx))
  pred_2 <- double(length = length(test_idx))
  Y_1 <- double(length = length(test_idx))
  Y_2 <- double(length = length(test_idx))

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
      metric,
      thining_gap = thining_gap,
      merge_groups = merge_groups,
      merged_group_name = merged_group_name,
      chain_combine = chain_combine
    )
  }

  # compute mean and sd of scores
  n_result_rows <- nrow(scores[[1]])
  n_cols <- ncol(scores[[1]])
  array_3d <- array(unlist(scores),
    dim = c(n_result_rows, n_cols, length(scores))
  )

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
    metric = NULL,
    num_cores = 1,
    thining_gap = 1,
    merge_groups = FALSE,
    merged_group_name = NULL,
    chain_models = NULL,
    chain_combine = "param_mean") {
  test_idx <- sort(test_idx)
  stopifnot("Not a ngme object." = inherits(ngme, "ngme"))
  repls <- attr(ngme, "fit")$replicate
  uni_repl <- unique(repls)

  scores <- sd_scores <- weight <- NULL
  n_scores <- 0
  pred_1 <- double(length = length(test_idx))
  pred_2 <- double(length = length(test_idx))
  Y_1 <- double(length = length(test_idx))
  Y_2 <- double(length = length(test_idx))

  for (i in seq_along(uni_repl)) {
    data_idx_rep <- ngme$replicates[[i]]$data_idx
    bool_train_idx <- data_idx_rep %in% train_idx # current rep has train
    bool_test_idx <- data_idx_rep %in% test_idx # current rep has test

    # skip this replicate if no target or train data
    if (sum(bool_train_idx) == 0 || sum(bool_test_idx) == 0) next
    n_scores <- n_scores + 1

    ngme_1rep <- ngme$replicates[[i]]
    ngme_chain_reps <- NULL
    if (chain_combine == "predictive_average" &&
        !is.null(chain_models) &&
        length(chain_models) > 1) {
      ngme_chain_reps <- lapply(chain_models, function(chain_fit) {
        chain_fit$replicates[[i]]
      })
    }
    result_1rep <- compute_err_1rep(
      ngme_1rep,
      bool_train_idx = bool_train_idx,
      bool_test_idx = bool_test_idx,
      N_sim = N_sim,
      n_gibbs_samples = n_gibbs_samples,
      seed = seed,
      keep_pred = keep_pred,
      parallel = parallel,
      metric = metric,
      num_cores = num_cores,
      thining_gap = thining_gap,
      merge_groups = merge_groups,
      merged_group_name = merged_group_name,
      ngme_chain_reps = ngme_chain_reps,
      chain_combine = chain_combine
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
##
# A custom metric function can be supplied to combine group-wise observations and predictions into a single quantity before scoring.
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
    metric = NULL,
    num_cores = 1, # Default to 1 core to avoid potential issues
    thining_gap = 1,
    merge_groups = FALSE,
    merged_group_name = NULL,
    ngme_chain_reps = NULL,
    chain_combine = "param_mean") {
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
  X_pred <- ngme_1rep$X[bool_test_idx, , drop = FALSE]
  noise_test_idx <- subset_noise(
    ngme_1rep$noise,
    sub_idx = bool_test_idx, compute_corr = TRUE
  )

  # Subset noise, X, Y in train location
  ngme_1rep$X <- ngme_1rep$X[bool_train_idx, , drop = FALSE]
  ngme_1rep$Y <- ngme_1rep$Y[bool_train_idx]

  ngme_1rep$noise <- subset_noise(
    ngme_1rep$noise,
    sub_idx = bool_train_idx, compute_corr = TRUE
  )

  # Subset A for test and train location
  A_preds <- list()
  for (i in seq_along(ngme_1rep$models)) {
    A_preds[[i]] <- ngme_1rep$models[[i]]$A[bool_test_idx, , drop = FALSE]
    ngme_1rep$models[[i]]$A <- ngme_1rep$models[[i]]$A[bool_train_idx, , drop = FALSE]
  }

  if (!is.null(ngme_chain_reps) && length(ngme_chain_reps) > 0) {
    ngme_chain_reps <- lapply(ngme_chain_reps, function(chain_rep) {
      chain_rep$X <- chain_rep$X[bool_train_idx, , drop = FALSE]
      chain_rep$Y <- chain_rep$Y[bool_train_idx]
      chain_rep$noise <- subset_noise(
        chain_rep$noise,
        sub_idx = bool_train_idx, compute_corr = TRUE
      )
      for (ii in seq_along(chain_rep$models)) {
        chain_rep$models[[ii]]$A <- chain_rep$models[[ii]]$A[bool_train_idx, , drop = FALSE]
      }
      chain_rep
    })
  }

  # A_pred_blcok <- [A1_pred .. An_pred]
  # extract A and cbind!
  A_pred_block <- Reduce(cbind, x = A_preds)

  if (is.null(seed)) seed <- as.integer(Sys.time())

  scores <- pred_1 <- pred_2 <- Y_1 <- Y_2 <- list()

  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    scores <- parallel::mclapply(1:N_sim, function(nn) {
      s <- compute_scores(
        ngme_1rep, n_gibbs_samples, n_burnin, seed + nn, A_pred_block, noise_test_idx, y_data, group_data, X_pred, metric, thining_gap, merge_groups, merged_group_name, ngme_chain_reps, chain_combine
      )
      s
    }, mc.cores = num_cores)
  } else {
    for (nn in 1:N_sim) {
      tryCatch(
        {
          scores[[nn]] <- compute_scores(
            ngme_1rep, n_gibbs_samples, n_burnin, seed + nn, A_pred_block,
            noise_test_idx, y_data, group_data, X_pred, metric, thining_gap, merge_groups, merged_group_name, ngme_chain_reps, chain_combine
          )
        },
        error = function(e) {
          warning(paste("Error in sequential computation:", e$message))
        }
      )
    }
  }

  # Check if we have any valid scores
  if (length(scores) == 0) {
    stop("Failed to compute any valid scores")
  }

  # compute extra pred_N?
  if (keep_pred) {
    tmp <- compute_pred_N(
      ngme_1rep, n_gibbs_samples, n_burnin, seed, A_pred_block, noise_test_idx, y_data, group_data, X_pred
    )
    pred_1 <- tmp$pred_N_1
    pred_2 <- tmp$pred_N_2
    Y_1 <- tmp$Y_N_1
    Y_2 <- tmp$Y_N_2
  }

  # Use the actual number of rows returned by scores instead of original group count
  # This handles the case where merge_groups=TRUE reduces 2 groups to 1 result
  n_result_rows <- nrow(scores[[1]])
  n_cols <- ncol(scores[[1]])

  # compute mean and sd
  array_3d <- array(unlist(scores),
    dim = c(n_result_rows, n_cols, length(scores))
  )

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


build_cv_aw_draws_matrix <- function(A_pred_block, Ws_block, n_obs, has_models) {
  if (!has_models) {
    return(matrix(0, nrow = n_obs, ncol = length(Ws_block)))
  }
  if (length(Ws_block) == 0) {
    return(matrix(0, nrow = n_obs, ncol = 0))
  }
  do.call(cbind, lapply(Ws_block, function(W) as.numeric(A_pred_block %*% W)))
}


collect_cv_chain_samples <- function(
    ngme_chain_reps,
    n_gibbs_samples,
    n_burnin,
    seed_int,
    X_pred) {
  Ws_block <- list()
  W2s_block <- list()
  chain_idx_1 <- integer()
  chain_idx_2 <- integer()
  fe_by_chain <- vector("list", length(ngme_chain_reps))

  for (ci in seq_along(ngme_chain_reps)) {
    chain_rep <- ngme_chain_reps[[ci]]
    fe_by_chain[[ci]] <- as.numeric(X_pred %*% chain_rep$feff)

    Ws_result <- sampling_cpp(
      chain_rep,
      n = n_gibbs_samples * 2,
      n_burnin = n_burnin,
      posterior = TRUE,
      seed = seed_int + ci - 1L
    )
    Ws <- Ws_result[["W"]]

    if (is.null(Ws) || length(Ws) == 0) {
      stop("sampling_cpp returned NULL or empty result for chain ", ci)
    }

    Ws_1 <- head(Ws, n_gibbs_samples)
    Ws_2 <- tail(Ws, n_gibbs_samples)
    Ws_block <- c(Ws_block, Ws_1)
    W2s_block <- c(W2s_block, Ws_2)
    chain_idx_1 <- c(chain_idx_1, rep(ci, length(Ws_1)))
    chain_idx_2 <- c(chain_idx_2, rep(ci, length(Ws_2)))
  }

  list(
    Ws_block = Ws_block,
    W2s_block = W2s_block,
    chain_idx_1 = chain_idx_1,
    chain_idx_2 = chain_idx_2,
    fe_by_chain = fe_by_chain
  )
}


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
    metric,
    thining_gap,
    merge_groups = FALSE,
    merged_group_name = NULL,
    ngme_chain_reps = NULL,
    chain_combine = "param_mean") {
  tryCatch(
    {
      seed_int <- as.integer(abs(seed) %% 2147483647)
      n_obs <- nrow(X_pred)
      has_models <- length(ngme_1rep$models) > 0
      use_chain_average <- chain_combine == "predictive_average" &&
        !is.null(ngme_chain_reps) &&
        length(ngme_chain_reps) > 1

      if (use_chain_average) {
        sampled <- collect_cv_chain_samples(
          ngme_chain_reps = ngme_chain_reps,
          n_gibbs_samples = n_gibbs_samples,
          n_burnin = n_burnin,
          seed_int = seed_int,
          X_pred = X_pred
        )

        AW_N_1 <- build_cv_aw_draws_matrix(
          A_pred_block,
          sampled$Ws_block,
          n_obs = n_obs,
          has_models = has_models
        )
        AW_N_2 <- build_cv_aw_draws_matrix(
          A_pred_block,
          sampled$W2s_block,
          n_obs = n_obs,
          has_models = has_models
        )

        fe_N_1 <- do.call(cbind, lapply(sampled$chain_idx_1, function(ci) sampled$fe_by_chain[[ci]]))
        fe_N_2 <- do.call(cbind, lapply(sampled$chain_idx_2, function(ci) sampled$fe_by_chain[[ci]]))

        n_draws <- ncol(fe_N_1)
        thin_idx <- seq(1, n_draws, by = thining_gap + 1)

        pred_N_1_thin <- fe_N_1[, thin_idx, drop = FALSE] + AW_N_1[, thin_idx, drop = FALSE]
        pred_N_2_thin <- fe_N_2[, thin_idx, drop = FALSE] + AW_N_2[, thin_idx, drop = FALSE]
      } else {
        Ws_result <- sampling_cpp(
          ngme_1rep,
          n = n_gibbs_samples * 2,
          n_burnin = n_burnin,
          posterior = TRUE,
          seed = seed_int
        )
        Ws <- Ws_result[["W"]]

        if (is.null(Ws) || length(Ws) == 0) {
          stop("sampling_cpp returned NULL or empty result")
        }

        Ws_block <- head(Ws, n_gibbs_samples)
        W2s_block <- tail(Ws, n_gibbs_samples)

        AW_N_1 <- build_cv_aw_draws_matrix(
          A_pred_block,
          Ws_block,
          n_obs = n_obs,
          has_models = has_models
        )
        AW_N_2 <- build_cv_aw_draws_matrix(
          A_pred_block,
          W2s_block,
          n_obs = n_obs,
          has_models = has_models
        )

        fe <- as.numeric(X_pred %*% ngme_1rep$feff)
        fe_N <- matrix(rep(fe, n_gibbs_samples), ncol = n_gibbs_samples, byrow = FALSE)
        thin_idx <- seq(1, n_gibbs_samples, by = thining_gap + 1)

        pred_N_1_thin <- fe_N[, thin_idx, drop = FALSE] + AW_N_1[, thin_idx, drop = FALSE]
        pred_N_2_thin <- fe_N[, thin_idx, drop = FALSE] + AW_N_2[, thin_idx, drop = FALSE]
      }

      n_thin <- ncol(pred_N_1_thin)
      seed_group_1 <- seed_int + seq_len(n_thin)
      seed_group_2 <- seed_int + seq_len(n_thin) + n_thin
      mn_N_1 <- sapply(seq_len(n_thin), function(x) simulate(test_noise, seed = seed_group_1[x])[[1]])
      mn_N_2 <- sapply(seq_len(n_thin), function(x) simulate(test_noise, seed = seed_group_2[x])[[1]])

      Y_N_1_thin <- pred_N_1_thin + mn_N_1
      Y_N_2_thin <- pred_N_2_thin + mn_N_2

      compute_score_given_pred(
        Y_N_1_thin, Y_N_2_thin,
        y_data,
        group_data,
        merge_groups = merge_groups,
        merged_group_name = merged_group_name,
        metric = metric
      )
    },
    error = function(e) {
      if (inherits(e, "invalid_metric_error")) {
        stop(e)
      }
      warning(paste("Error in compute_scores:", e$message))
      n_group <- length(levels(group_data))
      if (n_group == 0) n_group <- 1

      scores <- matrix(NA, nrow = n_group, ncol = 4)
      colnames(scores) <- c("MAE", "MSE", "neg.CRPS", "neg.sCRPS")
      if (n_group == 1) {
        rownames(scores) <- "all"
      } else {
        rownames(scores) <- levels(group_data)
      }
      scores
    }
  )
}



# helper function to generate predictions and simulate Y
compute_pred_N <- function(
    ngme_1rep, n_gibbs_samples, n_burnin, seed, A_pred_block, noise_test_idx, y_data, group_data, X_pred) {
  seed_int <- as.integer(abs(seed) %% 2147483647)

  Ws <- sampling_cpp(
    ngme_1rep,
    n = n_gibbs_samples * 2,
    n_burnin = n_burnin,
    posterior = TRUE,
    seed = seed_int
  )[["W"]]
  # )[["cond_W"]]

  Ws_block <- head(Ws, n_gibbs_samples)
  W2s_block <- tail(Ws, n_gibbs_samples)

  # Note: Ws_block is a list of N realizations of W of current replicate
  # Note: AW_N_1 is a matrix of n_test * N
  AW_N_1 <- if (length(ngme_1rep$models) == 0) {
    0
  } else {
    Reduce(cbind, sapply(Ws_block, function(W) A_pred_block %*% W))
  }
  AW_N_2 <- if (length(ngme_1rep$models) == 0) {
    0
  } else {
    Reduce(cbind, sapply(W2s_block, function(W) A_pred_block %*% W))
  }

  # sampling Y by, Y = X feff + (block_A %*% block_W) + eps
  # AW_N_1[[1]] is concat(A1 W1, A2 W2, ..)

  # generate fixed effect
  fe <- with(ngme_1rep, as.numeric(X_pred %*% feff))
  fe_N <- matrix(
    rep(fe, n_gibbs_samples),
    ncol = n_gibbs_samples,
    byrow = F
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
#' @param Y_N_1_thin posterior predictive draws (rows = observations, columns = samples)
#' @param Y_N_2_thin posterior predictive draws (rows = observations, columns = samples)
#' @param y_data a vector of length n_obs
#' @param group_data a vector of length n_obs
#' @param metric optional custom metric function used to combine group-wise values before scoring
compute_score_given_pred <- function(
    Y_N_1_thin, Y_N_2_thin,
    y_data, group_data,
    merge_groups = FALSE,
    merged_group_name = NULL,
    metric = NULL) {
  metric_data <- prepare_metric_data(metric, y_data, Y_N_1_thin, Y_N_2_thin, group_data)

  y_data <- metric_data$y
  group_data <- metric_data$group
  Y_N_1_thin <- metric_data$samples1
  Y_N_2_thin <- metric_data$samples2
  pred <- metric_data$pred

  # Now Y is of dim n_obs * N
  n_obs <- length(y_data)

  if (merge_groups && length(levels(group_data)) == 2) {
    # Handle merged groups for bivariate data (e.g., wind components)
    groups <- levels(group_data)
    group1_idx <- which(group_data == groups[1])
    group2_idx <- which(group_data == groups[2])

    # Check if we have matching observations for both groups
    if (length(group1_idx) != length(group2_idx)) {
      warning("Unequal number of observations for the two groups. Cannot merge.")
      # Fall through to original computation instead of setting merge_groups = FALSE
      # This ensures consistent return structure
    } else {
      # Organize data by groups
      pred_group1 <- pred[group1_idx]
      pred_group2 <- pred[group2_idx]
      y_data_group1 <- y_data[group1_idx]
      y_data_group2 <- y_data[group2_idx]

      # Compute vector-based metrics
      # MAE as Euclidean distance, MSE as squared Euclidean distance
      vec_mae <- sqrt((pred_group1 - y_data_group1)^2 + (pred_group2 - y_data_group2)^2)
      vec_mse <- (pred_group1 - y_data_group1)^2 + (pred_group2 - y_data_group2)^2

      # For CRPS and sCRPS, we need to work with the simulation arrays
      n_thin <- ncol(Y_N_1_thin)
      E_sim_data_thin <- E_sim_sim_thin <- double(length(group1_idx))

      for (i in 1:length(group1_idx)) {
        # Get simulation data for both components
        y_sim1_thin_comp1 <- as.numeric(Y_N_1_thin[group1_idx[i], ])
        y_sim2_thin_comp1 <- as.numeric(Y_N_2_thin[group1_idx[i], ])
        y_sim1_thin_comp2 <- as.numeric(Y_N_1_thin[group2_idx[i], ])
        y_sim2_thin_comp2 <- as.numeric(Y_N_2_thin[group2_idx[i], ])

        # Compute vector distances for CRPS
        obs_vec <- c(y_data_group1[i], y_data_group2[i])

        # E(||X_i - y_data||) where X_i ~ predictive distribution
        vec_dist_sim_data <- sqrt((y_sim1_thin_comp1 - obs_vec[1])^2 + (y_sim1_thin_comp2 - obs_vec[2])^2)
        E_sim_data_thin[i] <- mean(vec_dist_sim_data)

        # E(||X_i - Y_i||) where X_i, Y_i ~ predictive distribution
        vec_dist_sim_sim <- sqrt((y_sim1_thin_comp1 - y_sim2_thin_comp1)^2 + (y_sim1_thin_comp2 - y_sim2_thin_comp2)^2)
        E_sim_sim_thin[i] <- mean(vec_dist_sim_sim)
      }

      denom <- pmax(E_sim_sim_thin, .Machine$double.eps)
      CRPS_thin <- 0.5 * denom - E_sim_data_thin
      sCRPS_thin <- -E_sim_data_thin / denom - 0.5 * log(denom)

      # Create merged results
      scores <- data.frame(
        MAE = mean(vec_mae),
        MSE = mean(vec_mse),
        neg.CRPS = -mean(CRPS_thin),
        neg.sCRPS = -mean(sCRPS_thin)
      )

      # Use custom name if provided, otherwise use default "group1_group2" format
      if (!is.null(merged_group_name)) {
        rownames(scores) <- merged_group_name
      } else {
        rownames(scores) <- paste(groups, collapse = "_")
      }

      return(scores)
    }
  }

  # Original computation for non-merged groups
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
  pred_each_group <- split(pred, group_data)
  y_data_each_group <- split(y_data, group_data)

  denom <- pmax(E_sim_sim_thin, .Machine$double.eps)
  CRPS_thin <- split(0.5 * denom - E_sim_data_thin, group_data)
  sCRPS_thin <- split(
    -E_sim_data_thin / denom - 0.5 * log(denom),
    group_data
  )

  # Compute MAE and MSE within each group
  MAE <- MSE <- double(length(pred_each_group))
  for (j in seq_along(pred_each_group)) {
    MAE[[j]] <- mean(abs(pred_each_group[[j]] - y_data_each_group[[j]]))
    MSE[[j]] <- mean((pred_each_group[[j]] - y_data_each_group[[j]])^2)
  }

  scores <- data.frame(
    MAE = MAE,
    MSE = MSE,
    neg.CRPS = -sapply(CRPS_thin, mean), # mean over 1:n_obs_test within each group
    neg.sCRPS = -sapply(sCRPS_thin, mean) # same
  )

  scores
}


prepare_metric_data <- function(metric, y_data, Y_N_1_thin, Y_N_2_thin, group_data) {
  samples1 <- as.matrix(Y_N_1_thin)
  samples2 <- as.matrix(Y_N_2_thin)
  pred_default <- 0.5 * (rowMeans(samples1) + rowMeans(samples2))

  group_factor <- if (is.factor(group_data)) droplevels(group_data) else factor(group_data)

  if (is.null(metric)) {
    return(list(
      y = as.numeric(y_data),
      pred = pred_default,
      samples1 = samples1,
      samples2 = samples2,
      group = group_factor
    ))
  }

  if (!is.function(metric)) {
    stop("Custom metric must be provided as a function.")
  }

  groups <- levels(group_factor)
  idx_list <- split(seq_along(group_factor), group_factor)
  counts <- lengths(idx_list)
  if (length(unique(counts)) != 1) {
    abort_invalid_metric(
      paste0(
        "Custom metric requires the same number of observations for each group within a fold. ",
        "Consider supplying paired splits via `type = \"custom\"`. Available groups: ",
        paste(groups, collapse = ", "),
        "."
      )
    )
  }

  n_locations <- counts[[1]]
  n_thin <- ncol(samples1)

  metric_y <- numeric(n_locations)
  metric_samples1 <- matrix(NA_real_, n_locations, n_thin)
  metric_samples2 <- matrix(NA_real_, n_locations, n_thin)
  label <- NULL

  for (i in seq_len(n_locations)) {
    obs_vec <- setNames(vapply(groups, function(g) as.numeric(y_data[idx_list[[g]][i]]), numeric(1)), groups)
    sample1_mat <- do.call(rbind, lapply(groups, function(g) samples1[idx_list[[g]][i], , drop = FALSE]))
    rownames(sample1_mat) <- groups
    sample2_mat <- do.call(rbind, lapply(groups, function(g) samples2[idx_list[[g]][i], , drop = FALSE]))
    rownames(sample2_mat) <- groups

    metric_input <- list(
      y = obs_vec,
      samples1 = sample1_mat,
      samples2 = sample2_mat,
      group = groups
    )

    metric_output <- metric(metric_input)
    parsed <- parse_metric_output(
      metric_output,
      expected_length = n_thin,
      available_groups = groups
    )

    metric_y[i] <- parsed$value
    if (is.null(label) && !is.null(parsed$label) && parsed$label != "") {
      label <- parsed$label
    }

    samples1_values <- parsed$samples1
    if (is.null(samples1_values)) {
      samples1_values <- compute_metric_samples(metric, sample1_mat, groups)
    }
    if (length(samples1_values) != n_thin) {
      abort_invalid_metric(
        "`samples1` returned by the custom metric must match the number of posterior draws."
      )
    }

    samples2_values <- parsed$samples2
    if (is.null(samples2_values)) {
      samples2_values <- compute_metric_samples(metric, sample2_mat, groups)
    }
    if (length(samples2_values) != n_thin) {
      abort_invalid_metric(
        "`samples2` returned by the custom metric must match the number of posterior draws."
      )
    }

    metric_samples1[i, ] <- samples1_values
    metric_samples2[i, ] <- samples2_values
  }

  metric_pred <- 0.5 * (rowMeans(metric_samples1) + rowMeans(metric_samples2))
  if (is.null(label) || label == "") {
    label <- "metric"
  }

  list(
    y = metric_y,
    pred = metric_pred,
    samples1 = metric_samples1,
    samples2 = metric_samples2,
    group = factor(rep(label, n_locations), levels = label)
  )
}


parse_metric_output <- function(result, expected_length = NULL, available_groups = NULL) {
  label <- NULL
  samples1 <- NULL
  samples2 <- NULL
  value <- NULL

  if (is.list(result)) {
    if (!is.null(result$label)) {
      label <- as.character(result$label)[1]
    }

    if (!is.null(result$y)) {
      value <- result$y
    } else if (!is.null(result$value)) {
      value <- result$value
    } else if (!is.null(result$result)) {
      value <- result$result
    } else if (length(result) == 1 && is.numeric(result[[1]])) {
      value <- result[[1]]
    }

    if (is.null(value)) {
      abort_invalid_metric(
        build_metric_message(
          "Custom metric must return a numeric value via `y` (or `value`).",
          available_groups
        )
      )
    }

    if (!is.null(names(value)) && length(value) == 1 && is.null(label)) {
      label <- names(value)[1]
    }

    value <- as.numeric(value)
    if (length(value) != 1 || !is.finite(value)) {
      abort_invalid_metric(
        build_metric_message(
          "Custom metric must provide a finite numeric scalar for `y`.",
          available_groups
        )
      )
    }

    if (!is.null(result$samples1)) {
      samples1 <- as.numeric(result$samples1)
    }
    if (!is.null(result$samples2)) {
      samples2 <- as.numeric(result$samples2)
    }
  } else if (is.numeric(result)) {
    value <- as.numeric(result)
    if (length(value) != 1 || !is.finite(value)) {
      abort_invalid_metric(
        build_metric_message(
          "Custom metric must return a scalar numeric value.",
          available_groups
        )
      )
    }
    if (!is.null(names(result))) {
      label <- names(result)[1]
    }
  } else {
    abort_invalid_metric(
      build_metric_message(
        "Custom metric must return either a numeric scalar or a list.",
        available_groups
      )
    )
  }

  if (!is.null(expected_length)) {
    if (!is.null(samples1) && length(samples1) > 0 && length(samples1) != expected_length) {
      abort_invalid_metric(
        "`samples1` returned by the custom metric must match the number of posterior draws."
      )
    }
    if (!is.null(samples2) && length(samples2) > 0 && length(samples2) != expected_length) {
      abort_invalid_metric(
        "`samples2` returned by the custom metric must match the number of posterior draws."
      )
    }
  }

  if (!is.null(samples1) && any(!is.finite(samples1))) {
    abort_invalid_metric("`samples1` returned by the custom metric must be finite.")
  }
  if (!is.null(samples2) && any(!is.finite(samples2))) {
    abort_invalid_metric("`samples2` returned by the custom metric must be finite.")
  }

  list(value = value, label = label, samples1 = samples1, samples2 = samples2)
}


compute_metric_samples <- function(metric, sample_mat, groups) {
  if (ncol(sample_mat) == 0) {
    return(numeric(0))
  }

  apply(sample_mat, 2, function(col) {
    sample_input <- list(y = setNames(as.numeric(col), groups), group = groups)
    parsed <- parse_metric_output(
      metric(sample_input),
      available_groups = groups
    )
    parsed$value
  })
}


abort_invalid_metric <- function(message) {
  stop(
    structure(
      list(message = message, call = sys.call(-1)),
      class = c("invalid_metric_error", "error", "condition")
    )
  )
}


build_metric_message <- function(message, available_groups = NULL) {
  if (is.null(available_groups) || length(available_groups) == 0) {
    return(message)
  }
  paste0(
    message,
    " Available groups: ",
    paste(available_groups, collapse = ", "),
    "."
  )
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
    test_idx) {
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
    bool_test_idx <- data_idx_rep %in% test_idx # current rep has test

    ngme_1rep <- ngme$replicates[[i]]
    train_X[[i]] <- ngme_1rep$X[bool_train_idx, , drop = FALSE]
    train_Y[[i]] <- ngme_1rep$Y[bool_train_idx]
    train_noise[[i]] <- subset_noise(
      ngme_1rep$noise,
      sub_idx = bool_train_idx, compute_corr = TRUE
    )

    test_X[[i]] <- ngme_1rep$X[bool_test_idx, , drop = FALSE]
    test_Y[[i]] <- ngme_1rep$Y[bool_test_idx]
    test_group[[i]] <- ngme_1rep$group[bool_test_idx]
    test_noise[[i]] <- subset_noise(
      ngme_1rep$noise,
      sub_idx = bool_test_idx, compute_corr = TRUE
    )

    # Subset A for test and train location
    A_preds <- list()
    for (j in seq_len(n_latent)) {
      A_train[[j]][[i]] <- ngme_1rep$models[[j]]$A[bool_train_idx, , drop = FALSE]
      A_preds[[j]] <- ngme_1rep$models[[j]]$A[bool_test_idx, , drop = FALSE]
    }
    A_preds_block[[i]] <- Reduce(cbind, x = A_preds)
  }

  # Merge train data
  merged_rep$X <- Reduce(rbind, train_X)
  merged_rep$Y <- Reduce(c, train_Y)
  merged_rep$noise <- Reduce(merge_noise, train_noise)
  for (j in seq_len(n_latent)) {
    merged_rep$models[[j]]$A <- Reduce(rbind, A_train[[j]])
  }

  # organize test data
  test_Y <- Reduce(c, test_Y)
  test_group <- Reduce(c, test_group)
  test_X <- Reduce(rbind, test_X)
  test_noise <- Reduce(merge_noise, test_noise)
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

#' Create paired indices for bivariate cross-validation
#' Ensures that paired observations (e.g., u_wind and v_wind at same location)
#' are kept together in the same fold
#'
#' @param data data frame containing the observations
#' @param loc_col character vector of column names defining locations (e.g., c("lon", "lat") or c("x", "y"))
#' @param group character, name of the group column (e.g., "direction" for wind data)
#' @param k number of folds
#' @param seed random seed
#' @return list with test_idx and train_idx, each containing k lists of indices
#' @export
#' @examples
#' # Example with wind data
#' data_long <- data.frame(
#'   lon = rep(c(1, 2, 3, 4), 2),
#'   lat = rep(c(1, 1, 2, 2), 2),
#'   direction = rep(c("u_wind", "v_wind"), each = 4),
#'   wind = rnorm(8)
#' )
#' splits <- create_paired_cv_splits(data_long, c("lon", "lat"), "direction", k = 2)
create_paired_cv_splits <- function(data, loc_col, group, k, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Validate inputs
  stopifnot(
    "data must be a data frame" = is.data.frame(data),
    "loc_col must be character vector" = is.character(loc_col),
    "group must be character" = is.character(group) && length(group) == 1,
    "k must be positive integer" = k > 0,
    "location columns must exist in data" = all(loc_col %in% names(data)),
    "group column must exist in data" = group %in% names(data)
  )

  # Create location identifier by pasting location columns together
  data$location_id <- do.call(paste, c(data[loc_col], sep = "_"))

  # Get unique locations
  unique_locations <- unique(data$location_id)
  n_locations <- length(unique_locations)

  # Check that each location has observations for all groups
  groups_per_location <- aggregate(data[[group]],
    by = list(location = data$location_id),
    FUN = function(x) length(unique(x))
  )
  unique_groups <- unique(data[[group]])
  n_groups <- length(unique_groups)

  # Verify that all locations have the same number of groups
  if (!all(groups_per_location$x == n_groups)) {
    incomplete_locs <- groups_per_location$location[groups_per_location$x != n_groups]
    warning(paste(
      "Some locations don't have observations for all groups:",
      paste(head(incomplete_locs), collapse = ", ")
    ))
  }

  # Create folds based on unique locations
  if (n_locations < k) {
    warning(paste(
      "Number of unique locations (", n_locations,
      ") is less than k (", k, "). Some folds will be empty."
    ))
    k <- min(k, n_locations)
  }

  # Assign locations to folds
  if (k == 1) {
    # Special case: all locations go to fold 1
    location_folds <- rep(1, n_locations)
    names(location_folds) <- unique_locations
  } else {
    location_folds <- cut(sample(1:n_locations), breaks = k, label = FALSE)
    names(location_folds) <- unique_locations
  }

  # Create test and train indices
  test_idx <- list()
  train_idx <- list()

  for (i in 1:k) {
    # Get locations assigned to this fold
    test_locations <- unique_locations[location_folds == i]
    train_locations <- unique_locations[location_folds != i]

    # Get all observation indices for these locations (across all groups)
    test_idx[[i]] <- which(data$location_id %in% test_locations)
    train_idx[[i]] <- which(data$location_id %in% train_locations)

    # Sort indices to maintain order
    test_idx[[i]] <- sort(test_idx[[i]])
    train_idx[[i]] <- sort(train_idx[[i]])
  }

  # Add some diagnostic information
  attr(test_idx, "n_locations") <- n_locations
  attr(test_idx, "n_groups") <- n_groups
  attr(test_idx, "unique_groups") <- unique_groups
  attr(test_idx, "loc_col") <- loc_col
  attr(test_idx, "group_col") <- group

  return(list(test_idx = test_idx, train_idx = train_idx))
}
