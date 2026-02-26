#' Fit an additive linear mixed effect model over replicates
#'
#'  \code{ngme} function performs an analysis of non-gaussian additive models.
#'  It does the maximum likelihood estimation via stochastic gradient descent.
#'  The prediction of unknown location can be performed by leaving the response
#'  variable to be \code{NA}. The likelihood is specified by \code{family}.
#' The model estimation control can be setted in \code{control} using
#'  \code{control_opt()} function, see \code{?control_opt} for details.
#' See \code{ngme_model_types()} for available models.
#' @param formula formula
#' @param data    a dataframe or a list providing data
#'   (Only response variable can contain \code{NA} value,
#'    \code{NA} value in other columns will cause problem)
#' @param control_opt  control for optimizer. by default it is \code{control_opt()}. See \code{?control_opt} for details.
#' @param control_ngme control for ngme model. by default it is \code{control_ngme()}. See \code{?control_ngme} for details.
#' @param replicate factor, used for divide data into different replicates
#' @param group factor, used for bivariate model, indicating which group the observation belongs to
#' @param family likelihood type, same as measurement noise specification. It can be provided as:
#'   \enumerate{
#'     \item a string, e.g., \code{"normal"}, \code{"nig"}, \code{"t"}.
#'     \item an ngme noise object, e.g., \code{noise_normal()}, \code{noise_nig(mu = 0, sigma = 1, nu = 1)}, \code{noise_t(nu=5)}.
#'   }
#' @param prior_beta prior specification for fixed effects (`beta`), created by
#'   \code{prior_*()} or \code{priors(...)}.
#' @param start  starting ngme object (usually object from last fit)
#' @param moving_window number of iterations to average the estimation
#' @param debug  toggle debug mode
#'
#' @return random effects (for different replicate) + models(fixed effects, measuremnt noise, and latent process)
#' @export
#'
#' @examples
#' ngme(
#'   formula = Y ~ x1 + f(
#'     x2,
#'     model = "ar1",
#'     noise = noise_nig(),
#'     rho = 0.5
#'   ) + f(x1,
#'     model = "rw1",
#'     noise = noise_normal(),
#'   ),
#'   family = noise_normal(sd = 0.5),
#'   data = data.frame(Y = 1:5, x1 = 2:6, x2 = 3:7),
#'   control_opt = control_opt(
#'     estimation = FALSE
#'   )
#' )
ngme <- function(
    formula,
    data,
    family = "normal",
    control_opt = NULL,
    control_ngme = NULL,
    group = NULL,
    replicate = NULL,
    start = NULL,
    moving_window = 1, # return the average estimation of last .. iterations
    prior_beta = NULL,
    debug = FALSE) {
  # -------------  CHECK INPUT ---------------
  if (is.null(data)) {
    stop("Missing data.frame/list `data'. Leaving `data' empty might lead to\n\t\tuncontrolled behaviour, therefore is it required.")
  }
  if (!is.data.frame(data)) {
    stop("\n\tArgument `data' must be a data.frame.")
  }

  # configure control parameters
  if (is.null(control_ngme)) control_ngme <- control_ngme()
  if (is.null(control_opt)) control_opt <- control_opt()
  stopifnot(
    inherits(control_ngme, "control_ngme"),
    inherits(control_opt, "control_opt"),
    "data provide should be of the same length" =
      all(diff(sapply(data, length)) == 0)
  )
  control_ngme <- update_control_ngme(control_ngme, control_opt)

  group <- validate_rep_or_group(group, data)
  replicate <- validate_rep_or_group(replicate, data)

  # model fit information
  fit <- list(
    formula = formula,
    data = data,
    family = family,
    replicate = replicate,
    group = group,
    n_data = nrow(data)
  )
  if (debug) control_ngme$debug <- TRUE

  noise <- if (is.character(family)) {
    switch(family,
      "normal" = noise_normal(),
      "gaussian" = noise_normal(),
      "nig" = noise_nig(),
      "t" = noise_t(),
      "skew_t" = noise_skew_t(),
      stop("Unknown family!")
    )
  } else {
    family
  } # ngme noise object

  stopifnot(class(noise) == "ngme_noise")

  # parse the formula get a list of ngme_replicate
  ngme_model <- ngme_parse_formula(
    formula, data, control_ngme, noise, group, replicate, prior_beta,
    control_opt$standardize_fixed # convert fixed effects
  )

  if (contain_bv_model(ngme_model)) {
    stopifnot(
      "Please supply `group` argument in ngme() function, not in f()." =
        length(levels(group)) == 2
    )
  }
  attr(ngme_model, "fit") <- fit
  attr(ngme_model, "estimation_enabled") <- control_opt$estimation

  # Check if using bfgs for non-Gaussian model
  if (control_opt$sgd_method == "bfgs") {
    stopifnot(
      "Please use other optimizer for non-Gaussian model, BFGS is used for optimizing Gaussian model likelihood." =
        ngme_model$replicates[[1]]$all_gaussian &&
          noise$noise_type == "normal"
    )
  }

  ####### Use Last_fit ngme object to update Rcpp_list
  if (!is.null(start) && !inherits(start, "ngme")) {
    stop("start should be an ngme object.")
  }

  # update with start (list of ngmes)
  if (inherits(start, "ngme")) {
    for (i in seq_along(ngme_model$replicates)) {
      ngme_model$replicates[[i]] <- within(ngme_model$replicates[[i]], {
        # check if feff is the same, then overwrite the feff
        same_feff <- ncol(X) == ncol(start$replicates[[i]]$X)
        if (!is.null(colnames(X)) && !is.null(colnames(start$replicates[[i]]$X))) {
          same_feff <- same_feff && identical(colnames(X), colnames(start$replicates[[i]]$X))
        }
        if (same_feff) {
          start_feff <- start$replicates[[i]]$feff
          if (!is.null(names(start_feff)) &&
            !is.null(colnames(X)) &&
            all(colnames(X) %in% names(start_feff))) {
            start_feff <- start_feff[colnames(X)]
          }

          if (!ngme_model$replicates[[i]]$standardize) {
            feff <- as.numeric(start_feff)
          } else {
            # start_feff is on the raw scale after restoration; map it to the
            # current centered + SVD-standardized parameterization.
            feff_center <- ngme_beta_raw_to_center(start_feff, X_center, X_center_target)
            svd_idx <- if (!is.null(svd) && !is.null(svd$idx)) svd$idx else integer(0)
            if (length(svd_idx) > 0) {
              beta_fixed <- feff_center[svd_idx]
              beta_fixed <- as.numeric(svd$d * drop(t(svd$v) %*% beta_fixed))
              feff_center[svd_idx] <- beta_fixed
              feff <- feff_center
            } else {
              feff <- feff_center
            }
          }
          names(feff) <- colnames(X)
        }

        noise <- update_noise(noise, new_noise = start$replicates[[i]]$noise)
        for (j in seq_along(start$replicates[[i]]$models)) {
          prev_model_type <-
            start$replicates[[i]]$models[[j]]$model
          prev_model_ope <-
            start$replicates[[i]]$models[[j]]$operator

          # update parameter of K
          if (models[[j]]$model == "bv_matern" && prev_model_type == "bv_matern") {
            # update theta_K
            if (!models[[j]]$operator$fix_theta) {
              models[[j]]$theta_K <- models[[j]]$operator$theta_K <-
                c(0, prev_model_ope$theta_K)
            } else {
              models[[j]]$theta_K <- models[[j]]$operator$theta_K <-
                prev_model_ope$theta_K
            }

            # for printing
            models[[j]]$operator$first <- prev_model_ope$first
            models[[j]]$operator$second <- prev_model_ope$second
          } else {
            stopifnot(
              "Please make sure model type are the same" =
                models[[j]]$model == prev_model_type
            )
            # default case: seed new operator with previous estimates without
            # overwriting the new model configuration (e.g., fix_alpha changes)
            op_new <- models[[j]]$operator
            op_prev <- prev_model_ope

            # try name-based matching first so we only overwrite shared params
            if (!is.null(op_new$param_name) && !is.null(op_prev$param_name)) {
              theta_new <- op_new$theta_K
              names(theta_new) <- op_new$param_name

              theta_prev <- op_prev$theta_K
              names(theta_prev) <- op_prev$param_name

              common <- intersect(names(theta_new), names(theta_prev))
              if (length(common) > 0) {
                theta_new[common] <- theta_prev[common]
              }

              op_new$theta_K <- theta_new
              models[[j]]$theta_K <- theta_new
            } else if (length(op_new$theta_K) == length(op_prev$theta_K)) {
              # fall back to position-wise copy when dimensions align
              op_new$theta_K <- op_prev$theta_K
              models[[j]]$theta_K <- op_prev$theta_K
            } else {
              # keep op_new defaults when shapes differ (e.g., freeing alpha)
              models[[j]]$theta_K <- op_new$theta_K
            }

            op_new$K <- op_new$update_K(op_new$theta_K)

            models[[j]]$operator <- op_new
          }

          stopifnot(
            "length of W should be the same" =
              models[[j]]$W_size == length(start$replicates[[i]]$models[[j]]$W)
          )
          stopifnot(
            "length of V should be the same" =
              models[[j]]$V_size == length(start$replicates[[i]]$models[[j]]$noise$V)
          )
          # update the rest
          models[[j]]$W <- start$replicates[[i]]$models[[j]]$W
          models[[j]]$noise$V <- start$replicates[[i]]$models[[j]]$noise$V
          models[[j]]$noise <- update_noise(
            models[[j]]$noise,
            new_noise = start$replicates[[i]]$models[[j]]$noise
          )
        }
      })
    }
  }
  if (debug) {
    print(str(ngme_model$replicates[[1]]))
  }

  # configuration of controls

  # check all f has the same replicate
  ################# Run CPP ####################
  check_dim(ngme_model)
  if (control_opt$estimation) {
    cat("Starting estimation... \n")
    outputs <- tryCatch(
      estimate_cpp(ngme_model, control_opt),
      error = function(err) {
        stop("C++ estimation failed: ", conditionMessage(err), call. = FALSE)
      }
    )
    cat("\n")

    ################# Update the estimates ####################
    est_output <- mean_list(outputs)
    for (i in seq_along(ngme_model$replicates)) {
      ngme_model$replicates[[i]] <- update_ngme_est(ngme_model$replicates[[i]], est_output[[i]])
    }

    # return posterior samples of W and V
    cat("Starting posterior sampling... \n")
    for (i in seq_along(ngme_model$replicates)) {
      res <- tryCatch(
        sampling_cpp(
          ngme_model$replicates[[i]],
          n = control_ngme$n_post_samples,
          n_burnin = 1,
          posterior = TRUE,
          seed = control_opt$seed
        ),
        error = function(err) {
          stop(
            "C++ posterior sampling failed for replicate ", i, ": ",
            conditionMessage(err),
            call. = FALSE
          )
        }
      )

      df_V <- data.frame(res$V)
      colnames(df_V) <- paste0("sample_", 1:ncol(df_V))
      ngme_model$replicates[[i]]$post_V <- df_V

      df_W <- data.frame(res$W)
      colnames(df_W) <- paste0("sample_", 1:ncol(df_W))
      ngme_model$replicates[[i]]$post_W <- df_W

      sd_W <- mean(sapply(res$W, sd))
    }
    cat("Posterior sampling done! \n")
    cat("Average standard deviation of the posterior W: ", sd_W, "\n")

    # print R_hat
    # if (!is.null(attr(outputs, "R_hat"))) {
    #   cat("R_hat: \n")
    #   print(attr(outputs, "R_hat"))
    # }

    cat("Note:
      1. Use ngme_post_samples(..) to access the posterior samples.
      2. Use ngme_result(..) to access different latent models. \n")

    # mn_nu <- ngme_model$replicates[[1]]$noise$nu
    # if (length(mn_nu) > 1 && mn_nu > 100)
    #   cat("The parameter nu for measurement noise is too big, consider using family=normal instead. \n")

    if (isTRUE(control_opt$store_traj)) {
      # Transform trajectory
      traj_df_chains <- transform_traj(attr(outputs, "opt_traj"))
      # dispatch trajs to each latent and block
      idx <- 0
      for (i in seq_along(ngme_model$replicates[[1]]$models)) {
        n_params <- ngme_model$replicates[[1]]$models[[i]]$n_params
        lat_traj_chains <- list()
        for (j in seq_along(traj_df_chains)) {
          lat_traj_chains[[j]] <- traj_df_chains[[j]][idx + 1:n_params, ]
        }

        attr(ngme_model$replicates[[1]]$models[[i]], "lat_traj") <- lat_traj_chains
        idx <- idx + n_params
      }

      # measurement noise and feff
      block_traj <- list()
      n_feff <- length(ngme_model$replicates[[1]]$feff)
      n_chains <- length(traj_df_chains)
      for (j in seq_len(n_chains)) {
        block_traj[[j]] <- traj_df_chains[[j]][(idx + 1):ngme_model$replicates[[1]]$n_params, ]
      }

      n_block_params <- nrow(block_traj[[1]])
      # Map fixed-effects trajectories back to raw parameterization so traceplot
      # is comparable with ngme_result()/printed feff.
      if (n_feff > 0) {
        rep1 <- ngme_model$replicates[[1]]
        feff_idx <- (n_block_params - n_feff + 1):n_block_params
        feff_names <- names(rep1$feff)

        for (i in seq_along(block_traj)) {
          betas <- as.matrix(block_traj[[i]][feff_idx, , drop = FALSE])

          # 1) Undo SVD standardization on the standardized subset.
          if (isTRUE(rep1$standardize) && !is.null(rep1$svd) && length(rep1$svd$d) > 0) {
            svd <- rep1$svd
            svd_idx <- if (!is.null(svd$idx)) svd$idx else seq_len(length(svd$d))
            betas[svd_idx, ] <- svd$v %*% ngme_diag_vec(1 / svd$d) %*% betas[svd_idx, , drop = FALSE]
          }

          # 2) Undo demeaning shift (centered beta -> raw beta), per iteration.
          if (!is.null(rep1$X_center) && length(rep1$X_center) > 0) {
            for (k in seq_len(ncol(betas))) {
              beta_k <- betas[, k]
              names(beta_k) <- feff_names
              beta_k <- ngme_beta_center_to_raw(beta_k, rep1$X_center, rep1$X_center_target)
              betas[, k] <- as.numeric(beta_k)
            }
          }

          block_traj[[i]][feff_idx, ] <- betas
        }
      }

      attr(ngme_model$replicates[[1]], "block_traj") <- block_traj
      attr(outputs, "opt_traj") <- NULL
    } else {
      attr(outputs, "opt_traj") <- NULL
    }

    if (is.list(outputs) && length(outputs) > 1) {
      attr(ngme_model, "chain_outputs") <- outputs
    } else {
      attr(ngme_model, "chain_outputs") <- NULL
    }
  } else {
    # Estimation skipped: map fixed effects and X back to the raw covariate
    # scale so outputs are comparable to pre-demean runs.
    for (i in seq_along(ngme_model$replicates)) {
      ngme_replicate <- ngme_model$replicates[[i]]
      ngme_model$replicates[[i]] <- ngme_restore_fixed_effect_scale(ngme_replicate)
    }
  }
  ngme_model
}

# helper function
# get trajs from a list of estimates
get_trajs <- function(outputs) {
  ret <- list()
  for (i in seq_along(outputs)) {
    ret[[i]] <- list()
    ret[[i]]$block_traj <- attr(outputs[[i]], "trajectory")
    for (j in seq_along(outputs[[i]]$models)) {
      ret[[i]]$models[[j]] <- list()
      ret[[i]]$models[[j]] <- attr(outputs[[i]]$models[[j]], "trajectory")
    }
  }
  ret
}

# helper function to tranform the trajectory
# input: a list (n_chain) of all parameters
# return a list (n_chain) of all parameters transposed
transform_traj <- function(traj) {
  n_chain <- length(traj)
  dfs <- list()
  for (i in 1:n_chain) {
    df <- as.data.frame(traj[[i]])
    names(df) <- NULL
    dfs[[i]] <- df
  }
  dfs
}

ngme_diag_vec <- function(x) {
  x <- as.numeric(x)
  if (length(x) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }
  diag(x, nrow = length(x), ncol = length(x))
}

# use estimate result to update ngme object
update_ngme_est <- function(
    ngme_replicate, est_output) {
  # Fixed effects
  names(est_output$feff) <- names(ngme_replicate$feff)
  ngme_replicate$feff <- est_output$feff
  ngme_replicate$log_likelihood <- est_output$log_likelihood

  ngme_replicate <- ngme_restore_fixed_effect_scale(ngme_replicate)

  ngme_replicate$noise <- update_noise(ngme_replicate$noise, new_noise = est_output$noise)
  for (i in seq_along(ngme_replicate$models)) {
    # update theta_K and K
    theta_K <- est_output$models[[i]]$theta_K
    if ("param_name" %in% names(ngme_replicate$models[[i]]$operator)) {
      names(theta_K) <- ngme_replicate$models[[i]]$operator$param_name
    }
    ngme_replicate$models[[i]]$operator$theta_K <-
      ngme_replicate$models[[i]]$theta_K <- theta_K

    # new_K <- ngme_replicate$models[[i]]$operator$update_K(theta_K)
    # ngme_replicate$models[[i]]$operator$K <- ngme_as_sparse(new_K)

    # update W and noise
    ngme_replicate$models[[i]]$W <- est_output$models[[i]]$W
    ngme_replicate$models[[i]]$noise <- update_noise(
      ngme_replicate$models[[i]]$noise,
      new_noise = est_output$models[[i]]
    )

    # tedious special case
    if (ngme_replicate$models[[i]]$model == "tp") {
      n1 <- ngme_replicate$models[[i]]$operator$first$n_theta_K
      n2 <- ngme_replicate$models[[i]]$operator$second$n_theta_K
      ngme_replicate$models[[i]]$operator$first$theta_K <- ngme_replicate$models[[i]]$theta_K[1:n1]
      ngme_replicate$models[[i]]$operator$second$theta_K <- ngme_replicate$models[[i]]$theta_K[(n1 + 1):(n1 + n2)]

      # update output for tp-bv model
      lat <- ngme_replicate$models[[i]]
      if (lat$operator$second$model %in% c("bv", "bv2", "bv_matern")) {
        bv <- lat$operator$second
        # theta_K = (theta_K_tp_first, theta_K_bv)
        # notice that theta_K_bv is already updated

        # number of bv-specific parameters before sub-model params
        n0 <- switch(bv$model,
          "bv" = (if (!isTRUE(bv$fix_theta)) 1 else 0) + 1 + (if (isTRUE(bv$use_c_param)) 2 else 0),
          "bv2" = 3,
          "bv_matern" = 4 - isTRUE(bv$fix_theta),
          0
        )
        n1 <- bv$first$n_theta_K
        n2 <- bv$second$n_theta_K
        bv$first$theta_K <- bv$theta_K[(n0 + 1):(n0 + n1)]
        bv$second$theta_K <- bv$theta_K[(n0 + n1 + 1):(n0 + n1 + n2)]

        ngme_replicate$models[[i]]$operator$second <- bv
      }
    }

    if (ngme_replicate$models[[i]]$model %in% c("bv", "bv_matern")) {
      n1 <- ngme_replicate$models[[i]]$operator$first$n_theta_K
      n2 <- ngme_replicate$models[[i]]$operator$second$n_theta_K

      n_param_bv <- if (ngme_replicate$models[[i]]$model == "bv") {
        (if (!isTRUE(ngme_replicate$models[[i]]$operator$fix_theta)) 1 else 0) + 1 + (if (isTRUE(ngme_replicate$models[[i]]$operator$use_c_param)) 2 else 0)
      } else {
        4 - isTRUE(ngme_replicate$models[[i]]$operator$fix_theta)
      }
      ngme_replicate$models[[i]]$operator$first$theta_K <- ngme_replicate$models[[i]]$theta_K[(n_param_bv + 1):(n1 + n_param_bv)]
      ngme_replicate$models[[i]]$operator$second$theta_K <- ngme_replicate$models[[i]]$theta_K[(n1 + n_param_bv + 1):(n1 + n2 + n_param_bv)]
    }
  }
  ngme_replicate
}

#' Print an ngme model
#'
#' @param x ngme model object
#' @param ... ...
#'
#' @return a list (model specifications)
#' @export
print.ngme <- function(x, ...) {
  print(x$replicates[[1]])
  if (x$n_repls > 1) {
    cat("Number of global replicates is", x$n_repls, "\n")
  }
}

######
check_dim <- function(ngme_model) {
  for (ngme in ngme_model$replicates) {
    if (ncol(ngme$X) != length(ngme$feff)) {
      stop("The number of columns of X is not equal to the length of feff")
    }
    for (latent in ngme$models) {
      # ncol(A) = W_size
      if (ncol(latent[["A"]]) != latent$W_size) {
        stop("The number of columns of A is not equal to the W_size of the latent model")
      }

      stopifnot(
        nrow(latent$noise$B_sigma) == latent$V_size,
        nrow(latent$noise$B_mu) == latent$V_size
      )
    }
  }
}

ngme_center_non_intercept_cols <- function(X, active_rows = NULL, center_target = NULL) {
  if (ncol(X) == 0) {
    return(list(X = X, center = numeric(0), center_target = character(0)))
  }

  if (is.null(colnames(X))) {
    colnames(X) <- as.character(seq_len(ncol(X)))
  }

  center <- setNames(numeric(ncol(X)), colnames(X))
  center_target_vec <- setNames(rep(NA_character_, ncol(X)), colnames(X))
  has_intercept <- any(startsWith(colnames(X), "(Intercept)"))
  if (!has_intercept) {
    # Keep no-intercept models on their original scale:
    # centering without an intercept changes the model semantics.
    return(list(X = X, center = center, center_target = center_target_vec))
  }

  if (is.null(active_rows)) {
    active_rows <- rep(TRUE, nrow(X))
  } else if (length(active_rows) != nrow(X)) {
    stop("active_rows must have the same length as nrow(X).")
  }

  out <- X
  non_intercept <- !startsWith(colnames(out), "(Intercept)")
  if (any(non_intercept)) {
    active_idx <- which(active_rows)
    if (length(active_idx) > 0) {
      center_vals <- colMeans(out[active_idx, non_intercept, drop = FALSE])
      center[non_intercept] <- center_vals
      out[active_idx, non_intercept] <- sweep(
        out[active_idx, non_intercept, drop = FALSE],
        2,
        center_vals,
        FUN = "-"
      )
    }

    default_target <- center_target
    if (is.null(default_target)) {
      default_target <- if ("(Intercept)" %in% colnames(out)) "(Intercept)" else NA_character_
    }
    center_target_vec[non_intercept] <- default_target
  }

  list(X = out, center = center, center_target = center_target_vec)
}

ngme_guess_intercept_targets <- function(coef_name, intercept_names) {
  if (length(intercept_names) == 0) {
    return(character(0))
  }

  grouped <- setdiff(intercept_names, "(Intercept)")
  if (length(grouped) > 0) {
    suffixes <- sub("^\\(Intercept\\)", "", grouped)
    keep <- nchar(suffixes) > 0
    if (any(keep)) {
      grouped <- grouped[keep]
      suffixes <- suffixes[keep]
      match_idx <- which(endsWith(coef_name, suffixes))
      if (length(match_idx) > 0) {
        # If multiple suffixes match, prefer the most specific (longest).
        best <- match_idx[which.max(nchar(suffixes[match_idx]))]
        return(grouped[best])
      }
    }
  }

  if ("(Intercept)" %in% intercept_names) {
    return("(Intercept)")
  }

  # Legacy fallback: if no unambiguous mapping exists, use all intercepts.
  intercept_names
}

ngme_collect_intercept_shifts <- function(beta, center, center_target = NULL) {
  if (length(beta) == 0 || is.null(names(beta))) {
    return(numeric(0))
  }

  center_beta <- center[names(beta)]
  center_beta[is.na(center_beta)] <- 0
  nonzero <- which(abs(center_beta) > 0)
  if (length(nonzero) == 0) {
    return(numeric(0))
  }

  intercept_names <- names(beta)[startsWith(names(beta), "(Intercept)")]
  if (length(intercept_names) == 0) {
    return(numeric(0))
  }

  shifts <- numeric(0)
  for (idx in nonzero) {
    coef_name <- names(beta)[idx]
    delta <- center_beta[idx] * beta[idx]
    if (abs(delta) < 1e-12) {
      next
    }

    targets <- character(0)
    if (!is.null(center_target)) {
      explicit <- center_target[coef_name]
      explicit <- explicit[!is.na(explicit) & nzchar(explicit)]
      if (length(explicit) > 0) {
        targets <- explicit[[1]]
      }
    }
    if (length(targets) == 0) {
      targets <- ngme_guess_intercept_targets(coef_name, intercept_names)
    }
    if (length(targets) == 0) {
      next
    }

    for (target in targets) {
      if (!(target %in% names(shifts))) {
        shifts[target] <- 0
      }
      shifts[target] <- shifts[target] + delta
    }
  }

  shifts
}

ngme_shift_intercepts <- function(beta, shifts, direction = c("to_raw", "to_center")) {
  direction <- match.arg(direction)
  if (length(beta) == 0 || length(shifts) == 0 || is.null(names(beta))) {
    return(beta)
  }

  sign <- if (direction == "to_raw") -1 else 1
  for (target in names(shifts)) {
    if (target %in% names(beta)) {
      beta[target] <- beta[target] + sign * shifts[target]
    }
  }

  beta
}

ngme_beta_center_to_raw <- function(beta_center, center, center_target = NULL) {
  if (length(beta_center) == 0) {
    return(beta_center)
  }

  shifts <- ngme_collect_intercept_shifts(beta_center, center, center_target)
  ngme_shift_intercepts(beta_center, shifts, direction = "to_raw")
}

ngme_beta_raw_to_center <- function(beta_raw, center, center_target = NULL) {
  if (length(beta_raw) == 0) {
    return(beta_raw)
  }

  shifts <- ngme_collect_intercept_shifts(beta_raw, center, center_target)
  ngme_shift_intercepts(beta_raw, shifts, direction = "to_center")
}

ngme_check_isotropic_normal_beta_prior <- function(prior_beta, idx_svd, tol = 1e-10) {
  if (length(idx_svd) == 0) {
    return(list(ok = TRUE, reason = NULL))
  }
  if (length(prior_beta) < max(idx_svd)) {
    return(list(ok = FALSE, reason = "invalid prior length"))
  }

  priors_fixed <- prior_beta[idx_svd]
  is_normal_coef <- vapply(priors_fixed, function(p) {
    identical(p$type, "normal") && identical(p$target, "coef")
  }, logical(1))
  if (!all(is_normal_coef)) {
    return(list(ok = FALSE, reason = "requires target='coef' normal priors"))
  }

  params <- lapply(priors_fixed, function(p) p$param)
  valid_param <- vapply(params, function(x) {
    is.numeric(x) && length(x) >= 2 &&
      all(is.finite(as.numeric(x[1:2])))
  }, logical(1))
  if (!all(valid_param)) {
    return(list(ok = FALSE, reason = "invalid normal prior parameters"))
  }

  precision <- as.numeric(vapply(params, function(x) x[[2]], numeric(1)))
  if (any(precision <= 0)) {
    return(list(ok = FALSE, reason = "normal prior precision must be positive"))
  }

  ref_precision <- precision[[1]]
  if (!all(abs(precision - ref_precision) <= tol * max(1, abs(ref_precision)))) {
    return(list(ok = FALSE, reason = "requires same variance across standardized fixed effects"))
  }

  list(ok = TRUE, reason = NULL)
}

# Map isotropic normal priors from raw fixed-effect coefficients (X basis)
# to the internal SVD basis (U basis): beta* = D V^T beta.
# Only columns indexed by idx_svd (the SVD-standardized block) are transformed.
ngme_map_isotropic_normal_priors_to_svd <- function(prior_beta, svd, idx_svd, tol = 1e-10) {
  if (is.null(svd) || length(idx_svd) == 0 || length(prior_beta) < max(idx_svd)) {
    return(prior_beta)
  }

  priors_fixed <- prior_beta[idx_svd]
  is_normal_coef <- vapply(priors_fixed, function(p) {
    identical(p$type, "normal") && identical(p$target, "coef")
  }, logical(1))
  if (!all(is_normal_coef)) {
    return(prior_beta)
  }

  params <- lapply(priors_fixed, function(p) p$param)
  valid_param <- vapply(params, function(x) {
    is.numeric(x) && length(x) >= 2 &&
      all(is.finite(as.numeric(x[1:2])))
  }, logical(1))
  if (!all(valid_param)) {
    return(prior_beta)
  }

  mean_raw <- as.numeric(vapply(params, function(x) x[[1]], numeric(1)))
  precision_raw <- as.numeric(vapply(params, function(x) x[[2]], numeric(1)))
  if (any(precision_raw <= 0)) {
    return(prior_beta)
  }

  ref_precision <- precision_raw[[1]]
  if (!all(abs(precision_raw - ref_precision) <= tol * max(1, abs(ref_precision)))) {
    # Only isotropic (same variance) normal priors are mapped.
    return(prior_beta)
  }

  d_fixed <- svd$d
  mean_svd <- as.numeric(d_fixed * drop(t(svd$v) %*% mean_raw))
  precision_svd <- ref_precision / (d_fixed^2)

  for (i in seq_along(idx_svd)) {
    prior_beta[[idx_svd[[i]]]]$param <- c(mean = mean_svd[[i]], precision = precision_svd[[i]])
  }

  prior_beta
}

ngme_restore_fixed_effect_scale <- function(ngme_replicate) {
  center <- ngme_replicate$X_center
  center_target <- ngme_replicate$X_center_target
  if (is.null(center) || length(center) == 0 || ncol(ngme_replicate$X) == 0) {
    return(ngme_replicate)
  }

  X_centered <- ngme_replicate$X
  feff_centered <- ngme_replicate$feff
  has_intercept <- any(startsWith(names(feff_centered), "(Intercept)"))

  if (isTRUE(ngme_replicate$standardize) &&
    !is.null(ngme_replicate$svd) &&
    length(ngme_replicate$svd$d) > 0) {
    svd_idx <- ngme_replicate$svd$idx
    if (is.null(svd_idx)) {
      svd_idx <- seq_len(length(ngme_replicate$svd$d))
    }
    stopifnot(length(svd_idx) == length(ngme_replicate$svd$d))

    feff_centered[svd_idx] <- as.numeric(
      ngme_replicate$svd$v %*%
        ((1 / ngme_replicate$svd$d) * feff_centered[svd_idx])
    )

    nrow_one_repl <- nrow(ngme_replicate$X)
    u_fixed <- ngme_replicate$svd$u
    if (nrow(u_fixed) != nrow_one_repl) {
      idx <- ngme_replicate$data_idx
      if (is.null(idx) || length(idx) != nrow_one_repl) {
        stop(
          "Unable to restore standardized fixed effects: replicate row count ",
          "does not match svd$u and data_idx is unavailable or invalid."
        )
      }
      if (any(idx < 1) || any(idx > nrow(u_fixed))) {
        stop(
          "Unable to restore standardized fixed effects: data_idx is out of ",
          "bounds for svd$u."
        )
      }
      u_fixed <- u_fixed[idx, , drop = FALSE]
    }

    X_fixed <- u_fixed %*% ngme_diag_vec(ngme_replicate$svd$d) %*% t(ngme_replicate$svd$v)
    X_centered[, svd_idx] <- X_fixed
  }

  X_raw <- X_centered
  feff_raw <- feff_centered
  if (has_intercept) {
    center_x <- center[colnames(X_centered)]
    center_x[is.na(center_x)] <- 0

    nonzero <- abs(center_x) > 0
    if (any(nonzero)) {
      intercept_cols <- colnames(X_centered)[startsWith(colnames(X_centered), "(Intercept)")]
      for (col_name in names(center_x)[nonzero]) {
        shift <- center_x[[col_name]]

        target <- NA_character_
        if (!is.null(center_target)) {
          explicit <- center_target[col_name]
          explicit <- explicit[!is.na(explicit) & nzchar(explicit)]
          if (length(explicit) > 0) {
            target <- explicit[[1]]
          }
        }
        if (is.na(target) || !nzchar(target) || !(target %in% colnames(X_centered))) {
          guessed <- ngme_guess_intercept_targets(col_name, intercept_cols)
          if (length(guessed) == 1 && guessed %in% colnames(X_centered)) {
            target <- guessed
          }
        }

        if (!is.na(target) && nzchar(target) && target %in% colnames(X_centered)) {
          X_raw[, col_name] <- X_raw[, col_name] + shift * X_centered[, target]
        } else {
          X_raw[, col_name] <- X_raw[, col_name] + shift
        }
      }
    }

    feff_raw <- ngme_beta_center_to_raw(feff_centered, center, center_target)
  }
  names(feff_raw) <- names(ngme_replicate$feff)

  ngme_replicate$feff <- feff_raw
  ngme_replicate$X <- X_raw
  ngme_replicate
}

#' Parse the formula for ngme function
#'
#' @param fm formula
#' @param data data.frame
#' @param control_ngme control_ngme
#' @param noise noise
#' @param group group factor
#' @param replicate replicate vector
#' @param prior_beta prior specification for fixed effects
#' @param standardize logical, whether fixed effects are standardized
#'
#' @return a list (replicate) of ngme_replicate models
ngme_parse_formula <- function(
    fm,
    data,
    control_ngme,
    noise,
    group,
    replicate,
    prior_beta,
    standardize) {
  enclos_env <- list2env(as.list(parent.frame()), parent = parent.frame(2))
  global_env_first <- list2env(as.list(parent.frame(2)), parent = parent.frame())

  tf <- terms.formula(fm, specials = c("f", "fe"))
  terms <- attr(tf, "term.labels")
  intercept <- attr(tf, "intercept")

  # order of f terms in labels
  f_order <- attr(tf, "specials")$f - 1
  fe_order <- attr(tf, "specials")$fe - 1
  # construct plain formula without f
  # watch out! terms[-double(0)] -> character(0)
  fixf <- if (length(f_order) == 0) terms else terms[-f_order]
  fixf <- if (length(fe_order) == 0) fixf else fixf[-fe_order]

  # construct fixed effect with resposne ~ intercept + fixf
  response <- deparse(attr(tf, "variables")[[2]])
  plain_fm_str <- paste(response, "~", intercept, paste(c("", fixf), collapse = " + "))
  plain_fm <- formula(plain_fm_str)

  # eval the data
  ngme_response <- eval(stats::terms(fm)[[2]], envir = data, enclos = enclos_env)
  stopifnot("Have NA in your response variable" = all(!is.na(ngme_response)))
  X_full <- model.matrix(delete.response(terms(plain_fm)), as.data.frame(data))
  centered <- ngme_center_non_intercept_cols(X_full)
  X_full <- centered$X
  X_center <- centered$center
  X_center_target <- centered$center_target
  if (ncol(X_full) == 0) {
    idx_svd <- integer(0)
  } else {
    if (is.null(colnames(X_full))) {
      colnames(X_full) <- as.character(seq_len(ncol(X_full)))
    }
    idx_svd <- which(!startsWith(colnames(X_full), "(Intercept)"))
  }

  # SVD-standardize only non-intercept fixed effects.
  svd <- NULL
  if (length(idx_svd) >= 1) {
    svd <- svd(X_full[, idx_svd, drop = FALSE])
    if (any((svd$d) < 1e-10)) stop("The design matrix is not full rank.")
  } else {
    standardize <- FALSE
  }

  # adding fixed effect (fe() syntax used for bivariate model..)
  for (i in fe_order) {
    lang <- str2lang(terms[i])
    stopifnot(
      "Please provide which_group=<group_name>" = !is.null(lang$which),
      "Please make sure <group_name> is in group" = lang$which %in% levels(group)
    )
    mask <- group == lang$which
    X_sub <- model.matrix(as.formula(lang[[2]]), data)
    colnames(X_sub) <- paste0(colnames(X_sub), "_", lang$which)
    X_sub[!mask, ] <- 0
    intercept_sub <- colnames(X_sub)[startsWith(colnames(X_sub), "(Intercept)")]
    center_target_sub <- if (length(intercept_sub) == 1) intercept_sub else NA_character_
    centered_sub <- ngme_center_non_intercept_cols(
      X_sub,
      active_rows = mask,
      center_target = center_target_sub
    )
    X_sub <- centered_sub$X
    X_center <- c(X_center, centered_sub$center)
    X_center_target <- c(X_center_target, centered_sub$center_target)
    X_full <- cbind(X_full, X_sub)
  }

  # fe() terms are appended after the initial SVD-based check above.
  # Re-check rank here so aliased grouped columns fail early with a clear error.
  if (ncol(X_full) > 0) {
    qr_X <- qr(X_full)
    if (qr_X$rank < ncol(X_full)) {
      dep_cols <- colnames(X_full)[qr_X$pivot[(qr_X$rank + 1):ncol(X_full)]]
      stop(
        paste0(
          "The fixed-effects design matrix is not full rank after fe() expansion. ",
          "Please remove redundant terms. Non-identifiable columns: ",
          paste(dep_cols, collapse = ", ")
        )
      )
    }
  }

  prior_beta_compiled <- compile_beta_priors(prior_beta, colnames(X_full))
  if (standardize && !is.null(prior_beta)) {
    prior_check <- ngme_check_isotropic_normal_beta_prior(
      prior_beta = prior_beta_compiled,
      idx_svd = idx_svd
    )
    if (!isTRUE(prior_check$ok)) {
      warning(
        paste0(
          "control_opt(standardize_fixed = TRUE) was disabled because custom ",
          "prior_beta on standardized fixed effects is not isotropic normal (",
          prior_check$reason,
          ")."
        ),
        call. = FALSE
      )
      standardize <- FALSE
    }
  }

  if (standardize) {
    colnames(svd$u) <- colnames(X_full)[idx_svd]
    X_full[, idx_svd] <- svd$u
    prior_beta_compiled <- ngme_map_isotropic_normal_priors_to_svd(
      prior_beta = prior_beta_compiled,
      svd = svd,
      idx_svd = idx_svd
    )
    svd$idx <- idx_svd
  }

  # If the user supplies beta_init/feff on the original design scale, map it to
  # the standardized basis used internally (svd$u) while leaving any additional
  # fe() columns untouched.
  user_beta <- if (!is.null(control_ngme$beta_init)) control_ngme$beta_init else control_ngme$feff
  if (!is.null(user_beta)) {
    stopifnot(
      "The length of beta_init must equal the number of fixed-effect columns" =
        length(user_beta) == ncol(X_full)
    )

    # Internal estimation uses centered covariates; map user-provided beta from
    # raw scale to centered scale first.
    user_beta <- ngme_beta_raw_to_center(user_beta, X_center, X_center_target)
  }

  if (standardize && !is.null(user_beta)) {
    beta_fixed <- user_beta[idx_svd]
    beta_fixed <- as.numeric(svd$d * drop(t(svd$v) %*% beta_fixed))
    user_beta[idx_svd] <- beta_fixed
  }

  if (!is.null(user_beta)) {
    control_ngme$beta_init <- user_beta
    control_ngme$feff <- user_beta # legacy name retained for downstream code
  }

  ########## parse models terms
  pre_model <- list()
  all_gaussian <- noise$noise_type == "normal"
  idx_effect <- 1
  idx_field <- 1 # for setting names
  for (i in f_order) {
    str <- gsub("^f\\(", "ngme2::f(", terms[i])
    lang <- str2lang(str)
    if (is.null(lang$model)) stop("Please provide model=<model_name>.")
    # pass extra argument into f
    if (is.null(lang$data)) lang$data <- data
    if (is.null(lang$group)) lang$group <- group
    if (is.null(lang$name) && lang$model == "re") {
      lang$name <- paste0("effect", idx_effect)
      idx_effect <- idx_effect + 1
    }
    if (is.null(lang$name)) {
      lang$name <- paste0("field", idx_field)
      idx_field <- idx_field + 1
    }

    pre_model[[lang$name]] <- lang
  }

  levels <- levels(replicate)
  blocks_rep <- list() # of length n_repl

  # Validate mesh lists if any f() calls use them
  for (tmp in pre_model) {
    if (!is.null(tmp$mesh) && is.list(tmp$mesh) && !inherits(tmp$mesh, c("inla.mesh.1d", "inla.mesh", "fm_mesh_1d", "fm_mesh_2d", "metric_graph"))) {
      n_meshes <- length(tmp$mesh)
      n_repls <- length(levels)
      if (n_meshes < n_repls) {
        stop(paste("Insufficient meshes provided for field '", tmp$name, "'. ",
          "Expected ", n_repls, " meshes for ", n_repls, " replicates, ",
          "but only ", n_meshes, " meshes were provided.",
          sep = ""
        ))
      }
    }
  }

  noise_new <- update_noise(noise, n = length(ngme_response))
  if (noise_new$corr_measurement) {
    stopifnot(
      "Please make sure the len(index_corr) == observations" =
        length(ngme_response) == length(noise_new$index_corr)
    )
  }
  for (level in levels) {
    idx <- replicate == level
    data_idx <- which(idx) # record the original index

    Y <- ngme_response[idx]
    X <- X_full[idx, , drop = FALSE]
    svd_rep <- NULL
    if (standardize) {
      svd_rep <- list(
        u = svd$u[idx, , drop = FALSE],
        d = svd$d,
        v = svd$v,
        idx = svd$idx
      )
    }

    # re-evaluate each f model using idx
    models_rep <- list()
    for (tmp in pre_model) {
      tmp$subset <- idx

      # Evaluate mesh parameter to get actual value
      actual_mesh <- if (is.null(tmp$mesh)) NULL else eval(tmp$mesh, envir = data, enclos = global_env_first)

      # Handle mesh selection for different replicates
      if (
        tmp$model != "spacetime" &&
          !is.null(actual_mesh) &&
          is.list(actual_mesh) &&
          !inherits(actual_mesh, c("inla.mesh.1d", "inla.mesh", "fm_mesh_1d", "fm_mesh_2d", "metric_graph"))
      ) {
        # mesh is a list of meshes for different replicates
        mesh_list <- actual_mesh

        # Convert level to numeric index if needed
        replicate_idx <- which(levels == level)

        # Check if we have enough meshes for this replicate
        if (replicate_idx <= length(mesh_list)) {
          selected_mesh <- mesh_list[[replicate_idx]]
        } else {
          stop(paste("Not enough meshes provided for replicate", level, ". Expected at least", replicate_idx, "meshes, but only", length(mesh_list), "provided."))
        }

        # Replace the mesh parameter in the call with the selected mesh
        tmp$mesh <- selected_mesh

        # Force A matrix to be NULL so it gets rebuilt with the correct mesh
        tmp$A <- NULL
      }

      model_eval <- eval(tmp, envir = data, enclos = global_env_first)

      models_rep[[model_eval$name]] <- model_eval
      if (all(model_eval$noise$noise_type != "normal")) all_gaussian <- FALSE
    }
    # give initial value (whole dataset)
    lm.model <- stats::lm.fit(X_full, ngme_response)
    if (is.null(control_ngme$beta_init) && is.null(control_ngme$feff)) {
      control_ngme$beta_init <- control_ngme$feff <- lm.model$coeff
    }
    if (is.null(noise$theta_sigma)) noise$theta_sigma <- log(sd(lm.model$residuals))
    noise_rep <- subset_noise(noise_new, sub_idx = idx, compute_corr = FALSE)
    group_rep <- group[idx]

    # Re-order according to index_corr!
    # s.t. noise$index_corr=1,1,2,2,3,4,4,....

    # p_oder is the order after permutation
    # original is just 1 2 3, ...
    p_order <- seq_along(Y)
    if (noise$corr_measurement) {
      stopifnot(
        "The length of noise$index_corr should be the same as the number of observations" = length(noise_rep$index_corr) == sum(idx),
        "Now more than 2 locations are correlated in 1 replicate is not allowed" = !any(table(noise_rep$index_corr) > 2)
      )
      p_order <- order(noise_rep$index_corr)
      data_idx <- data_idx[p_order]
      X <- X[p_order, , drop = FALSE]
      Y <- Y[p_order]
      if (standardize) svd_rep$u <- svd_rep$u[p_order, , drop = FALSE]

      group_rep <- group_rep[p_order]
      for (j in seq_along(models_rep)) {
        models_rep[[j]]$A <- models_rep[[j]]$A[p_order, , drop = FALSE]
      }

      # update noise, consider index_corr
      noise_rep <- subset_noise(noise_rep, sub_idx = p_order, compute_corr = TRUE)
    }

    blocks_rep[[level]] <- ngme_replicate(
      data_idx = data_idx,
      Y = Y,
      X = X,
      prior_beta = prior_beta_compiled,
      group = group_rep,
      noise = noise_rep,
      models = models_rep,
      control_ngme = control_ngme,
      n_repl = length(levels),
      all_gaussian = all_gaussian,
      standardize = standardize,
      X_center = X_center,
      X_center_target = X_center_target,
      svd = svd_rep
    )
  }

  n_repls <- length(blocks_rep)
  n_params <- blocks_rep[[1]]$n_params

  structure(
    list(
      replicates   = blocks_rep,
      n_repls      = n_repls,
      n_params     = n_params,
      repls_ngme   = replicate,
      control_ngme = control_ngme
    ),
    class = "ngme"
  )
}

#' Helper function to compute the index_corr vector
#'
#' @param map used as location to compute distance, can be 1d (numerical) or 2d (data.frame)
#' @param eps threshold to determine if two points are close (if close, we consider them as the same point)
#'
#' @return the index_corr vector for ngme correlated measurement noise
#' @examples
#' x_coord <- c(1.11, 1.12, 2, 1.3, 1.3)
#' y_coord <- c(2.11, 2.11, 2, 3.3, 3.3)
#' coord <- data.frame(x_coord, y_coord)
#' compute_index_corr_from_map(map = coord, 0.1)
#' @export
compute_index_corr_from_map <- function(map, eps = 0.1) {
  if (is.null(map)) {
    return(NULL)
  }

  index_corr <- 1:length_map(map)
  if (length(index_corr) == 1) {
    return(index_corr)
  }
  for (i in 2:length_map(map)) {
    for (j in 1:(i - 1)) {
      # compute dist of i and j entry
      d <- dist(sub_map(map, c(i, j)))
      if (d < eps) {
        index_corr[j] <- index_corr[i]
      }
    }
  }

  as.integer(as.factor(index_corr))
}

# idx: integer vector, indicating which observations are correlated
compute_corr_index <- function(idx) {
  stopifnot(sum(round(idx - as.integer(idx))) < 1e8)
  n <- length(idx)
  rows <- cols <- 1:n

  has_correlation <- rep(FALSE, n)
  unique_idx <- unique(idx)
  count <- 1
  for (i in seq_along(unique_idx)) {
    idx_i <- which(idx == unique_idx[i])
    if (length(idx_i) == 1) next
    stopifnot(
      "Now we don't accept measurement noise over 2 places are correlated" = length(idx_i) == 2
    )
    rows[n + count] <- max(idx_i)
    cols[n + count] <- min(idx_i)
    count <- count + 1
    has_correlation[idx_i] <- TRUE
  }

  sort_idx <- order(cols, rows) # col_major order
  list(
    cor_rows = rows[sort_idx] - 1,
    cor_cols = cols[sort_idx] - 1,
    has_correlation = has_correlation,
    n_corr_pairs = sum(has_correlation) / 2
  )
}

#' Summary of ngme fit result
#' @param object an object of class \code{ngme}
#' @param name name of the latent model to be summarized (if NULL, will print all)
#' @param ... other arguments
#'
#' @return a list of summary
#' @export
summary.ngme <- function(
    object,
    name = NULL,
    ...) {
  stopifnot(inherits(object, "ngme"))

  result <- object

  if (!is.null(name)) {
    ngme_rep <- result$replicates[[1]]
    stopifnot(
      "Please provide the correct name of the model" =
        name %in% names(ngme_rep$models)
    )
    result <- ngme_rep$models[[name]]
  }

  result
}

#' Access the result of a ngme fitted model
#' @param ngme_object a ngme fitted model object
#' @param model latent model name to filter by (e.g., "my_ar1", "my_matern", "my_rw1", etc.), "data" for fixed effects and measurement noise, if NULL, return all models
#' @param transformed logical, if TRUE (default) return transformed parameters, if FALSE return raw parameters
#'
#' @return a list of parameters for the specified model, or all models if model is NULL
#' @export
#'
#' @examples
#' \dontrun{
#' # Fit a simple AR(1) model
#' Y <- 1:10
#' n_obs <- length(Y)
#' x1 <- runif(n_obs)
#' x2 <- rexp(n_obs)
#'
#' ngme_out <- ngme(
#'   Y ~ x1 + x2 + f(
#'     1:n_obs,
#'     name = "my_ar",
#'     model = "ar1",
#'     rho = 0.5,
#'     noise = noise_nig(mu = 2, sigma = 3, nu = 1)
#'   ),
#'   data = data.frame(x1 = x1, x2 = x2, Y = Y)
#' )
#'
#' # Get all model parameters (transformed)
#' all_params <- ngme_result(ngme_out)
#' # Returns: list(my_ar = list(rho = 0.5, mu = 2, sigma = 3, nu = 1),
#' #               data = list(fixed_effects = c(...), sigma = 0.5))
#'
#' # Get parameters for specific latent model
#' ar_params <- ngme_result(ngme_out, model = "my_ar")
#' # Returns: list(rho = 0.5, mu = 2, sigma = 3, nu = 1)
#'
#' # Get raw (untransformed) parameters
#' ar_raw <- ngme_result(ngme_out, model = "my_ar", transformed = FALSE)
#' # Returns: list(theta_rho = 1.099, mu = 2, sigma = 3, nu = 1)
#'
#' # Get fixed effects and measurement noise
#' data_params <- ngme_result(ngme_out, model = "data")
#' # Returns: list(fixed_effects = c(...), sigma = 0.5, ...)
#'
#' # For models with multiple latent processes
#' ngme_out2 <- ngme(
#'   Y ~ x1 + x2 + f(
#'     1:n_obs,
#'     name = "my_ar",
#'     model = "ar1",
#'     noise = noise_nig(mu = 2, sigma = 3, nu = 1)
#'   ) + f(
#'     1:n_obs,
#'     name = "my_ou",
#'     model = "ou",
#'     noise = noise_normal(sigma = 1)
#'   ),
#'   data = data.frame(x1 = x1, x2 = x2, Y = Y)
#' )
#'
#' # Get all models
#' all_models <- ngme_result(ngme_out2)
#' # Returns: list(my_ar = list(...), my_ou = list(...), data = list(...))
#'
#' # Get specific model
#' ou_params <- ngme_result(ngme_out2, model = "my_ou")
#' # Returns: list(theta_K1 = 0.5, sigma = 1)
#' }
#'
#' @seealso \code{\link{extract_parameters}} for the underlying function
ngme_result <- function(
    ngme_object,
    model = NULL,
    transformed = TRUE) {
  stopifnot(inherits(ngme_object, "ngme"))

  # Extract all parameters using extract_parameters
  all_params <- extract_parameters(ngme_object)

  # Choose transformed or raw parameters
  params <- if (transformed) all_params$transformed else all_params$raw

  # Return specific model or all models
  if (is.null(model)) {
    return(params)
  } else {
    stopifnot(
      "Please provide the correct model name" = model %in% names(params)
    )
    return(params[[model]])
  }
}
