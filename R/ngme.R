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
#' @param family likelihood type, same as measurement noise specification, 1. string 2. ngme noise obejct
#' @param start  starting ngme object (usually object from last fit)
#' @param moving_window number of iterations to average the estimation
#' @param debug  toggle debug mode
#'
#' @return random effects (for different replicate) + models(fixed effects, measuremnt noise, and latent process)
#' @export
#'
#' @examples
#' ngme(
#'  formula = Y ~ x1 + f(
#'    x2,
#'    model = "ar1",
#'    noise = noise_nig(),
#'    rho = 0.5
#'  ) + f(x1,
#'    model = "rw1",
#'    noise = noise_normal(),
#'  ),
#'  family = noise_normal(sd = 0.5),
#'  data = data.frame(Y = 1:5, x1 = 2:6, x2 = 3:7),
#'  control_opt = control_opt(
#'    estimation = FALSE
#'  )
#')
ngme <- function(
  formula,
  data,
  family        = "normal",
  control_opt   = NULL,
  control_ngme  = NULL,
  group         = NULL,
  replicate     = NULL,
  start         = NULL,
  moving_window = 1, # return the average estimation of last .. iterations
  debug         = FALSE
) {
   # -------------  CHECK INPUT ---------------
  if (is.null(data)) {
    stop("Missing data.frame/list `data'. Leaving `data' empty might lead to\n\t\tuncontrolled behaviour, therefore is it required.")
  }
  if (!is.data.frame(data)) {
    stop("\n\tArgument `data' must be a data.frame.")
  }

  # configure control parameters
  if (is.null(control_ngme)) control_ngme <- control_ngme()
  if (is.null(control_opt))  control_opt <- control_opt()
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

  noise <- if (is.character(family)) switch(family,
    "normal" = noise_normal(),
    "gaussian" = noise_normal(),
    "nig"    = noise_nig(),
    stop("Unknown family!")
  ) else family # ngme noise object

  stopifnot(class(noise) == "ngme_noise")

  # parse the formula get a list of ngme_replicate
  ngme_model <- ngme_parse_formula(
    formula, data, control_ngme, noise, group, replicate,
    control_opt$standardize_fixed # convert fixed effects
  )
  
  if (contain_bv_model(ngme_model)) {
    stopifnot("Please supply `group` argument in ngme() function, not in f()." =
      length(levels(group)) == 2)
  }
  attr(ngme_model, "fit") <- fit
  
  # Check if using bfgs for non-Gaussian model
  if (control_opt$sgd_method == "bfgs") {
    stopifnot(
      "Please use other optimizer for non-Gaussian model, BFGS is used for optimizing Gaussian model likelihood." =
      ngme_model$replicates[[1]]$all_gaussian &&
      noise$noise_type == "normal"
    )
  }

  ####### Use Last_fit ngme object to update Rcpp_list
  if (!is.null(start) && !inherits(start, "ngme"))
   stop("start should be an ngme object.")

  # update with start (list of ngmes)
  if (inherits(start, "ngme")) {
    for (i in seq_along(ngme_model$replicates)) {
      ngme_model$replicates[[i]] <- within(ngme_model$replicates[[i]], {
        # check if feff is the same, then overwrite the feff
        same_feff <- all(dim(X) == dim(start$replicates[[i]]$X))
        if (same_feff) {
          if (!ngme_model$replicates[[i]]$standardize)
            feff <- start$replicates[[i]]$feff
          else {
            # notice that the SVD effects
            svd <- start$replicates[[i]]$svd
            feff <- as.numeric(t(svd$v) %*% start$replicates[[i]]$feff) * svd$d
          }
        }
        
        noise <- update_noise(noise, new_noise = start$replicates[[i]]$noise)
        for (j in seq_along(start$replicates[[i]]$models)) {
          prev_model_type <- 
            start$replicates[[i]]$models[[j]]$model
          prev_model_ope <- 
            start$replicates[[i]]$models[[j]]$operator
          
          # update parameter of K
          if (models[[j]]$model == "bv_matern_nig" && prev_model_type == "bv_matern_normal") {
            # update theta_K
            models[[j]]$theta_K <- models[[j]]$operator$theta_K <- 
              c(0, prev_model_ope$theta_K)
            
            # for printing
            models[[j]]$operator$first <- prev_model_ope$first
            models[[j]]$operator$second <- prev_model_ope$second
          } else {
            stopifnot(
              "Please make sure model type are the same" =
              models[[j]]$model == prev_model_type)
            # default case
            # update operator representation
            models[[j]]$operator <- prev_model_ope

            models[[j]]$theta_K  <- models[[j]]$operator$theta_K
          }

stopifnot(
  "length of W should be the same" =
  models[[j]]$W_size == length(start$replicates[[i]]$models[[j]]$W))
stopifnot(
  "length of V should be the same" =
  models[[j]]$V_size == length(start$replicates[[i]]$models[[j]]$noise$V))
          # update the rest
          models[[j]]$W        <- start$replicates[[i]]$models[[j]]$W
          models[[j]]$noise$V  <- start$replicates[[i]]$models[[j]]$noise$V
          models[[j]]$noise    <- update_noise(
            models[[j]]$noise, new_noise = start$replicates[[i]]$models[[j]]$noise
          )
        }
      })
    }
  }
if (debug) {print(str(ngme_model$replicates[[1]]))}

# configuration of controls

# check all f has the same replicate
  ################# Run CPP ####################
  check_dim(ngme_model)
  if (control_opt$estimation) {
    cat("Starting estimation... \n")
    outputs <- estimate_cpp(ngme_model, control_opt)
    cat("\n")

  ################# Update the estimates ####################
    est_output <- mean_list(outputs)
    for (i in seq_along(ngme_model$replicates))
      ngme_model$replicates[[i]] <- update_ngme_est(ngme_model$replicates[[i]], est_output[[i]])

    # return posterior samples of W and V
    cat("Starting posterior sampling... \n")
    for (i in seq_along(ngme_model$replicates)) {
      res <- sampling_cpp(
        ngme_model$replicates[[i]],
        n = control_ngme$n_post_samples,
        n_burnin = 1,
        posterior = TRUE,
        seed = control_opt$seed
      )

      df_V <- data.frame(res$V)
      colnames(df_V) <- paste0("sample_", 1:ncol(df_V))
      ngme_model$replicates[[i]]$post_V <- df_V

      df_W <- data.frame(res$W)
      colnames(df_W) <- paste0("sample_", 1:ncol(df_W))
      ngme_model$replicates[[i]]$post_W <- df_W
    }
    cat("Posterior sampling done! \n")
    cat("Note:
      1. Use ngme_post_samples(..) to access the posterior samples.
      2. Use ngme_result(..) to access different latent models. \n"
    )

    # mn_nu <- ngme_model$replicates[[1]]$noise$nu
    # if (length(mn_nu) > 1 && mn_nu > 100)
    #   cat("The parameter nu for measurement noise is too big, consider using family=normal instead. \n")

    # Transform trajectory
    traj_df_chains <- transform_traj(attr(outputs, "opt_traj"))
    # dispatch trajs to each latent and block
      idx <- 0;
      for (i in seq_along(ngme_model$replicates[[1]]$models)) {
        n_params <- ngme_model$replicates[[1]]$models[[i]]$n_params
        lat_traj_chains = list()
        for (j in seq_along(traj_df_chains))
          lat_traj_chains[[j]] <- traj_df_chains[[j]][idx + 1:n_params, ]

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
      # update feff (if using svd)
      if (ngme_model$replicates[[1]]$standardize) {
        svd <- ngme_model$replicates[[1]]$svd
        # loop over num. of chains
        for (i in seq_along(block_traj)) {
          # last n_feff rows are fixed effects
          feff_idx <- (n_block_params - n_feff + 1):n_block_params
          betas = as.matrix(block_traj[[i]][feff_idx, ])
          block_traj[[i]][feff_idx,] = svd$v %*% diag(1/svd$d) %*% betas
        }
      }

      attr(ngme_model$replicates[[1]], "block_traj") <- block_traj
      attr(outputs, "opt_traj") <- NULL
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

# use estimate result to update ngme object
update_ngme_est <- function(
  ngme_replicate, est_output
) {
  # Fixed effects
  names(est_output$feff) <- names(ngme_replicate$feff)
  ngme_replicate$feff <- est_output$feff

  if (ngme_replicate$standardize) {
    # standardize feff (transform back)
    feff = as.numeric(ngme_replicate$svd$v %*% (1/ngme_replicate$svd$d * ngme_replicate$feff))
    names(feff) <- names(ngme_replicate$feff)
    ngme_replicate$feff <- feff

    # convert U back to UDV^t
    X <- ngme_replicate$svd$u %*% diag(ngme_replicate$svd$d) %*% t(ngme_replicate$svd$v)
    colnames(X) <- colnames(ngme_replicate$X)
    ngme_replicate$X <- X
    ngme_replicate$X
  }

  ngme_replicate$noise <- update_noise(ngme_replicate$noise, new_noise = est_output$noise)
  for (i in seq_along(ngme_replicate$models)) {
    # update theta_K and K
    theta_K <- est_output$models[[i]]$theta_K
    ngme_replicate$models[[i]]$operator$theta_K <-
      ngme_replicate$models[[i]]$theta_K <- theta_K

    new_K <- ngme_replicate$models[[i]]$operator$update_K(theta_K)
    ngme_replicate$models[[i]]$operator$K <- ngme_as_sparse(new_K)

    # update W and noise
    ngme_replicate$models[[i]]$W        <- est_output$models[[i]]$W
    ngme_replicate$models[[i]]$noise    <- update_noise(
      ngme_replicate$models[[i]]$noise, new_noise = est_output$models[[i]]
    )

    # tedious special case
    if (ngme_replicate$models[[i]]$model == "tp") {
      n1 <- ngme_replicate$models[[i]]$operator$first$n_theta_K
      n2 <- ngme_replicate$models[[i]]$operator$second$n_theta_K
      ngme_replicate$models[[i]]$operator$first$theta_K <- ngme_replicate$models[[i]]$theta_K[1:n1]
      ngme_replicate$models[[i]]$operator$second$theta_K <- ngme_replicate$models[[i]]$theta_K[(n1+1):(n1+n2)]

      # update output for tp-bv model
      lat <- ngme_replicate$models[[i]]
      if (lat$operator$second$model  %in% c("bv", "bv_normal")) {
        bv = lat$operator$second

        n1 = bv$first$n_theta_K
        n2 = bv$second$n_theta_K
        bv$first$theta_K = bv$theta_K[3:(n1+2)]
        bv$second$theta_K = bv$theta_K[(n1+3):(2+n1+n2)]

        ngme_replicate$models[[i]]$operator$second <- bv
      }
    }

    if (ngme_replicate$models[[i]]$model %in% c("bv", "bv_normal", "bv_matern_normal", "bv_matern_nig")) {
      n1 <- ngme_replicate$models[[i]]$operator$first$n_theta_K
      n2 <- ngme_replicate$models[[i]]$operator$second$n_theta_K

      n_param_bv <- switch(ngme_replicate$models[[i]]$model, 
        "bv" = 2,
        "bv_matern_nig" = 4,
        3
      )
      ngme_replicate$models[[i]]$operator$first$theta_K <- ngme_replicate$models[[i]]$theta_K[(n_param_bv+1) : (n1+n_param_bv)]
      ngme_replicate$models[[i]]$operator$second$theta_K <- ngme_replicate$models[[i]]$theta_K[(n1+n_param_bv+1) : (n1+n2+n_param_bv)]
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
  cat("Number of replicates is ", x$n_repls, "\n");
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

#' Parse the formula for ngme function
#'
#' @param fm formula
#' @param data data.frame
#' @param control_ngme control_ngme
#' @param noise noise
#' @param group group factor
#' @param replicate replicate vector
#'
#' @return a list (replicate) of ngme_replicate models
ngme_parse_formula <- function(
  fm,
  data,
  control_ngme,
  noise,
  group,
  replicate,
  standardize
) {
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

  # Do SVD if ncol > 1
  if (ncol(X_full) > 1) {
    svd <- svd(X_full)
    if (any((svd$d) < 1e-10)) stop("The design matrix is not full rank.")
  } else {
    standardize <- FALSE
  }

  if (standardize) {
    colnames(svd$u) <- colnames(X_full)
    X_full <- svd$u  # do regression wrt U
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
    X_full <- cbind(X_full, X_sub)
  }

  ########## parse models terms
  pre_model <- list(); all_gaussian <- noise$noise_type == "normal"
  idx_effect = 1; idx_field = 1; # for setting names
  for (i in f_order) {
    str <- gsub("^f\\(", "ngme2::f(", terms[i])
    lang <- str2lang(str)
    if (is.null(lang$model)) stop("Please provide model=<model_name>.")
    # pass extra argument into f
    if (is.null(lang$data)) lang$data <- data
    if (is.null(lang$group)) lang$group <- group
    if (is.null(lang$name) && lang$model == "re")
      {lang$name <- paste0("effect", idx_effect); idx_effect <- idx_effect + 1}
    if (is.null(lang$name))
      {lang$name <- paste0("field", idx_field); idx_field <- idx_field + 1}

    pre_model[[lang$name]] <- lang
  }

  levels <- levels(replicate)
  blocks_rep <- list() # of length n_repl

  noise_new <- update_noise(noise, n = length(ngme_response))
  if (noise_new$corr_measurement) {
      stopifnot("Please make sure the len(index_corr) == observations" =
        length(ngme_response) == length(noise_new$index_corr))
  }
  for (level in levels) {
    idx <- replicate == level
    data_idx <- which(idx) # record the original index

    Y <- ngme_response[idx]
    X <- X_full[idx, , drop = FALSE]

    # re-evaluate each f model using idx
    models_rep <- list();
    for (tmp in pre_model) {
      tmp$subset <- idx
      model_eval <- eval(tmp, envir = data, enclos = global_env_first)
      models_rep[[model_eval$name]] <- model_eval
      if (all(model_eval$noise$noise_type != "normal")) all_gaussian <- FALSE
    }
    # give initial value (whole dataset)
    lm.model <- stats::lm.fit(X_full, ngme_response)
    if (is.null(control_ngme$feff)) control_ngme$feff <- lm.model$coeff
    if (is.null(noise$theta_sigma)) noise$theta_sigma <- log(sd(lm.model$residuals))
    noise_rep <- subset_noise(noise_new, sub_idx=idx, compute_corr=FALSE)
    group_rep <- group[idx]

    # Re-order according to index_corr!
    # s.t. noise$index_corr=1,1,2,2,3,4,4,....
    
    # p_oder is the order after permutation
    # original is just 1 2 3, ...
    p_order <- seq_along(Y)
    if (noise$corr_measurement) {
      stopifnot(
        "The length of noise$index_corr should be the same as the number of observations"
          = length(noise_rep$index_corr) == sum(idx),
        "Now more than 2 locations are correlated in 1 replicate is not allowed"
          = !any(table(noise_rep$index_corr) > 2)
      )
      p_order <- order(noise_rep$index_corr)
      data_idx <- data_idx[p_order]
      X <- X[p_order, , drop = FALSE]
      Y <- Y[p_order]
      if (standardize) svd$u <- svd$u[p_order, , drop = FALSE]

      group_rep <- group_rep[p_order]
      for (j in seq_along(models_rep))
        models_rep[[j]]$A <- models_rep[[j]]$A[p_order, , drop = FALSE]

      # update noise, consider index_corr
      noise_rep <- subset_noise(noise_rep, sub_idx = p_order, compute_corr = TRUE)
    }

    blocks_rep[[level]] <- ngme_replicate(
      data_idx = data_idx,
      Y = Y,
      X = X,
      group = group_rep,
      noise = noise_rep,
      models = models_rep,
      control_ngme = control_ngme,
      n_repl = length(levels),
      all_gaussian = all_gaussian,
      standardize = standardize,
      svd = svd
    )
  }

  n_repls  <- length(blocks_rep)
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
#' coord = data.frame(x_coord, y_coord)
#' compute_index_corr_from_map(map = coord, 0.1)
#' @export
compute_index_corr_from_map <- function(map, eps=0.1) {
  if (is.null(map)) return(NULL)

  index_corr <- 1:length_map(map)
  if (length(index_corr) == 1) return(index_corr)
  for (i in 2:length_map(map)) {
    for (j in 1:(i-1)) {
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
    stopifnot("Now we don't accept measurement noise over 2 places are correlated"
      = length(idx_i) == 2)
    rows[n+count] <- max(idx_i)
    cols[n+count] <- min(idx_i)
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
  ...
) {
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

#' ngme fit result
#' @param ngme_object a ngme model
#' @param name name of the latent model to be summarized (if NULL, will print all)
#' @param replicate replicate number
#'
#' @return a list of summary
#' @export
ngme_result <- function(
  ngme_object,
  name = NULL,
  replicate = 1
) {
  stopifnot(inherits(ngme_object, "ngme"))
  result <- ngme_object$replicates[[replicate]]

  if (!is.null(name)) {
    names <- sapply(result$models, function(x) x$name)
    stopifnot(
      "Please provide only one name" = length(name) == 1,
      "Please provide the correct name of the model" =
      name %in% names
    )
    result <- result$models[[name]]
  }

  result
}