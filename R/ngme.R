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
#' @param group group factor, used for multivariate model
#' @param family likelihood type, same as measurement noise specification, 1. string 2. ngme noise obejct
#' @param start  starting ngme object (usually object from last fitting)
#' @param corr_measure estimate the correlation between measurement noise (only for bivariate model case)
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
#'    theta_K = 0.5
#'  ) + f(x1,
#'    model = "rw",
#'    order = 1,
#'    circular = TRUE,
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
  start         = NULL,
  corr_measure  = FALSE,
  debug         = FALSE
) {
   # -------------  CHECK INPUT ---------------
  if (is.null(data)) {
    stop("Missing data.frame/list `data'. Leaving `data' empty might lead to\n\t\tuncontrolled behaviour, therefore is it required.")
  }
  if (!is.data.frame(data)) {
    stop("\n\tArgument `data' must be a data.frame.")
  }

  if (is.null(control_ngme)) control_ngme <- control_ngme()
  if (is.null(control_opt))  control_opt <- control_opt()
  stopifnot(inherits(control_ngme, "control_ngme"))
  stopifnot(inherits(control_opt, "control_opt"))
  stopifnot("data provide should be of the same length" =
    all(diff(sapply(data, length)) == 0)
  )

  if (!is.null(group)) group <- as.factor(group)

  # model fitting information
  fitting <- list(
    formula = formula,
    data = data,
    family = family
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
  ngme_model <- ngme_parse_formula(formula, data, control_ngme, noise, group, corr_measure)
  attr(ngme_model, "fitting") <- fitting

  ####### Use Last_fit ngme object to update Rcpp_list
  stopifnot("start should be an ngme object"
    = inherits(start, "ngme") || is.null(start))

  # update with start (list of ngmes)
  if (inherits(start, "ngme")) {
    for (i in seq_along(ngme_model$replicates)) {
      ngme_model$replicates[[i]] <- within(ngme_model$replicates[[i]], {
        beta <- start[[i]]$beta
        noise <- update_noise(noise, new_noise = start[[i]]$noise)
        for (i in seq_along(start[[i]]$models)) {
          models[[i]]$theta_K  <- start[[i]]$models[[i]]$theta_K
          models[[i]]$W        <- start[[i]]$models[[i]]$W
          models[[i]]$noise    <- update_noise(
            models[[i]]$noise, new_noise = start[[i]]$models[[i]]$noise
          )
        }
      })
    }
  }
if (debug) {print(str(ngme_model$replicates[[1]]))}

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

    # return the mean of samples of W of posterior
    # cat("Starting posterior sampling... \nNote: Use ngme$models[[model_name]]$W  to access the posterior mean of process \n")
    # for (i in seq_along(ngme_model$replicates)) {
    #   ngme_replicate <- ngme_model$replicates[[i]]
    #   ngme_replicate$control_ngme$init_sample_W <- FALSE
    #   mean_post_W <- mean_list(
    #     sampling_cpp(ngme_replicate, control_ngme$post_samples_size, TRUE, control_opt$seed)[["W"]]
    #   )

    #   idx <- 1
    #   for (j in seq_along(ngme_replicate$models)) {
    #     ngme_replicate$models[[j]]$W <- mean_post_W[idx : (ngme_replicate$models[[j]]$W_size + idx - 1)]
    #     idx <- idx + ngme_replicate$models[[j]]$W_size
    #   }
    #   ngme_model$replicates[[i]] <- ngme_replicate
    # }
    cat("Posterior sampling done! \n")

    # transform trajectory
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
      # mn and beta
      block_traj <- list()
      for (j in seq_along(traj_df_chains))
        block_traj[[j]] <- traj_df_chains[[j]][(idx + 1):ngme_model$replicates[[1]]$n_params, ]
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

update_ngme_est <- function(
  ngme_replicate, est_output
) {
  ngme_replicate$beta <- est_output$beta
  ngme_replicate$noise <- update_noise(ngme_replicate$noise, new_noise = est_output$noise)
  for (i in seq_along(ngme_replicate$models)) {
    ngme_replicate$models[[i]]$operator$theta_K  <- ngme_replicate$models[[i]]$theta_K <- est_output$models[[i]]$theta_K
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
    }

    if (ngme_replicate$models[[i]]$model == "bv") {
      n1 <- ngme_replicate$models[[i]]$operator$first$n_theta_K
      n2 <- ngme_replicate$models[[i]]$operator$second$n_theta_K
      ngme_replicate$models[[i]]$operator$first$theta_K <- ngme_replicate$models[[i]]$theta_K[3:(n1+2)]
      ngme_replicate$models[[i]]$operator$second$theta_K <- ngme_replicate$models[[i]]$theta_K[(n1+3):(2+n1+n2)]
    }
  }
  ngme_replicate
}

#' @export
print.ngme <- function(x, ...) {
  print(x$replicates[[1]])
  cat("Number of replicates is ", x$n_repls, "\n");
}

######
check_dim <- function(ngme_model) {
  for (ngme in ngme_model$replicates) {
    if (ncol(ngme$X) != length(ngme$beta)) {
      stop("The number of columns of X is not equal to the length of beta")
    }
    for (latent in ngme$models) {
        if (latent$V_size != latent$noise$n_noise) {
          stop("The V_size of the latent model is not equal to the length of noise")
        }
        # ncol(A) = W_size
        if (ncol(latent$A) != latent$W_size) {
          stop("The number of columns of A is not equal to the W_size of the latent model")
        }
        if (!all(latent$operator$h == latent$noise$h)) {
          stop("The h of the latent model is not equal to the h of the noise")
        }

        stopifnot(nrow(latent$noise$B_sigma) == latent$noise$n_noise)
    }
  }
}

#' Parse the formula for ngme function
#'
#' @param fm Formula
#' @param data data.frame
#' @param control_ngme control_ngme
#' @param noise noise
#' @param group group factor
#' @param corr_measure logical, whether to estimate correlation
#'
#' @return a list (replicate) of ngme_replicate models
ngme_parse_formula <- function(
  fm,
  data,
  control_ngme,
  noise,
  group,
  corr_measure
) {
  enclos_env <- list2env(as.list(parent.frame()), parent = parent.frame(2))
  global_env_first <- list2env(as.list(parent.frame(2)), parent = parent.frame())

  tf <- terms.formula(fm, specials = c("f"))
  terms <- attr(tf, "term.labels")
  intercept <- attr(tf, "intercept")

  # order of f terms in labels
  spec_order <- attr(tf, "specials")$f - 1

  # construct plain formula without f
  # watch out! terms[-double(0)] -> character(0)
  fixf <- if (length(spec_order) == 0) terms else terms[-spec_order]

  response <- deparse(attr(tf, "variables")[[2]])
  plain_fm_str <- paste(response, "~", intercept, paste(c("", fixf), collapse = " + "))
  plain_fm <- formula(plain_fm_str)

  # eval the data
  ngme_response <- eval(stats::terms(fm)[[2]], envir = data, enclos = enclos_env)
  stopifnot("Have NA in your response variable" = all(!is.na(ngme_response)))
  X_full    <- model.matrix(delete.response(terms(plain_fm)), as.data.frame(data))

  ########## parse models terms
  pre_model <- list();
  idx_effect = 1; idx_field = 1; # for setting names
  for (i in spec_order) {
    str <- gsub("^f\\(", "ngme2::f(", terms[i])
    lang <- str2lang(str)

    # pass extra argument into f
    if (is.null(lang$data)) lang$data <- data
    if (is.null(lang$group)) lang$group <- group
    res <- eval(lang, envir = global_env_first)

    # give default name
    if (res$name == "field") {res$name <- paste0("field", idx_field); idx_field <- idx_field + 1}
    if (res$name == "effect") {res$name <- paste0("effect", idx_effect); idx_effect <- idx_effect + 1}
    pre_model[[res$name]] <- res
  }
  # splits f_eff, latent according to replicates
  repls <- if (length(pre_model) > 0)
      lapply(pre_model, function(x) x$replicate)
    else
      list(rep(1, length(ngme_response)))

  repl <- merge_repls(repls)
  uni_repl <- unique(repl)
  blocks_rep <- list() # of length n_repl
  for (i in seq_along(uni_repl)) {
    idx <- repl == uni_repl[[i]]
    # data
    Y <- ngme_response[idx]
    X <- X_full[idx, , drop = FALSE]

    # re-evaluate each f model using idx
    models_rep <- list();
    for (tmp in pre_model) {
      tmp$map <- sub_map(tmp$map, idx)
      tmp$replicate <- tmp$replicate[idx]
      tmp$eval = TRUE
      tmp$data <- data[idx, , drop = FALSE]
      model_eval <- eval(tmp, envir = data, enclos = global_env_first)
      models_rep[[model_eval$name]] <- model_eval
    }

    # give initial value (whole dataset)
    lm.model <- stats::lm.fit(X_full, ngme_response)
    if (is.null(control_ngme$beta)) control_ngme$beta <- lm.model$coeff
    if (is.null(noise$theta_sigma)) noise$theta_sigma <- log(sd(lm.model$residuals))

    noise_new <- update_noise(noise, n = length(Y))

    if (corr_measure) {
      stopifnot(
        "Please provide the group vector" = !is.null(group),
        "length of unique group should equal to 2" = length(levels(group)) == 2
      )

      bv_idx <- which(sapply(models_rep, function(x) x$model) == "bv")
      stopifnot("Please make sure there is 1 bivariate (bv) model" =
        length(bv_idx) == 1)
      bv_map <- models_rep[[bv_idx]]$map
      cov_rc <- cov_row_col(bv_map, group)
    }

    blocks_rep[[i]] <- ngme_replicate(
      Y = Y,
      X = X,
      noise = noise_new,
      models = models_rep,
      replicate = uni_repl[[i]],
      control_ngme = control_ngme,
      n_repl = length(uni_repl),
      corr_measure = corr_measure,
      cov_rows = if (corr_measure) cov_rc$cov_rows else NULL,
      cov_cols = if (corr_measure) cov_rc$cov_cols else NULL,
      mark_cov = if (corr_measure) cov_rc$mark_cov else NULL
    )
  }

  n_repls  <- length(blocks_rep)
  n_params <- blocks_rep[[1]]$n_params

  structure(
    list(
      replicates   = blocks_rep,
      n_repls      = n_repls,
      n_params     = n_params,
      repls_ngme   = repls,
      control_ngme = control_ngme
    ),
    class = "ngme"
  )
}

# bv_map: numeric or matrix of 2 col
# group: factor
# eps: consider correlation if |map1-map2| < eps
cov_row_col <- function(bv_map, group, eps = 0.05) {
  stopifnot(
    length_map(bv_map) == length(group)
  )
  group <- as.factor(group)

  n <- length(bv_map)
  idx1 <- group == levels(group)[1]
  idx2 <- group == levels(group)[2]

  # subset bv_map
  map1 <- sub_map(bv_map, idx1)
  map2 <- sub_map(bv_map, idx2)
  corr_map <- intersect(map1, map2)

  cov_cols <- cov_rows <- 1:n
  # cov_cols <- cov_rows <- double()
  for (i in seq_along(corr_map)) {
    rc <- which(bv_map == corr_map[i])
    cov_rows[n+i] <- max(rc)
    cov_cols[n+i] <- min(rc)
  }

  # sort cov_rows and cov_cols (col from small to big, then cov_rows)
  sort_idx <- order(cov_cols, cov_rows)
  cov_rows <- cov_rows[sort_idx]
  cov_cols <- cov_cols[sort_idx]

  list(
    cov_rows = cov_rows - 1,
    cov_cols = cov_cols - 1,
    mark_cov = 1:n %in% c(cov_rows, cov_cols)
  )
}