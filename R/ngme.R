#' Fit an additice linear mixed effect model
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
#' @param family likelihood type, same as measurement noise specification, 1. string 2. ngme noise obejct
#' @param start  starting ngme object (usually object from last fitting)
#' @param debug  toggle debug mode
#'
#' @return a list of outputs contains estimation of operator paramters, noise parameters
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
  start         = NULL,
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

  # parse the formula get a list of ngme_block
  parse_ngme <- ngme_parse_formula(formula, data, control_ngme, noise)
  list_ngmes <- parse_ngme$blocks_rep

  attr(list_ngmes, "fitting") <- fitting

  ####### Use Last_fit ngme object to update Rcpp_list
  stopifnot("start should be an ngme object"
    = inherits(start, "ngme") || is.null(start))

  # update with start (list of ngmes)
  if (inherits(start, "ngme")) {
    for (i in seq_along(list_ngmes)) {
      list_ngmes[[i]] <- within(list_ngmes[[i]], {
        beta <- start[[i]]$beta
        noise <- update_noise(noise, new_noise = start[[i]]$noise)
        for (i in seq_along(start[[i]]$latents)) {
          latents[[i]]$theta_K  <- start[[i]]$latents[[i]]$theta_K
          latents[[i]]$W        <- start[[i]]$latents[[i]]$W
          latents[[i]]$noise    <- update_noise(
            latents[[i]]$noise, new_noise = start[[i]]$latents[[i]]$noise
          )
        }
      })
    }
  }

if (debug) {print(str(list_ngmes[[1]]))}

# check all f has the same replicate
# otherwise change replicate to group="iid"
  ################# Run CPP ####################
  check_dim(list_ngmes)
  if (control_opt$estimation) {
    cat("Starting estimation... \n")
    outputs <- estimate_cpp(list_ngmes, control_opt)
    cat("\n")

  ################# Update the estimates ####################
    est_output <- mean_list(outputs)
    for (i in seq_along(list_ngmes))
      list_ngmes[[i]] <- update_ngme_est(list_ngmes[[i]], est_output[[i]])

    # return the mean of samples of W of posterior
    cat("Starting posterior sampling... \nNote: Use ngme$latents[[model_name]]$W  to access the posterior mean of process \n")
    for (i in seq_along(list_ngmes)) {
      ngme_block <- list_ngmes[[i]]
      ngme_block$control_ngme$init_sample_W <- FALSE
      mean_post_W <- mean_list(
        sampling_cpp(ngme_block, control_ngme$post_samples_size, TRUE, control_opt$seed)[["W"]]
      )

      idx <- 1
      for (j in seq_along(ngme_block$latents)) {
        ngme_block$latents[[j]]$W <- mean_post_W[idx : (ngme_block$latents[[j]]$W_size + idx - 1)]
        idx <- idx + ngme_block$latents[[j]]$W_size
      }
      list_ngmes[[i]] <- ngme_block
    }
    cat("Posterior sampling done! \n")

    # transform trajectory
    traj_df_chains <- transform_traj(attr(outputs, "opt_traj"))
    # dispatch trajs to each latent and block
      idx <- 0;
      for (i in seq_along(list_ngmes[[1]]$latents)) {
        n_params <- list_ngmes[[1]]$latents[[i]]$n_params
        lat_traj_chains = list()
        for (j in seq_along(traj_df_chains))
          lat_traj_chains[[j]] <- traj_df_chains[[j]][idx + 1:n_params, ]

        attr(list_ngmes[[1]]$latents[[i]], "lat_traj") <- lat_traj_chains
        idx <- idx + n_params
      }
      # mn and beta
      block_traj <- list()
      for (j in seq_along(traj_df_chains))
        block_traj[[j]] <- traj_df_chains[[j]][(idx + 1):list_ngmes[[1]]$n_params, ]
      attr(list_ngmes[[1]], "block_traj") <- block_traj
      attr(outputs, "opt_traj") <- NULL
  }

  # cat(paste("total time is", Sys.time() - time.start, " \n"))
  attr(list_ngmes, "control_opt") <- control_opt

  class(list_ngmes) <- "ngme_fit"
  list_ngmes
}

# create the general block model
ngme_block <- function(
  Y            = NULL,
  X            = NULL,
  noise        = noise_normal(),
  latents      = list(),
  control_ngme = list(),
  ...
) {
  # compute W_sizes and V_sizes
  W_sizes     = sum(unlist(lapply(latents, function(x) x[["n_rep"]] * x[["W_size"]])))   #W_sizes = sum(ncol_K)
  V_sizes     = sum(unlist(lapply(latents, function(x) x[["n_rep"]] * x[["V_size"]])))   #W_sizes = sum(nrow_K)
  n_la_params = sum(unlist(lapply(latents, function(x) x["n_params"])))
  n_feff <- ncol(X);
  n_params <- n_la_params + n_feff + noise$n_params

  latents_string <- rep(" ", 14) # padding of 14 spaces
  for (latent in latents)
    latents_string <- c(latents_string, latent$par_string)
  beta_str  <- if (ncol(X) > 0) paste0("  beta_", seq_along(beta)) else ""
  m_mu_str    <- paste0("    mu_", seq_along(noise$theta_mu))
  m_sigma_str <- paste0(" sigma_", seq_along(noise$theta_sigma))
  m_nu_str    <- "    nu_1"
  merr_str <- switch(noise$noise_type,
    normal  = m_sigma_str,
    nig     = c(m_mu_str, m_sigma_str, m_nu_str)
  )
  par_string <- do.call(paste0, as.list(c(latents_string, beta_str, merr_str)))

  structure(
    list(
      Y                 = Y,
      X                 = X,
      beta              = control_ngme$beta,
      latents           = latents,
      noise             = noise,
      control_ngme      = control_ngme,
      par_string        = par_string,
      W_sizes           = W_sizes,
      V_sizes           = V_sizes,
      n_merr            = noise$n_params,
      n_params          = n_params,
      n_la_params       = n_la_params,
      ...
    ),
    class = c("ngme", "list")
  )
}


#' Print ngme object
#'
#' @param x ngme object
#' @param ... ignored
#'
#' @return a list (noise specifications)
#' @export
print.ngme <- function(x, ...) {
  ngme <- x
  cat("*** Ngme object ***\n\n");

  cat("Fixed effects: \n");
  cat(paste("  ", ngme_format("beta", ngme$beta)));
  cat("\n\n")

  cat("Measurement noise: \n");
  print.ngme_noise(ngme$noise, padding = 2); cat("\n\n")

  cat("Latent models: \n");
  for (i in seq_along(ngme$latents)) {
    # cat("[["); cat(i); cat("]]")
    # cat("\""); cat(names(ngme$latents)[[i]]); cat("\"\n")
    cat("$"); cat(names(ngme$latents)[[i]]); cat("\n")
    print.ngme_model(ngme$latents[[i]], padding = 2)
  }
}

# helper function
# get trajs from a list of estimates
get_trajs <- function(outputs) {
  ret <- list()
  for (i in seq_along(outputs)) {
    ret[[i]] <- list()
    ret[[i]]$block_traj <- attr(outputs[[i]], "trajectory")
    for (j in seq_along(outputs[[i]]$latents)) {
      ret[[i]]$latents[[j]] <- list()
      ret[[i]]$latents[[j]] <- attr(outputs[[i]]$latents[[j]], "trajectory")
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
  ngme_block, est_output
) {
  ngme_block$beta <- est_output$beta
  ngme_block$noise <- update_noise(ngme_block$noise, new_noise = est_output$noise)
  for (i in seq_along(ngme_block$latents)) {
    ngme_block$latents[[i]]$theta_K  <- est_output$latents[[i]]$theta_K
    ngme_block$latents[[i]]$W        <- est_output$latents[[i]]$W
    ngme_block$latents[[i]]$noise    <- update_noise(
      ngme_block$latents[[i]]$noise, new_noise = est_output$latents[[i]]
    )
    if (ngme_block$latents[[i]]$model == "tp") {
      n1 <- ngme_block$latents[[i]]$left$n_theta_K
      n2 <- ngme_block$latents[[i]]$right$n_theta_K
      ngme_block$latents[[i]]$left$theta_K <- ngme_block$latents[[i]]$theta_K[1:n1]
      ngme_block$latents[[i]]$right$theta_K <- ngme_block$latents[[i]]$theta_K[(n1+1):(n1+n2)]
    }
  }
  ngme_block
}

#' @export
print.ngme_fit <- function(x, ...) {
  print(x[[1]])
  cat("\n");
  cat("Number of replicates is ", length(x), "\n");
}


######
check_dim <- function(list_ngmes) {
  for (ngme in list_ngmes) {
    if (ncol(ngme$X) != length(ngme$beta)) {
      stop("The number of columns of X is not equal to the length of beta")
    }
    for (latent in ngme$latents) {
        if (latent$model != "tp" && latent$W_size != ncol(latent$C)) {
         stop("The W_size of the latent model is not equal to the number of columns of K")
        }
        if (latent$model != "tp" && latent$V_size != nrow(latent$C)) {
          stop("The V_size of the latent model is not equal to the number of rows of K")
        }
        if (latent$V_size != latent$noise$n_noise) {
          stop("The V_size of the latent model is not equal to the length of noise")
        }
        # ncol(A) = W_size
        if (ncol(latent$A) != latent$W_size * latent$n_rep) {
          stop("The number of columns of A is not equal to the W_size of the latent model")
        }
        if (!all(latent$h == latent$noise$h) || length(latent$h) != latent$V_size) {
          stop("The h of the latent model is not equal to the h of the noise")
        }
    }
  }
}

# check dim. before run cpp
# for (latent in ngme_block$latents) {
#   if (latent$model != "tensor_prod") {
#     stopifnot("nrow(K) should be equal to length of noise, please check idx, replicate argument" =
#       nrow(latent$K) == latent$V_size)
#     stopifnot("ncol(K) should be equal to length of mesh, please check idx, replicate argument" =
#       ncol(latent$K) == latent$W_size)
#     # check given V > 0
#   }
# }

#' Parse the formula for ngme function
#'
#' @param fm Formula
#' @param data data.frame
#' @param control_ngme control_ngme
#' @param noise noise
#'
#' @return a list (replicate) of ngme_block models
ngme_parse_formula <- function(
  fm,
  data,
  control_ngme,
  noise
) {
  Fm <- Formula::Formula(fm)
  enclos_env <- list2env(as.list(parent.frame()), envir = parent.frame(2))
  if ((length(Fm)[[1]] > 1)) {
    # LHS more than 1 variable
    # now assume bivariate model
    y1 <- eval(formula(Fm, lhs=1)[[2]], envir=data, enclos=enclos_env)
    y2 <- eval(formula(Fm, lhs=2)[[2]], envir=data, enclos=enclos_env)

    stop("Bivariate model is not supported yet")
  } else {
    # adding special mark
    tf <- terms.formula(fm, specials = c("f"))
    terms <- attr(tf, "term.labels")
    intercept <- attr(tf, "intercept")

    # order of f terms in labels
    spec_order <- attr(tf, "specials")$f - 1

    # construct plain formula without f
    # watch out! terms[-double(0)] -> character(0)
    fixf <- if (length(spec_order) == 0) terms else terms[-spec_order]
    response <- as.character(attr(tf, "variables")[[2]])
    plain_fm_str <- paste(response, "~", intercept, paste(c("", fixf), collapse = " + "))
    plain_fm <- formula(plain_fm_str)

    # eval the data
    ngme_response <- eval(stats::terms(fm)[[2]], envir = data, enclos = enclos_env)
    stopifnot("Have NA in your response variable" = all(!is.na(ngme_response)))
    X_full    <- model.matrix(delete.response(terms(plain_fm)), as.data.frame(data))

    ########## parse latents terms
    pre_model <- list()
    for (i in spec_order) {
      if (!grepl("data *=", terms[i])) {
        # adding data=data if not specified
        str <- gsub("^f\\(", "ngme2::f(data=data,", terms[i])
      } else if (grepl("data *= *NULL", terms[i])) {
        # change data=NULL to data=data
        str <- gsub("^f\\(", "ngme2::f(", terms[i])
        str <- gsub("data *= *NULL", "data=data", str)
      } else {
        # keep data=sth.
        str <- gsub("^f\\(", "ngme2::f(", terms[i])
      }

      # add information of index_NA
      # adding 1 term for furthur use in f
      # str <- gsub("ngme2::f\\(", "ngme2::f(index_NA=index_NA,", str)
      res <- eval(parse(text = str), envir = data, enclos = enclos_env)

      # give default name
      if (res$name == "field") res$name <- paste0("field", length(pre_model) + 1)
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
      latents_rep <- list()
      for (latent in pre_model) {
        latent$map <- sub_locs(latent$map, idx)
        latent$replicate <- latent$replicate[idx]
        latent$eval = TRUE
        latents_rep[[latent$name]] <- eval(latent, envir = data, enclos = enclos_env)
      }
      # give initial value
      lm.model <- stats::lm.fit(X, Y)
      if (is.null(control_ngme$beta)) control_ngme$beta <- lm.model$coeff
      if (is.null(noise$theta_sigma)) noise$theta_sigma <- log(sd(lm.model$residuals))

      noise_new <- update_noise(noise, n = length(Y))

      blocks_rep[[i]] <- ngme_block(
        Y = Y,
        X = X,
        noise = noise_new,
        latents = latents_rep,
        replicate = uni_repl[[i]],
        control_ngme = control_ngme
      )
    }
  }

  list(
    replicate   = repls,
    blocks_rep  = blocks_rep
  )
}
