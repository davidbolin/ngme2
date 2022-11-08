#' Additive non-guassian model fitting
#'
#'  \code{ngme} function performs a analysis of non-gaussian additive models.
#'
#' @param formula formula
#' @param data    a dataframe or a list providing data
#'   (Only response variable can contain NA value,
#'    NA value in other columns will cause problem)
#' @param control control variables, see ?ngme.control
#' @param last_fit  can be ngme object from last fitting
#' @param family likelihood type, same as measurement noise specification, 1. string 2. ngme noise obejct
#' @param beta starting value for fixed effects
#' @param start  starting ngme object (usually object from last fitting)
#' @param seed  set the seed for pesudo random number generator
#'
#' @return a list of outputs contains estimation of operator paramters, noise parameters
#' @export
#'
#' @examples
#' ngme(
#'  formula = Y ~ x1 + f(
#'    index = x2,
#'    model = "ar1",
#'    noise = noise_nig(),
#'    theta_K = 0.5
#'  ) + f(
#'    model = model_rw1(1:5, circular = TRUE),
#'    noise = noise_normal(),
#'  ),
#'  family = noise_normal(sd = 0.5),
#'  data = data.frame(Y = 1:5, x1 = 2:6, x2 = 3:7),
#'  control = ngme_control(
#'    estimation = FALSE
#'  )
#')
ngme <- function(
  formula,
  data,
  control       = ngme_control(),
  family        = "normal",
  last_fit      = NULL,
  beta          = NULL,
  seed          = NULL,
  start         = NULL,
  debug         = FALSE
) {
  if (is.character(family))
    noise <- switch(family,
      "normal" = noise_normal(),
      "nig"    = noise_nig(),
      stop("Unknown family!")
    )
  else
    noise <- family # ngme noise object

  if (is.null(seed)) seed <- Sys.time()
  # -------------  CHECK INPUT ---------------
  if (is.null(formula)) {
    stop("Formula empty. See ?ngme\n")
  }

  if (is.null(data)) {
    stop("Missing data.frame/list `data'. Leaving `data' empty might lead to\n\t\tuncontrolled behaviour, therefore is it required.")
  }

  if (!is.data.frame(data) && !is.list(data)) {
    stop("\n\tArgument `data' must be a data.frame or a list.")
  }

  stopifnot(class(noise) == "ngme_noise")
  family_type <- noise$noise_type

  # 2. parse the formula
  time.start <- Sys.time()

  fm <- Formula::Formula(formula)
# Y1|Y2|Y3 ~ ..|..|..
  if (all(length(fm)==c(2,2))) { ######################### bivariate model
    lfm = formula(fm, lhs=1, rhs=1)
    rfm = formula(fm, lhs=2, rhs=2)
    ########## to-do

    # a list of B.theta.mu and B.theta.sigma and thetas...
  }
  else if (all(length(fm)==c(1,1))) {  ########################## univariate case
    fm <- formula(fm)

    # eval the response variable in data environment
    ngme_response <- eval(stats::terms(fm)[[2]], envir = data, enclos = parent.frame())
    data$ngme_response <- ngme_response # watch out! injection, for f to see

    # 1. extract f and eval  2. get the formula without f function
    res <- ngme_parse_formula(fm, data)
    latents_in <- res$latents_in
    plain_fm <- res$plain_fm

    # check if there is NA, and split data
    split_data <- parse_formula_NA(plain_fm, data)
      Y_data <- split_data$Y_data
      X_data <- split_data$X_data
      n_Y_data <- split_data$length
    ############### W_sizes is the dim of the block matrix
    W_sizes     = sum(unlist(lapply(latents_in, function(x) x["W_size"])))   #W_sizes = sum(ncol_K)
    V_sizes     = sum(unlist(lapply(latents_in, function(x) x["V_size"])))   #W_sizes = sum(nrow_K)
    n_la_params = sum(unlist(lapply(latents_in, function(x) x["n_params"])))

    n_feff <- ncol(X_data);
    if (family_type == "normal") {
      n_merr <- noise$n_theta_sigma
    } else if (family_type == "nig") {
      n_merr <- noise$n_theta_mu + noise$n_theta_sigma + noise$n_theta_V
    }

    # 3. prepare Rcpp_list for estimate
    lm.model <- stats::lm.fit(X_data, Y_data)
    if (is.null(beta)) beta <- lm.model$coeff
    n_params <- n_la_params + n_feff + n_merr

    noise <- update_noise(noise, n = n_Y_data)

    if (family_type == "normal" && is.null(noise$theta_sigma == 0))
      noise$theta_sigma <- sd(lm.model$residuals)

    ngme_block <- ngme.block_model(
      Y                 = Y_data,
      X                 = X_data,
      beta              = beta,
      W_sizes           = W_sizes,
      V_sizes           = V_sizes,
      n_la_params       = n_la_params,
      n_params          = n_params, # how many param to opt. in total
      latents           = latents_in,
      noise             = noise,
      seed              = seed,
      debug             = debug,
      control           = control
    )

  ####### Use Last_fit ngme object to update Rcpp_list
    stopifnot("start should be an ngme object"
      = inherits(start, "ngme") || is.null(start))

    if (inherits(start, "ngme")) {
      ngme_block <- within(ngme_block, {
        beta <- start$beta

        noise$theta_mu    <- start$noise$theta_mu
        noise$theta_sigma <- start$noise$theta_sigma
        noise$theta_V     <- start$noise$theta_V
        noise$V           <- start$noise$V

        # latents
        for (i in seq_along(latents_in)) {
          latents_in[[i]]$theta_K           <- start$latents[[i]][["theta_K"]]
          latents_in[[i]]$W                 <- start$latents[[i]][["W"]]

          latents_in[[i]]$noise$theta_mu    <- start$latents[[i]][["theta_mu"]]
          latents_in[[i]]$noise$theta_sigma <- start$latents[[i]][["theta_sigma"]]
          latents_in[[i]]$noise$theta_V     <- start$latents[[i]][["theta_V"]]
          latents_in[[i]]$noise$V           <- start$latents[[i]][["V"]]
        }
      })
    }

  } else {
    stop("unknown structure of formula")
  }

if (debug) print(str(ngme_block))

  ################# Run CPP ####################
  if (control$estimation) {
    cat("Starting estimation... \n")
    outputs <- estimate_cpp(ngme_block)
    cat("Estimation done! \n")

    # 1. update with estimates
    ngme_block <- update_ests(ngme_block, mean_list(outputs))
    # 2. get trajs
    attr(ngme_block, "trajectory") <- get_trajs(outputs)

  ################# doing prediction ####################
    if (split_data$contain_NA) {
      # form a linear predictor
      linear_predictor <- double(length(ngme_response))

      AW_pred <- 0; AW_data <- 0
      for (i in seq_along(latents_in)) {
        W <- ngme_block$latents[[i]]$W
        AW_pred <- AW_pred + drop(latents_in[[i]]$A_pred %*% W)
        AW_data <- AW_data + drop(latents_in[[i]]$A %*% W)
      }

      # fixed effects. watch out! Xb could be double(0)
      X_pred <- split_data$X_pred;
      Xb_pred <- drop(X_pred %*% ngme_block$beta)
      Xb_data <- drop(X_data %*% ngme_block$beta)

      # ngme_response[split_data$indeX_pred] <- if (length(Xb_pred) == 0) AW_pred else AW_pred + Xb_pred
      #
      linear_predictor[split_data$indeX_pred]   <- if (length(Xb_pred) == 0) AW_pred else AW_pred + Xb_pred
      linear_predictor[split_data$index_data] <- if (length(Xb_data) == 0) AW_data else AW_data + Xb_data

      attr(ngme_block, "prediction") <- list(
        linear_predictor  = linear_predictor,
        index_pred        = split_data$indeX_pred
      )
    }
  }
  # cat(paste("total time is", Sys.time() - time.start, " \n"))
  ngme_block
}

#' Print ngme object
#'
#' @param ngme ngme object
#'
#' @return a list (noise specifications)
#' @export
print.ngme <- function(ngme) {
  cat("*** Ngme object ***\n\n");

  cat("Fixed effects: \n");
  cat(paste("  ", ngme_format("beta", ngme$beta)));
  cat("\n\n")

  cat("Measurement noise: \n");
  print.ngme_noise(ngme$noise, padding = 2); cat("\n\n")

  cat("Latent models: \n");
  for (i in seq_along(ngme$latents)) {
    cat("[["); cat(i); cat("]]\n")
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

# helper function
# update the model using the mean of chains
update_ests <- function(ngme_block, est_output) {
  # helper function - update noise with est. values
  update_noise_with_est <- function(noise, noise_out) {
    if (noise_out$noise_type == "nig") {
      noise$theta_mu    <- noise_out$theta_mu
      noise$theta_sigma <- noise_out$theta_sigma
      noise$theta_V     <- noise_out$theta_V
      noise$V           <- noise_out$V
    } else if (noise_out$noise_type == "normal") {
      noise$theta_sigma <- noise_out$theta_sigma
    }
    noise
  }

  # helper function - update latent process with est. values
  update_latents_with_est <- function(latents, latents_out) {
    for (i in seq_along(latents_out)) {
      latents[[i]]$theta_K  <- latents_out[[i]]$theta_K
      latents[[i]]$W        <- latents_out[[i]]$W
      latents[[i]]$noise    <- update_noise_with_est(latents[[i]]$noise, latents_out[[i]])
    }
    latents
  }

  ngme_block$beta <- est_output$beta
  ngme_block$noise <- update_noise_with_est(ngme_block$noise, est_output$noise)
  ngme_block$latents <- update_latents_with_est(ngme_block$latents, est_output$latents)

  ngme_block
}

