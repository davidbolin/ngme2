#' Additive non-gaussian model fitting
#'
#'  \code{ngme} function performs a analysis of non-gaussian additive models.
#'
#' @param formula formula
#' @param data    a dataframe or a list providing data
#'   (Only response variable can contain NA value,
#'    NA value in other columns will cause problem)
#' @param control control variables, see ?ngme.control
#' @param last_fit  can be ngme object from last fitting
#' @param noise measurement noise specification
#' @param beta starting value for fixed effects; a vector for matern2D model
#' @param seed  set the seed for pesudo random number generator
#'
#' @return a list of outputs contains estimation of operator paramters, noise parameters
#' @export
#'
#' @examples
#' ngme(formula = (y1 | y2 ~ x1 + x2 + f(x1, model="SPDE", var="nig") + f(W) | x3 + f(X|I, model="ar1", var="nig")),
#'
ngme <- function(
  formula,
  data,
  control       = ngme.control(),
  noise         = ngme.noise.normal(),
  last_fit      = NULL,
  beta          = NULL,
  seed          = NULL,
  debug         = FALSE
) {
  if (is.null(seed)) seed <- Sys.time()
  # -------------  CHECK INPUT ---------------
  if (is.null(formula)) {
    stop("Formula empty. See ?ngme\n")
  }

  if (is.null(data)) {
    stop("Missing data.frame/list `data'. Leaving `data' empty might lead to\n\t\tuncontrolled behaviour, therefore it is required.")
  }

  if (!is.data.frame(data) && !is.list(data)) {
    stop("\n\tArgument `data' must be a data.frame or a list.")
  }

  stopifnot(class(noise) == "ngme_noise")
  family <- noise$noise_type

  # 2. parse the formula
  time.start <- Sys.time()

  fm <- Formula::Formula(formula)
# formula = (y1 | y2 ~ x1 + f(x1, model="SPDE2D", var="nig") | x2 + f(X2, model="SPDE2D", var="nig"))
  if (all(length(fm)==c(2,2))) { ######################### bivariate model
  # strucutre-wise the bivaraite model: list of lists
    lfm = formula(fm, lhs = 1, rhs = 1)
    rfm = formula(fm, lhs = 2, rhs = 2)
  
    ########## extract data for each field
    ngme_response1 <- eval(terms(lfm)[[2]], as.data.frame(data[[1]]))
    ngme_response2 <- eval(terms(rfm)[[2]], as.data.frame(data[[2]]))
    data$ngme_response <- c(ngme_response1, ngme_response2) ##use list instead?
    data$ngme_response1 <- ngme_response1
    data$ngme_response2 <- ngme_response2
    # # 1. extract f and eval  2. get the formula without f function
    #for 1st field
    res1 <- ngme.parse.formula(lfm, as.data.frame(data[[1]]))
    #for 2nd field
    res2 <- ngme.parse.formula(rfm, as.data.frame(data[[2]]))

    latents_in <- list(res1$latents_in, res2$latents_in) #latent model = SPDE2D
    plain_fm <- list(res1$plain.fm, res2$plain.fm)

    # check if there is NA, and split data
    split_data <- parse_formula_NA(plain_fm[[1]], as.data.frame(data[[1]]))
      Y_data <- split_data$Y_data
      X_data <- split_data$X_data
      n_Y_data <- split_data$length
    split_data2 <- parse_formula_NA(plain_fm[[2]], as.data.frame(data[[2]]))
      Y_data2 <- split_data2$Y_data
      X_data2 <- split_data2$X_data
      n_Y_data2 <- split_data2$length

#TODO check the correct sizing needed for bivariate model
    ############### W_sizes is the dim of the block matrix
    W_sizes     = sum(unlist(lapply(latents_in[[1]], function(x) x["W_size"])))   #W_sizes = sum(ncol_K)
    V_sizes     = sum(unlist(lapply(latents_in[[1]], function(x) x["V_size"])))   #W_sizes = sum(nrow_K)
    n_la_params = sum(unlist(lapply(latents_in[[1]], function(x) x["n_params"])))
    model.types = unlist(lapply(latents_in[[1]], function(x) x["model_type"]))
    var.types   = unlist(lapply(latents_in[[1]], function(x) x["var.type"]))
    # print(W_sizes)
    # print(V_sizes)
#TODO check the correct sizing
    n_feff <- ncol(X_data);
    if (family == "normal") {
      n_merr <- noise$n_theta_sigma
    } else if (family == "nig") {
      n_merr <- noise$n_theta_mu + noise$n_theta_sigma + noise$n_theta_V
    }
 # 3. prepare Rcpp_list for estimate
    lm.model <- lm.fit(X_data, Y_data)
    lm.model2 <- lm.fit(X_data2, Y_data2)
    if (is.null(beta[1])) beta[1] <- lm.model$coeff
    if (is.null(beta[2])) beta[2] <- lm.model2$coeff
#TODO check the correct sizing
    n_params <- n_la_params + n_feff + n_merr

    noise <- update.ngme.noise(noise, n = n_Y_data)
#TODO change the logical statement to check the both components of vector theta_sigma
    if (family == "normal" && is.null(noise$theta_sigma == 0)) #TODO what exactly it checks - being 0 or being NULL
      noise$theta_sigma <- c(sd(lm.model$residuals), sd(lm.model2$residuals))

    ngme_block <- ngme.block_model(
      Y                 = list(Y_data, Y_data2),
      X                 = list(X_data, X_data2),
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

    # a list of B.theta.mu and B.theta.sigma and thetas...
  }
  else if (all(length(fm)==c(1,1))) {  ########################## univariate case
    fm <- formula(fm)

    ngme_response <- eval(terms(fm)[[2]], data)
    data$ngme_response <- ngme_response # watch out! injection, for f to use

    # 1. extract f and eval  2. get the formula without f function
    res <- ngme.parse.formula(fm, data)
    latents_in <- res$latents_in
    plain_fm <- res$plain.fm

    # check if there is NA, and split data
    split_data <- parse_formula_NA(plain_fm, data)
      Y_data <- split_data$Y_data
      X_data <- split_data$X_data
      n_Y_data <- split_data$length

    ############### W_sizes is the dim of the block matrix
    W_sizes     = sum(unlist(lapply(latents_in, function(x) x["W_size"])))   #W_sizes = sum(ncol_K)
    V_sizes     = sum(unlist(lapply(latents_in, function(x) x["V_size"])))   #W_sizes = sum(nrow_K)
    n_la_params = sum(unlist(lapply(latents_in, function(x) x["n_params"])))
    model.types = unlist(lapply(latents_in, function(x) x["model_type"]))
    var.types   = unlist(lapply(latents_in, function(x) x["var.type"]))

    n_feff <- ncol(X_data);
    if (family == "normal") {
      n_merr <- noise$n_theta_sigma
    } else if (family == "nig") {
      n_merr <- noise$n_theta_mu + noise$n_theta_sigma + noise$n_theta_V
    }

    # 3. prepare Rcpp_list for estimate
    lm.model <- lm.fit(X_data, Y_data)
    if (is.null(beta)) beta <- lm.model$coeff
    n_params <- n_la_params + n_feff + n_merr

    noise <- update.ngme.noise(noise, n = n_Y_data)

    if (family == "normal" && is.null(noise$theta_sigma == 0))
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
    # stopifnot("last_fit should be an ngme object"
    #   = class(last_fit) == "ngme" || is.null(last_fit))

    # if (inherits(last_fit, "ngme")) {
    #   output <- last_fit$est_output
    #   # use last fit estimates

    #   # block noise
    #   noise$theta_mu    <- output$noise[["theta_mu"]]
    #   noise$theta_sigma <- output$noise[["theta_sigma"]]
    #   noise$theta_noise <- output$noise[["theta_V"]]
    #   noise$V           <- output$noise[["V"]]

    #   # latents
    #   for (i in seq_along(latents_in)) {
    #     last_fit_latent <- output$latent[[i]]
    #     latents_in[[i]]$operator$theta_K  <- last_fit_latent[["theta_K"]]
    #     latents_in[[i]]$W                 <- last_fit_latent[["W"]]

    #     latents_in[[i]]$noise$theta_mu    <- last_fit_latent[["theta_mu"]]
    #     latents_in[[i]]$noise$theta_sigma <- last_fit_latent[["theta_sigma"]]
    #     latents_in[[i]]$noise$theta_noise <- last_fit_latent[["theta_V"]]
    #     latents_in[[i]]$noise$V           <- last_fit_latent[["V"]]
    #   }
    # }

  } else {
    stop("unknown structure of formula")
  }

  ################# Run CPP ####################
  if (!control$estimation) {
    print("Start estimation by setting estimation = TRUE")
    return(ngme_block)
  }

if (debug) print(str(ngme_block))
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
    X_pred <- split_data$X_NA;
    Xb_pred <- drop(X_pred %*% ngme_block$beta)
    Xb_data <- drop(X_data %*% ngme_block$beta)

    # ngme_response[split_data$index_NA] <- if (length(Xb_pred) == 0) AW_pred else AW_pred + Xb_pred
    #
    linear_predictor[split_data$index_NA]   <- if (length(Xb_pred) == 0) AW_pred else AW_pred + Xb_pred
    linear_predictor[split_data$index_data] <- if (length(Xb_data) == 0) AW_data else AW_data + Xb_data

    attr(ngme_block, "prediction") <- list(
      linear_predictor  = linear_predictor,
      index_pred        = split_data$index_NA
    )
  }

  # cat(paste("total time is", Sys.time() - time.start, " \n"))
  ngme_block
}

clean_outputs <- function() {
  trajs <- list()
  # deal with multiple chain data
  for (i in seq_along(outputs)) { # what if 2d?
    # flat the lists
    beta_traj        <- unlist(attr(outputs[[i]], "trajectory")$beta)
    theta_mu_traj    <- unlist(attr(outputs[[i]], "trajectory")$theta_mu_traj)
    theta_sigma_traj <- unlist(attr(outputs[[i]], "trajectory")$theta_sigma_traj)
    theta_V_traj     <- attr(outputs[[i]], "trajectory")$theta_V_traj
    attr(outputs[[i]], "trajectory") <- NULL

    trajs_lat <- list()
    for (j in seq_along(outputs[[i]]$latents)) {
      la_K_traj     <- unlist(attr(outputs[[i]]$latents[[j]], "trajectory")$theta_K)
      la_mu_traj    <- unlist(attr(outputs[[i]]$latents[[j]], "trajectory")$theta_mu_traj)
      la_sigma_traj <- unlist(attr(outputs[[i]]$latents[[j]], "trajectory")$theta_sigma_traj)
      la_V_traj     <- attr(outputs[[i]]$latents[[j]], "trajectory")$theta_V_traj
      trajs_lat[[j]] <- list(
        theta_K       = la_K_traj,
        theta_mu      = la_mu_traj,
        theta_sigma   = la_sigma_traj,
        theta_V       = la_V_traj
      )
      attr(outputs[[i]]$latents[[j]], "trajectory") <- NULL
    }

    trajs[[i]] <- list(
      beta = beta_traj,
      theta_mu = theta_mu_traj,
      theta_sigma = theta_sigma_traj,
      theta_V = theta_V_traj,
      trajs_lat = trajs_lat
    )
  }
}

# the general block model
ngme.block_model <- function(
  Y           = NULL,
  X           = NULL,
  beta        = NULL,
  latents     = list(),
  noise       = list(),
  control     = list(),
  ...
) {
  structure(
    list(
      Y                 = Y,
      X                 = X,
      beta              = beta,
      latents           = latents,
      noise             = noise,
      control           = control,
      ...
    ),
    class = "ngme"
  )
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
  cat(paste("   ",  cat(ngme.format(ngme$beta)))); cat("\n\n")

  cat("Measurement noise: \n");
  print.ngme_noise(ngme$noise, padding = 2); cat("\n\n")

  cat("Latent models: \n");
  for (i in seq_along(ngme$latents)) {
    cat("[["); cat(i); cat("]]\n")
    print.ngme_model(ngme$latents[[i]], padding = 2)
  }
}

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

