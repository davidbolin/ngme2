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
#' @param noise measurement noise specification
#' @param beta starting value for fixed effects
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
  noise         = ngme.noise(),
  last_fit      = NULL,
  beta          = NULL,
  seed          = NULL
) {
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
  family <- noise$noise_type

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

    ############### n_meshs is the dim of the block matrix
    n_meshs     = sum(unlist(lapply(latents_in, function(x) x["n_mesh"])))
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
      n_meshs           = n_meshs,
      n_la_params       = n_la_params,
      n_params          = n_params, # how many param to opt. in total
      latents           = latents_in,
      noise             = noise,
      seed              = seed,
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

  cat("Starting estimation... \n")
  out <- estimate_cpp(ngme_block)
  cat("Estimation done! \n")

  # update the ngme_block using estimation.
  estimation <- out$estimation
    # 1. update fixed effects
    ngme_block$beta <- estimation$beta
    # 2. updates noise
    ngme_block$noise <- update_noise_with_est(ngme_block$noise, estimation$noise)
    # 3. update latents
    ngme_block$latents <- update_latents_with_est(ngme_block$latents, estimation$latents)

  attr(ngme_block, "opt_trajectory") <- out$opt_trajectory

  ################# doing prediction ####################
  if (split_data$contain_NA) {
    AW <- 0
    for (i in seq_along(latents_in)) {
      A_pred <- latents_in[[i]]$A_pred
      W <- out$est_output$latents[[i]]$W
      AW <- AW + drop(A_pred %*% W)
    }

    # fixed effects. watch out! Xb could be double(0)
    X_pred <- split_data$X_NA
    Xb <- drop(X_pred %*% out$est_output$fixed_effects)
    ngme_response[split_data$index_NA] <- if (length(Xb) == 0) AW else AW + Xb

    out$prediction <- list(
      linear_predictor  = ngme_response,
      index_pred        = split_data$index_NA
    )
  }

  # cat(paste("total time is", Sys.time() - time.start, " \n"))
  ngme_block
}


# helper function
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

update_latents_with_est <- function(latents, latents_out) {
  for (i in seq_along(latents_out)) {
    latents[[i]]$theta_K  <- latents_out[[i]]$theta_K
    latents[[i]]$W        <- latents_out[[i]]$W

    latents[[i]]$noise    <- update_noise_with_est(latents[[i]]$noise, latents_out[[i]])
  }

  latents
}


# the general block model
ngme.block_model <- function(
  Y       = NULL,
  X       = NULL,
  beta    = NULL,
  latents = list(),
  noise   = list(),
  control = list(),
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
  cat(paste("   ",  cat(format(ngme$beta)))); cat("\n")

  cat("Measurement noise: \n");
  print.ngme_noise(ngme$noise, padding = 2); cat("\n\n")

  cat("Latent models: \n");
  for (i in seq_along(ngme$latents)) {
    cat("[["); cat(i); cat("]]\n")
    print.ngme_model(ngme$latents[[i]], padding = 2)
  }
}