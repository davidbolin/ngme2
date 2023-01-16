#' Fit an additice linear mixed effect model
#'
#'  \code{ngme} function performs an analysis of non-gaussian additive models.
#'  It does the maximum likelihood estimation via stochastic gradient descent.
#'  The prediction of unknown location can be performed by leaving the response
#'  variable to be \code{NA}. The likelihood is specified by \code{family}.
#' The model estimation control can be setted in \code{control} using
#'  \code{ngme_control()} function, see \code{?ngme_control} for details.
#' See \code{ngme_model_types()} for available models.
#' @param formula formula
#' @param data    a dataframe or a list providing data
#'   (Only response variable can contain \code{NA} value,
#'    \code{NA} value in other columns will cause problem)
#' @param control control variables, see \code{?ngme.control}
#' @param family likelihood type, same as measurement noise specification, 1. string 2. ngme noise obejct
#' @param beta starting value for fixed effects
#' @param start  starting ngme object (usually object from last fitting)
#' @param seed  set the seed for pesudo random number generator
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
#'  ) + f(
#'    model = model_rw(1:5, order=1, circular = TRUE),
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
    index_NA <- is.na(ngme_response)
    data$index_NA <- index_NA # watch out! injection, for f to see

    # 1. extract f and eval  2. get the formula without f function
    res <-ngme_parse_formula(fm, data)
    latents_in <- res$latents_in
    plain_fm <- res$plain_fm
    # names(latents_in) <- sapply(latents_in, function(x) {x$name})

    # get Y and X
    Y_data    <- ngme_response[!index_NA]
    n_Y_data  <- length(Y_data)
    X_full    <- model.matrix(delete.response(terms(plain_fm)), as.data.frame(data))
    # if (length(X_full) == 0) X_full <- model.matrix(terms(plain_fm), data) # Y ~ 1 case
    X_data    <- X_full[!index_NA, , drop = FALSE]

    ############### W_sizes is the dim of the block matrix
    W_sizes     = sum(unlist(lapply(latents_in, function(x) x["W_size"])))   #W_sizes = sum(ncol_K)
    V_sizes     = sum(unlist(lapply(latents_in, function(x) x["V_size"])))   #W_sizes = sum(nrow_K)
    n_la_params = sum(unlist(lapply(latents_in, function(x) x["n_params"])))

    n_feff <- ncol(X_data);
    if (family_type == "normal") {
      n_merr <- noise$n_theta_sigma
    } else if (family_type == "nig") {
      n_merr <- noise$n_theta_mu + noise$n_theta_sigma + noise$n_nu
    }

    # 3. prepare Rcpp_list for estimate
    lm.model <- stats::lm.fit(X_data, Y_data)
    if (is.null(beta)) beta <- lm.model$coeff
    n_params <- n_la_params + n_feff + n_merr

    noise <- update_noise(noise, n = n_Y_data)

    if (family_type == "normal" && is.null(noise$theta_sigma == 0))
      noise$theta_sigma <- sd(lm.model$residuals)

    ngme_block <- ngme_block(
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
        noise <- update_noise(noise, new_noise = start$noise)
        for (i in seq_along(start$latents)) {
          latents[[i]]$theta_K  <- start$latents[[i]]$theta_K
          latents[[i]]$W        <- start$latents[[i]]$W
          latents[[i]]$noise    <- update_noise(
            latents[[i]]$noise, new_noise = start$latents[[i]]$noise
          )
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

  ################# Update the estimates ####################
    est_output <- mean_list(outputs)
    ngme_block$beta <- est_output$beta
    ngme_block$noise <- update_noise(ngme_block$noise, new_noise = est_output$noise)
    for (i in seq_along(ngme_block$latents)) {
      ngme_block$latents[[i]]$theta_K  <- est_output$latents[[i]]$theta_K
      ngme_block$latents[[i]]$W        <- est_output$latents[[i]]$W
      ngme_block$latents[[i]]$noise    <- update_noise(
        ngme_block$latents[[i]]$noise, new_noise = est_output$latents[[i]]
      )
    }
    # 2. get trajs
    attr(ngme_block, "trajectory") <- get_trajs(outputs)
  }

  ################# Prediction ####################
  if (any(data$index_NA)) {
    # posterior sampling
    # ngme_block <- sampling_cpp(ngme_block, 100, TRUE)

    # form a linear predictor
    lp <- double(length(ngme_response))

    AW_pred <- 0; AW_data <- 0
    for (i in seq_along(latents_in)) {
      W <- ngme_block$latents[[i]]$W
      lp[index_NA]  <- lp[index_NA] + drop(latents_in[[i]]$A_pred %*% W)
      lp[!index_NA] <- lp[!index_NA] + drop(latents_in[[i]]$A %*% W)
    }

    # fixed effects. watch out! beta could be double(0)
    fe <- if (length(ngme_block$beta) == 0) 0 else drop(X_full %*% ngme_block$beta)
    lp <- lp + fe

    attr(ngme_block, "prediction") <- list(
      fe        = fe,
      lp        = lp,
      index_NA  = index_NA
    )
  }

  # cat(paste("total time is", Sys.time() - time.start, " \n"))
  ngme_block
}

# create the general block model
ngme_block <- function(
  Y           = NULL,
  X           = NULL,
  beta        = NULL,
  noise       = noise_normal(),
  latents     = list(),
  control     = list(),
  debug       = FALSE,
  ...
) {

  latents_string <- rep(" ", 14) # padding of 14 spaces
  for (latent in latents)
    latents_string <- c(latents_string, latent$par_string)
  beta_str  <- if (length(beta) > 0) paste0("  beta_", seq_along(beta)) else ""
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
      beta              = beta,
      latents           = latents,
      noise             = noise,
      control           = control,
      n_merr            = noise$n_params,
      debug             = debug,
      par_string        = par_string,
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
