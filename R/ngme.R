#' Additive non-guassian model fitting
#'
#'  \code{ngme} function performs a analysis of non-gaussian additive models.
#'
#'
#' @param formula formula
#' @param data    a dataframe contains data
#' @param controls control variables
#' @param debug  debug option
#' @param last_fit  1. ngme object from last fitting
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
  controls      = ngme.control(),
  debug         = ngme.debug(),
  noise         = ngme.noise(),
  last_fit      = NULL,
  beta          = NULL,
  seed          = NULL
) {
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

  stopifnot(class(noise) == "noise")
  family <- noise$type

  if (is.null(seed)) seed <- Sys.time()

  # 2. parse the formula
  time.start <- Sys.time()

  fm = Formula::Formula(formula)
# Y1|Y2|Y3 ~ ..|..|..
  if (all(length(fm)==c(2,2))) { ######################### bivariate model
    lfm = formula(fm, lhs=1, rhs=1)
    rfm = formula(fm, lhs=2, rhs=2)
    ########## to-do

    # a list of B.theta.mu and B.theta.sigma and thetas...
  }
  else if (all(length(fm)==c(1,1))) {  ########################## univariate case
    fm = formula(fm)

    # 1. extract f and eval  2. get the formula without f function
    res = ngme.interpret.formula(fm, data)
    latents_in = res$latents_in
    plain.fm = res$plain.fm

    # eval part without f (fm = y ~ x1 + x2)
    Y = model.frame(plain.fm, data)[[1]]
    X = model.matrix(plain.fm, data) # design matrix

    ############### n_meshs is the dim of the block matrix
    n_meshs     = sum(unlist(lapply(latents_in, function(x) x["n_mesh"] )))
    n_la_params = sum(unlist(lapply(latents_in, function(x) x["n_la_params"] )))
    model.types = unlist(lapply(latents_in, function(x) x["model_type"] ))
    var.types   = unlist(lapply(latents_in, function(x) x["var.type"] ))

    n_feff = ncol(X);
    if (family == "normal") {
      n_merr = noise$n_theta_sigma
    } else if (family == "nig") {
      n_merr = noise$n_theta_mu + noise$n_theta_sigma + noise$n_theta_V
    }

    # 3. prepare in_list for estimate
    lm.model = lm.fit(X, Y)
    if (is.null(beta)) beta = lm.model$coeff
    n_params = n_la_params + n_feff + n_merr

    general_in <- list(
      seed             = seed,
      Y                = Y,
      X                = X,
      beta             = beta,
      n_meshs          = n_meshs,
      n_la_params      = n_la_params,
      n_params         = n_params # how many param to opt. in total
    )

    noise <- update.ngme.noise(noise, n = length(Y))

    if (family == "normal" && is.null(noise$theta_sigma == 0))
      noise$theta_sigma = sd(lm.model$residuals)

  ####### Use Last_fit ngme object to update in_list
    stopifnot("last_fit should be an ngme object"
      = class(last_fit) == "ngme" || is.null(last_fit))

    if (inherits(last_fit, "ngme")) {
      output <- last_fit$est_output
      # use last fit estimates

      # block noise
      noise$theta_mu    <- output$noise[["theta_mu"]]
      noise$theta_sigma <- output$noise[["theta_sigma"]]
      noise$theta_noise <- output$noise[["theta_V"]]
      noise$V           <- output$noise[["V"]]

      # latents
      for (i in seq_along(latents_in)) {
        last_fit_latent <- output$latent[[i]]
        latents_in[[i]]$operator$theta_K  <- last_fit_latent[["theta_K"]]
        latents_in[[i]]$W                 <- last_fit_latent[["W"]]

        latents_in[[i]]$noise$theta_mu    <- last_fit_latent[["theta_mu"]]
        latents_in[[i]]$noise$theta_sigma <- last_fit_latent[["theta_sigma"]]
        latents_in[[i]]$noise$theta_noise <- last_fit_latent[["theta_V"]]
        latents_in[[i]]$noise$V           <- last_fit_latent[["V"]]
      }
    }

    in_list <- list(general_in = general_in,
                    latents_in  = latents_in,
                    noise_in    = noise,
                    control_in  = controls,
                    debug       = debug,
                    seed        = seed)

  } else {
    stop("unknown structure of formula")
  }

# print
if (debug$debug) print(str(in_list))

  ################# Run CPP ####################
  if (debug$not_run) {
    print(str(in_list))
    stop("Start estimation by setting not_run = FALSE")
  }

  out <- estimate_cpp(in_list)

################# Construct Output ####################
    # out$input <- in_list

    # # fix_eff
    # out$n_fe     = ncol(X)

    # # m_err
    # out$family = general_in$family

    # # operator
    # out$n_la_params = unlist(lapply(latents_in, function(x) x["n_la_params"] ))

    # # process
    # out$n_latent = length(latents_in)
    # out$model.types = model.types
    # out$var.types   = var.types

  class(out) <- "ngme"
  print(paste("total time is", Sys.time() - time.start))

  out
}