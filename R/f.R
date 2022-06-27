#' f function for specifying a model
#'
#' @param x covariates
#' @param model 1. string: type of model, 2. ngme.spde object
#' @param noise 1. string: type of model, 2. ngme.noise object
#' @param debug debug variables
#' @param control
#' @param A
#' @param B.sigma
#' @param B.mu
#' @param theta.sigma
#' @param theta.mu
#' @param theta.K
#' @param theta.noise
#' @param start.V
#'
#' @return a list latent_in for constructing latent model, e.g. A, h, C, G,
#' which also has
#' 1. list operator_in for building operator,
#' 2. list var_in for variance component,
#' 3. list init_values of parameters
#'
#' @export
f <- function(
  x = NULL,
  model  = "ar1",
  noise = ngme.noise(),
  control = ngme.control.f(),
  debug  = FALSE,
  A = NULL,
  B.sigma = 1, # non-stationary case -> into matrix n_mesh * n_sigma
  B.mu = 1,
  theta.sigma = 0, # exp(0) = 1
  theta.mu = 0,
  # for these two you can specified inside ngme.noise and ngme.model function
  theta.K = NULL,
  theta.noise = NULL,
  start.V = NULL
) {
  ################## construct operator (n_ope, C, G, A, h) ##################
  if (is.character(model)) {  ######## string
    if (model=="ar1") {
      model.type = "ar1"
      ar1_in = ngme.ar1(x)
      A = ar1_in$A
      n = ncol(A)
      h = rep(1.0, n)

      operator_in = ar1_in$operator_in
    } else {
      stop("unknown model name")
    }
  }
  else if (inherits(model, "ngme.ar1")) {   ################ AR1 ##################
    model.type = "ar1"
    A = model$A
    n = ncol(A)
    h = rep(1.0, n)
    theta.K = model$alpha

    operator_in = model$operator_in
  }
  else if (inherits(model, "ngme.matern")) {  ########  stationary Matern ########
    model.type = "matern"
    if (is.null(A)) stop("Provide A matrix")
    n = length(x)
    h = rep(1, n)
    theta.K = model$kappa

    operator_in <- model$operator_in
  }
  else if (inherits(model, "ngme.spde")) { ######## nonstationary Matern ########
    model.type = "spde.matern"
    n = length(x)
    if (is.null(A)) stop("Provide A matrix")
    h = rep(1, n)
    theta.K = model$theta.kappa

    operator_in <- model$operator_in
  }
  else {
    stop("unknown model")
  }

  ################## construct noise (e.g. nig noise) ##################
  if (is.character(noise)) {
    if (noise=="nig")    noise = ngme.noise(type="nig")
    if (noise=="normal") noise = ngme.noise(type="normal")
  }
  if (is.null(theta.noise) && !(is.null(noise$theta.noise))) theta.noise = noise$theta.noise

  ################## construct Mu and Sigma ##################
  # turn B.sigma and B.mu into matrix
  n_mu <- length(theta.mu)
  n_sigma <- length(theta.sigma)

  B.mu <- matrix(B.mu, nrow=n, ncol=n_mu)
  B.sigma <- matrix(B.sigma, nrow=n, ncol=n_sigma)

  # construct latent_in
  latent_in <- list(
    model_type  = model.type,
    var_type    = noise$type,
    n_mesh      = n,        # !: make sure this is the second place
    A           = A,
    h           = h,

    # mu and sigma
    B_mu        = B.mu,
    B_sigma     = B.sigma,
    n_mu        = n_mu,
    n_sigma     = n_sigma,

    n_la_params = operator_in$n_params + noise$n_params + n_mu + n_sigma,

    # lists
    start       = list(
      theta_K     = theta.K,
      theta_mu    = theta.mu,
      theta_sigma = theta.sigma,
      theta_noise = theta.noise,
      V           = start.V
    ),
    operator_in   = operator_in,
    var_in        = noise,
    control_f     = control,
    debug         = debug
  )

  return (latent_in)
}

# what is a valid operator_in (for constructing K)
