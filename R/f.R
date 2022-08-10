#' f function for specifying a model
#'
#' @param ...
#' @param model     1. string: type of model, 2. ngme.spde object
#' @param noise     1. string: type of model, 2. ngme.noise object
#' @param replicates   Representing the replicates (can also be specified in each ngme model)
#' @param control      control variables for f model
#' @param A            A Matrix connecting observation and
#' @param B.sigma      Basis matrix for sigma
#' @param B.mu         Basis matrix for mu
#' @param theta.sigma  Starting value for theta.sigma
#' @param theta.mu     Starting value for theta.sigma
#' @param theta.K      Starting value for theta.sigma
#' @param theta.noise  Starting value for theta.sigma
#' @param start.V      Starting value for V
#' @param debug        Debug variables
#' @param data      if not specified or NULL, inherit the data from ngme function
#'
#' @return a list latent_in for constructing latent model, e.g. A, h, C, G,
#' which also has
#' 1. list operator_in for building operator,
#' 2. list var_in for variance component,
#' 3. list init_values of parameters
#'
#' @export
f <- function(
  ...,
  model  = "ar1",
  noise = ngme.noise(),
  control = ngme.control.f(),
  debug  = FALSE,
  A = NULL,
  B.sigma = 1, # non-stationary case -> into matrix n_mesh * n_sigma
  B.mu = 1,
  theta.mu = 0,
  theta.sigma = 0, # exp(0) = 1
  # for these two you can specified inside ngme.noise and ngme.model function
  replicates=NULL,
  theta.K = NULL,
  theta.noise = NULL,
  start.V = NULL,
  data = NULL
) {
  ################## Index and Replicates ##################

  # get the index from ddd
  ddd <- match.call(expand.dots = FALSE)$...
  if (length(ddd) > 1) {
    stop(paste("To many variables included in f():", paste(unlist(ddd), collapse=" ")))
  }
  numerical.or.symbol <- ddd[[1]]
  index <- eval(numerical.or.symbol, envir=data)

  # get the replicates
  if (is.null(replicates))
    replicates = rep(1, length(index))
  nrep <- length(unique(replicates))

  # no need for re-order
  # # re-order the values according to the replicates (to be block diagonal for C and G)
  # df <- data.frame(original.order=1:length(index), replicates=replicates, index=index)
  # df <- df[order(df$replicates), ]

  if(length(index)!= length(replicates)){
    stop("The index and the replicates should have the same length!")
  }

  # if (model=="ar1") {

  # }

  ################## construct Mu and Sigma #################
  # turn B.sigma and B.mu into matrix

  n_mu <- length(theta.mu)
  n_sigma <- length(theta.sigma)

  ################## construct operator (n_ope, C, G, A, h) ##################
  if (is.character(model)) {  ######## string
    if (model=="ar1") {
      model.type = "ar1"

      ar1_in = ngme.ar1(
        index=index,
        replicates=replicates
      )
      A = ar1_in$A
      n = ncol(A)

      h = rep(1.0, n)
      B.mu <- matrix(B.mu, nrow=n, ncol=n_mu)
      B.sigma <- matrix(B.sigma, nrow=n, ncol=n_sigma)

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

    B.mu <- matrix(B.mu, nrow=n, ncol=n_mu)
    B.sigma <- matrix(B.sigma, nrow=n, ncol=n_sigma)

    operator_in = model$operator_in
  }
  else if (inherits(model, "ngme.matern")) {  ########  stationary Matern ########
    # read in operator and modify
    operator_in = model$operator_in

    model.type = "matern"
    if (is.null(A)) stop("Provide A matrix")

    # compute nrep
    nrep <- ncol(A)/nrow(operator_in$C)

    if(!is.null(nrep)){
      operator_in$C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), operator_in$C)
      operator_in$G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), operator_in$G)
    }

    operator_in$C = ngme.as.sparse(operator_in$C)
    operator_in$G = ngme.as.sparse(operator_in$G)

    n = length(index)
    h = rep(1, n)
    theta.K = model$kappa

    B.mu <- matrix(B.mu, nrow=n, ncol=n_mu)
    B.sigma <- matrix(B.sigma, nrow=n, ncol=n_sigma)

    B.mu <- kronecker(matrix(1, ncol=1, nrow=nrep), B.mu)
    B.sigma <- kronecker(matrix(1, ncol=1, nrow=nrep), B.sigma)
    h <- rep(h, times=nrep)
  }
  else if (inherits(model, "ngme.spde")) { ######## nonstationary Matern ########
    model.type = "spde.matern"

    # compute nrep
    operator_in <- model$operator_in
    nrep <- ncol(A)/nrow(operator_in$C)

    if(!is.null(nrep)){
      operator_in$C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), operator_in$C)
      operator_in$G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), operator_in$G)
    }
    C = ngme.as.sparse(C)
    G = ngme.as.sparse(G)

    n = length(index)
    if (is.null(A)) stop("Provide A matrix")
    
    n <- model$operator_in$n
    h <- rep(1, n)
    theta.K = model$theta.kappa
    operator_in <- model$operator_in
  }
  else {
    stop("unknown model")
  }

  ################## construct noise (e.g. nig noise) ##################
  if (is.character(noise)) {
    if (noise=="nig")                         
      noise <- ngme.noise(type = "nig")
    else if (noise=="normal" || noise=="gaussian") 
      noise <- ngme.noise(type = "normal")
  }
  if (is.null(theta.noise) && !(is.null(noise$theta.noise))) theta.noise = noise$theta.noise

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
