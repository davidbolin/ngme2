#' Specifying a latent model
#'
#' Function used for defining of smooth and spatial terms 
#' within ngme model formulae. 
#' The function does not evaluate anything - 
#' it exists purely to help set up a model. 
#' (see ngme.models.types for available models).
#'
#' @param ...    index
#' @param model     1. string: type of model, 2. ngme.spde object
#' @param noise     1. string: type of model, 2. ngme.noise object
#' @param replicates   Representing the replicates 
#'  (can also be specified in each ngme model)
#' @param control      control variables for f model
#' @param A            A Matrix connecting observation and
#' @param theta_K      Unbounded parameter for K
#' @param debug        Debug variables
#' @param data      if not specified or NULL, then 
#' @param W          starting realization of the process
#'  inherit the data from ngme function
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
  replicates = NULL,
  model  = "ar1",
  noise = ngme.noise(),
  data = NULL,
  control = ngme.control.f(),
  debug  = FALSE,
  A = NULL,
  theta_K = NULL,
  W = NULL
) {
  ################## Index and Replicates ##################

  if (is.null(W) && control$fix_W) stop("Provide initial W to use fix_W")
  # get the index from ddd
  ddd <- match.call(expand.dots = FALSE)$...
  if (length(ddd) > 1) {
    stop(paste("To many variables included in f():", paste(unlist(ddd), collapse=" ")))
  }
  numerical.or.symbol <- ddd[[1]]
  index <- eval(numerical.or.symbol, envir = data)

  # get the replicates
  if (is.null(replicates))
    replicates = rep(1, length(index))
  nrep <- length(unique(replicates))

  # no need for re-order
  # # re-order the values according to the replicates (to be block diagonal for C and G)
  # df <- data.frame(original.order=1:length(index), replicates=replicates, index=index)
  # df <- df[order(df$replicates), ]

  stopifnot("The index and the replicates should have the same length!"
    = length(index) == length(replicates))

  ################## construct operator (n_ope, C, G, A, h) ##################
  nrep <- NULL
  n_mesh <- NULL

  if (is.character(model)) {  ######## string
    if (model == "ar1") {
      model_type <- "ar1"

      ar1_in = ngme.ar1(
        index = index,
        replicates = replicates
      )
      A <- ar1_in$A
      n_mesh <- ncol(A)

      h <- rep(1.0, n_mesh)
      operator_in = ar1_in$operator_in
    } else {
      stop("unknown model name")
    }
  }
  else if (inherits(model, "ngme.ar1")) {
    model_type = "ar1"

    n_mesh = ncol(model$A)
    h = rep(1.0, n_mesh)
    theta_K = model$alpha

    operator_in = model$operator_in
  }
  else if (inherits(model, "ngme.matern")) {  # stationary Matern
    # read in operator and modify
    operator_in = model$operator_in

    model_type = "matern"
    if (is.null(A)) stop("Provide A matrix")

    # replicates
    nrep <- ncol(A)/nrow(operator_in$C)
    operator_in$C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), operator_in$C)
    operator_in$G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), operator_in$G)
    operator_in$C <- ngme.as.sparse(operator_in$C)
    operator_in$G <- ngme.as.sparse(operator_in$G)

    n_mesh = length(index)
    h = rep(1, n_mesh)
    theta_K = model$kappa

    h <- rep(h, times=nrep)
  }
  else if (inherits(model, "ngme.spde")) { # nonstationary Matern
    model_type = "spde.matern"

    # compute nrep
    operator_in <- model$operator_in
    nrep <- ncol(A)/nrow(operator_in$C)
    operator_in$C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), operator_in$C)
    operator_in$G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), operator_in$G)
    operator_in$C <- ngme.as.sparse(operator_in$C)
    operator_in$G <- ngme.as.sparse(operator_in$G)

    n_mesh = length(index)
    if (is.null(A)) stop("Provide A matrix")

    h <- rep(1, n_mesh)
    theta_K = model$theta_Kappa
    operator_in <- model$operator_in
  }
  else {
    stop("unknown model")
  }

  ################## construct noise (e.g. nig noise) ##################
    # ?? check
    B_mu <- matrix(noise$B_mu, nrow = n_mesh, ncol = noise$n_theta_mu)
    B_sigma <- matrix(noise$B_sigma, nrow = n_mesh, ncol = noise$n_theta_sigma)

    # replicates
    if (is.integer(nrep)) {
      B_mu <- kronecker(matrix(1, ncol = 1, nrow = nrep), B_mu)
      B_sigma <- kronecker(matrix(1, ncol = 1, nrow = nrep), B_sigma)
    }

  # total params
  n_la_params = operator_in$n_params + noise$n_theta_mu + noise$n_theta_sigma + noise$n_theta_V

  # check initial value of W
  if (!is.null(W)) stopifnot(length(W) == n_mesh)

  # overwrites
  if (!is.null(theta_K)) operator_in$theta_K <- theta_K

  latent_in <- list(
    model_type  = model_type,
    noise_type  = noise$type,
    n_mesh      = n_mesh,        # !: make sure this is the second place
    A           = A,
    h           = h,
    n_la_params = n_la_params,
    W           = W,

    # mu and sigma
    B_mu          = B_mu,
    B_sigma       = B_sigma,
    n_theta_mu    = noise$n_theta_mu,
    n_theta_sigma = noise$n_theta_sigma,

    # lists
    operator      = operator_in,
    noise         = update.ngme.noise(noise, n = n_mesh),
    control_f     = control,
    debug         = debug
  )

  class(latent_in) <- "latent"
  return (latent_in)
}

# what is a valid operator_in (for constructing K)
