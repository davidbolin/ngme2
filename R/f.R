#' f function
#'
#' @param x covariates
#' @param model 1. string: type of model, 2. ngme.spde object
#' @param var 1. string: type of variance component
#' @param config control variables
#' @param debug debug variables
#'
#' @return a list latent_in for constructing latent model, e.g. A, h, C, G,
#' which also has
#' 1. list ope_in for building operator,
#' 2. list var_in for variance component,
#' 3. list init_values of parameters
#'
#' @export
f <- function(
  x = NULL,
  model  = "ar1",
  var    = "nig",
  control = control.f(),
  debug  = FALSE,
  A = NULL
  ) {

  # construct variance component
  build.var <- function(var) {
    if (var=="nig") {
      var_in = list(
        # V_type = "ind_IG",
        nu = control$init_var
      )
    }
    else if (var=="normal") {
      var_in = list(
        # V_type = "normal"
      )
    }
    else {
      stop("unknown var name")
    }
    var_in
  }

  # 1. construct operator (n_covar, n_ope, C, G, A, h)

  ################ SPDE ##################
  if (inherits(model, "ngme.spde")) {
    # # construct latent_in
    # latent_in <- list(
    #   model_type  = model.type,
    #   var_type    = var,
    #   n_mesh      = n,        # !: make sure this is the second place
    #   n_la_params = n_ope_params + 3,
    #   A           = A,
    #   h           = h,
    #
    #   # lists
    #   la_model_in   = la_model_in,
    #   operator_in   = operator_in,
    #   var_in        = var_in,
    #   control_f     = control,
    #   debug         = debug
    # )

    # 1. general
    n = length(x)
    model.type = "spde.matern"
    n_la_params = model$n_params + 3
    var = model$var # var from spde
    h = rep(1, n)

    if (is.null(A)) stop("Provide A matrix")

    # 2. ope
    operator_in <- model$operator_in

    # 3. var
    var_in <- build.var(var)

    special_in <- list()
  }
  ################ AR1 ##################
  else if (model=="ar1") {
    n_ope_params = 1

    # 1. for general latent model
    model.type = "ar1"
    n_la_params = n_ope_params + 3

    # A <- as(Matrix::Matrix(diag(n)), "dgCMatrix");
    A <- ngme.ts.make.A(x)
    n = ncol(A)
    # n = length(x) # length of covariates

    h <- rep(1.0, n)

    # for specific latent model
    special_in <- list()

    # 2. for operator K
    G <- Matrix::Matrix(diag(n));
      G <- as(G, "dgCMatrix");
    C <- Matrix::Matrix(0, n, n)
      C[seq(2, n*n, by=n+1)] <- -1
      C <- as(C, "dgCMatrix");

    operator_in <- list(
      kappa = control$init_operator,
      n_params=n_ope_params,
      C=C,
      G=G,
      use_num_dK=FALSE
    )

    var <- build.var(var)
  }
  ################ MATERN ##################
  else if (model=="matern1d") {
    model = "matern"
    n_params = 1

    mesh <- INLA::inla.mesh.1d(x)
    fem <- INLA::inla.mesh.1d.fem(mesh)
    n <- mesh$n

    G <- fem$g1
    C <- fem$c1
    A <- INLA::inla.spde.make.A(mesh, loc = x)
    h <- rep(1.0, n)

    operator_in  <- list(
      n_params=n_params,
      C=C,
      G=G,
      use_num_dK=FALSE
    )
    var <- build.var(var)

    special_in <- list()
  }
  else {
    stop("unknown model")
  }

  # construct latent_in
  latent_in <- list(
    model_type  = model.type,
    var_type    = var,
    n_mesh      = n,        # !: make sure this is the second place
    n_la_params = n_la_params,
    A           = A,
    h           = h,

    # lists
    special_in    = special_in, # for specific model input
    operator_in   = operator_in,
    var_in        = var_in,
    control_f     = control,
    debug         = debug
    )

  return (latent_in)
}


