#' f function
#'
#' @param x numerical vector, which is the covariate
#' @param model specify a model
#' @param var variance component type
#' @param config control variables
#' @param debug debug
#'
#' @return a list of objects
#' @export
f <- function(x = NULL,
              model  = "ar1",
              var    = "NIG",
              control = control.f(),
              debug  = FALSE
              ) {
  ## in ... is the name of the covariate  and possibly the location of the weights
  ## like f(x, weights)

  # construct operator
  if (model=="ar1") {
      n = length(x)

    # G
      G <- Matrix::Matrix(diag(n));
      G <- as(G, "dgCMatrix");

    # C
      C <- Matrix::Matrix(0, n, n)
      C[seq(2, n*n, by=n+1)] <- -1
      C <- as(C, "dgCMatrix");

    # A
      A <- as(Matrix::Matrix(diag(n)), "dgCMatrix");

    # h
      h <- rep(1.0, n)

    operator_in   = list(C=C, G=G, numerical_dK=FALSE)
  }
  else if (model=="matern1d") {
    model = "matern"

    mesh <- INLA::inla.mesh.1d(x)
    fem <- INLA::inla.mesh.1d.fem(mesh)
    n <- mesh$n
print(n)

    # G
    G <- fem$g1

    # C
    C <- fem$c1

    # A
    A <- INLA::inla.spde.make.A(mesh, loc = x)

    # h
    h <- rep(1.0, n)

    operator_in   = list(C=C, G=G, numerical_dK=FALSE)
  }
  else {
    stop("unknown model")
  }

  # construct variance component
  if (var=="nig" || var=="NIG") {
    var_in = list(type = "ind_IG")
  }
  else if (var=="normal") {
    var_in = list(type = "normal")
  }
  else {
    stop("unknown var name")
  }

  # construct latent_in
  la_in <- list(type          = model,
                var.type      = var,
                n_mesh        = n,        # !: make sure this is the second place
                A             = A,
                h             = h,
                opt_kappa     = control$opt_kappa,
                opt_mu        = control$opt_mu,
                opt_sigma     = control$opt_sigma,
                opt_var       = control$opt_var,
                numer_grad    = control$numer_grad,
                use_precond   = control$use_precond,
                eps           = control$eps,
                debug         = debug,
                operator_in   = operator_in,
                var_in        = var_in,
                init_value    = list(kappa     = control$init_kappa,
                                     mu        = control$init_mu,
                                     sigma     = control$init_sigma,
                                     nu        = control$init_nu)
                )

  return (la_in)
}
