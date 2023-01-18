######################### models #########################

#' Create a Matern SPDE model
#'
#' @param loc       numeric vector (1d) or matrix of column 2 (2d),
#'     location to make index
#' @param replicates replicates for the process
#' @param alpha     2 or 4, SPDE smoothness parameter
#' @param mesh      mesh argument
#' @param index_NA Logical vector, same as is.na(response var.)
#' @param kappa     parameterization for kappa^2 C + G, only for stationary
#' @param theta_kappa parameterization for non-stationary
#' @param B_kappa   bases for kappa
#' @param d         indicating the dimension of mesh (together with fem.mesh.matrices)
#' @param fem.mesh.matrices specify the FEM matrices
#' @param noise     1. string: type of model, 2. ngme.noise object
#' @param A         A Matrix connecting observation and mesh
#' @param A_pred    A Matrix connecting NA location and mesh
#' @param ... extra arguments in f()
#'
#' @return a list (n, C (diagonal), G, B.kappa) for constructing operator
#' @export
model_matern <- function(
  loc         = NULL,
  replicates  = NULL,
  alpha       = 2,
  kappa       = 1,
  theta_kappa = NULL,
  B_kappa     = NULL,
  mesh        = NULL,
  fem.mesh.matrices = NULL,
  d           = NULL,
  index_NA    = NULL,
  A           = NULL,
  A_pred      = NULL,
  noise       = noise_normal(),
  ...
) {
  # loc <- eval(substitute(loc), envir = data, enclos = parent.frame())

  if (is.null(mesh) && is.null(fem.mesh.matrices))
    stop("At least specify mesh or matrices")

  if (is.vector(mesh)) # mesh is 1d vector
    mesh <- INLA::inla.mesh.1d(loc = mesh)

  # deal with coords
  if ((is.data.frame(loc) || is.matrix(loc)) && ncol(loc) == 2) {
    loc <- as.matrix(loc)
    if (is.null(index_NA)) index_NA <- rep(FALSE, nrow(loc))
    if (is.null(A))      A <- INLA::inla.spde.make.A(mesh = mesh, loc = loc[!index_NA, ,drop = FALSE])
    if (is.null(A_pred)) A_pred <- INLA::inla.spde.make.A(mesh = mesh, loc = loc[index_NA, ,drop = FALSE])
  } else { # 1d case
    if (is.null(index_NA)) index_NA <- rep(FALSE, length(loc))
    if (is.null(A))      A <- INLA::inla.spde.make.A(mesh = mesh, loc = loc[!index_NA])
    if (is.null(A_pred) && any(index_NA)) A_pred <- INLA::inla.spde.make.A(mesh = mesh, loc = loc[index_NA])
  }


  if (alpha - round(alpha) != 0) {
    stop("alpha should be integer, now only 2 or 4")
  }

  stopifnot(alpha == 2 || alpha == 4)

  # use kappa as parameterization
  stopifnot("Don't overwrite kappa with NULL." = !is.null(kappa))
  stopifnot("kappa is only for stationary case." = length(kappa) == 1)
  stopifnot("kappa is greater than 0." = kappa > 0)
  if (is.null(theta_kappa)) theta_kappa <- log(kappa)

  if (is.null(B_kappa))
    B_kappa <- matrix(1, nrow = mesh$n, ncol = length(theta_kappa))

  # supply mesh
  if (!is.null(mesh)) {
    d <- get_inla_mesh_dimension(mesh)
    if (d == 1) {
      fem <- INLA::inla.mesh.1d.fem(mesh)
      C <- fem$c1
      G <- fem$g1
      h <- Matrix::diag(fem$c0)
      # h <- rowsum(fem$c0)
    } else {
      fem <- INLA::inla.mesh.fem(mesh, order = alpha)
      C <- fem$c0  # diag
      G <- fem$g1
      h <- Matrix::diag(fem$c0)
    }
  } else {
    C <- fem.mesh.matrices$C
    G <- fem.mesh.matrices$G
    h <- Matrix::diag(fem.mesh.matrices$c0)
  }

  if (!is.null(A)) {
    nrep <- ncol(A) / nrow(C)
    C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), C)
    G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), G)
    h <- rep(h, times = nrep)
  }

  if (noise$n_noise == 1) noise <- update_noise(noise, n = mesh$n)
  noise$h <- h
  # kappas <- drop(exp(theta_kappa %*% B_kappa))

  model <- ngme_model(
    model       = "matern",
    A           = A,
    A_pred      = A_pred,
    W_size      = mesh$n,
    V_size      = nrow(C),
    theta_K     = theta_kappa,
    alpha       = alpha,
    B_kappa     = B_kappa,
    C           = ngme_as_sparse(C),
    G           = ngme_as_sparse(G),
    # C           = asdgCMatrix_(C),
    # G           = asdgCMatrix_(G),
    # K           = kappas * kappas * C + G
    h           = h,
    noise       = noise,
    ...
  )
  model
}
