######################### models #########################

#' Create a Matern SPDE model
#'
#' @param loc       numeric vector (1d) or matrix of column 2 (2d),
#'     location to make index
#'     keep aligned with Y if make prediction with index_NA!!!
#' @param replicate replicate for the process
#' @param alpha     2 or 4, SPDE smoothness parameter
#' @param mesh      mesh argument
#' @param index_NA Logical vector, same as is.na(response var.)
#' @param kappa     parameterization for kappa^2 C + G, only for stationary
#' @param theta_kappa parameterization for non-stationary
#' @param B_kappa   bases for kappa
#' @param d         indicating the dimension of mesh (together with fem.mesh.matrices)
#' @param noise     1. string: type of model, 2. ngme.noise object
#' @param ... extra arguments in f()
#'
#' @return a list (n, C (diagonal), G, B.kappa) for constructing operator
#' @export
model_matern <- function(
  loc         = NULL,
  replicate   = NULL,
  alpha       = 2,
  kappa       = 1,
  theta_kappa = NULL,
  B_kappa     = NULL,
  mesh        = NULL,
  # fem.mesh.matrices = NULL,
  d           = NULL,
  index_NA    = NULL,
  # A           = NULL,
  # A_pred      = NULL,
  noise       = noise_normal(),
  ...
) {
  if (is.numeric(mesh)) # mesh is 1d vector
    mesh <- INLA::inla.mesh.1d(loc = mesh)

  if (is.null(loc)) {
    if (inherits(mesh, "inla.mesh.1d")) loc <- mesh$loc
    if (inherits(mesh, "inla.mesh")) loc <- as.matrix(mesh$loc[, 1:2])
  }
stopifnot("loc is NULL" = !is.null(loc))
  n_loc <- if(is.null(dim(loc))) length(loc) else nrow(loc)
  # loc <- eval(substitute(loc), envir = data, enclos = parent.frame())

  if (is.null(mesh)) stop("Please provide mesh!")

  # deal with coords
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
  }
  # else {
  #   C <- fem.mesh.matrices$C
  #   G <- fem.mesh.matrices$G
  #   h <- Matrix::diag(fem.mesh.matrices$c0)
  # }

  tmp <- ngme_make_A(
    mesh = mesh,
    map = loc,
    n_map = n_loc,
    idx_NA = index_NA,
    replicate = replicate
  )
  A <- tmp$A; A_pred <- tmp$A_pred

  if (!is.null(A)) {
    nrep <- ncol(A) / ncol(C)
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
    W_size      = nrow(C),
    V_size      = ncol(C),
    theta_K     = theta_kappa,
    alpha       = alpha,
    B_kappa     = B_kappa,
    C           = ngme_as_sparse(C),
    G           = ngme_as_sparse(G),
    # K           = kappas * kappas * C + G
    h           = h,
    noise       = noise,
    mesh        = mesh,
    map         = loc,
    n_map       = n_loc,
    replicate   = replicate,
    ...
  )
  model
}
