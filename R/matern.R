######################### models #########################

#' Create a Matern SPDE model
#'
#' @param alpha
#' @param mesh mesh argument
#' @param fem.mesh.matrices specify the FEM matrices
#' @param d indicating the dimension of mesh (together with fem.mesh.matrices)
#' @param kappa # parameterization for k^2 C + G, only for stationary
#' @param theta_kappa
#' @param B_kappa bases for kappa
#'
#' @return a list (n, C (diagonal), G, B.kappa) for constructing operator
#' @export
#'
#' @examples
model_matern <- function(
  index       = NULL,
  replicates  = NULL,
  alpha       = 2,
  kappa       = 1,
  theta_kappa = NULL,
  B_kappa     = NULL,
  mesh        = NULL,
  fem.mesh.matrices = NULL,
  d           = NULL,
  A           = NULL,
  noise       = NULL,
  ...
) {
  if (is.null(mesh) && is.null(fem.mesh.matrices))
    stop("At least specify mesh or matrices")

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
    } else {
      fem <- INLA::inla.mesh.fem(mesh, order = alpha)
      C <- fem$c0  # diag
      G <- fem$g1
    }
  } else {
    C <- fem.mesh.matrices$C
    G <- fem.mesh.matrices$G
  }
  # h <- diag(C)
  h <- rep(1, mesh$n)

  if (!is.null(A)) {
    nrep <- ncol(A) / nrow(C)
    C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), C)
    G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), G)
    h <- rep(h, times = nrep)
  }

  if (noise$n_noise == 1) noise <- update_noise(noise, n = mesh$n)

  model <- ngme_model(
    model       = "matern",
    A           = A,
    W_size      = mesh$n,
    V_size      = nrow(C),
    theta_K     = theta_kappa,
    alpha       = alpha,
    B_kappa     = B_kappa,
    C           = ngme_as_sparse(C),
    G           = ngme_as_sparse(G),
    h           = h,
    noise       = noise,
    ...
  )
  model
}
