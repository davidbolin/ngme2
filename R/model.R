# function for specify ngme.model basic structure
ngme.model <- function(
  model,
  n_mesh      = NULL,
  theta_K     = NULL,
  fix_theta_K = FALSE,
  W           = NULL,
  fix_W       = NULL,
  A           = NULL,
  A_pred      = NULL,
  noise_type  = NULL,
  noise       = list(),
  control     = ngme.control.f(),
  ...
) {
  stopifnot(is.character(model))

  structure(
    list(
      model         = model,
      n_mesh        = n_mesh,
      theta_K       = theta_K,
      n_theta_K     = length(theta_K),
      A             = A,
      A_pred        = A_pred,
      fix_theta_K   = fix_theta_K,
      noise_type    = noise_type,
      noise         = noise,
      W             = W,
      fix_W         = fix_W,
      control       = control,
      ...
    ),
    class = "ngme_model"
  )
}

#' ngme ar1 model specification
#' 
#' Generating C, G and A given index and replicates
#'
#' @param index index for the process
#' @param replicates replicates for the process
#' @param alpha initial value for alpha
#' @param range range for the mesh
#'
#' @return a list of specification of model
#' @export
#'
#' @examples
ngme.ar1 <- function(
  index,
  replicates = NULL,
  alpha = 0.5,
  range = c(1, max(index)),

  index_pred = NULL,
  use_num_dK = FALSE,
  ...
) {
  # overwirte the default
  if (!is.null(list(...)$theta_K)) alpha <- list(...)$theta_K
  if (is.null(replicates)) replicates <- rep(1, length(index))

  unique_rep <- unique(replicates)
  nrep <- length(unique_rep)

  # e.g. index = c(1:200, 1:100)
  #      replicates = c(rep(1, 200), rep(2, 100))
  #      n =200 in this case

  n <- range[2] - range[1] + 1

  # construct G
    G <- Matrix::Matrix(diag(n));

  # construct C
    C <- Matrix::Matrix(0, n, n)
    C[seq(2, n*n, by = n+1)] <- -1

    if(!is.null(nrep)) {
      C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), C)
      G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), G)
    }

  model <- ngme.model(
    model       = "ar1",
    theta_K     = alpha,
    n_mesh      = n,
    A           = ngme.ts.make.A(loc = index, replicates = replicates, range = range),
    A_pred      = ngme.ts.make.A(index_pred, replicates = replicates, range = range),
    h           = rep(1.0, n),
    C = ngme.as.sparse(C),
    G = ngme.as.sparse(G)
  )

  model
}

#' Create a Matern SPDE model
#'
#' @param alpha
#' @param mesh mesh argument
#' @param fem.mesh.matrices specify the FEM matrices
#' @param d indicating the dimension of mesh (together with fem.mesh.matrices)
#' @param theta_kappa 
#' @param B_kappa bases for kappa
#'
#' @return a list (n, C (diagonal), G, B.kappa) for constructing operator
#' @export
#'
#' @examples
ngme.matern <- function(
  index = NULL,
  alpha = 2,
  mesh = NULL,
  replicates = NULL,
  fem.mesh.matrices = NULL,
  d = NULL,
  theta_kappa = 0,
  A = NULL,        # watch out! Can also specify in f function, not used for now
  B_kappa = NULL
) {
  if (is.null(mesh) && is.null(fem.mesh.matrices)) 
    stop("At least specify mesh or matrices")

  if (alpha - round(alpha) != 0) {
    stop("alpha should be integer, now only 2 or 4")
  }

  stopifnot(alpha == 2 || alpha == 4)

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

  if (!is.null(A)) {
    nrep <- ncol(A) / nrow(C)
    C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), C)
    G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), G)
  }

  spde <- list(
    # general
    A = A,
    n_params = length(theta_kappa),
    operator = list(
      theta_kappa = theta_kappa,
      alpha = alpha,
      B_kappa = B_kappa,
      n_params = length(theta_kappa),
      n_mesh = mesh$n,
      C = ngme.as.sparse(C),
      G = ngme.as.sparse(G),
      use_num_dK = FALSE
    )
  )

  # create precision matrix
  class(spde) <- "ngme.matern"
  spde
}

# rw1, rw2
# nodes = 100 (inla.group)

# ?inla.spde.make.A
# inla.spde.make.A(
#   index
#  replicates=)

# ngme.rw1 <- function(
#   index,
#   replicates = NULL,

# ) {

# }