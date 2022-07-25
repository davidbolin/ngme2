#' Create Matern spde model object with nonstationary kappa
#'
#' @param alpha
#' @param mesh mesh argument
#' @param fem.mesh.matrices specify the FEM matrices
#' @param d indicating the dimension of mesh (together with fem.mesh.matrices)
#' @param B.kappa bases for kappa
#' @param theta.init
#'
#' @return a list (n, C (diagonal), G, B.kappa) for constructing operator
#' @export
#'
#' @examples

ngme.spde.matern <- function(
    alpha = 2,
    mesh = NULL,
    replicates = NULL,
    fem.mesh.matrices = NULL,
    d = NULL,
    B.kappa = matrix(1, 1, 2),
    theta.kappa = c(0, 0)
    )
{
  if (is.null(mesh) && is.null(fem.mesh.matrices)) stop("At least specify mesh or matrices")
  if (alpha - round(alpha) != 0) {
    stop("alpha should be integer")
  }

  # supply mesh
  if (!is.null(mesh)) {
    n <- mesh$n
    d <- get_inla_mesh_dimension(mesh)
    if (d == 1) {
      fem <- INLA::inla.mesh.1d.fem(mesh)
      C <- fem$c1
      G <- fem$g1
    } else {
      fem <- INLA::inla.mesh.fem(mesh, order = alpha)
      C <- fem$c0 # diag
      G <- fem$g1
    }

    if (is.null(replicates)) replicates <- rep(1, n)
    nrep <- length(unique(replicates))
    if(!is.null(nrep)){
      C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), C)
      G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), G)
    }

    n <- mesh$n
    spde.spec <- list(
      # general
      n_params = ncol(B.kappa),
      theta.kappa = theta.kappa,

      # spde
      operator_in = list(
        alpha = alpha,
        B.kappa = B.kappa,
        n_params = length(theta.kappa),
        n = n,
        C = ngme.as.sparse(C),
        G = ngme.as.sparse(G),
        use_num_dK = FALSE
      )
    )
  } else {
    stop("please supply mesh")
  }

  # create precision matrix
  class(spde.spec) <- "ngme.spde"

  return(spde.spec)
}


#' Create Matern spde model object with stationary kappa
#'
#' @param alpha
#' @param mesh
#' @param fem.mesh.matrices
#' @param d
#' @param kappa
#'
#' @return
#' @export
#'
#' @examples
ngme.matern <- function(
    mesh = NULL,
    replicates = NULL,
    alpha = 2,
    fem.mesh.matrices = NULL,
    d = NULL,
    kappa = 1,
    use_num_dK = FALSE
)
{
  if (is.null(mesh) && is.null(fem.mesh.matrices)) stop("At least specify mesh or matrices")
  if (alpha - round(alpha) != 0) {
    stop("alpha should be integer")
  }

  # supply mesh
  if (!is.null(mesh)) {
    n <- mesh$n
    d <- get_inla_mesh_dimension(mesh)
    if (d == 1) {
      fem <- INLA::inla.mesh.1d.fem(mesh)
      C <- fem$c1
      G <- fem$g1
    } else {
      fem <- INLA::inla.mesh.fem(mesh, order = alpha)
      C <- fem$c0 # diag
      G <- fem$g1
    }

    if (is.null(replicates)) replicates <- rep(1, n)
    nrep <- length(unique(replicates))
    if(!is.null(nrep)){
      C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), C)
      G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), G)
    }

    n <- mesh$n
    spde.spec <- list(
      # general
      n_params = 1,
      kappa = kappa,

      # spde
      operator_in = list(
        alpha = alpha,
        n_params = 1,
        n = n,
        C = ngme.as.sparse(C),
        G = ngme.as.sparse(G),
        use_num_dK = use_num_dK
      )
    )
    class(spde.spec) <- "ngme.matern"
  }

  spde.spec
}


#' @name get_inla_mesh_dimension
#' @title Get the dimension of an INLA mesh
#' @description Get the dimension of an INLA mesh
#' @param inla_mesh An INLA mesh
#' @return The dimension of an INLA mesh.
#' @noRd
#'
get_inla_mesh_dimension <- function(inla_mesh) {
  cond1 <- inherits(inla_mesh, "inla.mesh.1d")
  cond2 <- inherits(inla_mesh, "inla.mesh")
  stopifnot(cond1 || cond2)
  if (inla_mesh$manifold == "R1") {
    d <- 1
  } else if (inla_mesh$manifold == "R2") {
    d <- 2
  } else {
    stop("The mesh should be from a flat manifold.")
  }
  return(d)
}


# # test
# mesh_size = 5
# range = 0.2
# sigma = 1
# nu=0.8
# kappa = sqrt(8*nu)/range
# mesh_grid = inla.mesh.lattice(x=seq(from=0,to=1, length=mesh_size),y=seq(from=0,to=1,length=mesh_size))
# mesh_grid_2d = inla.mesh.create(lattice = mesh_grid, extend = FALSE, refine = FALSE)
# mesh_2d = mesh_grid_2d
# obs_coords = mesh_2d$loc[,c(1,2)]
# nob = nrow(obs_coords)
# fem = inla.mesh.fem(mesh_2d, order = 2)
# C = fem$c0
# G = fem$g1
# A = inla.spde.make.A(mesh_2d, loc = obs_coords)
#
# o = ngme.spde.matern(alpha=2, mesh=mesh_2d)
# o$Q(1, 1)
#
# o$Q(1,1)
