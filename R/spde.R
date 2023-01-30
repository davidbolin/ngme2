# #### deprecated, delete soon

# #' Create a Matern SPDE model
# #'
# #' @param alpha
# #' @param mesh mesh argument
# #' @param fem.mesh.matrices specify the FEM matrices
# #' @param d indicating the dimension of mesh (together with fem.mesh.matrices)
# #' @param theta.kappa
# #' @param B.kappa bases for kappa
# #'
# #' @return a list (n, C (diagonal), G, B.kappa) for constructing operator
# #' @export
# #'
# #' @examples

# ngme.nonstation.matern <- function(
#     alpha = 2,
#     mesh = NULL,
#     replicate = NULL,
#     fem.mesh.matrices = NULL,
#     d = NULL,
#     B.kappa = matrix(1, 1, 2),
#     theta.kappa = c(0, 0)
#     )
# {
#   if (is.null(mesh) && is.null(fem.mesh.matrices)) stop("At least specify mesh or matrices")
#   if (alpha - round(alpha) != 0) {
#     stop("alpha should be integer")
#   }

#   # supply mesh
#   if (!is.null(mesh)) {
#     n <- mesh$n
#     d <- get_inla_mesh_dimension(mesh)
#     if (d == 1) {
#       fem <- INLA::inla.mesh.1d.fem(mesh)
#       C <- fem$c1
#       G <- fem$g1
#     } else {
#       fem <- INLA::inla.mesh.fem(mesh, order = alpha)
#       C <- fem$c0 # diag
#       G <- fem$g1
#     }

#     n <- mesh$n
#     spde.spec <- list(
#       # general
#       n_params = ncol(B.kappa),
#       theta.kappa = theta.kappa,

#       # spde
#       operator_in = list(
#         alpha = alpha,
#         B.kappa = B.kappa,
#         n_params = length(theta.kappa),
#         n = n,
#         C = ngme_as_sparse(C),
#         G = ngme_as_sparse(G),
#         use_num_dK = FALSE
#       )
#     )
#   }

#   # create precision matrix
#   class(spde.spec) <- "ngme.spde"

#   return(spde.spec)
# }


# #' Create Matern spde model object with stationary kappa
# #'
# #' @param alpha
# #' @param mesh
# #' @param fem.mesh.matrices
# #' @param d
# #' @param kappa
# #'
# #' @return
# #' @export
# #'
# #' @examples
# ngme.station.matern <- function(
#   mesh = NULL,
#   nrep = 1,
#   alpha = 2,
#   fem.mesh.matrices = NULL,
#   d = NULL,
#   kappa = 1,
#   use_num_dK = FALSE
# )
# {
#   if (is.null(mesh) && is.null(fem.mesh.matrices)) stop("At least specify mesh or matrices")
#   if (alpha - round(alpha) != 0) {
#     stop("alpha should be integer")
#   }

#   # supply mesh
#   if (!is.null(mesh)) {
#     n <- mesh$n
#     d <- get_inla_mesh_dimension(mesh)
#     if (d == 1) {
#       fem <- INLA::inla.mesh.1d.fem(mesh)
#       C <- fem$c1
#       G <- fem$g1
#     } else {
#       fem <- INLA::inla.mesh.fem(mesh, order = alpha)
#       C <- fem$c0 # diag
#       G <- fem$g1
#     }

#     n <- mesh$n
#     spde.spec <- list(
#       # general
#       n_params = 1,
#       kappa = kappa,
#       # spde
#       operator_in = list(
#         alpha = alpha,
#         n_params = 1,
#         n = n,
#         C = C,
#         G = G,
#         use_num_dK = use_num_dK
#       )
#     )
#     class(spde.spec) <- "ngme.matern"
#   }

#   spde.spec
# }


