# define the interface for Matern2d latent model

# y1 | y2 ~ f(model="matern2d", ..)

#' Create a Matern2D SPDE model
#'
#' @param loc list of length 2, each element is numeric vector (1d) or matrix of column 2 (2d)
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
model_matern2d <- function(
  loc = list(NULL, NULL),
  names = list("matern1", "matern2"),
  mesh = list(NULL, NULL),
  theta_kappa = list(NULL, NULL),
  kappa = list(1, 1),
  rho = 0.5, theta = 1,

  alpha = 2,
  replicate = NULL,
  index_NA = NULL,
  noise = noise_normal(),
  ...
) {
  # check input
  stopifnot(
    length(loc) == 2,
    length(mesh) == 2,
    length(theta_kappa) == 2,
    length(kappa) == 2,
    length(names) == 2
  )

  # check loc / mesh dimension
  # to-do

  # build C and G
  # to-do

  # build A
  # to-do

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