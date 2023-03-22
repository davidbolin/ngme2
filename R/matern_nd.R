# define the interface for Matern2d latent model

# y1 | y2 ~ f(model="matern2d", ..)

#' Create a Matern2D SPDE model
#'
#' @param map list of length 2 or 3, each element is numeric vector (1d) or matrix of column 2 (2d)
#' @param replicate replicate for the process
#' @param alpha     2 or 4, SPDE smoothness paramete
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
model_matern_nd <- function(
  map = list(NULL, NULL),
  names = list("matern1", "matern2"),
  theta_kappa = list(NULL, NULL),
  kappa = list(1, 1),
  rho = 0.5, theta = 1,
  mesh = NULL,
  alpha = 2,
  replicate = NULL,
  index_NA = NULL,
  noise = list(noise_normal(),noise_normal()),
  ...
) {
  # check input
  n_dim <- length(map[[1]])
  stopifnot(
    length(mesh) == n_dim,
    length(names) == n_dim,
    length(theta_kappa) == n_dim,
    length(kappa) == n_dim
  )
  stopifnot(n_dim == 2 || n_dim == 3)

  if (n_dim == 2) {
    Q  <- matrix(c(1, 0, 0, 1), 2, 2)
    Dl <- matrix(c(1, 0, 0, 1), 2, 2)
    D  <- Q %*% Dl
  }
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

# Y <- c(1:5, 1:5)
# group <- c(rep(1, 5), rep(2, 5))
# Y ~ Fixeff(~1+x, group=1) + matern_nd(..., group=group)
