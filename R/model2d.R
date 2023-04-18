# define the interface for a bivariate latent model

#' Create a bivariate SPDE model
#'
#' @param map list of length 2 or 3, each element is numeric vector (1d) or matrix of column 2 (2d)
#' @param replicate replicate for the process
#' @param alpha     2 or 4, SPDE smoothness paramete
#' @param mesh      mesh argument
#' @param index_NA Logical vector, same as is.na(response var.)
#' @param kappa     parameterization for kappa^2 C + G, only for stationary
#' @param theta_kappa parameterization for non-stationary
#' @param B_kappa   bases for kappa
#' @param noise     1. string: type of model, 2. ngme.noise object
#' @param ... extra arguments in f()
#'
#' @return a list (n, C (diagonal), G, B.kappa) for constructing operator
#' @export
bv <- function(
  m1, m2,
  rho = 0.5, theta = 1,
  group = c(1, 2),
  ...
) {
  stopifnot(
    inherits(m1, "ngme_model"),
    inherits(m2, "ngme_model"),
    length_map(m1$map) == length_map(m2$map)
  )

  model <- ngme_model(
    model       = "bv",
    A           = A,
    A_pred      = A_pred,
    W_size      = nrow(C),
    V_size      = ncol(C),
    theta_K     = theta_kappa,
    alpha       = alpha,
    B_kappa     = B_kappa,
    C           = ngme_as_sparse(C),
    G           = ngme_as_sparse(G),
    rho         = rho,
    theta       = theta,
    h           = h,
    noise       = noise,
    mesh        = mesh,
    map         = c(m1$map, m2$map),
    n_map       = m1$n_map * 2,
    replicate   = replicate,
    ...
  )
  model
}

