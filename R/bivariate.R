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
model_bv <- function(
  m1, m2,
  theta = 1, rho = 1,
  group = c(1, 2),
  replicate = NULL,
  noise = noise_normal(),
  share_param = FALSE,
  ...
) {
  stopifnot(
    inherits(m1, "ngme_model"),
    inherits(m2, "ngme_model"),
    length_map(m1$map) == length_map(m2$map),
    m1$W_size == m2$W_size,
    m1$V_size == m1$W_size, m2$V_size == m2$W_size    # square K
  )

  replicate <- rep(1, length_map(m1$map) * 2)
  A <- Matrix::bdiag(m1$A, m2$A)

  model <- ngme_model(
    model       = "bv",
    m1          = m1,
    m2          = m2,
    A           = A,
    A_pred      = NULL,
    W_size      = m1$W_size*2,
    V_size      = m1$V_size*2,
    theta_K     = c(theta, rho, m1$theta_K, m2$theta_K),
    rho         = rho,
    theta       = theta,
    h           = c(m1$h, m2$h),
    noise       = noise,
    # mesh        = mesh,
    map         = c(m1$map, m2$map),
    n_map       = m1$n_map * 2,
    replicate   = replicate,
    share_param = share_param,
    ...
  )
  model
}

