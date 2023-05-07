#' ngme Matern SPDE model specification
#'
#' Generating C, G and A given index and replicate
#'
#' @param map integer vector, time index for the AR(1) process
#' @param mesh mesh for build the SPDE model
#' @param replicate replicate for the process
#' @param alpha 2 or 4, SPDE smoothness parameter
#' @param theta_K initial value for theta_K, kappa = exp(B_K %*% theta_K)
#' @param B_K bases for theta_K
#'
#' @param ... extra arguments
#'
#' @return ngme_operator object
#' @export
matern <- function(
  map,
  mesh,
  replicate = rep(1, length_map(map)),
  alpha = 2,
  theta_K = 0,
  B_K = NULL,
  ...
) {
  stopifnot(alpha == 2 || alpha == 4)

  n <- mesh$n; nrep <- length(unique(replicate))

  if (is.null(B_K) && length(theta_K) == 1)
    B_K <- matrix(1, nrow = mesh$n, ncol = 1)
  else if (is.null(B_K) && length(theta_K) > 1)
    stop("Please provide B_K for non-stationary case.")

  d <- get_inla_mesh_dimension(mesh)
  if (d == 1) {
    fem <- INLA::inla.mesh.1d.fem(mesh)
    C <- fem$c1
    G <- fem$g1
    h <- Matrix::diag(fem$c0)
  } else {
    fem <- INLA::inla.mesh.fem(mesh, order = alpha)
    C <- fem$c0  # diag
    G <- fem$g1
    h <- Matrix::diag(fem$c0)
  }

  # for replicate
  if (nrep > 1) {
    G <- Matrix::kronecker(diag(nrep), G)
    C <- Matrix::kronecker(diag(nrep), C)
    h <- rep(h, nrep)
  }

  kappas <- as.numeric(exp(B_K %*% theta_K))
  if (length(theta_K) == 1) {
    # stationary
    K <- kappas[1]**alpha * C + G
  } else {
    # non-stationary
    K <- if (alpha == 2) diag(kappas) %*% C %*% diag(kappas)  + G
    else diag(kappas) %*% C %*% diag(kappas) %*% C %*% diag(kappas) + G
  }

  A <- if (!is.null(map)) INLA::inla.spde.make.A(mesh = mesh, loc = map, repl=replicate) else NULL

  ngme_operator(
    map = map,
    mesh = mesh,
    alpha = alpha,
    model = "matern",
    theta_K = theta_K,
    B_K = B_K,
    C = ngme_as_sparse(C),
    G = ngme_as_sparse(G),
    K = ngme_as_sparse(K),
    h = h,
    A = A,
    symmetric = TRUE,
    zero_trace = FALSE
  )
}

