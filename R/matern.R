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
  # build C, G, K, A
  n <- mesh$n; nrep <- length(unique(replicate))
  stopifnot(alpha == 2 || alpha == 4)

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

  kappas <- as.numeric(exp(B_K %*% theta_K))
  if (length(theta_K) == 1) {
    # stationary
    K <- kappas[1]**alpha * C + G
  }
  else {
    # non-stationary
    K <- if (alpha == 2) diag(kappas) %*% C %*% diag(kappas)  + G
    else diag(kappas) %*% C %*% diag(kappas) %*% C %*% diag(kappas) + G
  }

  ngme_operator(
    alpha = alpha,
    model = "matern",
    theta_K = theta_K,
    B_K = B_K,
    C = ngme_as_sparse(C),
    G = ngme_as_sparse(G),
    K = ngme_as_sparse(K),
    h = h,
    A = INLA::inla.spde.make.A(mesh = mesh, loc = map),
    symmetric = TRUE,
    zero_trace = FALSE
  )
}


#' Create a Matern SPDE model
#'
#' @param map       numeric vector (1d) or matrix of column 2 (2d),
#'     location to make index
#'     keep aligned with Y if make prediction with index_NA!!!
#' @param replicate replicate for the process
#' @param alpha     2 or 4, SPDE smoothness parameter
#' @param mesh      mesh argument
#' @param index_NA Logical vector, same as is.na(response var.)
#' @param theta_kappa parameterization for non-stationary
#' @param B_kappa   bases for kappa
#' @param d         indicating the dimension of mesh (together with fem.mesh.matrices)
#' @param noise     1. string: type of model, 2. ngme.noise object
#' @param ... extra arguments in f()
#'
#' @return a list (n, C (diagonal), G, B.kappa) for constructing operator
#' @export
model_matern <- function(
  map,
  replicate   = NULL,
  alpha       = 2,
  theta_kappa = 0,
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
  # if (is.null(loc)) {
  #   if (inherits(mesh, "inla.mesh.1d")) loc <- mesh$loc
  #   if (inherits(mesh, "inla.mesh")) loc <- as.matrix(mesh$loc[, 1:2])
  # }
  loc <- map

  if (is.numeric(mesh)) # mesh is 1d vector
    mesh <- INLA::inla.mesh.1d(loc = mesh)

  n_loc <- if(is.null(dim(loc))) length(loc) else nrow(loc)
  # loc <- eval(substitute(loc), envir = data, enclos = parent.frame())

  if (is.null(mesh)) stop("Please provide mesh!")

  if (is.null(replicate)) replicate <- rep(1, n_loc)

  # deal with coords
  if (alpha - round(alpha) != 0) {
    stop("alpha should be integer, now only 2 or 4")
  }

  stopifnot(alpha == 2 || alpha == 4)

  if (is.null(B_kappa) && length(theta_kappa) == 1)
    B_kappa <- matrix(1, nrow = mesh$n, ncol = 1)
  else if (is.null(B_kappa) && length(theta_kappa) > 1)
    stop("Please provide B_kappa for non-stationary case.")

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

  tmp <- ngme_make_A(
    mesh = mesh,
    map = loc,
    n_map = n_loc,
    idx_NA = index_NA,
    replicate = replicate
  )
  A <- tmp$A; A_pred <- tmp$A_pred

  nrep <- ncol(A) / ncol(C)

  if (noise$n_noise == 1) noise <- update_noise(noise, n = mesh$n)
  noise$h <- h
  if (is.null(noise$V)) noise$V <- h
  kappas <- drop(exp(B_kappa %*% theta_kappa))

  model <- ngme_model(
    model       = "matern",
    A           = A,
    A_pred      = A_pred,
    W_size      = nrow(C),
    V_size      = ncol(C),
    theta_K     = theta_kappa,
    alpha       = alpha,
    B_K         = B_kappa,
    C           = ngme_as_sparse(C),
    G           = ngme_as_sparse(G),
    K           = ngme_as_sparse(diag(kappas) %*% C %*% diag(kappas)  + G),
    h           = h,
    noise       = noise,
    mesh        = mesh,
    map         = loc,
    n_map       = n_loc,
    replicate   = replicate,
    n_rep       = nrep,
    # group       = NULL,
    ...
  )
  model
}