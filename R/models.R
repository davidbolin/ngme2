# -------- ngme operators --------

#' ngme iid model specification
#'
#' @param map integer vector, time index for the AR(1) process
#' @param replicate replicate for the process
#'
#' @param ... extra arguments
#'
#' @return ngme_operator object
#' @export
#'
iid <- function(
  map,
  replicate = rep(1, length_map(map)),
  ...
) {
  if (inherits(map, "formula")) map <- model.matrix(map)[, -1]
  replicate <- as.integer(as.factor(replicate))
  K <- ngme_as_sparse(Matrix::Diagonal(length(map)))

  ngme_operator(
    map = map,
    mesh = NULL,
    n_rep = length(unique(replicate)),
    model = "iid",
    theta_K = double(0),
    K = K,
    h = rep(1, length(map)),
    A = K,
    symmetric = TRUE,
    zero_trace = FALSE
  )
}

#' ngme AR(1) model specification
#'
#' Generating C, G and A given index and replicate
#'
#' @param map integer vector, time index for the AR(1) process
#' @param mesh mesh for build the model
#' @param replicate replicate for the process
#' @param theta_K initial value for theta_K, if want to specify, can use e.g. ar1_a2th(0.5)
#'
#' @param ... extra arguments
#'
#' @return ngme_operator object
#' @export
#'
#' @examples
#' ar1(c(1:3, 1:3), replicate = c(1,1,1,2,2,2))
ar1 <- function(
  map,
  mesh      = NULL,
  replicate = rep(1, length_map(map)),
  theta_K   = 0,
  ...
) {
  if (inherits(map, "formula")) map <- model.matrix(map)[, -1]

  if (is.null(mesh)) mesh <- ngme_build_mesh(map)
  replicate <- as.integer(as.factor(replicate))

  stopifnot("The map should be integers."
      = all(map == round(map)))

  stopifnot("length of map and replicate should be the same." = length(map) == length(replicate))

  n <- mesh$n; nrep <- length(unique(replicate))

  h <- c(diff(mesh$loc), 1)
  G <- Matrix::Diagonal(n);
  C <- Matrix::sparseMatrix(j=1:(n-1), i=2:n, x=-1, dims=c(n,n))
  stopifnot("The mesh should be 1d and has gap 1." = all(h == 1))

  # for replicate
  if (nrep > 1) {
    G <- Matrix::kronecker(diag(nrep), G)
    C <- Matrix::kronecker(diag(nrep), C)
    h <- rep(h, nrep)
  }

  stopifnot("The length of theta_K should be 1." = length(theta_K) == 1)
  alpha <- ar1_th2a(theta_K)

  A <- if (!is.null(map)) INLA::inla.spde.make.A(mesh = mesh, loc = map, repl=replicate) else NULL

  ngme_operator(
    map = map,
    mesh = mesh,
    n_rep = nrep,
    model = "ar1",
    theta_K = theta_K,
    C = ngme_as_sparse(C),
    G = ngme_as_sparse(G),
    K = alpha * C + G,
    h = h,
    A = A,
    symmetric = FALSE,
    zero_trace = TRUE
  )
}

#' ngme random walk model of order 1
#'
#' generate K matrix of size (n-1) x n (non-cyclic case), where n is size of map
#'
#' @param map        numerical vector, covariates to build index for the process
#' @param replicate replicate for the process
#' @param cyclic  whether the mesh is circular, i.e. the first one is connected to the last
#'   if it is circular, we will treat the 1st location and the last location as neigbour, with distance of average distance.
#' @param mesh      mesh for the process, if not specified, will use inla.mesh.1d(loc = map)
#' @param ...       additional arguments
#'
#' @return a list
#' @export
#'
#' @examples
#' r1 <- rw1(1:7, cyclic = TRUE); r1$K
rw1 <- function(
  map,
  mesh      = NULL,
  replicate = rep(1, length_map(map)),
  cyclic    = FALSE,
  ...
) {
  if (inherits(map, "formula")) map <- model.matrix(map)[, -1]
  replicate <- as.integer(as.factor(replicate))
  stopifnot("length of map and replicate should be the same." = length(map) == length(replicate))

  if (is.null(mesh)) mesh <- ngme_build_mesh(map)
  n <- mesh$n; nrep <- length(unique(replicate))

  x <- map
  h <- diff(mesh$loc);
  n <- mesh$n
  if (!cyclic) {
    C <- Matrix::sparseMatrix(i = 1:n, j=1:n, x=1, dims=c(n,n))
    G <- Matrix::sparseMatrix(i = 2:n, j=1:(n-1), x=-1, dims=c(n,n))
    h <- c(0.01, h) # assume first point fixed to 0
  } else {
    stopifnot(length(x) >= 4)
    C <- Matrix::Diagonal(n)
    G <- Matrix::sparseMatrix(i = 1:n, j=c(2:n, 1), x=-1, dims=c(n,n))
    h <- c(h, mean(h))
  }

  # for replicate
  G <- Matrix::kronecker(diag(nrep), G)
  C <- Matrix::kronecker(diag(nrep), C)
  h <- rep(h, nrep)

  ngme_operator(
    map = map,
    mesh = mesh,
    n_rep = nrep,
    model = "rw1",
    theta_K = double(0),
    K = ngme_as_sparse(C + G),
    h = h,
    A = INLA::inla.spde.make.A(mesh = mesh, loc = map, repl=replicate),
    symmetric = FALSE,
    zero_trace = FALSE
  )
}

#' ngme random walk model of order 2
#'
#' generate K matrix of size (n-2) x n (non-cyclic case), where n is size of map
#'
#' @param map        numerical vector, covariates to build index for the process
#' @param replicate replicate for the process
#' @param cyclic  whether the mesh is circular, i.e. the first one is connected to the last
#'   if it is circular, we will treat the 1st location and the last location as neigbour, with distance of average distance.
#' @param mesh      mesh for the process, if not specified, will use inla.mesh.1d(loc = map)
#' @param ...       additional arguments
#'
#' @return a list
#' @export
#'
#' @examples
#' r2 <- rw2(1:7); r2$K
rw2 <- function(
  map,
  mesh      = map,
  replicate = rep(1, length_map(map)),
  cyclic    = FALSE,
  ...
) {
  if (inherits(map, "formula")) map <- model.matrix(map)[, -1]
  replicate <- as.integer(as.factor(replicate))
  if (is.null(mesh)) mesh <- ngme_build_mesh(map)

  stopifnot("length of map and replicate should be the same." = length(map) == length(replicate))

  x <- map
  n <- mesh$n; nrep <- length(unique(replicate))
  stopifnot("mesh too small" = n >= 3)

  h <- diff(mesh$loc)
  if (!cyclic) {
    C <- Matrix::sparseMatrix(i = 3:n, j=2:(n-1), x=-2, dims=c(n,n))
    G <- Matrix::sparseMatrix(i = c(1:n, 3:n), j=c(1:n, 1:(n-2)), x=1, dims=c(n,n))
    h <- c(0.01, h)
  } else {
    C <- Matrix::sparseMatrix(i = 1:n, j=c(2:n,1), x=-2, dims=c(n,n))
    G <- Matrix::sparseMatrix(i = rep(1:n,2), j=c(1:n, 3:n, 1, 2), x=1, dims=c(n,n))
    h <- c(h, mean(h))
  }

  # for replicate
  G <- Matrix::kronecker(diag(nrep), G)
  C <- Matrix::kronecker(diag(nrep), C)
  h <- rep(h, nrep)

  ngme_operator(
    map = map,
    mesh = mesh,
    n_rep = nrep,
    model = "rw2",
    theta_K = double(0),
    K = ngme_as_sparse(C + G),
    h = h,
    A = INLA::inla.spde.make.A(mesh = mesh, loc = map, repl=replicate),
    symmetric = FALSE,
    zero_trace = FALSE
  )
}

#' ngme Ornstein–Uhlenbeck process specification
#'
#' @param map integer vector, time index for the AR(1) process
#' @param replicate replicate for the process
#' @param mesh mesh for build the model
#' @param theta_K initial value for theta_K, kappa = exp(B_K %*% theta_K)
#' @param B_K bases for theta_K
#' @param ... extra arguments
#'
#' @return ngme_operator object
#' @export
ou <- function(
  map,
  replicate = rep(1, length_map(map)),
  mesh      = map,
  theta_K   = 0,
  B_K       = NULL,
  ...
) {
  if (inherits(map, "formula")) map <- model.matrix(map)[, -1]
  replicate <- as.integer(as.factor(replicate))
  n <- mesh$n; nrep <- length(unique(replicate))

  stopifnot("length of map and replicate should be the same." = length(map) == length(replicate))

  if (is.null(mesh)) mesh <- ngme_build_mesh(map)
  h <- diff(mesh$loc); h <- c(h, mean(h))

  if (is.null(B_K)) B_K <- matrix(1, nrow = length_map(map), ncol = 1)
  stopifnot("B_theta is a matrix" = is.matrix(B_K))
  stopifnot("ncol(B_K) == length(theta_K)"
    = ncol(B_K) == length(theta_K))

  G <- Matrix::bandSparse(n=n,m=n,k=c(-1,0),diagonals=cbind(-rep(1,n), rep(1,n)))
  C <- Ce <- Matrix::bandSparse(n=n,m=n,k=c(-1,0),diagonals=cbind(0.5*c(h[-1],0), 0.5*h))
  Ci = Matrix::sparseMatrix(i=1:n,j=1:n,x=1/h,dims = c(n,n))

  # for replicate
  if (nrep > 1) {
    G <- Matrix::kronecker(diag(nrep), G)
    C <- Matrix::kronecker(diag(nrep), C)
    h <- rep(h, nrep)
  }

  kappas <- exp(as.numeric(B_K %*% theta_K))
  K <- Matrix::Diagonal(x=kappas) %*% C + G

  ngme_operator(
    map         = map,
    mesh        = mesh,
    n_rep       = nrep,
    model       = "ou",
    B_K         = B_K,
    theta_K     = theta_K,
    C           = ngme_as_sparse(C),
    G           = ngme_as_sparse(G),
    K           = ngme_as_sparse(K),
    h           = h,
    A           = INLA::inla.spde.make.A(mesh = mesh, loc = map, repl=replicate),
    symmetric   = FALSE,
    zero_trace  = TRUE
  )
}

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
  if (inherits(map, "formula")) map <- model.matrix(map)[, -1]
  replicate <- as.integer(as.factor(replicate))
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
    n_rep = nrep,
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

#' ngme random effect model
#'
#' @param map integer vector, time index for the AR(1) process
#' @param replicate replicate for the process
#' @param alpha 2 or 4, SPDE smoothness parameter
#' @param theta_K initial value for theta_K
#'
#' @param ... extra arguments
#'
#' @return ngme_operator object
#' @export
re <- function(
  map,
  replicate = rep(1, length_map(map)),
  theta_K = NULL,
  ...
) {
  if (inherits(map, "formula")) map <- model.matrix(map)
    else map <- as.matrix(map)

  replicate <- as.integer(as.factor(replicate))
  nrep <- length(unique(replicate))

  B_K <- map
  n_reff <- ncol(B_K) # number of random effects
  n_theta_K <- sum(1:n_reff) # number of theta_K
  h <- rep(1, n_reff)

  # provide initial value for theta_K
  if (!is.null(theta_K)) {
    stopifnot(length(theta_K) == n_theta_K)
  } else {
    theta_K <- rep(0, n_theta_K)
  }

  # build K
  K <- diag(n_reff); diag(K) <- exp(theta_K[1:n_reff])
  if (n_reff > 1)
    K[lower.tri(K)] <- theta_K[(n_reff+1):n_theta_K]

  ngme_operator(
    map = map,
    mesh = NULL,
    n_rep = nrep,
    model = "re",
    theta_K = theta_K,
    K = ngme_as_sparse(K),
    h = h,
    A = ngme_as_sparse(B_K),
    symmetric = FALSE,
    zero_trace = FALSE
  )
}
