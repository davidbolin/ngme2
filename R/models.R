# -------- ngme operators --------

#' ngme iid model specification
#'
#' @param map index vector
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
#' @param rho the correlation parameter (between -1 and 1)
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
  mesh      = INLA::inla.mesh.1d(map),
  replicate = rep(1, length_map(map)),
  rho       = 0,
  ...
) {
  # check map
  if (inherits(map, "formula")) map <- model.matrix(map)[, -1]
  stopifnot("The map should be integers." = all(map == round(map)))

  # check replicate
  replicate <- as.integer(as.factor(replicate))
  stopifnot("length of map and replicate should be the same." = length(map) == length(replicate))

  # check mesh
  stopifnot("Mesh should be inla.mesh.1d." = inherits(mesh, c("inla.mesh.1d")))

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

  theta_K <- ar1_a2th(rho)
  stopifnot("The length of rho(theta_K) should be 1." = length(theta_K) == 1)

  A <- if (!is.null(map)) INLA::inla.spde.make.A(mesh = mesh, loc = map, repl=replicate) else NULL

  ngme_operator(
    map = map,
    mesh = mesh,
    n_rep = nrep,
    model = "ar1",
    theta_K = theta_K,
    C = ngme_as_sparse(C),
    G = ngme_as_sparse(G),
    K = rho * C + G,
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
  mesh      = INLA::inla.mesh.1d(map),
  replicate = rep(1, length_map(map)),
  cyclic    = FALSE,
  ...
) {
  if (inherits(map, "formula")) map <- model.matrix(map)[, -1]
  replicate <- as.integer(as.factor(replicate))
  stopifnot("length of map and replicate should be the same." = length(map) == length(replicate))

  stopifnot("Mesh should be inla.mesh.1d." = inherits(mesh, c("inla.mesh.1d")))
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
#' @param map  numerical vector, covariates to build index for the process
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
  mesh      = INLA::inla.mesh.1d(map),
  replicate = rep(1, length_map(map)),
  cyclic    = FALSE,
  ...
) {
  if (inherits(map, "formula")) map <- model.matrix(map)[, -1]
  replicate <- as.integer(as.factor(replicate))
  if (is.null(mesh)) mesh <- ngme_build_mesh(map)

  stopifnot("length of map and replicate should be the same." = length(map) == length(replicate))

  x <- map
  stopifnot("Mesh should be inla.mesh.1d." = inherits(mesh, c("inla.mesh.1d")))
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
#' @param map numerical vector, covariates to build index for the process
#' @param replicate replicate for the process
#' @param mesh mesh for build the model
#' @param theta_K initial value for theta_K, kappa = exp(B_K * theta_K)
#' @param B_K bases for theta_K
#' @param ... extra arguments
#'
#' @return ngme_operator object
#' @export
ou <- function(
  map,
  mesh      = INLA::inla.mesh.1d(map),
  replicate = rep(1, length_map(map)),
  theta_K   = 0,
  B_K       = NULL,
  ...
) {
  if (inherits(map, "formula")) map <- model.matrix(map)[, -1]
  replicate <- as.integer(as.factor(replicate))

  stopifnot("Mesh should be inla.mesh.1d." = inherits(mesh, c("inla.mesh.1d")))
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
#' @param map  numerical vector, covariates to build index for the process
#' @param mesh mesh for build the SPDE model
#' @param replicate replicate for the process
#' @param alpha 2 or 4, SPDE smoothness parameter
#' @param theta_K initial value for theta_K, kappa = exp(B_K * theta_K)
#' @param B_K bases for theta_K
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
  A <- INLA::inla.spde.make.A(mesh = mesh, loc = map, repl=replicate)

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
#' @param map numerical vector, covariates to build index for the process (can be formula, provided data)
#' @param replicate replicate for the process
#' @param alpha 2 or 4, SPDE smoothness parameter
#' @param theta_K initial value for theta_K (build covariance matrix)
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

# ----  For computing precision matrix of multivariate model
# p: dimension
# cor_mat: controls the correlation (only look at upper.tri part)
D_l <- function(p, cor_mat) {
  stopifnot(
    "cor_mat should be of dim p*p" =
      is.matrix(cor_mat) &&
      ncol(cor_mat) == p &&
      nrow(cor_mat) == p
  )
  D_l <- diag(p)
  D_l[upper.tri(D_l)] <- cor_mat[upper.tri(cor_mat)]
  D_l <- t(D_l)
  # compute k(j)
  k <- double(p); k[1] <- 1
  for (j in 2:p) {
    # print(D_l[j, 1:(j-1)])
    k[j] <- sqrt(1 + sum(D_l[j, 1:(j-1)] ^ 2))
  }
  D_l <- solve(D_l, diag(k))
  D_l
}

dependence_matrix <- function(p, cor_mat, zeta=NULL, Q=NULL) {
  Q_2d <- function(zeta) {
    Q <- matrix(0, nrow = 2, ncol = 2)
    Q[1, 1] <- cos(zeta)
    Q[2, 2] <- cos(zeta)
    Q[1, 2] <- -sin(zeta)
    Q[2, 1] <- sin(zeta)
    Q
  }
  stopifnot(
    p-round(p)==0, p > 1,
    "Please provide zeta (p <= 3) or Q matrix, see ?precision_matrix_multivariate"
      = !is.null(zeta) | !is.null(Q)
  )

  # compute D_l
  D_l <- D_l(p, cor_mat)
  # compute Q
  if (p == 2) {
    stopifnot("Length of zeta should be 1 for p=2 case"
      = length(zeta) == 1)
    Q <- Q_2d(zeta)
  } else if (p == 3) {
    stopifnot("Length of zeta should be 3 for p=3 case"
      = length(zeta) == 3)

    Q_3x <- Matrix::bdiag(Q_2d(zeta[1]), 1)
    Q_3z <- Matrix::bdiag(1, Q_2d(zeta[3]))
    Q_3y <- matrix(0, nrow = 3, ncol = 3)
    Q_3y[c(1, 3, 7, 9)] <- Q_2d(zeta[2])
    Q <- Q_3x %*% Q_3y %*% Q_3z
  } else {
    if (is.null(Q)) stop("Please provide Q (p*p) for p > 3 case.")
  }

  Q %*% D_l
}

#' compute the precision matrix for multivariate model
#'
#' @param p dimension, should be integer and greater than 1
#' @param model an ngme_operator object
#' @param cor_mat matrix of dim p*p, controls the correlation (only look at upper.tri part)
#' @param zeta parameter for Q matrix (length of 1 when p=2, length of 3 when p=3)
#' @param Q orthogonal matrix of dim p*p (provide when p > 3)
#'
#' @return the precision matrix of the multivariate model
#' @details The general model is defined as $D diag(L_1, ..., L_p) x = M$. D is the dependence matrix, it is paramterized by $D = Q(zeta) * D_l(cor_mat)$, where $Q$ is the orthogonal matrix, and $D_l$ is matrix controls the cross-correlation.
#' See the section 2.2 of Bolin and Wallin (2020) for exact parameterization of Dependence matrix.
#' @references
#' Bolin, D. and Wallin, J. (2020), Multivariate type G Matérn stochastic partial differential equation random fields. J. R. Stat. Soc. B, 82: 215-239. https://doi.org/10.1111/rssb.12351
#' @export
#' @examples
#' rho_mat <- matrix(0, nrow = 3, ncol = 3)
#' rho_mat[1, 2] <- 0.4; rho_mat[1, 3] <- -0.5; rho_mat[2,3] <- 0.8
#' precision_matrix_multivariate(3, ar1(1:5), rho_mat, zeta=c(1,2,3))
precision_matrix_multivariate <- function(p, model, cor_mat, zeta=NULL, Q=NULL) {
  D <- dependence_matrix(p, cor_mat, zeta, Q)
  bigD <- kronecker(D, Matrix::Diagonal(nrow(model$K)))
  K <- bigD %*% Matrix::bdiag(rep(list(model$K), p))
  Matrix::t(K) %*% K
}

