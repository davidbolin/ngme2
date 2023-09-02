# -------- ngme operators --------

#' ngme iid model specification
#'
#' @param n integer, number of iid model
#' @param ... ignore
#'
#' @return ngme_operator object
#' @export
#'
iid <- function(
  n, ...
) {
  K <- ngme_as_sparse(Matrix::Diagonal(n))

  ngme_operator(
    mesh = fmesher::fm_mesh_1d(loc = 1:n),
    model = "iid",
    theta_K = double(0),
    K = K,
    h = rep(1, n),
    A = K,
    symmetric = TRUE,
    zero_trace = FALSE
  )
}

#' ngme AR(1) model specification
#'
#' @param mesh integer vector or inla.mesh.1d object, index to build the mesh
#' @param rho the correlation parameter (between -1 and 1)
#' @param ... ignore
#'
#' @return ngme_operator object
#' @export
#'
#' @examples
#' ar1(c(1:3, 1:3))
ar1 <- function(
  mesh, rho = 0, ...
) {
  stopifnot("rho should be between -1 and 1" = rho >= -1 && rho <= 1)

  mesh <- ngme_build_mesh(mesh)
  n <- mesh$n

  h <- c(diff(mesh$loc), 1)
  G <- Matrix::Diagonal(n); G[1, 1] = sqrt(1-rho**2)
  C <- Matrix::sparseMatrix(j=1:(n-1), i=2:n, x=-1, dims=c(n,n))
  stopifnot("The mesh should be 1d and has gap 1." = all(h == 1))

  theta_K <- ar1_a2th(rho)
  stopifnot("The length of rho(theta_K) should be 1." = length(theta_K) == 1)

  ngme_operator(
    mesh = mesh,
    model = "ar1",
    theta_K = theta_K,
    C = ngme_as_sparse(C),
    G = ngme_as_sparse(G),
    K = rho * C + G,
    h = h,
    symmetric = FALSE,
    zero_trace = TRUE
  )
}

#' ngme random walk model of order 1
#'
#' @param mesh numerical vector or inla.mesh.1d object, index to build the mesh
#' @param cyclic  whether the mesh is circular, i.e. the first one is connected to the last
#'   if it is circular, we will treat the 1st location and the last location as neigbour, with distance of average distance.
#' @param mesh mesh for build the model
#' @param ...       additional arguments
#'
#' @return a list
#' @export
#'
#' @examples
#' r1 <- rw1(1:7, cyclic = TRUE); r1$K
rw1 <- function(
  mesh,
  cyclic    = FALSE,
  ...
) {
  mesh <- ngme_build_mesh(mesh)

  n <- mesh$n
  h <- diff(mesh$loc);
  if (!cyclic) {
    C <- Matrix::sparseMatrix(i = 1:n, j=1:n, x=1, dims=c(n,n))
    G <- Matrix::sparseMatrix(i = 2:n, j=1:(n-1), x=-1, dims=c(n,n))
    h <- c(0.01, h) # assume first point fixed to 0
  } else {
    stopifnot("Too less data point" = length(h) >= 3)
    C <- Matrix::Diagonal(n)
    G <- Matrix::sparseMatrix(i = 1:n, j=c(2:n, 1), x=-1, dims=c(n,n))
    h <- c(h, mean(h))
  }

  ngme_operator(
    mesh = mesh,
    model = "rw1",
    theta_K = double(0),
    K = ngme_as_sparse(C + G),
    h = h,
    symmetric = FALSE,
    zero_trace = FALSE
  )
}

#' ngme random walk model of order 2
#'
#' generate K matrix of size (n-2) x n (non-cyclic case), where n is size of map
#'
#' @param mesh numerical vector or inla.mesh.1d object, index to build the mesh
#' @param mesh mesh for build the model
#' @param cyclic  whether the mesh is circular, i.e. the first one is connected to the last
#'   if it is circular, we will treat the 1st location and the last location as neigbour, with distance of average distance.
#' @param ...       additional arguments
#'
#' @return a list
#' @export
#'
#' @examples
#' r2 <- rw2(1:7); r2$K
rw2 <- function(
  mesh,
  cyclic = FALSE,
  ...
) {
  mesh <- ngme_build_mesh(mesh)

  stopifnot("Mesh should be inla.mesh.1d." = inherits(mesh, c("inla.mesh.1d")))
  n <- mesh$n
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

  ngme_operator(
    mesh = mesh,
    model = "rw2",
    theta_K = double(0),
    K = ngme_as_sparse(C + G),
    h = h,
    symmetric = FALSE,
    zero_trace = FALSE
  )
}

#' ngme Ornstein–Uhlenbeck process specification
#'
#' @param mesh numerical vector or inla.mesh.1d object, index to build the mesh
#' @param mesh mesh for build the model
#' @param theta_K initial value for theta_K, kappa = exp(B_K * theta_K)
#' @param B_K bases for theta_K
#' @param ... ignore
#'
#' @return ngme_operator object
#' @export
ou <- function(
  mesh,
  theta_K   = 0,
  B_K       = NULL,
  ...
) {
  mesh <- ngme_build_mesh(mesh)
  n <- mesh$n

  h <- diff(mesh$loc); h <- c(h, mean(h))

  if (is.null(B_K)) B_K <- matrix(1, nrow = length_map(mesh$loc), ncol = 1)
  stopifnot("B_theta is a matrix" = is.matrix(B_K))
  stopifnot("ncol(B_K) == length(theta_K)"
    = ncol(B_K) == length(theta_K))

  G <- Matrix::bandSparse(n=n,m=n,k=c(-1,0),diagonals=cbind(-rep(1,n), rep(1,n)))
  C <- Ce <- Matrix::bandSparse(n=n,m=n,k=c(-1,0),diagonals=cbind(0.5*c(h[-1],0), 0.5*h))
  Ci <- Matrix::sparseMatrix(i=1:n,j=1:n,x=1/h,dims = c(n,n))

  kappas <- exp(as.numeric(B_K %*% theta_K))
  K <- Matrix::Diagonal(x=kappas) %*% C + G

  ngme_operator(
    mesh        = mesh,
    model       = "ou",
    B_K         = B_K,
    theta_K     = theta_K,
    C           = ngme_as_sparse(C),
    G           = ngme_as_sparse(G),
    K           = ngme_as_sparse(K),
    h           = h,
    symmetric   = FALSE,
    zero_trace  = TRUE
  )
}

#' ngme Matern SPDE model specification
#'
#' @param mesh an inla.mesh.2d object, mesh for build the SPDE model
#' @param alpha 2 or 4, SPDE smoothness parameter
#' @param theta_K initial value for theta_K, kappa = exp(B_K * theta_K)
#' @param B_K bases for theta_K, ignore if use the stationary model
#' @param ... ignore
#'
#' @return ngme_operator object
#' @export
matern <- function(
  mesh,
  alpha = 2,
  theta_K = 0,
  B_K = NULL,
  ...
) {
  mesh <- ngme_build_mesh(mesh)
  stopifnot("alpha should be 2 or 4" = alpha == 2 || alpha == 4)

  n <- mesh$n
  if (is.null(B_K) && length(theta_K) == 1)
    B_K <- matrix(1, nrow = mesh$n, ncol = 1)
  else if (is.null(B_K) && length(theta_K) > 1)
    stop("Please provide B_K for non-stationary case.")

  d <- get_inla_mesh_dimension(mesh)
  if (d == 1) {
    fem <- fmesher::fm_mesh_1d.fem(mesh)
    C <- fem$c1
    G <- fem$g1
    h <- Matrix::diag(fem$c0)
  } else {
    # fem <- INLA::inla.mesh.fem(mesh, order = alpha)
    fem <- fmesher::fm_fem(mesh, order = alpha)
    C <- fem$c0  # diag
    G <- fem$g1
    h <- Matrix::diag(fem$c0)
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

  ngme_operator(
    mesh = mesh,
    alpha = alpha,
    model = "matern",
    theta_K = theta_K,
    B_K = B_K,
    C = ngme_as_sparse(C),
    G = ngme_as_sparse(G),
    K = ngme_as_sparse(K),
    h = h,
    symmetric = TRUE,
    zero_trace = FALSE
  )
}

#' ngme random effect model
#'
#' @param map numerical vector, covariates to build index for the process (can be formula, provided data)
#' @param theta_K initial value for theta_K (build covariance matrix)
#' @param ... ignore
#'
#' @return ngme_operator object
#' @export
re <- function(
  map,
  theta_K = NULL,
  ...
) {
  if (inherits(map, "formula")) map <- model.matrix(map)
    else map <- as.matrix(map)

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
  K <- Matrix::Matrix(K)

  ngme_operator(
    mesh = NULL,
    model = "re",
    theta_K = theta_K,
    K = ngme_as_sparse(K),
    h = h,
    B_K = B_K,
    symmetric = FALSE,
    zero_trace = FALSE
  )
}

# ----  For computing precision matrix of multivariate model

# p: dimension
# cor_mat: controls the correlation (only look at upper.tri part)
D_l <- function(p, rho) {
  stopifnot(
    "rho should be of length p(p-1)/2" =
      length(rho) == p*(p-1)/2
  )
  D_l <- diag(p)
  D_l[lower.tri(D_l)] <- rho
  # compute k(j)
  k <- double(p); k[1] <- 1
  for (j in 2:p) {
    # print(D_l[j, 1:(j-1)])
    k[j] <- sqrt(1+sum(D_l[j, 1:(j-1)] ^ 2))
  }
  D_l <- solve(D_l, diag(k))
  D_l
}

dependence_matrix <- function(p, cor_mat, theta=NULL, Q=NULL) {
  Q_2d <- function(theta) {
    Q <- matrix(0, nrow = 2, ncol = 2)
    Q[1, 1] <- cos(theta)
    Q[2, 2] <- cos(theta)
    Q[1, 2] <- -sin(theta)
    Q[2, 1] <- sin(theta)
    Q
  }
  stopifnot(
    p-round(p)==0, p > 1,
    "Please provide theta (p <= 3) or Q matrix, see ?precision_matrix_multivariate"
      = !is.null(theta) | !is.null(Q)
  )

  # compute D_l
  D_l <- D_l(p, rho)
  # compute Q
  if (p == 2) {
    stopifnot("Length of theta should be 1 for p=2 case"
      = length(theta) == 1)
    Q <- Q_2d(theta)
  } else if (p == 3) {
    stopifnot("Length of theta should be 3 for p=3 case"
      = length(theta) == 3)

    Q_3x <- Matrix::bdiag(Q_2d(theta[1]), 1)
    Q_3z <- Matrix::bdiag(1, Q_2d(theta[3]))
    Q_3y <- diag(3)
    Q_3y[c(1, 3, 7, 9)] <- Q_2d(theta[2])
    Q <- Q_3x %*% Q_3y %*% Q_3z
  } else {
    if (is.null(Q)) stop("Please provide Q (p*p) for p > 3 case.")
  }

  Q %*% D_l
}

#' Compute the precision matrix for multivariate model
#'
#' @param p dimension, should be integer and greater than 1
#' @param operator_list a list of ngme_operator object (length should be p)
#' @param rho vector with the p(p-1)/2 correlation parameters rho_11, rho_21, rho_22, ... rho_p1, rho_p2, ... rho_p(p-1)
#' @param theta parameter for Q matrix (length of 1 when p=2, length of 3 when p=3)
#' @param Q orthogonal matrix of dim p*p (provide when p > 3)
#' @param scale A vector of length p with constants to multiply each operator matrix with
#'
#' @return the precision matrix of the multivariate model
#' @details The general model is defined as $D diag(L_1, ..., L_p) x = M$. D is the dependence matrix, it is paramterized by $D = Q(theta) * D_l(cor_mat)$, where $Q$ is the orthogonal matrix, and $D_l$ is matrix controls the cross-correlation.
#' See the section 2.2 of Bolin and Wallin (2020) for exact parameterization of Dependence matrix.
#' @references
#' Bolin, D. and Wallin, J. (2020), Multivariate type G Matérn stochastic partial differential equation random fields. J. R. Stat. Soc. B, 82: 215-239. https://doi.org/10.1111/rssb.12351
#' @export
#' @examples
#' rho_mat <- matrix(0, nrow = 3, ncol = 3)
#' rho_mat[1, 2] <- 0.4; rho_mat[1, 3] <- -0.5; rho_mat[2,3] <- 0.8
#' operator_list <- list(ar1(1:5, rho=0.4), ar1(1:5, rho=0.5), ar1(1:5, rho=0.6))
#' precision_matrix_multivariate(3, operator_list, rho_mat, theta=c(1,2,3))
precision_matrix_multivariate <- function(p,
                                          operator_list,
                                          rho,
                                          theta=NULL,
                                          Q=NULL,
                                          scale = NULL) {
  stopifnot("Please provide a list of p models" =  length(operator_list) == p)
  stopifnot(
    "rho should be of length p(p-1)/2" =
      length(rho) == p*(p-1)/2
  )

  D <- dependence_matrix(p, rho, theta, Q)
  bigD <- kronecker(D, Matrix::Diagonal(length(operator_list[[1]]$h)))

  if(is.null(scale)){
    scale <- rep(1,p)
  }
  K_list <- list()
  for(i in 1:p){
    K_list[[i]] <- scale[i]*operator_list[[i]]$K
  }
  K <- bigD %*% Matrix::bdiag(K_list)

  # mass lumping version
  Cinv <- rep(1 / operator_list[[1]]$h, p)
  Matrix::t(K) %*% Matrix::Diagonal(x = Cinv) %*% K
}


#' Compute the precision matrix for multivariate spde Matern model
#'
#' @param p dimension, should be integer and greater than 1
#' @param mesh an inla.mesh.2d object, mesh for build the SPDE model
#' @param alpha 2 or 4, SPDE smoothness parameter
#' @param theta_K_list a list (length is p) of theta_K
#' @param B_K_list a list (length is p) of B_K (non-stationary case)
#' @param variance_list If provided, it should be a vector of length p, where the
#' kth element corresponds to a desired variance of the kth field. The kth operator
#' is then scaled by a constant c so that this variance is achieved in the stationary case
#' (default no scaling)
#' @param rho vector with the p(p-1)/2 correlation parameters rho_11, rho_21, rho_22, ... rho_p1, rho_p2, ... rho_p(p-1)
#' @param theta parameter for Q matrix (length of 1 when p=2, length of 3 when p=3)
#' @param Q orthogonal matrix of dim p*p (provide when p > 3)
#'
#' @return the precision matrix of the multivariate model
#' @details The general model is defined as $D diag(L_1, ..., L_p) x = M$. D is the dependence matrix, it is paramterized by $D = Q(theta) * D_l(cor_mat)$, where $Q$ is the orthogonal matrix, and $D_l$ is matrix controls the cross-correlation.
#' See the section 2.2 of Bolin and Wallin (2020) for exact parameterization of Dependence matrix.
#' @references
#' Bolin, D. and Wallin, J. (2020), Multivariate type G Matérn stochastic partial differential equation random fields. J. R. Stat. Soc. B, 82: 215-239. https://doi.org/10.1111/rssb.12351
#' @export
#' @examples
#' library(INLA)
#' library(fields)
#' # Define mesh
#' x <- seq(from=0,to=1,length.out = 40)
#' mesh <- inla.mesh.create(lattice = inla.mesh.lattice(x,x), extend = FALSE)
#' # Set parameters
#' rho <- c(-0.5, 0.5,-0.25) #correlation parameters
#' log_kappa <- list(2,2,2) #log(kappa)
#' variances <- list(1,1,1) #set marginal variances to 1
#' alpha <- list(2,2,2) #smoothness parameters
#' # Compute precision
#' Q <- precision_matrix_multivariate_spde(p, mesh = mesh, rho = rho,
#'                                        alpha = alpha, theta_K_list = log_kappa,
#'                                        variance_list = variances)
#'
# Plot the cross covariances
#' A <- as.vector(inla.spde.make.A(mesh,loc = matrix(c(0.5,0.5),1,2)))
#' Sigma <- as.vector(solve(Q,c(A,rep(0,2*mesh$n))))
#' r11 <- Sigma[1:mesh$n]
#' r12 <- Sigma[(mesh$n+1):(2*mesh$n)]
#' r13 <- Sigma[(2*mesh$n+1):(3*mesh$n)]
#' Sigma <- as.vector(solve(Q,c(rep(0,mesh$n),A,rep(0,mesh$n))))
#' r22 <- Sigma[(mesh$n+1):(2*mesh$n)]
#' r23 <- Sigma[(2*mesh$n+1):(3*mesh$n)]
#' Sigma <- as.vector(solve(Q,v <- c(rep(0,2*mesh$n),A)))
#' r33 <- Sigma[(2*mesh$n+1):(3*mesh$n)]
#'
#' proj <- inla.mesh.projector(mesh)
#'
#' par(mfrow=c(3,3))
#' image.plot(inla.mesh.project(proj,r11), main = "Cov(X_1(s0),X_1(s)")
#' plot.new()
#' plot.new()
#' image.plot(inla.mesh.project(proj,r12), main = "Cov(X_1(s0),X_2(s)")
#' image.plot(inla.mesh.project(proj,r22), main = "Cov(X_2(s0),X_2(s)")
#' plot.new()
#' image.plot(inla.mesh.project(proj,r13), main = "Cov(X_1(s0),X_3(s)")
#' image.plot(inla.mesh.project(proj,r23), main = "Cov(X_2(s0),X_3(s)")
#' image.plot(inla.mesh.project(proj,r33), main = "Cov(X_3(s0),X_3(s)")

precision_matrix_multivariate_spde <- function(
  p,
  mesh,
  rho,
  alpha_list = NULL,
  theta_K_list = NULL,
  variance_list = NULL,
  B_K_list = NULL,
  theta = NULL,
  Q = NULL
) {
  if (is.null(B_K_list)) {
    B_K_list <- lapply(1:p, function(i) NULL)
  } else {
    stopifnot("Please provide a list of B_K matrices for each matern model"
      = length(B_K_list) == p)
  }

  if (is.null(theta_K_list)) {
    theta_K_list <- lapply(1:p, function(i) 0)
  } else {
    stopifnot("Please provide a list of p parameters"
      = length(theta_K_list) == p)
  }

  if (is.null(alpha_list)) {
    alpha_list <- lapply(1:p, function(i) 2)
  } else {
    stopifnot("Please provide a list of p alpha parameters"
              = length(alpha_list) == p)
  }

  if (is.null(variance_list)) {
    variance_list <- lapply(1:p, function(i) NULL)
    c <- NULL
  } else {
    stopifnot("Please provide a list of p variance parameters"
              = length(variance_list) == p)
    if(mesh$manifold %in% c("R2", "S2") ) {
      d = 2
    } else if(mesh$manifold == "R1") {
      d = 1
    }
    scale <- rep(0,p)
    for(i in 1:p){
      scale[i] <- sqrt(gamma(alpha_list[[i]] - d/2)/(variance_list[[i]]*(4*pi)^(d/2)*exp(theta_K_list[[i]][1])^(2*(alpha_list[[i]]-d/2))*gamma(alpha_list[[i]])))
    }
  }

  operator_list <- NULL
  for (i in 1:p) {
    operator_list[[i]] <- matern(mesh, alpha_list[[i]], theta_K_list[[i]], B_K_list[[i]])
  }
  if(is.null(theta)) {
    theta <- rep(0,(p-1)*p/2)
  }
  precision_matrix_multivariate(p, operator_list, rho, theta, Q, scale = scale)
}
