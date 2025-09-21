# -------- ngme operators --------

#' ngme iid model specification
#'
#' @param map integer or factor, index to build the mesh
#' @param ... ignore
#'
#' @return ngme_operator object
#' @export
#'
iid <- function(
  map, ...
) {
  n <- length(levels(as.factor(map)))
  K <- ngme_as_sparse(Matrix::Diagonal(n))

  ngme_operator(
    matrices = list(K),
    mesh = fmesher::fm_mesh_1d(loc = 1:n),
    model = "iid",
    theta_K = double(0),
    update_K = function(theta_K) {K},
    K = K,
    h = rep(1, n),
    A = K,
    symmetric = TRUE,
    zero_trace = FALSE,
    param_name = character(0),
    generic_type = "generic"
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
  G <- ngme_as_sparse(G)
  C <- ngme_as_sparse(C)
  stopifnot("The mesh should be 1d and has gap 1." = all(h == 1))

  theta_K <- ar1_a2th(rho)
  stopifnot("The length of rho(theta_K) should be 1." = length(theta_K) == 1)

  update_K <- function(theta_K) {ar1_th2a(theta_K) * C + G}

  # Create a generic model internally
  g <- name2fun("tanh", inv=TRUE)

  ngme_operator(
    mesh = mesh,
    model = "ar1",
    theta_K = c(rho=g(rho)),
    trans = c(rho="tanh"),
    matrices = list(C, G),
    h = h,
    C = C,
    G = G,
    update_K = update_K,
    K = ngme_as_sparse(update_K(theta_K)),
    symmetric = FALSE,
    zero_trace = FALSE,
    param_name = "rho",
    param_trans = list(ar1_th2a),
    # using the generic structure to build the operator
    generic_type = "generic"
  )
}

#' Helper function to compute AR(p) autocovariance
#' 
#' @param rho vector of AR coefficients
#' @param p AR order
#' @return vector of autocovariances from lag 0 to lag p-1
#' @keywords internal
ARcov <- function(rho, p) {
  # Solve Yule-Walker equations for autocovariances
  # gamma(0), gamma(1), ..., gamma(p-1)
  
  # Create the coefficient matrix A
  A <- matrix(0, p, p)
  for (i in 1:p) {
    for (j in 1:p) {
      lag <- abs(i - j)
      if (lag == 0) {
        A[i, j] <- 1
      } else {
        A[i, j] <- -rho[lag]
      }
    }
  }
  
  # Right hand side vector b
  b <- c(1, rep(0, p-1))
  
  # Solve A * gamma = b
  gamma <- solve(A, b)
  
  return(gamma)
}

#' ngme AR(p) model specification
#'
#' @param mesh integer vector or inla.mesh.1d object, index to build the mesh
#' @param rho vector of AR coefficients of length p (should satisfy stationarity conditions)
#' @param order integer, the AR order p. If provided and rho is not specified, rho will be initialized as rep(0, p)
#' @param ... ignore
#'
#' @return ngme_operator object
#' @export
#'
#' @examples
#' # AR(2) model with specified coefficients
#' ar(c(1:10), rho = c(0.5, -0.3))
#' # AR(3) model with specified coefficients 
#' ar(c(1:10), rho = c(0.4, 0.2, -0.1))
#' # AR(2) model with default coefficients (zeros)
#' ar(c(1:10), order = 2)
#' # AR(3) model with specified order and coefficients
#' ar(c(1:10), order = 3, rho = c(0.4, 0.2, -0.1))
ar <- function(
  mesh, rho = NULL, order = NULL, ...
) {
  # Determine p and rho based on input parameters
  if (is.null(order) && is.null(rho)) {
    stop("Either 'rho' or 'order' must be provided")
  } else if (is.null(order) && !is.null(rho)) {
    # Infer order from rho length
    p <- length(rho)
  } else if (!is.null(order) && is.null(rho)) {
    # Use specified order and initialize rho as zeros
    p <- order
    rho <- rep(0, p)
  } else {
    # Both order and rho provided - check consistency
    p <- order
    stopifnot("Length of rho must match the specified order" = length(rho) == p)
  }
  
  stopifnot("p should be >= 1" = p >= 1)
  
  # Check stationarity conditions for AR(p)
  # For now, we'll do a simple check for AR(2)
  # if (p == 2) {
  #   stopifnot("AR(2) stationarity conditions violated" = 
  #             rho[1] + rho[2] < 1 && 
  #             rho[2] - rho[1] < 1 && 
  #             abs(rho[2]) < 1)
  # } else if (p > 2) {
  #   # For higher order AR, we check if all roots of characteristic polynomial 
  #   # are outside unit circle (simplified check)
  #   char_poly <- c(1, -rho)
  #   roots <- polyroot(char_poly)
  #   if (any(abs(roots) <= 1.001)) {  # small tolerance for numerical errors
  #     warning("AR(p) model may not be stationary")
  #   }
  # }

  mesh <- ngme_build_mesh(mesh)
  n <- mesh$n
  h <- c(diff(mesh$loc), 1)
  
  stopifnot("The mesh should be 1d and has gap 1." = all(h[1:(n-1)] == 1))
  stopifnot("n should be >= p" = n >= p)

  # Construct C1, ..., Cp and G matrices
  Cs <- vector("list", p)
  for (j in 1:p) {
    # Cj: for lag j, put -1 at (t, t-j) for t = (j+1):n
    # The first p rows should be all 0, so we only fill from row (p+1) to n
    Cs[[j]] <- Matrix::sparseMatrix(
      i = (p+1):n,
      j = (p+1-j):(n-j),
      x = rep(-1, n-p),
      dims = c(n, n)
    )
  }
  # G: identity matrix
  G <- Matrix::Diagonal(n)
  G <- ngme_as_sparse(G)

  # Create transformation functions
  rho_trans <- rep(list(identity), p)
  names(rho_trans) <- paste0("rho", 1:p)
  trans <- setNames(rep("tanh", p), paste0("rho", 1:p))

  # Initial theta_K values (no transformation for now)
  theta_K <- rho
  names(theta_K) <- paste0("rho", 1:p)

  # update_K: sum_j rho[j] * Cj + G
  update_K <- function(theta_K) {
    K <- G
    for (j in 1:p) {
      K <- K + theta_K[j] * Cs[[j]]
    }
    K
  }

  K <- ngme_as_sparse(update_K(theta_K))

  ngme_operator(
    mesh = mesh,
    model = paste0("ar", p),
    theta_K = theta_K,
    trans = trans,
    matrices = c(Cs, list(G)),
    h = h,
    update_K = update_K,
    K = K,
    symmetric = FALSE,
    zero_trace = FALSE,
    param_name = paste0("rho", 1:p),
    param_trans = rho_trans,
    generic_type = "generic"
  )
}

#' Random Walk Model of Order 1 (RW1)
#'
#' Constructs a first-order random walk model for spatial or temporal processes.
#' The RW1 model assumes that first-order differences \eqn{\Delta W_i = W_i - W_{i-1}}
#' are independent and identically distributed Gaussian variables.
#'
#' @details
#' The RW1 model is defined by the precision matrix \eqn{K} that penalizes 
#' first-order differences. The model has different structures depending on the
#' constraints:
#' 
#' **Non-cyclic, constrained (default)**: The first row of \eqn{K} enforces the 
#' constraint \eqn{\sum_{i=1}^n h_i W_i = 0}, where \eqn{h_i} are the mesh weights.
#' 
#' **Non-cyclic, unconstrained**: The first element is fixed at 0 (\eqn{W_1 = 0}).
#' 
#' **Cyclic**: Treats the domain as circular, connecting the first and last locations
#' as neighbors. The constraint \eqn{\sum_{i=1}^n h_i W_i = 0} is always enforced.
#' The precision matrix is expanded to size \eqn{(n+1) \times (n+1)} to handle
#' the constraint properly.
#'
#' @param mesh numerical vector or inla.mesh.1d object, locations to build the mesh.
#'   For numerical vectors, assumes equally spaced locations.
#' @param cyclic logical, whether the mesh is circular. If TRUE, the first and last
#'   locations are treated as neighbors. Cannot be FALSE when constr = TRUE.
#' @param constr logical, whether to enforce the sum-to-zero constraint 
#'   \eqn{\sum_{i=1}^n h_i W_i = 0}. If FALSE, fixes the first element \eqn{W_1 = 0}.
#'   Must be TRUE for cyclic models.
#' @param ... additional arguments (currently unused)
#'
#' @return An `ngme_operator` object containing the precision matrix and related
#'   components for the RW1 model.
#' @export
#'
#' @examples
#' # Non-cyclic constrained RW1 (default)
#' rw1_default <- rw1(1:5)
#' print(rw1_default$K)
#' 
#' # Non-cyclic unconstrained RW1 (fixes first element to 0)
#' rw1_fixed <- rw1(1:5, constr = FALSE)
#' print(rw1_fixed$K)
#' 
#' # Cyclic RW1 (connects first and last locations)
#' rw1_cyclic <- rw1(1:5, cyclic = TRUE)
#' print(rw1_cyclic$K)
#' 
#' # Using with unequally spaced locations
#' locations <- c(0, 1, 3, 6, 10)
#' rw1_unequal <- rw1(locations)
#' print(rw1_unequal$K)
rw1 <- function(
  mesh,
  cyclic    = FALSE,
  constr    = TRUE,
  ...
) {
  mesh <- ngme_build_mesh(mesh)

  n <- mesh$n
  h_left <- c(0, diff(mesh$loc))
  h_right <- c(diff(mesh$loc), 0)
  h <- 0.5 * (h_left + h_right)

  if (!cyclic) {
    C <- Matrix::sparseMatrix(i = 1:n, j=1:n, x=1, dims=c(n,n))
    G <- Matrix::sparseMatrix(i = 2:n, j=1:(n-1), x=-1, dims=c(n,n))
    K_mat <- C + G
    if (constr) K_mat[1, ] <- h else K_mat[1, 1] <- 1
    K = ngme_as_sparse(K_mat)
  } else {
    stopifnot("constraint cannot be false for cyclic case" = constr)
    stopifnot("Too less data point" = length(h) >= 3)
    C <- Matrix::Diagonal(n)
    G <- Matrix::sparseMatrix(i = 1:n, j=c(2:n, 1), x=-1, dims=c(n,n))
    K_mat <- C + G
    # expand K to be (0 1 1 K)
    K_mat <- rbind(0, K_mat)
    K_mat <- cbind(1, K_mat)
    K_mat[1, -1] <- h; 
    K_mat[1, 1] <- 0
    K = ngme_as_sparse(K_mat)
    h <- c(1, h)
  }

  ngme_operator(
    mesh = mesh,
    cyclic = cyclic,
    matrices = list(K),
    model = "rw1",
    generic_type = "generic",
    theta_K = double(0),
    update_K = function(theta_K) {K},
    K = K,
    h = h,
    symmetric = FALSE,
    zero_trace = FALSE,
    param_name = character(0)
  )
}

#' Random Walk Model of Order 2 (RW2)
#'
#' Constructs a second-order random walk model for spatial or temporal processes.
#' The RW2 model assumes that second-order differences 
#' \eqn{\Delta^2 W_i = W_i - 2W_{i-1} + W_{i-2}} are independent and identically 
#' distributed Gaussian variables.
#'
#' @details
#' The RW2 model is defined by the precision matrix \eqn{K} that penalizes 
#' second-order differences, making it smoother than RW1. The model enforces
#' constraints to ensure identifiability:
#' 
#' **Non-cyclic (default)**: The first two rows of \eqn{K} enforce constraints:
#' the first row implements \eqn{\sum_{i=1}^n h_i W_i = 0} (sum-to-zero), 
#' and the second row implements \eqn{\sum_{i=1}^n h_i \cdot i \cdot W_i = 0} 
#' (removes linear trend). The remaining rows penalize second-order differences.
#' 
#' **Cyclic**: Treats the domain as circular, connecting the first and last locations
#' as neighbors. The constraint \eqn{\sum_{i=1}^n h_i W_i = 0} is enforced.
#' The precision matrix is expanded to size \eqn{(n+1) \times (n+1)} to handle
#' the constraint and circular structure properly.
#'
#' @param mesh numerical vector or inla.mesh.1d object, locations to build the mesh.
#'   For numerical vectors, assumes equally spaced locations. Must have at least 3 locations.
#' @param constr logical, whether to enforce the sum-to-zero constraint \eqn{\sum_{i=1}^n h_i W_i = 0}
#'   and the linear trend constraint \eqn{\sum_{i=1}^n (\sum_{j=1}^i h_j - h_1) \cdot W_i = 0}.
#'   If FALSE, fixes the first and second elements \eqn{W_1 = 0} and \eqn{W_2 = 0}.
#' @param cyclic logical, whether the mesh is circular. If TRUE, the first and last
#'   locations are treated as neighbors with second-order differences computed
#'   across the boundary.
#' @param ... additional arguments (currently unused)
#'
#' @return An `ngme_operator` object containing the precision matrix and related
#'   components for the RW2 model.
#' @export
#'
#' @examples
#' # Non-cyclic RW2 with constraints (default)
#' rw2_default <- rw2(1:6)
#' print(rw2_default$K)
#' 
#' # Cyclic RW2 (connects first and last locations)
#' rw2_cyclic <- rw2(1:6, cyclic = TRUE)
#' print(rw2_cyclic$K)
#' 
#' # Using with unequally spaced locations
#' locations <- c(0, 1, 3, 6, 10, 15)
#' rw2_unequal <- rw2(locations)
#' print(rw2_unequal$K)
rw2 <- function(
  mesh,
  cyclic = FALSE,
  constr = TRUE,
  ...
) {
  mesh <- ngme_build_mesh(mesh)

  stopifnot("Mesh should be inla.mesh.1d." = inherits(mesh, c("inla.mesh.1d")))
  n <- mesh$n
  h_left <- c(0, diff(mesh$loc))
  h_right <- c(diff(mesh$loc), 0)
  h <- 0.5 * (h_left + h_right)
  stopifnot("mesh too small" = n >= 3)

  if (!cyclic) {
    C <- Matrix::sparseMatrix(i = 3:n, j=2:(n-1), x=-2, dims=c(n,n))
    G <- Matrix::sparseMatrix(i = c(1:n, 3:n), j=c(1:n, 1:(n-2)), x=1, dims=c(n,n))
    K <- as.matrix(C + G)
    if (constr) {
      K[1, ] <- h
      K[2, ] <- cumsum(h) - h[1]
    } else {
      K[1, 1] <- 1
      K[2, 2] <- 1
    }
    K <- ngme_as_sparse(K)
  } else {
    C <- Matrix::sparseMatrix(i = 1:n, j=c(2:n,1), x=-2, dims=c(n,n))
    G <- Matrix::sparseMatrix(i = rep(1:n,2), j=c(1:n, 3:n, 1, 2), x=1, dims=c(n,n))
    K <- as.matrix(C + G)
    K <- rbind(0, K)
    K <- cbind(1, K)
    K[1, -1] <- h;
    K[1, 1] <- 0
    K = ngme_as_sparse(K)
    h <- c(1, h)
    K <- ngme_as_sparse(K)
  }

  ngme_operator(
    mesh = mesh,
    model = "rw2",
    generic_type = "generic",
    cyclic = cyclic,
    matrices = list(K),
    theta_K = double(0),
    update_K = function(theta_K) {K},
    K = K,
    h = h,
    symmetric = FALSE,
    zero_trace = FALSE,
    param_name = character(0)
  )
}

#' ngme Ornstein–Uhlenbeck process specification
#'
#' @param mesh numerical vector or inla.mesh.1d object, index to build the mesh
#' @param mesh mesh for build the model
#' @param theta_K initial value for theta_K, kappa = exp(B_theta_K * theta_K)
#' @param B_theta_K bases for theta_K
#' @param ... ignore
#'
#' @return ngme_operator object
#' @export
ou <- function(
  mesh,
  theta_K   = 0,
  B_theta_K = NULL,
  ...
) {
  mesh <- ngme_build_mesh(mesh)
  n <- mesh$n

  # h <- diff(mesh$loc); h <- c(h, mean(h))
  h_left <- c(0, diff(mesh$loc))
  h_right <- c(diff(mesh$loc), 0)
  h <- 0.5 * (h_left + h_right)

  if (is.null(B_theta_K)) B_theta_K <- matrix(1, nrow = length_map(mesh$loc), ncol = 1)
  stopifnot("B_theta_K is a matrix" = is.matrix(B_theta_K))
  stopifnot("ncol(B_theta_K) == length(theta_K)"
    = ncol(B_theta_K) == length(theta_K))

  G <- Matrix::bandSparse(n=n,m=n,k=c(-1,0),diagonals=cbind(-rep(1,n), rep(1,n)))
  C <- Ce <- Matrix::bandSparse(n=n,m=n,k=c(-1,0),diagonals=cbind(0.5*c(h[-1],0), 0.5*h))
  Ci <- Matrix::sparseMatrix(i=1:n,j=1:n,x=1/h,dims = c(n,n))

  kappas <- exp(as.numeric(B_theta_K %*% theta_K))
  K <- Matrix::Diagonal(x=kappas) %*% C + G

  update_K <- function(theta_K) {
    kappas <- exp(as.numeric(B_theta_K %*% theta_K))
    Matrix::Diagonal(x=kappas) %*% C + G
  }

  C <- ngme_as_sparse(C)
  G <- ngme_as_sparse(G)
  K <- ngme_as_sparse(K)

  generic_ns(
    model = "ou",
    C = C,
    G = G,
    param_name  = paste("theta_K", seq_len(length(theta_K)), sep = ""),
    param_trans = rep(list(identity), length(theta_K)),

    matrices = list(C, G),
    theta_K = list(theta_K = theta_K),
    trans = list(theta_K = "exp"),
    B_theta_K = list(theta_K = B_theta_K),
    position = list(
      c(1, 2),
      c(3)
    ),
    h = h,
    mesh = mesh
  )
}

#' ngme Matern SPDE model specification
#'
#' @param mesh an fmesher::fm_mesh_2d object, mesh for build the SPDE model
#' @param alpha 2 or 4, SPDE smoothness parameter
#' @param kappa the range parameter, for the stationary model, it is the range parameter
#' @param theta_K kappa = exp(B_theta_K * theta_K), for the non-stationary model, it is the range parameter
#' @param B_theta_K bases for theta_K, ignore if use the stationary model
#' @param ... ignore
#'
#' @return ngme_operator object
#' @export
matern <- function(
  mesh,
  alpha = 2,
  kappa = NULL,
  theta_K = NULL,
  B_theta_K = NULL,
  ...
) {
  mesh <- ngme_build_mesh(mesh)
  stopifnot("alpha should be 2 or 4" = alpha == 2 || alpha == 4)

  if (is.null(kappa)) {
    if (is.null(theta_K)) {
      theta_K <- 0
    }
  } else {
    stopifnot("kappa should be a single value" = length(kappa) == 1)
    theta_K <- log(kappa)
  }

  if (inherits(mesh, "metric_graph")) {
    if (is.null(mesh$mesh$C)) {
      mesh$compute_fem()
    }
    C <- mesh$mesh$C
    G <- mesh$mesh$G
    h <- mesh$mesh$weights
    d <- 0
  } else {
    d <- get_inla_mesh_dimension(mesh)
    if (d == 1) {
      fem <- fmesher::fm_fem(mesh, order = alpha)
      C <- fem$c0
      G <- fem$g1
      h <- Matrix::diag(fem$c0)
    } else if (d == 2) {
      # including S2 mesh
      fem <- fmesher::fm_fem(mesh, order = alpha)
      C <- fem$c0  # diag
      G <- fem$g1
      h <- Matrix::diag(fem$c0)
    }
  }
  Cinv <- C;
  diag(Cinv) <- 1 / Matrix::diag(C)

  if (is.null(theta_K)) {
    if (inherits(mesh, "metric_graph"))
      theta_K = 0
    else {
      # loc <- if (d == 1) mesh$loc else mesh$loc[, c(1,2)]
      # dist_mat <- as.matrix(dist(loc))
      # max_dist <- max(dist_mat[lower.tri(dist_mat)])
      # range = max(4*min(h), max_dist)
      # theta_K = log(sqrt(8*3/2)/(0.2*range))
      theta_K = 0
    }
  }

  mesh_n <- length(h)
  if (is.null(B_theta_K) && length(theta_K) == 1)
    B_theta_K <- matrix(1, nrow = mesh_n, ncol = 1)
  else if (is.null(B_theta_K) && length(theta_K) > 1)
    stop("Please provide B_theta_K for non-stationary case.")

  stationary <- is_stationary(B_theta_K)
  update_K <- function(theta_K) {
    kappas <- as.numeric(exp(B_theta_K %*% theta_K))
    if (stationary) {
      if (alpha == 2) {
        kappas[1]^2 * C + G
      } else {
        Cinv <- C;
        diag(Cinv) <- 1 / Matrix::diag(C)
        (G + kappas[1]^2 * C) %*% Cinv %*% (G + kappas[1]^2 * C)
        # same as kappa^4 C + 2 kappa^2 G + G Cinv G
      }
    } else {
      GpKCK <- diag(kappas) %*% C %*% diag(kappas) + G
      if (alpha == 2)
        GpKCK
      else {
        Cinv <- C;
        diag(Cinv) <- 1 / Matrix::diag(C)
        (GpKCK) %*% Cinv %*% GpKCK
      }
    }
  }
  K <- update_K(theta_K)

  C <- ngme_as_sparse(C)
  G <- ngme_as_sparse(G)
  Cinv <- ngme_as_sparse(Cinv)

  if (stationary && alpha == 2) {
    matrices = list(C, G)
    theta_K = c(theta_K = theta_K)
    trans = c(theta_K = "exp2")
    position = NULL
  } else if (stationary && alpha == 4) {
    matrices = list(C, 2*G, G %*% Cinv %*% G)
    theta_K = c(theta_K=theta_K)
    trans = list(theta_K=c("exp4", "exp2", "null"))
    position = NULL
  } else if (alpha == 2) {
    theta_K = list(theta_K=theta_K)
    trans = list(theta_K=c("exp"))
    B_theta_K = list(theta_K = B_theta_K)
    matrices = list(C, G, Cinv)
    position = list(
      c(1, 2, 1), 
      c(3) 
    )
  } else if (alpha == 4) {
    theta_K = list(theta_K=theta_K)
    trans = list(theta_K=c("exp"))
    B_theta_K = list(theta_K = B_theta_K)
    matrices = list(C, G, Cinv)
    position = list(
      c(1, 2, 1, 4, 1, 2, 1), 
      c(1, 2, 1, 4, 3), 
      c(3, 4, 1, 2, 1), 
      c(3, 4, 3)
    )
  }

  if (stationary) g <- generic else g <- generic_ns

  g(
    model = "matern",
    mesh = mesh,
    alpha = alpha,
    theta_K = theta_K,
    trans = trans,
    matrices = matrices,
    position = position,
    B_theta_K = B_theta_K,
    C = C,
    G = G,
    # update_K = update_K,
    h = h,
    stationary = stationary,
    param_name =
      if (stationary) "kappa"
      else paste("theta_K", seq_len(length(theta_K)), sep = " "),
    param_trans =
      if (stationary) list(exp)
      else rep(list(identity), length(theta_K))
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

  B_theta_K <- map
  n_reff <- ncol(B_theta_K) # number of random effects
  n_theta_K <- sum(1:n_reff) # number of theta_K
  h <- rep(1, n_reff)

  # provide initial value for theta_K
  if (!is.null(theta_K)) {
    stopifnot(length(theta_K) == n_theta_K)
  } else {
    theta_K <- rep(0, n_theta_K)
  }

  # build K
  update_K <- function(theta_K) {
    K <- diag(n_reff); diag(K) <- exp(theta_K[1:n_reff])
    if (n_reff > 1)
      K[lower.tri(K)] <- theta_K[(n_reff+1):n_theta_K]
    Matrix::Matrix(K)
  }
  K <- update_K(theta_K)

  ngme_operator(
    mesh = NULL,
    model = "re",
    theta_K = theta_K,
    K = ngme_as_sparse(K),
    update_K = update_K,
    h = h,
    B_theta_K = B_theta_K,
    symmetric = FALSE,
    zero_trace = FALSE,
    param_name = paste("re_coeff", seq_len(length(theta_K)), sep = ""),
    param_trans = rep(list(identity), length(theta_K))
  )
}

# ----  For computing precision matrix of multivariate model

# p: dimension
# rho: a list of rho parameter, see paper D_l
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

# return Q * D_l
dependence_matrix <- function(p, rho, theta=NULL, Q=NULL) {
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
#' rho <- c(-0.5, 0.5,-0.25) #correlation parameters
#' operator_list <- list(ar1(1:5, rho=0.4), ar1(1:5, rho=0.5), ar1(1:5, rho=0.6))
#' precision_matrix_multivariate(3, operator_list, rho, theta=c(1,2,3))
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
#' @param mesh an fmesher::fm_mesh_2d object, mesh for build the SPDE model
#' @param alpha_list a list of SPDE smoothness parameter
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
#' library(fmesher)
#' library(fields)
#' # Define mesh
#' x <- seq(from=0,to=1,length.out = 40)
#' mesh <-  fm_rcdt_2d_inla(lattice = fm_lattice_2d(x,x), extend = FALSE)
#' # Set parameters
#' p <- 3 #number of fields
#' rho <- c(-0.5, 0.5,-0.25) #correlation parameters
#' log_kappa <- list(2,2,2) #log(kappa)
#' variances <- list(1,1,1) #set marginal variances to 1
#' alpha <- list(2,2,2) #smoothness parameters
#' # Compute precision
#' Q <- precision_matrix_multivariate_spde(p, mesh = mesh, rho = rho,
#'                                        alpha = alpha, theta_K_list = log_kappa,
#'                                        variance_list = variances)
#' # Plot the cross covariances
#' A <- as.vector(fm_basis(mesh,loc = matrix(c(0.5,0.5),1,2)))
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
#' proj <- fm_evaluator(mesh)
#'
#' par(mfrow=c(3,3))
#' image.plot(fm_evaluate(proj,r11), main = "Cov(X_1(s0),X_1(s)")
#' plot.new()
#' plot.new()
#' image.plot(fm_evaluate(proj,r12), main = "Cov(X_1(s0),X_2(s)")
#' image.plot(fm_evaluate(proj,r22), main = "Cov(X_2(s0),X_2(s)")
#' plot.new()
#' image.plot(fm_evaluate(proj,r13), main = "Cov(X_1(s0),X_3(s)")
#' image.plot(fm_evaluate(proj,r23), main = "Cov(X_2(s0),X_3(s)")
#' image.plot(fm_evaluate(proj,r33), main = "Cov(X_3(s0),X_3(s)")

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
