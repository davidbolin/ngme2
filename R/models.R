# -------- ngme operators --------

#' ngme iid model specification
#'
#' @param mesh integer or factor, index to build the mesh
#'
#' @return ngme_operator object
#' @export
#'
iid <- function(mesh = NULL) {
  if (is.null(mesh) || (is.list(mesh) && !inherits(mesh, "inla.mesh.1d"))) {
    return(structure(list(model = "iid", args = as.list(environment())), class = "ngme_operator_def"))
  }

  mesh <- ngme_build_mesh(mesh, "iid")
  n <- mesh$n
  K <- ngme_as_sparse(Matrix::Diagonal(n))

  ngme_operator(
    matrices = list(K),
    mesh = mesh,
    model = "iid",
    theta_K = double(0),
    update_K = function(theta_K) {
      K
    },
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
#'
#' @return ngme_operator object
#' @export
#'
ar1 <- function(mesh = NULL, rho = 0) {
  if (is.null(mesh) || (is.list(mesh) && !inherits(mesh, "inla.mesh.1d"))) {
    return(structure(list(model = "ar1", args = as.list(environment())), class = "ngme_operator_def"))
  }
  stopifnot("rho should be between -1 and 1" = rho >= -1 && rho <= 1)

  mesh <- ngme_build_mesh(mesh)
  n <- mesh$n
  h <- c(diff(mesh$loc), 1)
  G <- Matrix::Diagonal(n)
  G[1, 1] <- sqrt(1 - rho**2)
  C <- Matrix::sparseMatrix(j = 1:(n - 1), i = 2:n, x = -1, dims = c(n, n))
  G <- ngme_as_sparse(G)
  C <- ngme_as_sparse(C)
  stopifnot("The mesh should be 1d and has gap 1." = all(h == 1))

  theta_K <- ar1_a2th(rho)
  stopifnot("The length of rho(theta_K) should be 1." = length(theta_K) == 1)

  update_K <- function(theta_K) {
    ar1_th2a(theta_K) * C + G
  }

  # Create a generic model internally
  g <- name2fun("tanh", inv = TRUE)

  ngme_operator(
    mesh = mesh,
    model = "ar1",
    theta_K = c(rho = g(rho)),
    trans = c(rho = "tanh"),
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

#' ngme ARMA(p, q) model specification
#'
#' General ARMA(p, q) latent operator under AZW: K W = Z eps,
#' where K encodes the AR part and Z encodes the MA part.
#'
#' @param mesh Integer vector or `inla.mesh.1d` object.
#' @param ar_order Integer p (0, 1, or 2).
#' @param ma_order Integer q (0, 1, or 2).
#' @param ar Numeric length-p vector of AR coefficients in (-1, 1). Defaults to zeros.
#' @param ma Numeric length-q vector of MA coefficients in (-1, 1). Defaults to zeros.
#' @param fix_ar Logical or length-p logical vector to fix AR coefficients.
#' @param fix_ma Logical or length-q logical vector to fix MA coefficients.
#'
#' @details
#' - Supports p, q <= 2. The user provides standard AR/MA coefficients (ar, ma).
#'   For estimation, these are transformed in R into an unbounded parameter vector
#'   `theta_K` via PACF-style mappings so optimizers (e.g., Adam/SGD) operate on
#'   unconstrained variables while K and Z remain valid:
#'   - AR(1): phi1 = tanh(raw1/2). AR(2): pi1=tanh(raw1/2), pi2=tanh(raw2/2),
#'     phi2 = pi2, phi1 = pi1 * (1 - pi2).
#'   - MA(1): theta1 = tanh(raw1/2). MA(2): psi1=tanh(raw1/2), psi2=tanh(raw2/2),
#'     theta2 = psi2, theta1 = psi1 * (1 - psi2).
#' - The returned operator includes both K (AR part) and Z (MA part). During
#'   estimation the backend updates Z synchronously every time `theta_K` changes.
#'
#' @return An `ngme_operator` describing the ARMA(p, q) latent model.
#'
#' @examples
#' mesh <- 1:10
#' op <- arma(mesh, ar_order = 2, ma_order = 2, ar = c(0.5, 0.3), ma = c(0.6, 0.4))
#'
#' @export
arma <- function(
    mesh = NULL,
    ar_order = 1,
    ma_order = 1,
    ar = NULL, # AR coefficients (friendly scale)
    ma = NULL, # MA coefficients (friendly scale)
    fix_ar = FALSE,
    fix_ma = FALSE) {
  if (is.null(mesh) || (is.list(mesh) && !inherits(mesh, "inla.mesh.1d"))) {
    return(structure(list(model = "arma", args = as.list(environment())), class = "ngme_operator_def"))
  }
  p <- as.integer(ar_order)
  q <- as.integer(ma_order)
  stopifnot("ar_order >= 0" = p >= 0, "ma_order >= 0" = q >= 0)

  if (is.null(ar)) ar <- rep(0, p) else stopifnot(length(ar) == p)
  if (is.null(ma)) ma <- rep(0, q) else stopifnot(length(ma) == q)

  mesh <- ngme_build_mesh(mesh)
  n <- mesh$n
  h <- c(diff(mesh$loc), 1)
  stopifnot("The mesh should be 1d and has gap 1." = all(h[1:(n - 1)] == 1))

  # Build AR bases for placeholder K printing (Cj and G)
  Cs <- vector("list", max(1, p))
  if (p > 0) {
    for (j in 1:p) {
      Cs[[j]] <- Matrix::sparseMatrix(i = (p + 1):n, j = (p + 1 - j):(n - j), x = -1, dims = c(n, n))
    }
  } else {
    Cs[[1]] <- Matrix::sparseMatrix(i = 1, j = 1, x = 0, dims = c(n, n))
  }
  G <- Matrix::Diagonal(n)
  Cs <- lapply(Cs, ngme_as_sparse)
  G <- ngme_as_sparse(G)

  # Helper mappings (friendly <-> raw via AR1 transforms)
  clamp_unit <- function(x) pmax(-1 + 1e-8, pmin(1 - 1e-8, x))
  to_raw_vec <- function(coef, order) {
    if (order == 0) {
      return(numeric(0))
    }
    coef <- clamp_unit(coef)
    if (order == 1) {
      return(ar1_a2th(coef[1]))
    } else if (order == 2) {
      u2 <- coef[2]
      u1 <- coef[1] / (1 - u2)
      # clamp derived PACF before transform
      u1 <- clamp_unit(u1)
      u2 <- clamp_unit(u2)
      return(c(ar1_a2th(u1), ar1_a2th(u2)))
    } else {
      stop("Only order <= 2 supported")
    }
  }
  from_raw_vec <- function(raw, order) {
    if (order == 0) {
      return(numeric(0))
    }
    if (order == 1) {
      v1 <- ar1_th2a(raw[1])
      return(v1)
    } else if (order == 2) {
      t1 <- ar1_th2a(raw[1])
      t2 <- ar1_th2a(raw[2])
      v2 <- t2
      v1 <- t1 * (1 - t2)
      return(c(v1, v2))
    } else {
      stop("Only order <= 2 supported")
    }
  }

  # Build raw theta_K from friendly ar/ma
  theta_ar_raw <- to_raw_vec(ar, p)
  theta_ma_raw <- to_raw_vec(ma, q)
  theta_K <- c(theta_ar_raw, theta_ma_raw)
  names(theta_K) <- c(if (p > 0) paste0("ar", 1:p), if (q > 0) paste0("ma", 1:q))

  # placeholder K using friendly AR coefficients (printing)
  K0 <- G
  if (p > 0) {
    for (j in 1:p) K0 <- K0 + ar[j] * Cs[[j]]
  }
  K0 <- ngme_as_sparse(K0)

  # Attach spec for MA part (names for MA params now 'ma1..ma_q')
  Z_spec <- list(type = "ma", order = q, param = paste0("ma", 1:q))

  # Build MA shift matrices L^k = -Cs[[k]] (since Cs has -1 on sub-diagonal positions)
  Ls <- vector("list", max(1, q))
  if (q > 0) {
    for (j in 1:q) Ls[[j]] <- -Cs[[j]]
  } else {
    Ls[[1]] <- Matrix::sparseMatrix(i = 1, j = 1, x = 0, dims = c(n, n))
  }

  # placeholder Z using friendly MA coefficients: Z = I + sum theta_j L^j
  Z0 <- G
  if (q > 0) {
    for (j in 1:q) Z0 <- Z0 + ma[j] * Ls[[j]]
  }
  Z0 <- ngme_as_sparse(Z0)

  # Build per-parameter transforms mapping raw -> user-friendly proxy.
  # Note: for order 2, transforms map to PACF components (tanh(raw/2)),
  # because ar1 depends on both raw1 and raw2 and cannot be reconstructed
  # element-wise. The print method shows true AR/MA; traceplot will display PACF
  # unless we add grouped transforms (can be added on request).
  make_trans_list <- function(p, q) {
    trans <- list()
    if (p > 0) trans <- c(trans, rep(list(ar1_th2a), p))
    if (q > 0) trans <- c(trans, rep(list(ar1_th2a), q))
    trans
  }

  ngme_operator(
    mesh = mesh,
    model = "arma",
    theta_K = theta_K,
    # param_space omitted: R passes raw theta_K now
    matrices = c(Cs, list(G)),
    h = h,
    K = K0,
    Z = Z0,
    symmetric = FALSE,
    zero_trace = FALSE,
    param_name = names(theta_K),
    param_trans = make_trans_list(p, q),
    Z_spec = Z_spec,
    fix_rho = rep(as.logical(fix_ar), p),
    fix_phi = rep(as.logical(fix_ma), q),
    p = p, q = q,
    generic_type = "none"
  )
}

#' Convenience wrapper for ARMA(1,1)
#' @export
arma11 <- function(
    mesh = NULL,
    ar = 0,
    ma = 0,
    fix_ar = FALSE,
    fix_ma = FALSE) {
  arma(mesh,
    ar_order = 1,
    ma_order = 1,
    ar = ar,
    ma = ma,
    fix_ar = fix_ar,
    fix_ma = fix_ma
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
  b <- c(1, rep(0, p - 1))

  # Solve A * gamma = b
  gamma <- solve(A, b)

  return(gamma)
}

#' ngme AR(p) model specification
#'
#' @param mesh integer vector or inla.mesh.1d object, index to build the mesh
#' @param rho vector of AR coefficients of length p (should satisfy stationarity conditions)
#' @param order integer, the AR order p. If provided and rho is not specified, rho will be initialized as rep(0, p)
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
    mesh = NULL, rho = NULL, order = NULL) {
  if (is.null(mesh) || (is.list(mesh) && !inherits(mesh, "inla.mesh.1d"))) {
    return(structure(list(model = "ar", args = as.list(environment())), class = "ngme_operator_def"))
  }
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

  stopifnot("The mesh should be 1d and has gap 1." = all(h[1:(n - 1)] == 1))
  stopifnot("n should be >= p" = n >= p)

  # Construct C1, ..., Cp and G matrices
  Cs <- vector("list", p)
  for (j in 1:p) {
    # Cj: for lag j, put -1 at (t, t-j) for t = (j+1):n
    # The first p rows should be all 0, so we only fill from row (p+1) to n
    Cs[[j]] <- Matrix::sparseMatrix(
      i = (p + 1):n,
      j = (p + 1 - j):(n - j),
      x = rep(-1, n - p),
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
    mesh = NULL,
    cyclic = FALSE,
    constr = TRUE) {
  if (is.null(mesh) || (is.list(mesh) && !inherits(mesh, "inla.mesh.1d"))) {
    return(structure(list(model = "rw1", args = as.list(environment())), class = "ngme_operator_def"))
  }
  mesh <- ngme_build_mesh(mesh)

  n <- mesh$n
  h_left <- c(0, diff(mesh$loc))
  h_right <- c(diff(mesh$loc), 0)
  h <- 0.5 * (h_left + h_right)

  if (!cyclic) {
    C <- Matrix::sparseMatrix(i = 1:n, j = 1:n, x = 1, dims = c(n, n))
    G <- Matrix::sparseMatrix(i = 2:n, j = 1:(n - 1), x = -1, dims = c(n, n))
    K_mat <- C + G
    if (constr) K_mat[1, ] <- h else K_mat[1, 1] <- 1
    K <- ngme_as_sparse(K_mat)
  } else {
    stopifnot("constraint cannot be false for cyclic case" = constr)
    stopifnot("Too less data point" = length(h) >= 3)
    C <- Matrix::Diagonal(n)
    G <- Matrix::sparseMatrix(i = 1:n, j = c(2:n, 1), x = -1, dims = c(n, n))
    K_mat <- C + G
    # expand K to be (0 1 1 K)
    K_mat <- rbind(0, K_mat)
    K_mat <- cbind(1, K_mat)
    K_mat[1, -1] <- h
    K_mat[1, 1] <- 0
    K <- ngme_as_sparse(K_mat)
    h <- c(1, h)
  }

  ngme_operator(
    mesh = mesh,
    cyclic = cyclic,
    matrices = list(K),
    model = "rw1",
    generic_type = "generic",
    theta_K = double(0),
    update_K = function(theta_K) {
      K
    },
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
    mesh = NULL,
    cyclic = FALSE,
    constr = TRUE) {
  if (is.null(mesh) || (is.list(mesh) && !inherits(mesh, "inla.mesh.1d"))) {
    return(structure(list(model = "rw2", args = as.list(environment())), class = "ngme_operator_def"))
  }
  mesh <- ngme_build_mesh(mesh)

  stopifnot("Mesh should be inla.mesh.1d." = inherits(mesh, c("inla.mesh.1d")))
  n <- mesh$n
  h_left <- c(0, diff(mesh$loc))
  h_right <- c(diff(mesh$loc), 0)
  h <- 0.5 * (h_left + h_right)
  stopifnot("mesh too small" = n >= 3)

  if (!cyclic) {
    C <- Matrix::sparseMatrix(i = 3:n, j = 2:(n - 1), x = -2, dims = c(n, n))
    G <- Matrix::sparseMatrix(i = c(1:n, 3:n), j = c(1:n, 1:(n - 2)), x = 1, dims = c(n, n))
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
    C <- Matrix::sparseMatrix(i = 1:n, j = c(2:n, 1), x = -2, dims = c(n, n))
    G <- Matrix::sparseMatrix(i = rep(1:n, 2), j = c(1:n, 3:n, 1, 2), x = 1, dims = c(n, n))
    K <- as.matrix(C + G)
    K <- rbind(0, K)
    K <- cbind(1, K)
    K[1, -1] <- h
    K[1, 1] <- 0
    K <- ngme_as_sparse(K)
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
    update_K = function(theta_K) {
      K
    },
    K = K,
    h = h,
    symmetric = FALSE,
    zero_trace = FALSE,
    param_name = character(0)
  )
}

#' Ornstein-Uhlenbeck Process Model
#'
#' Implements the exact discrete-time representation of the continuous
#' Ornstein-Uhlenbeck (OU) process using a matrix formulation
#' \eqn{\mathbf{K}\mathbf{X} = \boldsymbol{\delta}}.
#'
#' @details
#'
#' The OU process is defined by the stochastic differential equation (SDE):
#'
#' \deqn{dX_t = -\theta X_t dt + \sigma dW_t}
#'
#' where \eqn{\theta > 0} is the mean-reversion rate (stiffness) and
#' \eqn{\sigma} is the volatility (diffusion coefficient).
#'
#' ## Exact Discretization
#'
#' Unlike Euler-Maruyama approximations, the OU process has an exact solution
#' between time points \eqn{t_i} and \eqn{t_{i+1}} with step size \eqn{\Delta t_i}:
#'
#' \deqn{X_{i+1} = X_i e^{-\theta \Delta t_i} + \epsilon_i}
#'
#' This is an AR(1) form with autoregressive coefficient:
#'
#' \deqn{\rho_i = e^{-\theta \Delta t_i}}
#'
#' and Gaussian noise:
#'
#' \deqn{\epsilon_i \sim \mathcal{N}\left(0, \frac{\sigma^2}{2\theta}(1 - \rho_i^2)\right)}
#'
#' ## Matrix Representation
#'
#' The precision matrix \eqn{\mathbf{K}} for the linear system
#' \eqn{\mathbf{K}\mathbf{X} = \boldsymbol{\delta}} is constructed as:
#'
#' \deqn{\mathbf{K} = \begin{bmatrix}
#' \sqrt{1-\rho_1^2} & 0 & 0 & \cdots & 0 \\
#' -\rho_1 & 1 & 0 & \cdots & 0 \\
#' 0 & -\rho_2 & 1 & \cdots & 0 \\
#' \vdots & \vdots & \ddots & \ddots & \vdots \\
#' 0 & 0 & \cdots & -\rho_{n-1} & 1
#' \end{bmatrix}}
#'
#' where \eqn{\rho_i = \exp(-\theta \Delta t_i)} with \eqn{\Delta t_i = t_{i+1} - t_i}.
#'
#' **Why the \eqn{\sqrt{1-\rho_1^2}} term?** The first row ensures that \eqn{X_1}
#' is drawn from the stationary distribution \eqn{\mathcal{N}(0, \sigma^2/(2\theta))}.
#' This scaling factor matches the variance of the first component of
#' \eqn{\boldsymbol{\delta}} with the transition noise variance in subsequent steps.
#'
#' ## Non-uniform Mesh
#'
#' For non-uniform meshes, each \eqn{\rho_i} is computed locally based on the
#' varying time step \eqn{\Delta t_i}, allowing the model to naturally handle
#' irregularly spaced observations.
#'
#' ## The \eqn{h} vector
#'
#' The \eqn{h} vector represents integration weights for the mass matrix.
#' For the OU model, it is set as:
#'
#' \deqn{h = [\Delta t_1, \Delta t_2, \ldots, \Delta t_{n-1}, \Delta t_{n-1}]}
#'
#' where the last element is duplicated. This ensures proper weighting in the
#' finite element formulation.
#'
#' ## The \eqn{\boldsymbol{\delta}} vector
#'
#' In the equation \eqn{\mathbf{K}\mathbf{X} = \boldsymbol{\delta}}, the vector
#' \eqn{\boldsymbol{\delta}} contains independent identically distributed (i.i.d.)
#' random variables with unit variance. The precision matrix \eqn{\mathbf{K}}
#' is constructed such that the resulting process \eqn{\mathbf{X}} has the
#' correct OU covariance structure.
#'
#' @param mesh numerical vector or `inla.mesh.1d` object giving the ordered
#'   time locations. Must be strictly increasing.
#' @param theta positive mean-reversion rate (stiffness parameter). Internally
#'   stored on the log scale so optimization remains unconstrained.
#'
#' @return An `ngme_operator` object containing:
#' \itemize{
#'   \item `K`: the precision matrix
#'   \item `h`: integration weights
#'   \item `theta_K`: log-transformed theta for unconstrained optimization
#'   \item `mesh`: the mesh object
#' }
#'
#' @references
#' Uhlenbeck, G. E., & Ornstein, L. S. (1930). On the theory of the Brownian motion.
#' Physical review, 36(5), 823.
#'
#' @export
#' @examples
#' # Uniform mesh
#' mesh_uniform <- seq(0, 10, by = 0.5)
#' ou_uniform <- ou(mesh = mesh_uniform, theta = 0.5)
#' print(ou_uniform)
#'
#' # Non-uniform mesh
#' mesh_nonuniform <- c(0, 1, 2.5, 3.5, 5, 8, 11)
#' ou_nonunif <- ou(mesh = mesh_nonuniform, theta = 0.8)
#' print(ou_nonunif)
ou <- function(
    mesh = NULL,
    theta = 1) {
  if (is.null(mesh) || (is.list(mesh) && !inherits(mesh, "inla.mesh.1d"))) {
    return(structure(list(model = "ou", args = as.list(environment())), class = "ngme_operator_def"))
  }

  mesh <- ngme_build_mesh(mesh)
  stopifnot("Mesh should be inla.mesh.1d." = inherits(mesh, c("inla.mesh.1d")))

  n <- mesh$n
  stopifnot("mesh should have at least two locations" = n >= 2)

  dt <- diff(mesh$loc)
  stopifnot("mesh locations must be strictly increasing" = all(dt > 0))

  theta_raw <- log(theta)
  stopifnot("theta must be a single positive value" = length(theta) == 1 && is.finite(theta) && theta > 0)

  build_K <- function(theta_val) {
    rho <- exp(-theta_val * dt)
    stopifnot("All rho values must lie in (0, 1)" = all(rho > 0 & rho < 1))

    diag_main <- c(sqrt(pmax(0, 1 - rho[1]^2)), rep(1, n - 1))
    sub_diag <- -rho

    Matrix::bandSparse(
      n = n,
      k = c(0, -1),
      diagonals = list(diag_main, sub_diag)
    )
  }

  update_K <- function(theta_K) {
    ngme_as_sparse(build_K(exp(theta_K)))
  }

  K <- update_K(theta_raw)
  h <- c(dt, utils::tail(dt, 1))

  ngme_operator(
    mesh = mesh,
    model = "ou",
    theta_K = c(theta = theta_raw),
    trans = c(theta = "exp"),
    update_K = update_K,
    K = K,
    h = h,
    symmetric = FALSE,
    zero_trace = FALSE,
    param_name = "theta",
    param_trans = list(exp),
    generic_type = "none"
  )
}

#' ngme Matern SPDE model specification
#'
#' @param mesh an fmesher::fm_mesh_2d object, mesh for build the SPDE model
#' @param alpha SPDE smoothness parameter (alpha = 2beta) 2 or 4 for integer case, otherwise for fractional case
#' @param kappa the range parameter, for the stationary model, it is the range parameter
#' @param theta_kappa kappa = exp(B_kappa * theta_kappa), for the non-stationary model, it is the range parameter
#' @param fix_alpha if FALSE, enable fractional modeling
#' @param rational_order order of rational approximation for fractional operators (default: 2)
#' @param B_kappa bases for theta_kappa, ignore if use the stationary model
#'
#' @return ngme_operator object
#' @export
matern <- function(
    mesh = NULL,
    alpha = 2,
    fix_alpha = TRUE,
    kappa = NULL,
    theta_kappa = NULL,
    rational_order = 2,
    B_kappa = NULL) {
  if (is.null(mesh) || (is.list(mesh) && !inherits(mesh, c("inla.mesh.1d", "inla.mesh", "fm_mesh_1d", "fm_mesh_2d", "metric_graph")))) {
    return(structure(list(model = "matern", args = as.list(environment())), class = "ngme_operator_def"))
  }
  mesh <- ngme_build_mesh(mesh)
  if (fix_alpha && alpha != 2 && alpha != 4) {
    if (!requireNamespace("rSPDE", quietly = TRUE)) {
      stop("For fixed alpha values not equal to 2 or 4, the 'rSPDE' package is required. Please install it with: install.packages('rSPDE')")
    }
  }

  # Get mesh dimensions first (needed for B_kappa initialization)
  if (inherits(mesh, "metric_graph")) {
    if (is.null(mesh$mesh$C)) mesh$compute_fem()
    C <- mesh$mesh$C
    G <- mesh$mesh$G
    h <- mesh$mesh$weights
    d <- 0
    # C is not mass lumped..
    C <- Matrix::Diagonal(nrow(C), Matrix::rowSums(C))
  } else {
    d <- get_inla_mesh_dimension(mesh)
    fem <- fmesher::fm_fem(mesh, order = alpha)
    C <- fem$c0
    G <- fem$g1
    h <- Matrix::diag(fem$c0)
  }
  Ci <- Matrix::Diagonal(x = 1 / Matrix::diag(C))

  mesh_n <- length(h)

  # Optimized initialization logic for kappa, theta_kappa, and B_kappa
  if (is.null(kappa)) {
    if (is.null(theta_kappa)) {
      # Case 1: kappa = theta_kappa = NULL
      if (is.null(B_kappa)) {
        # Stationary case: set theta_kappa=0, B_kappa as n_mesh*1 matrix of ones
        theta_kappa <- 0
        B_kappa <- matrix(1, nrow = mesh_n, ncol = 1)
      } else {
        # Non-stationary case: B_kappa is given, auto-set theta_kappa
        theta_kappa <- rep(0, ncol(B_kappa))
      }
    } else {
      # Case 2: theta_kappa is provided, kappa is NULL
      if (is.null(B_kappa)) {
        if (length(theta_kappa) == 1) {
          # Stationary case
          B_kappa <- matrix(1, nrow = mesh_n, ncol = 1)
        } else {
          stop("Please provide B_kappa for non-stationary case.")
        }
      }
    }
  } else {
    # Case 3: kappa is provided
    stopifnot("kappa should be a single value" = length(kappa) == 1)
    theta_kappa <- log(kappa)
    if (is.null(B_kappa)) {
      B_kappa <- matrix(1, nrow = mesh_n, ncol = 1)
    }
  }

  stationary <- is_stationary(B_kappa)
  # Build parameter vector theta_K
  # theta_K = c(alpha, theta_kappa) if alpha is not fixed, otherwise theta_K = theta_kappa
  L <- d / 2
  alpha_clamped <- max(L + 1e-6, min(4 - 1e-6, alpha))
  if (fix_alpha) {
    # Do not include alpha in theta when fixed; only pass kappa params
    theta_K <- if (stationary) c(theta_kappa) else c(as.numeric(theta_kappa))
  } else {
    eta_alpha <- qlogis((alpha_clamped - L) / (4 - L))
    theta_K <- if (stationary) c(eta_alpha, theta_kappa) else c(eta_alpha, as.numeric(theta_kappa))
  }

  C <- ngme_as_sparse(C)
  G <- ngme_as_sparse(G)
  Ci <- ngme_as_sparse(Ci)
  B_K <- if (stationary) matrix(1, nrow = mesh_n, ncol = 1) else B_kappa

  # Placeholder K for printing; real K is built in C++
  update_K <- function(theta_K) {
    if (stationary) {
      K <- exp(if (fix_alpha) theta_K[1] else theta_K[2])^2 * C + G
    } else {
      kappas <- exp(as.numeric(B_K %*% (if (fix_alpha) theta_K else theta_K[-1])))
      K <- Matrix::Diagonal(x = kappas^2) %*% C + G
    }
    ngme_as_sparse(K)
  }

  # Build operator args and conditionally include roots and param meta
  op_args <- list(
    mesh = mesh,
    model = "matern",
    update_K = update_K,
    K = update_K(theta_K),
    h = h,
    theta_K = theta_K,
    generic_type = "none",
    C = C,
    Ci = Ci,
    G = G,
    B_K = B_K,
    alpha = alpha,
    fix_alpha = fix_alpha,
    spatial_dim = d,
    symmetric = TRUE,
    stationary = stationary,
    rational_order = rational_order,
    param_name = NULL,
    param_trans = NULL
  )
  # Parameter names/transforms for traceplot
  if (fix_alpha) {
    if (stationary) {
      op_args$param_name <- c("kappa")
      op_args$param_trans <- list(exp)
    } else {
      p <- ncol(B_kappa)
      op_args$param_name <- paste0("theta_kappa", seq_len(p))
      op_args$param_trans <- rep(list(identity), p)
    }
  } else {
    L <- d / 2
    alpha_trans <- function(eta) L + (4 - L) * (1 / (1 + exp(-eta)))
    if (stationary) {
      op_args$param_name <- c("alpha", "kappa")
      op_args$param_trans <- list(alpha_trans, exp)
    } else {
      p <- ncol(B_kappa)
      op_args$param_name <- c("alpha", paste0("theta_kappa", seq_len(p)))
      op_args$param_trans <- c(list(alpha_trans), rep(list(identity), p))
    }
  }
  do.call(ngme_operator, op_args)
}

#' ngme random effect model
#'
#' @param map numerical vector, covariates to build index for the process (can be formula, provided data)
#' @param theta_K initial value for theta_K (build covariance matrix)
#'
#' @return ngme_operator object
#' @export
re <- function(
    map,
    theta_K = NULL,
    ...) {
  if (missing(map) || is.null(map)) {
    return(structure(list(model = "re", args = list(theta_K = theta_K, ...)), class = "ngme_operator_def"))
  }
  if (inherits(map, "formula")) {
    map <- model.matrix(map)
  } else {
    map <- as.matrix(map)
  }

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
    K <- diag(n_reff)
    diag(K) <- exp(theta_K[1:n_reff])
    if (n_reff > 1) {
      K[lower.tri(K)] <- theta_K[(n_reff + 1):n_theta_K]
    }
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
      length(rho) == p * (p - 1) / 2
  )
  D_l <- diag(p)
  D_l[lower.tri(D_l)] <- rho
  # compute k(j)
  k <- double(p)
  k[1] <- 1
  for (j in 2:p) {
    # print(D_l[j, 1:(j-1)])
    k[j] <- sqrt(1 + sum(D_l[j, 1:(j - 1)]^2))
  }
  D_l <- solve(D_l, diag(k))
  D_l
}

# return Q * D_l
dependence_matrix <- function(p, rho, theta = NULL, Q = NULL) {
  Q_2d <- function(theta) {
    Q <- matrix(0, nrow = 2, ncol = 2)
    Q[1, 1] <- cos(theta)
    Q[2, 2] <- cos(theta)
    Q[1, 2] <- -sin(theta)
    Q[2, 1] <- sin(theta)
    Q
  }
  stopifnot(
    p - round(p) == 0, p > 1,
    "Please provide theta (p <= 3) or Q matrix, see ?precision_matrix_multivariate" = !is.null(theta) | !is.null(Q)
  )

  # compute D_l
  D_l <- D_l(p, rho)
  # compute Q
  if (p == 2) {
    stopifnot(
      "Length of theta should be 1 for p=2 case" = length(theta) == 1
    )
    Q <- Q_2d(theta)
  } else if (p == 3) {
    stopifnot(
      "Length of theta should be 3 for p=3 case" = length(theta) == 3
    )

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
#' rho <- c(-0.5, 0.5, -0.25) # correlation parameters
#' operator_list <- list(ar1(1:5, rho = 0.4), ar1(1:5, rho = 0.5), ar1(1:5, rho = 0.6))
#' precision_matrix_multivariate(3, operator_list, rho, theta = c(1, 2, 3))
precision_matrix_multivariate <- function(p,
                                          operator_list,
                                          rho,
                                          theta = NULL,
                                          Q = NULL,
                                          scale = NULL) {
  stopifnot("Please provide a list of p models" = length(operator_list) == p)
  stopifnot(
    "rho should be of length p(p-1)/2" =
      length(rho) == p * (p - 1) / 2
  )

  D <- dependence_matrix(p, rho, theta, Q)
  bigD <- kronecker(D, Matrix::Diagonal(length(operator_list[[1]]$h)))

  if (is.null(scale)) {
    scale <- rep(1, p)
  }
  K_list <- list()
  for (i in 1:p) {
    K_list[[i]] <- scale[i] * operator_list[[i]]$K
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
#' x <- seq(from = 0, to = 1, length.out = 40)
#' mesh <- fm_rcdt_2d_inla(lattice = fm_lattice_2d(x, x), extend = FALSE)
#' # Set parameters
#' p <- 3 # number of fields
#' rho <- c(-0.5, 0.5, -0.25) # correlation parameters
#' log_kappa <- list(2, 2, 2) # log(kappa)
#' variances <- list(1, 1, 1) # set marginal variances to 1
#' alpha <- list(2, 2, 2) # smoothness parameters
#' # Compute precision
#' Q <- precision_matrix_multivariate_spde(p,
#'   mesh = mesh, rho = rho,
#'   alpha = alpha, theta_K_list = log_kappa,
#'   variance_list = variances
#' )
#' # Plot the cross covariances
#' A <- as.vector(fm_basis(mesh, loc = matrix(c(0.5, 0.5), 1, 2)))
#' Sigma <- as.vector(solve(Q, c(A, rep(0, 2 * mesh$n))))
#' r11 <- Sigma[1:mesh$n]
#' r12 <- Sigma[(mesh$n + 1):(2 * mesh$n)]
#' r13 <- Sigma[(2 * mesh$n + 1):(3 * mesh$n)]
#' Sigma <- as.vector(solve(Q, c(rep(0, mesh$n), A, rep(0, mesh$n))))
#' r22 <- Sigma[(mesh$n + 1):(2 * mesh$n)]
#' r23 <- Sigma[(2 * mesh$n + 1):(3 * mesh$n)]
#' Sigma <- as.vector(solve(Q, v <- c(rep(0, 2 * mesh$n), A)))
#' r33 <- Sigma[(2 * mesh$n + 1):(3 * mesh$n)]
#'
#' proj <- fm_evaluator(mesh)
#'
#' par(mfrow = c(3, 3))
#' image.plot(fm_evaluate(proj, r11), main = "Cov(X_1(s0),X_1(s)")
#' plot.new()
#' plot.new()
#' image.plot(fm_evaluate(proj, r12), main = "Cov(X_1(s0),X_2(s)")
#' image.plot(fm_evaluate(proj, r22), main = "Cov(X_2(s0),X_2(s)")
#' plot.new()
#' image.plot(fm_evaluate(proj, r13), main = "Cov(X_1(s0),X_3(s)")
#' image.plot(fm_evaluate(proj, r23), main = "Cov(X_2(s0),X_3(s)")
#' image.plot(fm_evaluate(proj, r33), main = "Cov(X_3(s0),X_3(s)")
precision_matrix_multivariate_spde <- function(
    p,
    mesh,
    rho,
    alpha_list = NULL,
    theta_K_list = NULL,
    variance_list = NULL,
    B_K_list = NULL,
    theta = NULL,
    Q = NULL) {
  if (is.null(B_K_list)) {
    B_K_list <- lapply(1:p, function(i) NULL)
  } else {
    stopifnot(
      "Please provide a list of B_K matrices for each matern model" = length(B_K_list) == p
    )
  }

  if (is.null(theta_K_list)) {
    theta_K_list <- lapply(1:p, function(i) 0)
  } else {
    stopifnot(
      "Please provide a list of p parameters" = length(theta_K_list) == p
    )
  }

  if (is.null(alpha_list)) {
    alpha_list <- lapply(1:p, function(i) 2)
  } else {
    stopifnot(
      "Please provide a list of p alpha parameters" = length(alpha_list) == p
    )
  }

  if (is.null(variance_list)) {
    variance_list <- lapply(1:p, function(i) NULL)
    c <- NULL
  } else {
    stopifnot(
      "Please provide a list of p variance parameters" = length(variance_list) == p
    )
    if (mesh$manifold %in% c("R2", "S2")) {
      d <- 2
    } else if (mesh$manifold == "R1") {
      d <- 1
    }
    scale <- rep(0, p)
    for (i in 1:p) {
      scale[i] <- sqrt(gamma(alpha_list[[i]] - d / 2) / (variance_list[[i]] * (4 * pi)^(d / 2) * exp(theta_K_list[[i]][1])^(2 * (alpha_list[[i]] - d / 2)) * gamma(alpha_list[[i]])))
    }
  }

  operator_list <- NULL
  for (i in 1:p) {
    operator_list[[i]] <- matern(mesh, alpha_list[[i]], theta_K_list[[i]], B_K_list[[i]])
  }
  if (is.null(theta)) {
    theta <- rep(0, (p - 1) * p / 2)
  }
  precision_matrix_multivariate(p, operator_list, rho, theta, Q, scale = scale)
}
