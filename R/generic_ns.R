#' Non-stationary precision matrix operator with custom matrix combinations
#'
#' Creates a flexible non-stationary precision matrix K by allowing both linear combinations 
#' and matrix multiplications with parameter-dependent diagonal matrices. This enables modeling
#' spatially varying parameters through basis expansions.
#'
#' The key distinction from `generic` is that `generic_ns`:
#' 1. Requires explicit position specification for matrix combinations
#' 2. Converts parameters to diagonal matrices for multiplication with fixed matrices
#' 3. Allows for basis expansions via B_theta_K
#'
#' @param theta_K List of parameter vectors (in real scale)
#' @param matrices List of fixed matrices to be used in the model
#' @param position List of vectors specifying matrix combinations (required)
#' @param h The integration weights vector
#' @param B_theta_K List of basis matrices for parameters, defaults to matrices of 1s if not provided
#' @param trans List of transformations for each parameter (one transformation per parameter)
#' @param mesh The mesh object
#' @param zero_trace Whether the trace of K should be zero
#' @param ... Additional arguments (ignored)
#' 
#' @details 
#' The `generic_ns` operator constructs K through the following steps:
#' 
#' 1. Create diagonal matrices for each parameter using basis expansions: D_param = diag(B_param * theta_param)
#' 2. Combine parameter matrices with fixed matrices according to position specifications
#' 3. Sum all the resulting combinations
#' 
#' The `position` parameter defines how matrices are combined:
#' - Each element of the list represents a term in the sum
#' - Each vector contains indices of matrices to multiply (parameter matrices first, then fixed matrices)
#' - For example: `position = list(c(1, 3), c(2, 4))` with one parameter means:
#'   Multiply the parameter diagonal matrix (index 1) with fixed matrix 1 (index 3),
#'   then add parameter diagonal matrix 2 (index 2) multiplied by fixed matrix 2 (index 4)
#' 
#' Spatially-varying parameters can be created using basis matrices in B_theta_K.
#' 
#' @examples
#' \dontrun{
#' # Simple example with one parameter
#' n <- 5
#' A <- matrix(1, n, n)
#' B <- matrix(2, n, n)
#' alpha <- 0.4
#' 
#' # Create a simple generic_ns model: D_alpha * A + B
#' model <- generic_ns(
#'   theta_K = list(alpha = c(alpha)),
#'   matrices = list(A, B),
#'   position = list(c(1, 2), c(3)),  # D_alpha * A + B
#'   h = rep(1, n)
#' )
#' 
#' # AR1 model representation: rho * C + G
#' ar1_obj <- ar1(1:10, rho = 0.5)
#' g <- name2fun("tanh", inv = TRUE)
#' generic_ar1 <- generic_ns(
#'   theta_K = list(x = c(g(0.5))),  # Transform from real scale
#'   trans = list(x = "tanh"),        # Apply tanh transformation
#'   matrices = list(ar1_obj$C, ar1_obj$G),
#'   position = list(c(1, 2), c(3)),  # D_rho * C + G
#'   h = ar1_obj$h,
#'   mesh = 1:10
#' )
#' 
#' # Spatially varying Matern model
#' # Create a simple 1D mesh
#' mesh <- fmesher::fm_mesh_1d(seq(0, 1, length.out = 10))
#' 
#' # Create basis for spatially-varying kappa
#' n <- 10
#' B_kappa <- matrix(0, n, 2)
#' B_kappa[1:(n/2), 1] <- 1  # First half of the domain
#' B_kappa[(n/2 + 1):n, 2] <- 1  # Second half of the domain
#' 
#' # Create standard Matern components
#' matern_model <- matern(mesh)
#' 
#' # Create model with space-varying kappa
#' ns_model <- generic_ns(
#'   theta_K = list(kappa = c(log(1), log(2))),  # Different values for each region
#'   matrices = list(matern_model$C, matern_model$G),
#'   B_theta_K = list(kappa = B_kappa),
#'   trans = list(kappa = "exp2"),  # kappa^2 transformation for C
#'   h = matern_model$h,
#'   position = list(c(1, 2), c(3)),  # D_kappa * C + G
#'   mesh = mesh
#' )
#' }
#'
#' @export
generic_ns <- function(
    theta_K,
    matrices,
    position,
    h,
    model = "generic_ns",
    B_theta_K = NULL,
    trans = NULL,
    mesh = NULL,
    zero_trace = FALSE,
    param_name = NULL,
    param_trans = NULL,
    ...
) {
  if (!is.null(mesh)) mesh <- ngme_build_mesh(mesh)
  
  # Validate inputs
  stopifnot(
    "theta_K should be a list" = is.list(theta_K),
    "matrices should be a list" = is.list(matrices),
    "matrices should not be NULL" = !is.null(matrices),
    "matrices should not be a list of NULL" = all(!sapply(matrices, is.null)),
    "all matrices should have the same dimensions" = all(sapply(matrices, function(m) nrow(m) == nrow(matrices[[1]])) &
                                                       sapply(matrices, function(m) ncol(m) == ncol(matrices[[1]]))),
    "h should be of the same dimension as the matrices" = length(h) == nrow(matrices[[1]]),
    "position is required" = !is.null(position) && is.list(position)
  )
  
  # Initialize transformations if not provided
  if (is.null(trans)) {
    trans <- list()
    for (theta_K_name in names(theta_K)) {
      trans[[theta_K_name]] <- "identity"
    }
  }
  
  # Validate transformation list
  stopifnot(
    "trans should be a list" = is.list(trans),
    "trans should have names matching theta_K" = all(names(trans) %in% names(theta_K))
  )
  
  # Initialize B_theta_K with matrices of 1s if not provided
  if (is.null(B_theta_K)) {
    B_theta_K <- list()
    n <- nrow(matrices[[1]])
    for (theta_K_name in names(theta_K)) {
      # Create a matrix of 1s with appropriate dimensions
      # Use a single column matrix since most parameters are scalar
      B_theta_K[[theta_K_name]] <- matrix(1, nrow = n, ncol = 1)
    }
  }
  
  # Validate B_theta_K
  stopifnot(
    "B_theta_K should be a list" = is.list(B_theta_K),
    "B_theta_K should have names matching theta_K" = all(names(B_theta_K) %in% names(theta_K))
  )
  
  # Process basis matrices for parameters
  param_matrices <- list()
  for (theta_K_name in names(B_theta_K)) {
    B <- B_theta_K[[theta_K_name]]
    theta <- theta_K[[theta_K_name]]
    
    # Validate dimensions
    if (length(theta) > 1 && ncol(B) != length(theta)) {
      stop(sprintf("B_theta_K[[%s]] should have ncol matching length of theta_K[[%s]]", theta_K_name, theta_K_name))
    }
    
    # Apply transformation to create parameter matrix
    trans_type <- trans[[theta_K_name]]
    transformed_value <- name2fun(trans_type)(B %*% theta)
    
    if (trans_type %in% c("identity", "exp", "exp2", "exp4", "sqrt", "square", "log", "tanh")) {
      # For element-wise transformations, create diagonal matrix
      param_matrices[[theta_K_name]] <- Matrix::Diagonal(x = as.numeric(transformed_value))
    } else {
      # For other transformations (like tanh for correlations), create the matrix directly
      param_matrices[[theta_K_name]] <- transformed_value
    }
  }
  
  # Combine parameters first, then fixed matrices
  all_matrices <- c(param_matrices, matrices)
  
  # Create flat parameter vector and keep track of mappings
  flat_theta_K <- numeric(0)
  param_map <- list()
  idx <- 1
  
  for (theta_K_name in names(theta_K)) {
    theta_K_len <- length(theta_K[[theta_K_name]])
    param_map[[theta_K_name]] <- idx:(idx + theta_K_len - 1)
    flat_theta_K <- c(flat_theta_K, theta_K[[theta_K_name]])
    idx <- idx + theta_K_len
  }
  names(flat_theta_K) <- rep(names(theta_K), sapply(theta_K, length))
  
  # Extract parameter transformation functions
  if (is.null(param_trans)) {
    param_trans <- list()
    for (theta_K_name in names(trans)) {
      param_trans[[theta_K_name]] <- name2fun(trans[[theta_K_name]])
    }
  }
  if (is.null(param_name)) {
    param_name <- names(theta_K)
  }
  
  # Helper function to compute K (theta_K is a flat vector)
  update_K <- function(theta_K) {
    # Initialize K as zero matrix
    n <- nrow(matrices[[1]])
    K <- Matrix::Matrix(0, n, n)
    
    # Create parameter matrices using basis expansions
    param_matrices <- list()
    all_matrices <- list()
    
    # First, create parameter matrices
    for (theta_K_name in names(param_map)) {
      indices <- param_map[[theta_K_name]]
      
      # Get basis matrix
      B <- if (!is.null(B_theta_K[[theta_K_name]])) {
        B_theta_K[[theta_K_name]]
      } else {
        # Default to matrix of ones (1 column)
        Matrix::Matrix(1, nrow = n, ncol = 1)
      }
      
      # Extract parameter values from theta_K
      param_values <- theta_K[indices]
      
      # Apply basis expansion
      expanded_values <- as.vector(B %*% param_values)
      
      # Apply transformation
      trans_type <- if (!is.null(trans[[theta_K_name]])) trans[[theta_K_name]] else "identity"
      transformed_values <- name2fun(trans_type)(expanded_values)
      
      # Create diagonal matrix from expanded values
      param_matrices[[theta_K_name]] <- Matrix::Diagonal(x = as.numeric(transformed_values))
    }
    
    # Create combined list of matrices (parameters first, then fixed matrices)
    all_matrices <- c(param_matrices, matrices)
    
    # Apply matrix combinations according to position
    for (pos in position) {
      # Start with identity matrix
      term <- Matrix::Diagonal(n)
      
      # Multiply matrices in the combination
      for (idx in pos) {
        if (idx > 0 && idx <= length(all_matrices)) {
          term <- term %*% all_matrices[[idx]]
        } else {
          stop(sprintf("Position index %d is out of bounds (max: %d)", 
                      idx, length(all_matrices)))
        }
      }
      
      # Add to K
      K <- K + term
    }
    
    return(K)
  }
  
  # Create operator
  ngme_operator(
    model = model, 
    mesh = mesh,
    K = ngme_as_sparse(update_K(theta_K = flat_theta_K)),
    generic_type = "generic_ns",
    h = h,
    theta_K = flat_theta_K,
    trans = trans,  # Include the original transforms
    B_theta_K = B_theta_K,
    param_map = param_map,
    position = position,
    matrices = sapply(matrices, function(m) ngme_as_sparse(m)),  # Store original matrices, not all_matrices
    update_K_func = update_K,
    symmetric = all(sapply(all_matrices, Matrix::isSymmetric)),
    zero_trace = zero_trace,
    param_name = param_name,
    param_trans = param_trans,
    ... 
  )
} 