#' Generic precision matrix operator
#'
#' Creates a flexible precision matrix K by allowing linear combinations of base matrices
#' with parameter-dependent coefficients. The model follows the form:
#' K = sum_i f_i(theta) * matrices_i, where f_i depends on the parameter transformation.
#'
#' This function can be used to represent various models like AR1, Matern, and RW1 in a unified framework.
#'
#' @param theta_K The parameter vector (in real scale). If missing, initialized as double(0)
#' @param matrices A list of matrices to be combined with parameter-dependent coefficients
#' @param h The integration weights vector
#' @param trans Transformation for each parameter. Must be a list where each element
#'   is a character vector defining transformations for each matrix.
#'   For example, `trans=list(kappa=c("exp2", "identity", "null"))` means parameter 'kappa'
#'   affects matrix 1 with exp2 transformation, matrix 2 with identity, and doesn't affect matrix 3.
#'   Available transformations: "identity", "exp", "exp2", "exp4", "sqrt", "square", "log", "tanh"
#' @param mesh The mesh object
#' @param zero_trace Whether the trace of K should be zero
#' @param ... Additional arguments (ignored)
#' 
#' @details 
#' The `generic` function provides a flexible way to construct precision matrices by combining 
#' existing matrices with transformed parameters. Here's how it works:
#' 
#' 1. Each parameter in `theta_K` can affect one or more matrices in the `matrices` list
#' 2. The `trans` list defines how each parameter transforms before multiplying each matrix
#' 3. For each matrix, the coefficients from all parameters are multiplied together
#' 4. The final K is the sum of all transformed matrices
#' 
#' Common use cases:
#' 
#' - **AR1 model**: K = rho * C + G (where rho is between -1 and 1)
#' - **Matern (alpha=2)**: K = kappa^2 * C + G (where kappa > 0)
#' - **Matern (alpha=4)**: K = kappa^4 * C + 2*kappa^2 * G + G * Cinv * G
#' 
#' @examples
#' \dontrun{
#' # AR1 model with rho = 0.5
#' ar1_obj <- ar1(1:10, rho = 0.5)
#' g <- name2fun("tanh", inv = TRUE)
#' generic_ar1 <- generic(
#'   theta_K = c(rho = g(0.5)),  # Transform from real scale
#'   trans = list(rho = c("tanh", "null")),  # Apply tanh to first matrix
#'   matrices = list(ar1_obj$C, ar1_obj$G),
#'   h = ar1_obj$h
#' )
#' 
#' # Matern model with kappa = 2
#' mesh <- fmesher::fm_mesh_1d(seq(0, 1, length.out = 10))
#' matern_obj <- matern(mesh)
#' log_kappa <- log(2)
#' generic_matern <- generic(
#'   theta_K = c(kappa = log_kappa),
#'   trans = list(kappa = c("exp2", "null")),  # exp(2*kappa) for C, nothing for G
#'   matrices = list(matern_obj$C, matern_obj$G),
#'   h = matern_obj$h
#' )
#' 
#' # Matern model with alpha = 4 and kappa = 2
#' mesh <- fmesher::fm_mesh_2d(cbind(x = runif(20), y = runif(20)))
#' matern_obj <- matern(mesh, alpha = 4)
#' C <- matern_obj$C
#' G <- matern_obj$G
#' Cinv <- C; diag(Cinv) <- 1 / Matrix::diag(C)
#' 
#' generic_matern4 <- generic(
#'   theta_K = c(kappa = log(2)),
#'   trans = list(kappa = c("exp4", "exp2", "null")),
#'   matrices = list(C, 2*G, G %*% Cinv %*% G),
#'   h = matern_obj$h
#' )
#' }
#'
#' @export
generic <- function(
    theta_K,
    matrices,
    h,
    trans = NULL,
    mesh = NULL,
    zero_trace = FALSE,
    ...) {
  if (!is.null(mesh)) mesh <- ngme_build_mesh(mesh)

  # If theta_K is not provided, initialize it as double(0)
  if (missing(theta_K)) {
    theta_K <- double(0)
  }

  # Validate matrices
  stopifnot(
    "matrices should not be NULL" = !is.null(matrices),
    "matrices should be a list" = is.list(matrices),
    "matrices should not be a list of NULL" = all(!sapply(matrices, is.null)),
    "all matrices should have the same dimensions" = all(sapply(matrices, function(m) nrow(m) == nrow(matrices[[1]])) &
      sapply(matrices, function(m) ncol(m) == ncol(matrices[[1]]))),
    "h should be of the same dimension as the matrices" = length(h) == nrow(matrices[[1]])
  )

  # Convert old format to new format if needed
  if (is.null(trans) && length(theta_K) > 0) {
    # If trans is NULL, create default trans with only 1s for the identity matrices
    trans <- list()
    for (i in seq_along(theta_K)) {
      trans[[i]] <- rep("null", length(matrices))
      trans[[i]][i] <- "identity"
    }
    # If no parameters provided, just use default behavior (all matrices with coefficient 1)
  } else if (is.character(trans)) {
    # Convert old format to new format
    old_trans <- trans
    trans <- list()
    for (i in seq_along(theta_K)) {
      param_name <- names(theta_K)[i]
      trans[[param_name]] <- rep("null", length(matrices))
      trans[[param_name]][i] <- old_trans[i]
    }
  }

  # Validate the structure of trans
  if (length(theta_K) > 0) {
    # Only validate trans if we have parameters
    for (param_name in names(trans)) {
      if (length(trans[[param_name]]) != length(matrices)) {
        stop(sprintf(
          "For parameter '%s', the transformation vector length (%d) must match the number of matrices (%d)",
          param_name, length(trans[[param_name]]), length(matrices)
        ))
      }
    }
  }

  # Helper function to compute the coefficient of each matrix
  update_K <- function(theta_K) {
    coef <- compute_coef_new(theta_K, trans, length(matrices))

    K <- 0
    for (i in seq_along(matrices)) {
      K <- K + coef[i] * matrices[[i]]
    }
    return(K)
  }

  ngme_operator(
    model = "generic",
    mesh = mesh,
    K = ngme_as_sparse(update_K(theta_K)),
    generic = TRUE,
    h = h,
    theta_K = theta_K,
    trans = trans,
    update_K = update_K,
    matrices = sapply(matrices, ngme_as_sparse),
    symmetric = all(sapply(matrices, Matrix::isSymmetric)),
    zero_trace = zero_trace,
    param_name = names(theta_K),
    param_trans = get_param_trans(trans)
  )
}

# Helper function to get param_trans from trans list
# For each parameter, if only one transformation type is used (excluding null),
# use that transformation function, otherwise use identity
get_param_trans <- function(trans) {
  if (is.null(trans) || length(trans) == 0) {
    return(NULL)
  }
  
  result <- list()
  
  for (param_name in names(trans)) {
    # Get all transformations for this parameter (excluding "null")
    trans_types <- unique(trans[[param_name]][trans[[param_name]] != "null"])
    
    if (length(trans_types) == 1) {
      # If only one transformation type is used, use that
      result[[param_name]] <- name2fun(trans_types)
    } else {
      # If multiple or no transformation types are used, use identity
      result[[param_name]] <- name2fun("identity")
    }
  }
  
  return(result)
}

# Calculate coefficients for each matrix based on parameters and transformations
compute_coef_new <- function(theta_K, trans, n_matrices) {
  coef <- rep(1, n_matrices)

  # If no parameters, return all ones (default behavior)
  if (length(theta_K) == 0) {
    return(coef)
  }

  for (param_name in names(theta_K)) {
    if (param_name %in% names(trans)) {
      for (i in seq_along(trans[[param_name]])) {
        trans_type <- trans[[param_name]][i]
        if (trans_type != "null") {
          # Apply the transformation function to the parameter value
          transformed_value <- name2fun(trans_type)(theta_K[param_name])
          # Multiply the coefficient by the transformed parameter value
          coef[i] <- coef[i] * transformed_value
        }
      }
    }
  }

  return(coef)
}

#' Convert transformation name to function
#'
#' This function converts a transformation name to the corresponding function.
#' Convertion is from original scale to real scale by default.
#' Available transformations:
#' \itemize{
#'   \item \code{exp4}: $\exp(4x)$, inverse is $\log(x)/4$
#'   \item \code{exp2}: $\exp(2x)$, inverse is $\log(x)/2$
#'   \item \code{tanh}: Hyperbolic tangent transformation used for AR1 parameter, uses ar1_th2a and ar1_a2th
#'   \item \code{identity}: Identity function, no transformation
#'   \item \code{exp}: $\exp(x)$, inverse is $\log(x)$
#'   \item \code{sqrt}: $\sqrt(x)$, inverse is $x^2$
#'   \item \code{square}: $x^2$, inverse is $\sqrt(x)$
#'   \item \code{log}: $\log(x)$, inverse is $\exp(x)$
#' }
#'
#' @param trans Character string specifying the transformation type from original scale to real scale
#' @param inv Logical, if TRUE returns the transformation function from real scale to original scale
#'
#' @return A function that applies the specified transformation
#' @export
name2fun <- function(trans, inv = FALSE) {
  # Define transformation functions in a list
  transformations <- list(
    exp4 = list(
      forward = function(x) exp(4 * x),
      inverse = function(x) log(x) / 4
    ),
    exp2 = list(
      forward = function(x) exp(2 * x),
      inverse = function(x) log(x) / 2
    ),
    tanh = list(
      forward = ar1_th2a,  # restrict the range of rho to (-1, 1)
      inverse = ar1_a2th
    ),
    identity = list(
      forward = function(x) x,
      inverse = function(x) x
    ),
    exp = list(
      forward = function(x) exp(x),
      inverse = function(x) log(x)
    ),
    sqrt = list(
      forward = function(x) sqrt(x),
      inverse = function(x) x^2
    ),
    square = list(
      forward = function(x) x^2,
      inverse = function(x) sqrt(x)
    ),
    log = list(
      forward = function(x) log(x),
      inverse = function(x) exp(x)
    )
  )
  
  # Check if transformation exists
  if (!trans %in% names(transformations)) {
    stop("Unknown transformation: ", trans)
  }
  
  # Return the appropriate function
  if (!inv) {
    return(transformations[[trans]]$forward)
  } else {
    return(transformations[[trans]]$inverse)
  }
}
