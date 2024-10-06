# #' ngme general K operator
# #'
# #' K = sum_i f_i(theta_i) * matrices_i
# #' f_i depends on the parameter transformation type
# #'
# #' @param theta_K the parameter vector
# #' @param trans_type the type of parameter transformation
# #' @param matrices the matrices
# #' @param h the h vector
# #' @param zero_trace whether the trace of K is zero
# #' @param mesh the mesh
# #' @param ... ignore
# #'
# #' @return ngme_operator object
# #' @export
# general <- function(
#   theta_K, 
#   trans_type,
#   matrices,
#   h,
#   zero_trace = FALSE,
#   mesh = NULL,
#   ...
# )  {
#   theta <- param_trans_fun(theta_K, trans_type)
#   stopifnot(
#     "length of theta and matrices should be the same" = length(matrices) == length(theta),
#     "matrices should be of the same dimension" = 
#       all(sapply(matrices, nrow) == nrow(matrices[[1]])) &&
#       all(sapply(matrices, ncol) == ncol(matrices[[1]])),
#     "h should be of the same dimension as the matrices" = length(h) == nrow(matrices[[1]])
#   )
#   degree <- length(theta)
  
#   update_K <- function(theta_K) {
#     theta <- param_trans_fun(theta_K, trans_type)
#     K <- 0
#     for (i in 1:length(theta)) {
#       K <- K + theta[i] * matrices[[i]]
#     }
#     return(K)
#   }

#   K <- update_K(theta_K)
#   symmetric <- all(sapply(matrices, Matrix::isSymmetric))

#   ngme_operator(
#     model = "general",
#     mesh = mesh,
#     K = K,
#     h = h,
#     degree = degree,
#     theta_K = theta_K,
#     update_K = update_K,
#     matrices = matrices,
#     symmetric = symmetric,
#     zero_trace = zero_trace,
#     trans_type = trans_type,
#     param_name = paste("theta_K", seq_len(length(theta_K)), sep = " "),
#     param_trans = rep(list(identity), length(theta_K))
#   )
# }

# pre-defined parameter transformation functions
param_trans_fun <- function(theta_K, name) {
  th2a <- ar1_th2a

  if (name == "matern") {
    stopifnot(length(theta_K) == 1)
    return(c(exp(2*theta_K), 1))
  } else if (name == "ar1 x matern") {
    stopifnot(length(theta_K) == 2)
    # (aC1 + G1) x (kappa^2 C2 + G2)
    # matrices = (C1C2, C1G2, G1C2, G1G2)
    # theta = (a kappa^2, a, kappa^2, 1)
    # theta = (th2a(theta1) * exp(2theta2), th2a(theta1), exp(2theta2), 1)
    return(c(
      th2a(theta_K[1]) * exp(2*theta_K[2]),
      th2a(theta_K[1]),
      exp(2*theta_K[2]),
      1
    ))
  } else {
    stop("Unknown parameter transformation name")
  }
}


#' ngme general K operator
#'
#' K = sum_i f_i(theta_i) * matrices_i
#' f_i depends on the parameter transformation type
#'
#' @param theta_K the parameter vector
#' @param trans_type the type of parameter transformation
#' @param matrices the matrices
#' @param h the h vector
#' @param interact the interaction type, "multiply" or "kronecker"
#' @param theta_K2 the second parameter vector
#' @param trans2 the second parameter transformation type
#' @param matrices2 the second matrices
#' @param h2 the second h vector
#' @param zero_trace whether the trace of K is zero
#' @param mesh the mesh
#' @param ... ignore
general <- function(
  theta_K,
  trans,
  matrices,
  h,
  interact = NULL,
  theta_K2 = NULL,
  trans2 = NULL,
  matrices2 = NULL,
  h2 = NULL,
  zero_trace = FALSE,
  mesh = NULL,
  ...
) {
  stopifnot(
    "length of theta_K and trans should be the same" = length(theta_K) == length(trans),
    "theta_K and trans should be named" = !is.null(names(theta_K)) && !is.null(names(trans)),
    "theta_K and trans should have same names" = all(names(theta_K) == names(trans)),
    "h should be of the same dimension as the matrices" = length(h) == nrow(matrices[[1]])
  )

  if (is.null(interact)) {
    idx <- list()
    for (i in seq_len(length(theta_K))) {
      name <- names(theta_K)[i]
      idx[[name]] <- rep(FALSE, length(matrices))
      idx[[name]][[i]] <- TRUE
    }
    update_K <- function(theta_K) {
      coef <- compute_coef(theta_K, idx, trans)
      K <- 0
      for (i in 1:length(coef)) {
        K <- K + coef[i] * matrices[[i]]
      }
      return(K)
    }
  } else {
    param <- names(theta_K); param2 <- names(theta_K2)
    trans <- c(trans, trans2)
    
    if (length(param) < length(matrices)) {
      param <- c(param, rep(1, length(matrices) - length(param)))
    }
    if (length(param2) < length(matrices2)) {
      param2 <- c(param2, rep(1, length(matrices2) - length(param2)))
    }
    
    expanded_param <- expand.grid(param, param2)
    idx <- list()
    for (p in c(param, param2)) {
      # which rows contains p
      idx[[p]] <- apply(expanded_param, 1, function(x) any(x == p))
    }
    # expand matrices and combine with interaction
    matrices <- expand_matrices(matrices, matrices2, interact)

    update_K <- function(theta_K) {
      coef <- compute_coef(theta_K, idx, trans)
      K <- 0
      for (i in 1:length(coef)) {
        K <- K + coef[i] * matrices[[i]]
      }
      return(K)
    }
    
    theta_K <- c(theta_K, theta_K2); trans <- c(trans, trans2)
    for (name in names(theta_K)) {
      theta_K[name] <- name2fun(trans[[name]], inv=TRUE)(theta_K[name])
    }

    if (interact == "kronecker") h <- h %x% h2 else h <- h * h2
  }

  # remove `1` from idx
  idx[["1"]] <- NULL
  idx_mat <- as.matrix(as.data.frame(idx))
  
  ngme_operator(
    model = "general",
    mesh = mesh,
    K = update_K(theta_K),
    h = h,
    idx_mat = idx_mat,
    theta_K = theta_K,
    trans = trans,
    update_K = update_K,
    matrices = matrices,
    symmetric = all(sapply(matrices, Matrix::isSymmetric)),
    zero_trace = zero_trace,
    param_name = names(theta_K),
    param_trans = lapply(trans, name2fun)
  )
}

compute_coef <- function(theta_K, idx, trans) {
  idx <- as.data.frame(idx)
  nrow <- nrow(idx)
  coef <- rep(1, nrow)
  for (i in 1:nrow) {
    p <- apply(idx, 2, function(x) x[i])
    for (name in names(theta_K)) {
      if (p[name]) {
        coef[i] <- coef[i] * name2fun(trans[[name]], inv=FALSE)(theta_K[name])
      }
    }
  }
  return(coef)
}

name2fun <- function(trans, inv=FALSE) {
  if (!inv) {
    if (trans == "exp2") {
      return(function(x) exp(2*x))
    } else if (trans == "tanh") {
      return(ar1_th2a)
    } else if (trans == "identity") {
      return(function(x) x)
    } else {
      stop("Unknown transformation")
    }
  } else {
    if (trans == "exp2") {
      return(function(x) log(x) / 2)
    } else if (trans == "tanh") {
      return(ar1_a2th)
    } else if (trans == "identity") {
      return(function(x) x)
    } else {
      stop("Unknown transformation")
    }
  }
}

expand_matrices <- function(matrices, matrices2, interact) {
  result <- list()
  idx <- 1
  
  for (j in seq_along(matrices2)) {
    for (i in seq_along(matrices)) {
      # Multiply the matrices and store the result
      if (interact == "multiply") {
        result[[idx]] <- matrices[[i]] %*% matrices2[[j]]
      } else if (interact == "kronecker") {
        result[[idx]] <- kronecker(matrices[[i]], matrices2[[j]])
      }

      # Optionally, name the result for clarity
      names(result)[idx] <- paste0("m1_", i, "_m2_", j)
      # Increment the index counter
      idx <- idx + 1
    }
  }
  return(result)
}