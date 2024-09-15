# implementation of the general K operator
# K = sum_i f_i(theta_i) * matrices_i

general_ope <- function(
  theta_K, 
  theta_trans,
  matrices,
  h,
  zero_trace = FALSE,
  mesh = NULL,
  ...
)  {
  theta <- param_trans_fun(theta_K, theta_trans)
  stopifnot(
    "length of theta and matrices should be the same" = length(matrices) == length(theta),
    "matrices should be of the same dimension" = 
      all(sapply(matrices, nrow) == nrow(matrices[[1]])) &&
      all(sapply(matrices, ncol) == ncol(matrices[[1]])),
    "h should be of the same dimension as the matrices" = length(h) == nrow(matrices[[1]])
  )
  degree <- length(theta)
  
  update_K <- function(theta_K) {
    theta <- param_trans_fun(theta_K, theta_trans)
    K <- 0
    for (i in 1:length(theta)) {
      K <- K + theta[i] * matrices[[i]]
    }
    return(K)
  }

  K <- update_K(theta_K)
  symmetric <- all(sapply(matrices, Matrix::isSymmetric))

  ngme_operator(
    model = "general",
    mesh = mesh,
    K = K,
    h = h,
    degree = degree,
    theta_K = theta_K,
    update_K = update_K,
    matrices = matrices,
    symmetric = symmetric,
    zero_trace = zero_trace,
    theta_trans = theta_trans,
    param_name = paste("theta_K", seq_len(length(theta_K)), sep = " "),
    param_trans = rep(list(identity), length(theta_K))
  )
}

# pre-defined parameter transformation functions
param_trans_fun <- function(theta, name) {
  if (name == "matern") {
    return(c(exp(2*theta), 1))
  } else {
    stop("Unknown parameter transformation name")
  }
}