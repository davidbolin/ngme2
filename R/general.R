# implementation of the general K operator
# K = sum_i f_i(theta_i) * matrices_i

general_ope <- function(
  theta_K, 
  trans_type,
  matrices,
  h,
  zero_trace = FALSE,
  mesh = NULL,
  ...
)  {
  theta <- param_trans_fun(theta_K, trans_type)
  stopifnot(
    "length of theta and matrices should be the same" = length(matrices) == length(theta),
    "matrices should be of the same dimension" = 
      all(sapply(matrices, nrow) == nrow(matrices[[1]])) &&
      all(sapply(matrices, ncol) == ncol(matrices[[1]])),
    "h should be of the same dimension as the matrices" = length(h) == nrow(matrices[[1]])
  )
  degree <- length(theta)
  
  update_K <- function(theta_K) {
    theta <- param_trans_fun(theta_K, trans_type)
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
    trans_type = trans_type,
    param_name = paste("theta_K", seq_len(length(theta_K)), sep = " "),
    param_trans = rep(list(identity), length(theta_K))
  )
}

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