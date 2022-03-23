# library(Matrix)

ngme.model.ar1 <- function(term) {
  x = eval(term)

  n <- length(x)
  # construct A
  A <- as(Matrix(diag(n_obs)), "dgCMatrix");

  # construct G
  G <- Matrix(diag(n_obs));
  G <- as(G, "dgCMatrix");

  # construct C
  C <- Matrix(0, n_obs, n_obs)
  C[seq(2, n_obs*n_obs, by=n_obs+1)] <- -1
  C <- as(C, "dgCMatrix");

  var_in <- list(type="ind_IG", v_init=list(nu=1))
  ope_in <- list(C=C, G=G, a_init=0.5)

  ar_in <- list(type="ar1",
                A=A,
                n_reg=n,
                opt_kappa     = T,
                opt_mu        = T,
                opt_sigma     = T,
                opt_var       = T,
                operator_in   = ope_in,
                var_in        = var_in
                )

  return (ar_in)
}


# # return a function of alpha
# # construct K(a)
# K <- function(a){
#   K_temp <- diag(n)
#   K_temp[seq(2, n*n, by=n+1)] <- -a
#   K_temp
# }
#
# # construct R
# R <- NULL
#
# # m(theta, v)
# m <- function(a) {
#   rep(0, n)
# }
#
# # construct the derivative of K
# dK <- function(a) {
#   K_temp <- matrix(data=0, nrow=n, ncol=n)
#   K_temp[seq(2, n*n, by=n+1)] <- -1
#   K_temp
# }
#
# d2K <- function(a) {
#   diag(n) * 0
# }
