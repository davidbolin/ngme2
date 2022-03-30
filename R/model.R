# usage: ngme.model.ar1(1:3, var="NIG)

# parameter:
#   x - numeric
#   var - string
# return:
#   list of stuff

ngme.model.ar1 <- function(x,
                           var="NIG",
                           init,
                           config)
  {
  # x is a numeric vector - covariates
  n <- length(x)
  # construct A
  A <- as(Matrix(diag(n)), "dgCMatrix");

  # construct G
  G <- Matrix(diag(n));
  G <- as(G, "dgCMatrix");

  # construct C

  C <- Matrix(0, n, n)
  C[seq(2, n*n, by=n+1)] <- -1
  C <- as(C, "dgCMatrix");

  var_in = list()
  if (var=="NIG") {
    var_in <- list(type="ind_IG", v_init=list(nu=1))
  }

  ope_in <- list(C=C, G=G, a_init=0.5, use_numerical=F)

  ar_in <- list(type="ar1",
                A=A,
                n_reg=n,
                operator_in   = ope_in,
                var_in        = var_in,
                opt_kappa     = TRUE,
                opt_mu        = TRUE,
                opt_sigma     = TRUE,
                opt_var       = TRUE
                )

  return (ar_in)
}
