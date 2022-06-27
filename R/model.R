# usage: ngme.model.ar1(1:3, var="NIG)

# parameter:
#   x - numeric
#   var - string
# return:
#   list of stuff

ngme.ar1 <- function(
  x,
  alpha=0.5,
  use_num_dK = FALSE
) {
  # construct A
  A <- ngme.ts.make.A(x)

  x <- unique(x); n <- length(x)

  # construct G
    G <- Matrix::Matrix(diag(n));
    G <- as(G, "dgCMatrix");

  # construct C
    C <- Matrix::Matrix(0, n, n)
    C[seq(2, n*n, by=n+1)] <- -1
    C <- as(C, "dgCMatrix");

  ar1_in <- list(
    A     = A,
    alpha = alpha,
    operator_in = list(
      n_params    = 1,
      C           = C,
      G           = G,
      use_num_dK  = use_num_dK
    )
  )

  class(ar1_in) <- "ngme.ar1"
  return (ar1_in)
}
