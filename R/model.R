# usage: ngme.model.ar1(1:3, var="NIG)

# parameter:
#   x - numeric
#   var - string
# return:
#   list of stuff

ngme.ar1 <- function(
  index,
  replicates,
  alpha=0.5,
  use_num_dK = FALSE,
  nrep=NULL
) {
  # construct A
  unique_rep = unique(replicates)
  nrep = length(unique_rep)

  A <- ngme.ts.make.A(index, replicates)

  # e.g. index = c(1:200, 1:100)
  #      replicates = c(rep(1, 200), rep(2, 100))
  #      n =200 in this case
  n <- length(unique(index)) 

  # construct G
    G <- Matrix::Matrix(diag(n));

  # construct C
    C <- Matrix::Matrix(0, n, n)
    C[seq(2, n*n, by=n+1)] <- -1

    if(!is.null(nrep)){
      C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), C)
      G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), G)
    }

  ar1_in <- list(
    A     = A,
    alpha = alpha,
    operator_in = list(
      n_params    = 1,
      C           = ngme.as.sparse(C),
      G           = ngme.as.sparse(G),
      use_num_dK  = use_num_dK
    )
  )

  class(ar1_in) <- "ngme.ar1"
  return (ar1_in)
}
