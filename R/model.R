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

  if (!is.null(replicates)) {
    rep_tmp <- unique_rep[1]
    index.tmp <- index[replicates==rep_tmp]
    A = ngme.ts.make.A(index.tmp, index)

    if(nrep > 1){
      for (j in 2:nrep){
        index.tmp <- index[replicates==unique_rep[j]]
        Atmp <- ngme.ts.make.A(index.tmp, index)
        A <- Matrix::bdiag(A, Atmp)
      }
    }
  }

  x <- unique(index); n <- length(x)

  # construct G
    G <- Matrix::Matrix(diag(n));

  # construct C
    C <- Matrix::Matrix(0, n, n)
    C[seq(2, n*n, by=n+1)] <- -1

    if(!is.null(nrep)){
      C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), C)
      G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), G)
    }

  # G <- as(G, "dgCMatrix");

# convert C and G
C <- as(C, "dgCMatrix");

G = as(G, "dgTMatrix")
idx <- which(G@i <= G@j)
G = Matrix::sparseMatrix(i=G@i[idx], j = G@j[idx], x= G@x[idx],
                         symmetric=TRUE, index1 = FALSE)
G = as(G, "dgCMatrix")


print("C = ")
print(C)
print("G = ")
print(G)
print("A = ")
print(A)

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
