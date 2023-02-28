  load_all()
  n <<- 10
  z1 = arima.sim(n, model = list(ar = 0.5), sd = 0.5)
  z2 = arima.sim(n, model = list(ar = 0.5), sd = 0.5)

  out2 <- ngme(
    z ~ f(c(1:n ,1:n), model="rw1", replicate=rep(1:2,each=n), noise=noise_normal()),
    data = list(z = c(z1, z2)),
    control = ngme_control(
      estimation = T,
      iterations = 20,
      n_parallel_chain = 1
    ),
    debug = TRUE
  )


check()

# rw1

n = 5


# order = 1
C <- Matrix::sparseMatrix(i = 1:(n-1), j=2:n, x=1, dims=c(n-1,n)); C
G <- Matrix::sparseMatrix(i = 1:(n-1), j=1:(n-1), x=-1, dims=c(n-1,n)); G

# order = 2
a <- model_rw(1:10000000, order=1)

C <- Matrix::sparseMatrix(i = 1:(n-2), j=2:(n-1), x=-2, dims=c(n-2,n)); C
G <- Matrix::sparseMatrix(i = rep(1:(n-2),2), j=c(1:(n-2), 3:n), x=1, dims=c(n-2,n)); G
C+G

# order = 1, circular
model_rw(1:5, circular=T, order=1)$K
n=5
C <- Matrix::Diagonal(n); C
G <- Matrix::sparseMatrix(i = 1:n, j=c(2:n, 1), x=-1, dims=c(n,n)); G
C+G


# order = 2, circular
model_rw(1:5, circular=T, order=2)$K
C <- Matrix::sparseMatrix(i = 1:n, j=c(2:n,1), x=-2, dims=c(n,n))
G <- Matrix::sparseMatrix(i = rep(1:n,2), j=c(1:n, 3:n, 1, 2), x=1, dims=c(n,n))
C+G
