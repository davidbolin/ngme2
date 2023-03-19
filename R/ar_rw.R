# new models should have
# model_name <- function(
# index, replicate, index_NA, noise, data
# )

# ar1 and rw

#' ngme AR(1) model specification
#'
#' Generating C, G and A given index and replicate
#'
#' @param map integer vector, time index for the AR(1) process
#' @param replicate replicate for the process
#' @param index_NA Logical vector, same as is.na(response var.)
#' @param alpha initial value for alpha
#'
#' @param noise noise, can be specified in f()
#' @param data data, can be specified in f(), ngme()
#' @param ... extra arguments in f()
#'
#' @return a list of specification of model
#' @export
#'
#' @examples
#' model_ar1(c(1:3, 1:3), replicate = c(1,1,1,2,2,2))
#' f(xx, model = "ar1", data=list(xx = c(2,4,5)), noise=noise_nig())
model_ar1 <- ar1 <- function(
  map,   # time index
  replicate   = NULL,
  index_NA    = NULL,
  data        = NULL,
  noise       = noise_normal(),
  alpha       = 0.5,
  ...
) {
  # capture symbol in index
  index <- eval(substitute(map), envir = data, enclos = parent.frame())
  stopifnot("The index should be integers." = all(index == round(index)))
  if (is.null(replicate)) replicate <- rep(1, length(index))
  if (is.null(index_NA)) index_NA <- rep(FALSE, length(index))

  stopifnot("Make sure length(idx)==length(replicate)" = length(index) == length(replicate))

  replicate <- if (!is.null(list(...)$replicate)) list(...)$replicate
    else rep(1, length(map))

  # e.g. index      = 1 2 3 1 2 3 4
  #      replicate  = 1 1 1 2 2 2 2
  mesh <- INLA::inla.mesh.1d(min(index):max(index))
  n <- mesh$n
  nrep <- length(unique(replicate))
  nrep <- 1

  # construct G
  G <- Matrix::Diagonal(n);
  C <- Matrix::sparseMatrix(j=1:(n-1), i=2:n, x=-1, dims=c(n,n))

  # update noise with length n
  if (noise$n_noise == 1) noise <- update_noise(noise, n = n)
  # make A and A_pred
  tmp <- ngme_make_A(
    mesh = mesh,
    map = index,
    n_map = length(index),
    idx_NA = index_NA,
    replicate = replicate
  )
  A <- tmp$A; A_pred <- tmp$A_pred

  stopifnot(nrep == ncol(A) / ncol(C))

  # remove duplicate symbol in ... (e.g. theta_K)
  args <- within(list(...), {
    mesh        = mesh
    model       = "ar1"
    theta_K     = if (exists("theta_K")) ar1_a2th(theta_K) else ar1_a2th(alpha)
    W_size      = n
    V_size      = n
    A           = A
    A_pred      = A_pred
    h           = noise$h
    C           = ngme_as_sparse(C)
    G           = ngme_as_sparse(G)
    K           = alpha * C + G
    noise       = noise
    map         = index
    n_map       = length(index)
    replicate   = replicate
    replicate   = replicate
    n_rep       = nrep
  })
  do.call(ngme_model, args)
}

#' ngme model - random walk of order 1
#'
#' Generating C, G and A given index and replicate
#' size of C and G is (n-1) * n, size of V is n-1
#'
#' @param map        numerical vector, covariates to build index for the process
#' @param order     1 or 2, order of random walk model
#' @param replicate replicate for the process
#' @param data      specifed or inherit from ngme formula
#' @param circular  whether the mesh is circular, i.e. the first one is connected to the last
#'   if it is circular, we will treat the 1st location and the last location same.
#' @param index_NA  Logical vector, same as is.na(response variable)
#' @param noise     1. string: type of model, 2. ngme.noise object
#'  (can also be specified in each ngme model)
#' @param ...       additional arguments
#'
#' @return a list
#' @export
#'
#' @examples
#' r1 <- model_rw(1:7, order = 1, circular = TRUE); r1$C + r1$G
#' r2 <- model_rw(1:7, order = 1); r2$C + r2$G
model_rw <- rw <- function(
  map,
  order       = 1,
  replicate  = NULL,
  data        = NULL,
  circular    = FALSE,
  index_NA    = NULL,
  noise       = noise_normal(),
  # extra A matrix
  ...
) {
  stopifnot(order == 1 || order == 2)
# capture symbol in index
  x <- eval(substitute(map), envir = data, enclos = parent.frame())
  stopifnot("length of index at least >= 3" = length(x) >= 3)
  if (is.null(replicate)) replicate <- rep(1, length(x))
  if (is.null(index_NA)) index_NA <- rep(FALSE, length(x))
  stopifnot("Make sure length(x)==length(replicate)" = length(x) == length(replicate))

  # create mesh using index
  # sorted_index <- sort(index, index.return = TRUE)
  # h <- diff(sorted_index$x)
  # permutation matrix, same as A <- diag(length(index))[sorted_index$ix, ]

  # create A in one way
  # A <- Matrix::sparseMatrix(seq_along(index)[!index_NA], sorted_index$ix, x = 1)
  # if (any(index_NA))
  #   A_pred <- Matrix::sparseMatrix(seq_along(index)[index_NA], sorted_index$ix, x = 1)
  # else
  #   A_pred <- NULL

  # K is of size (n-1) * n
  # then mesh is of size (n-1)
  # lose 1 location (first) for non-circular rw1 case

  if (!circular) {
    # regular mesh
    if (order == 1) {
      mesh <- INLA::inla.mesh.1d(loc = unique(x)[-1])
      n <- mesh$n + 1
      C <- Matrix::sparseMatrix(i = 1:(n-1), j=2:n, x=1, dims=c(n-1,n))
      G <- Matrix::sparseMatrix(i = 1:(n-1), j=1:(n-1), x=-1, dims=c(n-1,n))
    } else if (order == 2) {
      mesh <- INLA::inla.mesh.1d(loc = unique(x)[-c(1,2)])
      n <- mesh$n + 2
      stopifnot(n >= 2)
      C <- Matrix::sparseMatrix(i = 1:(n-2), j=2:(n-1), x=-2, dims=c(n-2,n))
      G <- Matrix::sparseMatrix(i = rep(1:(n-2),2), j=c(1:(n-2), 3:n), x=1, dims=c(n-2,n))
    }
  } else {
    # circular mesh (remove last element)
      stopifnot(length(x) >= 4)
    mesh <- INLA::inla.mesh.1d(loc = x)
    n <- mesh$n
    if (order == 1) {
      C <- Matrix::Diagonal(n)
      G <- Matrix::sparseMatrix(i = 1:n, j=c(2:n, 1), x=-1, dims=c(n,n))
    } else if (order == 2) {
      C <- Matrix::sparseMatrix(i = 1:n, j=c(2:n,1), x=-2, dims=c(n,n))
      G <- Matrix::sparseMatrix(i = rep(1:n,2), j=c(1:n, 3:n, 1, 2), x=1, dims=c(n,n))
    }
  }
  noise$h <- diff(mesh$loc)
  noise$h <- c(noise$h, mean(noise$h))

  if (is.null(index_NA)) index_NA <- rep(FALSE, length(x))
  if (is.null(replicate)) replicate <- rep(1, length(x))

  mesh_W <- INLA::inla.mesh.1d(loc=x)
  tmp <- ngme_make_A(
    mesh = mesh_W,
    map = x,
    n_map = length(x),
    idx_NA = index_NA,
    replicate = replicate
  )
  A <- tmp$A; A_pred <- tmp$A_pred

  nrep <- ncol(A) / ncol(C)
  stopifnot(nrep == as.integer(nrep))

  # update noise with length n
  if (noise$n_noise == 1) noise <- update_noise(noise, n = mesh$n)

  args <- within(list(...), {
    model       = "rw" # if (order == 1) "rw1" else "rw2"
    theta_K     = 1
    fix_theta_K = TRUE
    W_size      = ncol(C) # mesh$n
    V_size      = mesh$n
    A           = ngme_as_sparse(A)
    A_pred      = A_pred
    C           = ngme_as_sparse(C)
    G           = ngme_as_sparse(G)
    K           = C + G
    noise       = noise
    h           = noise$h
    mesh        = mesh_W
    map         = x
    n_map       = length(x)
    replicate   = replicate
    n_rep       = nrep
    rw_order    = order
  })

  do.call(ngme_model, args)
}


#' ngme model - Ornstein–Uhlenbeck process
#'
#' @param map     numerical vector, covariates to build index for the process
#' @param replicate replicate for the process
#' @param data      specifed or inherit from ngme formula
#'   if it is circular, we will treat the 1st location and the last location same.
#' @param index_NA  Logical vector, same as is.na(response variable)
#' @param noise     1. string: type of model, 2. ngme.noise object
#'  (can also be specified in each ngme model)
#' @param B_theta_K   Basis matrix for theta, by default it is a matrix of 1
#' @param theta_K     theta parameter, will be exp(B_theta * theta)
#' @param ...       additional arguments
#'
#' @return a list
#' @export
#'
#' @examples
#' r2 <- model_ou(1:7, theta_K = 2); r2$K
#'
model_ou <- ou <- function(
  map,
  B_theta_K   = NULL,
  theta_K     = 0,
  replicate   = NULL,
  data        = NULL,
  index_NA    = NULL,
  noise       = noise_normal(),
  # extra A matrix
  ...
) {
# capture symbol in index
  x <- eval(substitute(map), envir = data, enclos = parent.frame())
  stopifnot("length of index at least >= 3" = length(x) >= 3)
  if (is.null(replicate)) replicate <- rep(1, length(x))
  if (is.null(index_NA)) index_NA <- rep(FALSE, length(x))
  stopifnot("Make sure length(x)==length(replicate)" = length(x) == length(replicate))

  if (is.null(B_theta_K)) B_theta_K <- matrix(1, nrow = length(x), ncol = 1)
  stopifnot("B_theta is a matrix" = is.matrix(B_theta_K))
  stopifnot("ncol(B_theta_K) == length(theta_K)"
    = ncol(B_theta_K) == length(theta_K))

  mesh <- INLA::inla.mesh.1d(loc=x)

  # h <- diff(x) ?
  h <- diff(mesh$loc)
  noise$h <- h <- c(h, mean(h))
  n <- length(h)
  G <- Matrix::bandSparse(n=n,m=n,k=c(-1,0),diagonals=cbind(-rep(1,n), rep(1,n)))
  C <- Ce <- Matrix::bandSparse(n=n,m=n,k=c(-1,0),diagonals=cbind(0.5*c(h[-1],0), 0.5*h))
  Ci = Matrix::sparseMatrix(i=1:n,j=1:n,x=1/h,dims = c(n,n))

  if (is.null(index_NA)) index_NA <- rep(FALSE, n)
  if (is.null(replicate)) replicate <- rep(1, n)

  tmp <- ngme_make_A(
    mesh = mesh,
    map = x,
    n_map = length(x),
    idx_NA = index_NA,
    replicate = replicate
  )
  A <- tmp$A; A_pred <- tmp$A_pred

  nrep <- ncol(A) / ncol(C)
  stopifnot(nrep == as.integer(nrep))

  # update noise with length n
  if (noise$n_noise == 1) noise <- update_noise(noise, n = n)

kappas <- exp(as.numeric(B_theta_K %*% theta_K))
  args <- within(list(...), {
    model       = "ou"
    B_K         = B_theta_K
    theta_K     = theta_K
    W_size      = ncol(C) # mesh$n
    V_size      = mesh$n
    A           = ngme_as_sparse(A)
    A_pred      = A_pred
    C           = ngme_as_sparse(C)
    G           = ngme_as_sparse(G)
    K           = Matrix::Diagonal(x=kappas) %*% C + G
    noise       = noise
    h           = noise$h
    mesh        = mesh
    map         = x
    n_map       = length(x)
    replicate   = replicate
    n_rep       = nrep
  })

  do.call(ngme_model, args)
}
