# new models should have
# model_name <- function(
# index, replicate, index_NA, noise, data
# )

# ar1 and rw

#' ngme AR(1) model specification
#'
#' Generating C, G and A given index and replicate
#'
#' @param x integer vector, time index for the AR(1) process
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
model_ar1 <- function(
  x,   # time index
  replicate   = NULL,
  index_NA    = NULL,
  data        = NULL,
  noise       = noise_normal(),
  alpha       = 0.5,
  ...
) {
  # capture symbol in index
  index <- eval(substitute(x), envir = data, enclos = parent.frame())
  stopifnot("The index should be integers." = all(index == round(index)))

  if (is.null(index_NA)) index_NA <- rep(FALSE, length(index))
  if (is.null(replicate)) replicate <- rep(1, length(index))

  # e.g. index = c(1:200, 1:100)
  #      replicate = c(rep(1, 200), rep(2, 100))
  #      n =200 in this case

  mesh <- INLA::inla.mesh.1d(x)
  n <- mesh$n

  # construct G
  G <- Matrix::Matrix(diag(n));

  # construct C
  C <- Matrix::Matrix(0, n, n)
  C[seq(2, n*n, by = n+1)] <- -1

  # update noise with length n
  if (noise$n_noise == 1) noise <- update_noise(noise, n = n)

  # make A and A_pred
  tmp <- ngme_make_A(
    mesh = mesh,
    map = x,
    n_map = length(x),
    idx_NA = index_NA,
    replicate = replicate
  )
  A <- tmp$A; A_pred <- tmp$A_pred

  if (is.null(A)) stop("A is NULL")
  nrep <- ncol(A) / ncol(C)
  stopifnot(nrep == as.integer(nrep))
  C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), C)
  G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), G)
  noise$h <- rep(noise$h, times = nrep)

  # remove duplicate symbol in ... (e.g. theta_K)
  args <- within(list(...), {
    mesh        = mesh
    model       = "ar1"
    theta_K     = if (exists("theta_K")) ar1_a2th(theta_K) else ar1_a2th(alpha)
    W_size      = n
    V_size      = n
    A           = A
    A_pred      = A_pred
    # A           = ngme_ts_make_A(loc = index[!index_NA], replicate = replicate, range = range)
    # A_pred      = ngme_ts_make_A(loc = index[index_NA], replicate = replicate, range = range)
    h           = noise$h
    C           = ngme_as_sparse(C)
    G           = ngme_as_sparse(G)
    K           = alpha * C + G
    noise       = noise
    map         = x
    n_map       = length(x)
  })

  do.call(ngme_model, args)
}

#' ngme model - random walk of order 1
#'
#' Generating C, G and A given index and replicate
#' size of C and G is (n-1) * n, size of V is n-1
#'
#' @param x         numerical vector, covariates to build index for the process
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
model_rw <- function(
  x,
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
  x <- eval(substitute(x), envir = data, enclos = parent.frame())

  stopifnot("length of index at least > 3" = length(x) > 3)
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
      n_mesh <- mesh$n
      G <- -Matrix::Matrix(diag(n_mesh));
      G <- cbind(G, rep(0, n_mesh))
      C <- Matrix::Matrix(0, n_mesh, n_mesh)
      C[seq(n_mesh+1, n_mesh*n_mesh, by = n_mesh+1)] <- 1
      C <- cbind(C, c(rep(0, n_mesh-1), 1))
    } else if (order == 2) {
      mesh <- INLA::inla.mesh.1d(loc = unique(x)[-c(1,2)])
      n_mesh <- mesh$n
      stopifnot(n_mesh >= 2)
      C <- Matrix::Matrix(0, n_mesh, n_mesh+2)
      G <- Matrix::Matrix(diag(n_mesh));
      G <- cbind(G, rep(0, n_mesh), rep(0, n_mesh))
      G[seq(n_mesh+1, n_mesh*(n_mesh+2), by = n_mesh+1)] <- -2
      G[seq(2*n_mesh+1, n_mesh*(n_mesh+2), by = n_mesh+1)] <- 1
    }
  } else {
    # circular mesh (remove last element)
      stopifnot(length(x) >= 4)
    mesh <- INLA::inla.mesh.1d(loc = x)
    n_mesh <- mesh$n
    if (order == 1) {
      G <- Matrix::Matrix(diag(n_mesh))
      C <- Matrix::Matrix(0, n_mesh, n_mesh)
      C[seq(n_mesh+1, n_mesh*n_mesh, by = n_mesh+1)] <- -1
      C[n_mesh, 1] <- -1
    } else if (order == 2) {
      C <- Matrix::Matrix(diag(n_mesh))
      G <- Matrix::Matrix(0, n_mesh, n_mesh)
      G[seq(n_mesh+1, n_mesh*n_mesh, by = n_mesh+1)] <- -2
      G[seq(2*n_mesh+1, n_mesh*n_mesh, by = n_mesh+1)] <- 1
      G[n_mesh] <- -2
      G[c(n_mesh-1, 2*n_mesh)] <- 1
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
  if (is.null(A)) stop("A is NULL")
  nrep <- ncol(A) / ncol(C)
  stopifnot(nrep == as.integer(nrep))
  C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), C)
  G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), G)
  noise$h <- rep(noise$h, times = nrep)

  # update noise with length n
  if (noise$n_noise == 1) noise <- update_noise(noise, n = n_mesh)

  args <- within(list(...), {
    model       = if (order == 1) "rw1" else "rw2"
    theta_K     = 1
    fix_theta_K = TRUE
    W_size      = ncol(C) # mesh$n
    V_size      = n_mesh
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
  })

  do.call(ngme_model, args)
}
