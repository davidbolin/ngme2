# new models should have
# model_name <- function(
# index, replicates, index_NA, noise, data
# )

# ar1 and rw

#' ngme AR(1) model specification
#'
#' Generating C, G and A given index and replicates
#'
#' @param x integer vector, time index for the AR(1) process
#' @param replicates replicates for the process
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
#' model_ar1(c(1:3, 1:3), replicates = c(1,1,1,2,2,2))
#' f(xx, model = "ar1", data=list(xx = c(2,4,5)), noise=noise_nig())
model_ar1 <- function(
  x,   # time index
  replicates  = NULL,
  index_NA    = NULL,
  data        = NULL,
  noise       = noise_normal(),
  alpha       = 0.5,
  ...
) {
  # capture symbol in index
  index <- eval(substitute(x), envir = data, enclos = parent.frame())
  stopifnot("The index should be integers." = all(index == round(index)))
  range <- c(min(index), max(index))

  if (is.null(index_NA)) index_NA <- rep(FALSE, length(index))

  if (is.null(replicates)) replicates <- rep(1, length(index))
  unique_rep <- unique(replicates)
  nrep <- length(unique_rep)

  # e.g. index = c(1:200, 1:100)
  #      replicates = c(rep(1, 200), rep(2, 100))
  #      n =200 in this case
  n <- range[2] - range[1] + 1

  # construct G
    G <- Matrix::Matrix(diag(n));

  # construct C
    C <- Matrix::Matrix(0, n, n)
    C[seq(2, n*n, by = n+1)] <- -1

    if(!is.null(nrep)) {
      C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), C)
      G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), G)
    }

  # update noise with length n
  if (noise$n_noise == 1) noise <- update_noise(noise, n = n)

  # remove duplicate symbol in ... (e.g. theta_K)
  args <- within(list(...), {
    model       = "ar1"
    theta_K     = if (exists("theta_K")) ar1_a2th(theta_K) else ar1_a2th(alpha)
    W_size      = n
    V_size      = n
    A           = ngme_ts_make_A(loc = index[!index_NA], replicates = replicates, range = range)
    A_pred      = ngme_ts_make_A(loc = index[index_NA], replicates = replicates, range = range)
    h           = noise$h
    C           = ngme_as_sparse(C)
    G           = ngme_as_sparse(G)
    K           = alpha * C + G
    noise       = noise
  })

  do.call(ngme_model, args)
}

#' ngme model - random walk of order 1
#'
#' Generating C, G and A given index and replicates
#' size of C and G is (n-1) * n, size of V is n-1
#'
#' @param x         numerical vector, covariates to build index for the process
#' @param order     1 or 2, order of random walk model
#' @param replicates replicates for the process
#' @param data      specifed or inherit from ngme formula
#' @param circular  whether the mesh is circular, i.e. the first one is connected to the last
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
  replicates  = NULL,
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
  if (is.null(index_NA)) index_NA <- rep(FALSE, length(x))

  if (is.null(replicates)) replicates <- rep(1, length(x))
  unique_rep <- unique(replicates)
  nrep <- length(unique_rep)

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

  # create A in INLA way
  mesh <- INLA::inla.mesh.1d(loc = x)
  # p_matrix <- as(mesh$idx$loc, "pMatrix") # bug! p_matrix

  A <- INLA::inla.spde.make.A(mesh = mesh, loc = x[!index_NA])
  A_pred <- if (!any(index_NA)) NULL else
   INLA::inla.spde.make.A(mesh = mesh, loc = x[index_NA])

  noise$h <- diff(mesh$loc)
# browser()
  n <- mesh$n
  if (order == 1) {
    # k is of size n-1 * n
    n <- n - 1
    if (!circular) {
      G <- -Matrix::Matrix(diag(n));
      G <- cbind(G, rep(0, n))
      C <- Matrix::Matrix(0, n, n)
      C[seq(n+1, n*n, by = n+1)] <- 1
      C <- cbind(C, c(rep(0, n-1), 1))
    } else {
      G <- Matrix::Matrix(diag(n))
      C <- Matrix::Matrix(0, n, n)
      C[seq(n+1, n*n, by = n+1)] <- -1
      C[n, 1] <- -1
    }
  } else if (order == 2) {
    n <- n - 2
    if (!circular) { # n * n+2
      stopifnot(n >= 2)
      C <- Matrix::Matrix(0, n, n+2)
      G <- Matrix::Matrix(diag(n));
      G <- cbind(G, rep(0, n), rep(0, n))
      G[seq(n+1, n*(n+2), by = n+1)] <- -2
      G[seq(2*n+1, n*(n+2), by = n+1)] <- 1
    } else {
      stopifnot(n >= 3)
      C <- Matrix::Matrix(diag(n))
      G <- Matrix::Matrix(0, n, n)
      G[seq(n+1, n*n, by = n+1)] <- -2
      G[seq(2*n+1, n*n, by = n+1)] <- 1
      G[n] <- -2
      G[c(n-1, 2*n)] <- 1
    }
    noise$h <- noise$h[-1]
  }

  if(!is.null(nrep)) {
    C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), C)
    G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), G)
    A <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), A)
  }

  # update noise with length n
  if (noise$n_noise == 1) noise <- update_noise(noise, n = n)

  args <- within(list(...), {
    model       = if (order == 1) "rw1" else "rw2"
    theta_K     = 1
    fix_theta_K = TRUE
    W_size      = ncol(C) # mesh$n
    V_size      = n
    A           = ngme_as_sparse(A)
    A_pred      = A_pred
    C           = ngme_as_sparse(C)
    G           = ngme_as_sparse(G)
    K           = C + G
    noise       = noise
    h           = noise$h
  })

  do.call(ngme_model, args)
}
