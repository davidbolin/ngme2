# ar1 and rw1

#' ngme ar1 model specification
#'
#' Generating C, G and A given index and replicates
#'
#' @param index index for the process
#' @param replicates replicates for the process
#' @param alpha initial value for alpha
#' @param range range for the mesh
#'
#' @param noise noise, can be specified in f()
#' @param data data, can be specified in f(), ngme()
#'
#' @return a list of specification of model
#' @export
#'
#' @examples
#' ar1(index = c(1:3, 1:3), replicates = c(1,1,1,2,2,2))
#' f(index = xx, model = "ar1", data=list(xx = c(2,4,5)), noise=noise_nig())
model_ar1 <- function(
  index,
  replicates  = NULL,
  alpha       = 0.5,
  range       = c(1, max(index)),

  index_pred  = NULL,
  use_num_dK  = FALSE,
  data        = NULL,
  noise       = noise_normal(),
  ...
) {
  # capture symbol in index
  index <- eval(substitute(index), envir = data)

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
    A           = ngme_ts_make_A(loc = index, replicates = replicates, range = range)
    A_pred      = ngme_ts_make_A(index_pred, replicates = replicates, range = range)
    h           = rep(1.0, n)
    C           = ngme_as_sparse(C)
    G           = ngme_as_sparse(G)
    noise       = noise
  })

  do.call(ngme_model, args)
}

#' ngme model - random walk of order 1
#'
#' Generating C, G and A given index and replicates
#' size of C and G is (n-1) * n, size of V is n-1
#'
#' @param index index for the process
#' @param replicates replicates for the process
#' @param mesh inla.1d.mesh
#' @param n_points or num of points, evenly spaced mesh
#' @return a list
#' @export
#'
#' @examples
#' r1 <- model_rw1(1:7, circular = T); r1$C + r1$G
#' r2 <- model_rw1(1:7); r2$C + r2$G
model_rw1 <- function(
  index,
  replicates = NULL,
  circular = FALSE,
  n_points = NULL,
  noise = noise_normal(),
  # extra A matrix
  ...
) {
  stopifnot(length(index) > 2)
  # create mesh using index
  sorted_index <- sort(index, index.return = TRUE)
  h <- diff(sorted_index$x)
  # permutation matrix, same as A <- diag(length(index))[sorted_index$ix, ]
  A <- Matrix::sparseMatrix(seq_along(index), sorted_index$ix, x = 1)

  n <- length(index) - 1
  if (!circular) {
    # construct G
      G <- Matrix::Matrix(diag(n));
      G <- cbind(G, rep(0, n))

    # construct C
      C <- Matrix::Matrix(0, n, n)
      C[seq(n+1, n*n, by = n+1)] <- -1
      C <- cbind(C, c(rep(0, n-1), -1))
  } else {
    G <- Matrix::Matrix(diag(n))
    C <- Matrix::Matrix(0, n, n)
      C[seq(n+1, n*n, by = n+1)] <- -1
      C[n, 1] <- -1
  }

  # update noise with length n
  if (noise$n_noise == 1) noise <- update_noise(noise, n = n)

  args <- within(list(...), {
    model       = "rw1"
    theta_K     = 1
    fix_theta_K = TRUE
    W_size      = ncol(C) # n + 1
    V_size      = n
    A           = A
    # A_pred      = ngme.ts.make.A(index_pred, replicates = replicates, range = range),
    h           = h
    C           = ngme_as_sparse(C)
    G           = ngme_as_sparse(G)
    noise       = noise
  })

  do.call(ngme_model, args)
}

##### duplicate of rw1 except C and G

#' ngme model - random walk of order 2
#'
#' Generating C, G and A given index and replicates
#' size of C and G is (n-1) * n, size of V is n-1
#'
#' @param index index for the process
#' @param replicates replicates for the process
#' @param mesh inla.1d.mesh
#' @param n_points or num of points, evenly spaced mesh
#' @return a list
#' @export
#'
#' @examples
#' r2 <- model_rw2(1:7, circular = T); r2$C + r2$G
#' r3 <- model_rw2(1:7); r3$C + r3$G
model_rw2 <- function(
  index,
  replicates = NULL,
  circular = FALSE,
  n_points = NULL,
  noise = noise_normal(),
  # extra A matrix
  ...
) {
  stopifnot(length(index) > 3)
  # create mesh using index
  sorted_index <- sort(index, index.return = TRUE)
  h <- diff(sorted_index$x)
  # permutation matrix, same as A <- diag(length(index))[sorted_index$ix, ]
  A <- Matrix::sparseMatrix(seq_along(index), sorted_index$ix, x = 1)

  n <- length(index) - 2
  if (!circular) { # n * n+2
    C <- Matrix::Matrix(0, n, n+2)
    G <- Matrix::Matrix(diag(n));
    G <- cbind(G, rep(0, n), rep(0, n))
    G[seq(n+1, n*(n+2), by = n+1)] <- -2
    G[seq(2*n+1, n*(n+2), by = n+1)] <- 1
  } else {
    C <- Matrix::Matrix(0, n, n)
    G <- Matrix::Matrix(diag(n));
    G[seq(n+1, n*n, by = n+1)] <- -2
    G[seq(2*n+1, n*n, by = n+1)] <- 1
    G[n] <- -2
    G[c(n-1, 2*n)] <- 1
  }

  # update noise with length n
  if (noise$n_noise == 1) noise <- update_noise(noise, n = n)

  args <- within(list(...), {
    model       = "rw1"
    theta_K     = 1
    fix_theta_K = TRUE
    W_size      = ncol(C) # n + 1
    V_size      = n
    A           = A
    # A_pred      = ngme.ts.make.A(index_pred, replicates = replicates, range = range),
    h           = h
    C           = ngme_as_sparse(C)
    G           = ngme_as_sparse(G)
    noise       = noise
  })

  do.call(ngme_model, args)
}