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
  noise       = NULL,
  ...
) {
# print(as.list(match.call())[-1])
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
    theta_K     = if (is.null(theta_K)) ar1_a2th(alpha) else ar1_a2th(theta_K)
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
model_rw1 <- function(
  index,
  replicates = NULL,
  circular = FALSE,
  n_points = NULL,
  noise = NULL,
  # extra A matrix
  ...
) {
  # create mesh using index
  sorted_index <- sort(index, index.return = TRUE)
  h <- diff(sorted_index$x)
  # permutation matrix, same as A <- diag(length(index))[sorted_index$ix, ]
  A <- Matrix::sparseMatrix(seq_along(index), sorted_index$ix, x = 1)

  n <- length(index) - 1
  # construct G
    G <- Matrix::Matrix(diag(n));
    G <- cbind(G, rep(0, n))

  # construct C
    C <- Matrix::Matrix(0, n, n)
    C[seq(n+1, n*n, by = n+1)] <- -1
    C <- cbind(C, c(rep(0, n-1), -1))

  # update noise with length n
  if (noise$n_noise == 1) noise <- update_noise(noise, n = n)

  args <- within(list(...), {
    model       = "rw1"
    theta_K     = 1
    fix_theta_K = TRUE
    W_size      = n + 1
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

#' ngme model - random walk of order 2
#'
#' Generating C, G and A given index and replicates
#' size of C and G is (n-1) * n, size of V is n-1
#'
#' @param index index for the process
#' @param replicates replicates for the process
#' @param mesh inla.1d.mesh
#' @return a list
#' @export
#'
#' @examples
model_rw2 <- function(
  index,
  replicates = NULL,
  n_points = NULL,
  # extra A matrix
  ...
) {
  # create mesh using index
  sorted_index <- sort(index, index.return = TRUE)
  h <- diff(sorted_index$x)
  # permutation matrix, same as A <- diag(length(index))[sorted_index$ix, ]
  A <- Matrix::sparseMatrix(seq_along(index), sorted_index$ix, x = 1)

  n <- length(index) - 1
  # construct G
    G <- Matrix::Matrix(diag(n));
    G <- cbind(G, rep(0, n))

  # construct C
    C <- Matrix::Matrix(0, n, n)
    C[seq(n+1, n*n, by = n+1)] <- -1
    C <- cbind(C, c(rep(0, n-1), -1))

  # update noise with length n
  if (noise$n_noise == 1) noise <- update_noise(noise, n = n)

  args <- within(list(...), {
    model       = "rw2"
    theta_K     = 1
    fix_theta_K = TRUE
    W_size      = n + 1
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
