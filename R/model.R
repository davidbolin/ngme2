#' ngme model specification
#' 
#' Generating C, G and A given index and replicates
#'
#' @param index index for the process
#' @param replicates replicates for the process
#' @param alpha initial value for alpha
#' @param range range for the mesh
#'
#' @return a list of specification of model
#' @export
#'
#' @examples
ngme.ar1 <- function(
  index,
  replicates = NULL,
  alpha = 0.5,
  range = c(1, max(index)),

  index_pred = NULL,
  use_num_dK = FALSE,
  ...
) {
  # overwirte the default
  if (!is.null(list(...)$theta_K)) alpha <- list(...)$theta_K

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

  ar1_in <- list(
    A       = ngme.ts.make.A(loc = index, replicates = replicates, range = range),
    A_pred  = ngme.ts.make.A(index_pred, replicates = replicates, range = range),
    operator_in = list(
      n_params    = 1,
      theta_K     = alpha,
      C           = ngme.as.sparse(C),
      G           = ngme.as.sparse(G),
      use_num_dK  = use_num_dK
    )
  )

  class(ar1_in) <- "ngme.ar1"
  ar1_in
}


