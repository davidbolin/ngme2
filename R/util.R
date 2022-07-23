#' Make observation matrix for time series
#'
#' @param loc integers (after sorting, no gaps > 1)
#'
#' @return A matrix (length(loc) * length(unique(loc)))
#' @export
#'
#' @examples
#' ngme.ts.make.A(c(11, 13, 12, 12))
#'

# ngme.ts.make.A <- function(loc) {
#   n_loc = length(loc)
#   n_range = max(loc)-min(loc)+1
#   if (any((diff(sort(loc))) > 1)) stop("no gaps allowed")
#
#   A = matrix(0, nrow=n_loc, ncol=n_range)
#   for (i in 1:n_loc) {
#     A[i, loc[i]-min(loc)+1] = 1
#   }
# #Do an na.rm
#   as(A, "dgCMatrix")
# }

ngme.ts.make.A <- function(loc, replicates) {
  n_loc = length(loc)
  n_range = max(loc)-min(loc)+1
  nrep = 1
  if(!is.null(replicates)){
    unique_rep = unique(replicates)
    nrep = length(unique_rep)
  }
  A = matrix(0, nrow=n_loc, ncol=n_range * nrep)
  for (i in 1:n_loc) {
    ncol_rep <- which(unique_rep == replicates[i])
    A[i, (ncol_rep-1)*n_range + loc[i]-min(loc)+1] = 1
  }
  as(A, "dgCMatrix")
}

#' NGME starting point for block model
#'
#' @param sigma.eps mesurement noise
#' @param fixed.effects fixed effects
#' @param W initial W for the block model
#'
#' @return
#' @export
#'
#' @examples
ngme.start <- function(
  fixed.effects = NULL,
  sigma.eps     = NULL,
  W = NULL
) {

  start = list(
    fixed.effects     = fixed.effects,
    mesurement.noise  = sigma.eps,
    block.W           = W
  )

  return (start)
}


#' Convert sparse matrix into dgCMatrix
#'
#' @param G sparse matrix
#'
#' @return
#' @export
#'
#' @examples
ngme.as.sparse <- function(G) {
  
  if (!inherits(G, "dgCMatrix")) {
    tryCatch(
      expr = {
        G <- as(G, "dgCMatrix")
      },
      error = {
        # first dgT, then convert to dgC
        G = as(G, "dgTMatrix")
        # idx <- which(G@i <= G@j)
        # G = Matrix::sparseMatrix(
        #   i = G@i[idx], 
        #   j = G@j[idx], 
        #   x= G@x[idx],
        #   index1 = FALSE
        # )
        G = as(G, "dgCMatrix")
      }
    )
  }
  G
}