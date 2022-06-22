#' Make observation matrix for time series
#'
#' @param loc integers (after sorting, no gaps > 1)
#'
#' @return A matrix (length(loc) * length(unique(loc)))
#' @export
#'
#' @examples
#' ngme.ts.make.A(c(11, 13, 12, 12))
ngme.ts.make.A <- function(loc) {
  n_loc = length(loc)
  n_range = max(loc)-min(loc)+1
  if (any((diff(sort(loc))) > 1)) stop("no gaps allowed")

  A = matrix(0, nrow=n_loc, ncol=n_range)
  for (i in 1:n_loc) {
    A[i, loc[i]-min(loc)+1] = 1
  }

  as(A, "dgCMatrix")
}


# noise.spec <- function(noise="nig") {
#
# }
