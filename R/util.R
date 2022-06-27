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

ngme.start <- function(
  last.fit      = NULL,
  sigma.eps     = 1,
  fixed.effects = NULL,
  latent.model  = NULL
) {
  if (inherits(last.fit, "ngme")) {
    start=NULL
  } else if (!is.null(last.fit)) {
    stop("last.fit should be an output of ngme object")
  } else {
    start=NULL
  }

  start
}
# example
#   List of 3
#  $ mesurement.noise: num 0.495
#  $ fixed.effects   : num [1:3] -3.01 -1.03 2.02
#  $ latent.model    :List of 1
#   ..$ :List of 3
#   .. ..$ W        : num [1:2000] 7.136 3.349 0.836 -1.486 20.241 ...
#   .. ..$ V        : num [1:2000] 1.948 0.705 1.318 0.382 4.53 ...
#   .. ..$ estimates:List of 4
#   .. .. ..$ alpha      : num 0.491
#   .. .. ..$ theta.mu   : num 1.83
#   .. .. ..$ theta.sigma: num 1.13
#   .. .. ..$ theta.noise: num 1

