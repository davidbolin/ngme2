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

# example
#   List of 3
#  $ mesurement.noise: num 0.495
#  $ fixed.effects   : num [1:3] -3.01 -1.03 2.02
#  $ block.W   : num [1:3] -3.01 -1.03 2.02
#  $ latent.model    :List of 1
#   ..$ :List of 3
#   .. ..$ W        : num [1:2000] 7.136 3.349 0.836 -1.486 20.241 ...
#   .. ..$ V        : num [1:2000] 1.948 0.705 1.318 0.382 4.53 ...
#   .. ..$ estimates:List of 4
#   .. .. ..$ alpha      : num 0.491
#   .. .. ..$ theta.mu   : num 1.83
#   .. .. ..$ theta.sigma: num 1.13
#   .. .. ..$ theta.noise: num 1


# ngme.as.sparse() {
#
# }


# ngme.simulate() {
# }

