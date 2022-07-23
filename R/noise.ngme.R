#' ngme noise specification
#'
#' @param type        type of noise
#' @param theta.noise parameter of noise
#'
#' @return
#' @export
#'
#' @examples
#' 

ngme.noise <- function(
  type="nig",
  theta.noise=1
  ) {

  list(
    type=type,
    theta.noise=theta.noise,
    n_params = length(theta.noise)
  )
}
