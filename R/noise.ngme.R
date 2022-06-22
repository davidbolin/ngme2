#' ngme noise specification
#'
#' @param type
#' @param theta.noise
#'
#' @return
#' @export
#'
#' @examples
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
