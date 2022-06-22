#' ngme noise specification
#'
#' @param noise noise type
#' @param init.param initial value
#'
#' @return
#' @export
#'
#' @examples
ngme.noise.spec <- function(
  noise="nig",
  init.param=1
  ) {

  list(
    noise=noise,
    init.param=init.param
  )
}
