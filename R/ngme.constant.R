#' show ngme noise types
#'
#' @return
#' @export
#'
ngme.noise.types <- function() {
    c("normal", "nig")
}

#' Show ngme model types
#'
#' @return a list of available models
#' @export
#'
ngme.model.types <- function() {
    c("ar1", "matern", "matern_ns", "rw1", "rw2", "matern2D")
}
