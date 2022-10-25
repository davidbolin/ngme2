#' show ngme noise types
#'
#' @return
#' @export
#'
ngme_noise_types <- function() {
    c("normal", "nig")
}

#' Show ngme mdoel types
#'
#' @return a list of available models
#' @export
ngme_model_types <- function() {
    c("ar1", "matern", "matern_ns", "rw1", "rw2")
}
