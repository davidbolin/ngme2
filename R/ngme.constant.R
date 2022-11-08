#' show ngme noise types
#'
#' @return available types for noise
#' @export
ngme_noise_types <- function() {
    c("normal", "nig")
}

#' Show ngme mdoel types
#'
#' @return available types for models
#' @export
ngme_model_types <- function() {
    c("ar1", "matern", "rw1", "rw2")
}
