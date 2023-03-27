#' Show ngme noise types
#'
#' @return available types for noise
#' @export
ngme_noise_types <- function() {
    c("normal", "nig", "gal", "normal_nig")
}

#' Show ngme mdoel types for K (operator)
#'
#' @return available types for models
#' @export
ngme_model_types <- function() {
    c("ar1", "matern", "rw", "ou", "tp", "iid")
}

#' Show ngme random effects types
#'
#' @return available types for models
#' @export
ngme_randeff_types <- function() {
    c("normal", "nig")
}

