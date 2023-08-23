#' Show ngme noise types
#'
#' @return available types for noise
#' @export
ngme_noise_types <- function() {
    c("normal", "nig", "gal", "normal_nig")
}

#' Show ngme model types
#'
#' @return available types for models
#' @export
ngme_model_types <- function() {
    c("ar1", "matern", "rw", "ou", "tp", "iid", "re", "bv")
}

#' Show ngme priors
#'
#' @return available types of priors
#' @export
ngme_prior_types <- function() {
    c("flat", "normal", "pc.sd")
}

