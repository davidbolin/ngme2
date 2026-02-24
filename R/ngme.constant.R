#' Show ngme noise types
#'
#' @return available types for noise
#' @export
ngme_noise_types <- function() {
    c("normal", "nig", "gal", "t", "skew_t", "normal_nig")
}

#' Show ngme model types
#'
#' @return available types for models
#' @export
ngme_model_types <- function() {
    c(
        "iid", "re",
        "ar1", "ar", "arma",
        "rw1", "rw2",
        "ou", "matern",
        "tp",
        "bv", "bv2", "bv_matern",
        "spacetime",
        "generic", "generic_ns"
    )
}

#' Show ngme priors
#'
#' @return available types of priors
#' @export
ngme_prior_types <- function() {
    c("none", "normal", "pc.sd", "half.cauchy", "inv.exponential")
}
