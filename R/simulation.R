#' Simulation of stochastic process
#'
#' @param model
#' @param B.mu
#' @param B.sigma
#' @param theta.mu
#' @param theta.sigma
#' @param noise
#'
#' @return
#' @export
#'
#' @examples
ngme.simulate <- function(
    model = model,
    noise = ngme.noise(),
    n = 100
) {
    # create nig noise
    if (noise$type == "nig") {
        stopifnot(noise$theta_sigma)
    } else if (noise$type == "normal") {
        stopifnot("non-stationary normal not yet" = noise$theta_sigma == 1)
        rnorm(n) * exp(theta_sigma)
    }

    # create operator structure
    
}
