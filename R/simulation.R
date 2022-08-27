#' simulate one latent model with noise
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
ngme.simulate.latent <- function(
    latent = NULL,
    noise  = ngme.noise(),
    seed   = 1
) {
    n <- latent$n_mesh
    e <- ngme.simulate.noise(noise, n = n)

    # create operator structure
    if (latent$model_type == "ar1") {
        alpha <- latent$operator$theta_K
        W <- Reduce(function(x, y) {y + alpha*x}, e, accumulate = T)
    } else {
        stop("not implement yet")
    }

    W
}

ngme.simulate <- function() {
    UseMethod("ngme.simulate")
}

#' simulate noise
#'
#' @param model
#'
#' @return
#' @export
#'
#' @examples
ngme.simulate.noise <- function(noise, n, seed = NULL) {
    if (is.null(seed)) seed <- as.numeric(Sys.time())

    noise <- update.ngme.noise(noise, n = n)
    # create nig noise
    if (noise$type == "nig") {
        mu <- drop(noise$B_mu %*% noise$theta_mu)
        sigma <- drop(exp(noise$B_sigma %*% noise$theta_sigma))
        eta <- noise$theta_V
        V <- ngme2::rig(n, eta, eta, seed = seed)
        e <- mu * (-1 + V) + sigma * sqrt(V) * rnorm(n)
    } else if (noise$type == "normal") {
        sigma <- drop(exp(noise$B_sigma %*% noise$theta_sigma))
        e <- rnorm(n, sd = sigma)
    }

    e
}