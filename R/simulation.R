#' Generic function for generating
#' the realization of latent or noise
#'
#' @param x ngme latent or noise object
#' @param ...  extra argument
#'
#' @return a realization of latent or noise
#' @export
#'
#' @examples
#' n_obs <- 10
#' ar1 <- f(1:n_obs, model = "ar1", theta_K = 0.4)
#' ngme.simulate(
#'   ar1,
#'   noise = ngme.noise(
#'     theta_mu = 2,
#'     theta_sigma = 0,
#'     theta_V = 2
#'   ),
#'   seed=NULL
#' )
#'
ngme.simulate <- function(x, ...) {
    UseMethod("ngme.simulate")
}

#' Simulate latent model with noise
#'
#' @param latent latent model
#' @param noise noise object
#' @param seed seed
#'
#' @return a realization of latent model
#' @export
#'
#' @examples
#' n_obs <- 10
#' ar1 <- f(1:n_obs, model = "ar1", theta_K = 0.4)
#' ngme.simulate(
#'   ar1,
#'   noise = ngme.noise(
#'     theta_mu = 2,
#'     theta_sigma = 0,
#'     theta_V = 2
#'   ),
#'   seed=NULL
#' )
ngme.simulate.latent <- function(
    latent,
    noise  = ngme.noise(),
    seed   = NULL
) {
    n <- latent$n_mesh
    e <- ngme.simulate.noise(noise, n = n, seed = seed)

    # create operator structure
    if (latent$model_type == "ar1") {
        alpha <- latent$operator$theta_K
        W <- Reduce(function(x, y) {y + alpha*x}, e, accumulate = T)
    } else {
        stop("not implement yet")
    }

    W
}

#' Simulate ngme noise
#'
#' @param noise noise object
#' @param n number of realization
#' @param seed seed
#'
#' @return a realization of ngme noise
#' @export
#'
#' @examples
#' ngme.simulate(
#'   noise = ngme.noise(
#'     theta_mu = 2,
#'     theta_sigma = 0,
#'     theta_V = 2
#'   ),
#'   n = 10,
#'   seed=NULL
#' )
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
