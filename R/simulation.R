#' Generic function for generating
#' the realization of latent or noise
#'
#' @param x ngme latent or noise object
#' @param ...  extra argument
#'
#' @return a realization and a noise object
#' @export
#'
#' @examples
#' n_obs <- 10
#' ar1 <- ngme.simulate.process(
#'    f(1:10, model = "ar1", theta_K = 0.5, noise = ngme.noise.nig()),
#'    seed = 10
#' )
#'
ngme.simulate <- function(x, ...) {
    UseMethod("ngme.simulate")
}

#' Simulate latent process with noise
#'
#' @param model latent process specified by f() function
#' @param noise noise object
#' @param seed seed
#'
#' @return a realization of latent model
#' @export
#'
#' @examples
#' n_obs <- 10
#' ngme.simulate(
#'   f(1:n_obs, model = "ar1", theta_K = 0.4),
#'   noise = ngme.noise(
#'     theta_mu = 2,
#'     theta_sigma = 0,
#'     theta_V = 2
#'   ),
#'   seed=NULL
#' )$realization
ngme.simulate.process <- function(
    model,
    seed   = NULL
) {
    n <- model$n_mesh
    noise_obj <- ngme.simulate.noise(model$noise, n = n, seed = seed)
    noise_real <- noise_obj$realization

    # create operator structure
    if (model$model_type == "ar1") {
        alpha <- model$operator$theta_K
        W <- Reduce(function(x, y) {y + alpha * x}, noise_real, accumulate = T)
    } else if (model$model_type == "matern") {
        # to-do
    } else {
        stop("not implement yet")
    }

    list(
        realization = W,
        noise = noise_obj$noise
    )
}

#' Simulate ngme noise
#'
#' @param noise ngme noise object
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
#' )$realization
ngme.simulate.noise <- function(noise, n, seed = NULL) {
    if (is.null(seed)) seed <- as.numeric(Sys.time())

    noise <- update.ngme.noise(noise, n = n)
    # create nig noise
    if (noise$type == "nig") {
        mu <- drop(noise$B_mu %*% noise$theta_mu)
        sigma <- drop(exp(noise$B_sigma %*% noise$theta_sigma))
        eta <- noise$theta_V
        V <- ngme2::rig(n, eta, eta, seed = seed);
        noise$V <- V
        e <- mu * (-1 + V) + sigma * sqrt(V) * rnorm(n)
    } else if (noise$type == "normal") {
        sigma <- drop(exp(noise$B_sigma %*% noise$theta_sigma))
        e <- rnorm(n, sd = sigma)
    }

    list(
        realization = e,
        noise = noise
    )
}
