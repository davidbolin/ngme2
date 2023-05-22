#' Simulate latent process with noise
#'
#' @param object  ngme model specified by f() function
#' @param nsim ignored
#' @param seed seed
#' @param ... ignored
#'
#' @return a realization of latent model
#' @export
#'
#' @examples
#' B_sigma <- matrix(1:10 / 10, nrow=10)
#' simulate(noise_nig(n = 10, B_sigma = B_sigma))
#' simulate(noise_normal(theta_sigma = 1.5, B_sigma = B_sigma))
#' simulate(model_ar1(1:10, theta_K = 0.4, noise = noise_nig()))
simulate.ngme_model <- function(
    object,
    nsim   = 1,
    seed   = NULL,
    ...
) {
    model <- object
    sim_noise <- simulate.ngme_noise(model$noise, seed = seed)
    W <- as.numeric(solve(model$operator$K, sim_noise))

    # attach noise attributes
    attr(W, "noise") <- attr(sim_noise, "noise")
    W
}

#' Simulate ngme noise
#'
#' @param object ngme noise object
#' @param nsim ignored
#' @param seed seed
#' @param ... can take n_noise = 3 (for measurement noise)
#'
#' @return a realization of ngme noise
#' @export
#'
#' @examples
#' simulate(noise_normal(sd = 5, n = 10))
simulate.ngme_noise <- function(
    object,
    nsim   = 1,
    seed   = NULL,
    ...
) {
    noise <- object
    stopifnot(
        is.numeric(noise$h),
        noise$n_noise == length(noise$h)
    )

    if (is.null(seed)) seed <- as.numeric(Sys.time())

    if (length(noise$noise_type) == 2) {
        # bivariate noise
        e <- c(simulate.ngme_noise(noise$bv_noises[[1]], seed = seed),
            simulate.ngme_noise(noise$bv_noises[[2]], seed = seed))
    } else {
        n <- noise$n_noise
        if (noise$noise_type == "nig") {
            nu <- noise$nu
            V <- ngme2::rig(n, a=nu, b=nu*(noise$h)^2, seed = seed)
            mu <- drop(noise$B_mu %*% noise$theta_mu)
            sigma <- drop(exp(noise$B_sigma %*% noise$theta_sigma))
            noise$V <- V
            e <- mu * (V - noise$h) + sigma * sqrt(V) * rnorm(n)
        } else if (noise$noise_type == "normal") {
            sigma <- drop(exp(noise$B_sigma %*% noise$theta_sigma))
            e <- rnorm(n, sd = sigma * sqrt(noise$h))
        } else if (noise$noise_type == "gal") {
            nu <- noise$nu
            V <- rgamma(n, shape = noise$h * nu, rate = nu)
            mu <- drop(noise$B_mu %*% noise$theta_mu)
            sigma <- drop(exp(noise$B_sigma %*% noise$theta_sigma))
            noise$V <- V
            e <- mu * (V - noise$h) + sigma * sqrt(V) * rnorm(n)
        }
    }

    attr(e, "noise") <- noise
    e
}
