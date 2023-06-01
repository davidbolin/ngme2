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
    noise <- model$noise

    # simulate noise
    h <- model$operator$h
    mu    <- as.numeric(noise$B_mu %*% noise$theta_mu)
    sigma <- as.numeric(exp(noise$B_sigma %*% noise$theta_sigma))
    nu <- noise$nu
    n <- length(mu)

    if (length(noise$noise_type) == 2) {
        # bivariate noise
        e1 <- simulate_noise(noise$noise_type[[1]],
          head(h, n/2), head(mu, n/2), head(sigma, n/2), nu[[1]], rnorm(1), noise$single_V)
        if (noise$share_V) e2 <- e1
        else e2 <- simulate_noise(noise$noise_type[[2]],
            tail(h, n/2), tail(mu, n/2), tail(sigma, n/2), nu[[2]], rnorm(1), noise$single_V)
        e <- c(e1, e2);
        attr(e, "V") <- c(attr(e1, "V"), attr(e2, "V"))
    } else {
        e <- simulate_noise(noise$noise_type, h, mu, sigma, nu, rnorm(1), noise$single_V)
    }
    W <- as.numeric(solve(model$operator$K, e))

    # attach noise attributes
    attr(W, "noise") <- noise
    attr(W, "V") <- attr(e, "V")
    W
}

simulate_noise <- function(
    noise_type, h_vec, mu_vec, sigma_vec, nu, seed, single_V
) {
    stopifnot(
        length(mu_vec) == length(sigma_vec),
        length(mu_vec) == length(h_vec)
    )
    n <- length(mu_vec)

    if (noise_type == "normal") {
        V <- h_vec
    } else if (noise_type == "nig" || noise_type == "normal_nig") {
        V <- if (single_V) h_vec * ngme2::rig(1, a=nu, b=nu, seed=seed)
            else ngme2::rig(n, a=nu, b=nu * (h_vec)^2, seed = seed)
    } else if (noise_type == "gal") {
        V <- if (single_V) h_vec * rgamma(1, nu, nu)
            else rgamma(n, shape = h_vec * nu, rate = nu)
    } else {
        stop("unknown noise to simulate")
    }

    e <- mu_vec * (V - h_vec) + sigma_vec * sqrt(V) * rnorm(n)
    attr(e, "V") <- V
    e
}