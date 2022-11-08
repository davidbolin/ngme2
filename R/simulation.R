#' Simulate latent process with noise
#'
#' @param model latent process specified by f() function
#' @param seed seed
#'
#' @return a realization of latent model
#' @export
#'
#' @examples
#' B_sigma <- matrix(1:10 / 10, nrow=10)
#' simulate(noise_nig(n = 10, B_sigma =  B_sigma))
#' simulate(noise_normal(theta_sigma = 1.5, B_sigma = B_sigma))
#' simulate(f(1:10, model = "ar1", theta_K = 0.4, noise = noise_nig()))
simulate.ngme_model <- function(
    object,
    nsim   = 1,
    seed   = NULL,
    ...
) {
    model <- object
    n <- model$W_size
    sim_noise <- simulate.ngme_noise(model$noise, nsim = n, seed = seed)

    # create operator structure
    if (model$model == "ar1") {
        alpha <- ar1_th2a(model$theta_K)
        # for loop
        W <- Reduce(function(x, y) {y + alpha * x}, sim_noise, accumulate = T)
    } else if (model$model == "matern") {
        # K_a %*% W = noise
        W <- with(model, {
            C.inv <- as(Matrix::diag(1 / Matrix::diag(C)), "sparseMatrix")
            kappas <- drop(exp(B_kappa %*% theta_K))
            Kappa <- diag(kappas)
            # build K_a
            if (alpha == 2) {
                K_a <- (Kappa %*% C %*% Kappa + G)
            } else if (alpha == 4) {
                K_a <- (Kappa %*% C %*% Kappa + G) %*% C.inv %*% (Kappa %*% C %*% Kappa + G)
            }
            drop(solve(K_a, sim_noise))
        })
    } else {
        stop("not implement yet")
    }

    # attach noise attributes
    attr(W, "noise") <- attr(sim_noise, "noise")
    W
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
#' simulate(noise_normal(sd = 5, n = 10))
simulate.ngme_noise <- function(
    object,
    seed   = NULL,
    ...
) {
    noise <- object

    if (is.null(seed)) seed <- as.numeric(Sys.time())

    n <- noise$n_noise
    # create nig noise
    if (noise$noise_type == "nig") {
        mu <- drop(noise$B_mu %*% noise$theta_mu)
        sigma <- drop(exp(noise$B_sigma %*% noise$theta_sigma))
        eta <- noise$theta_V
        V <- ngme2::rig(n, eta, eta, seed = seed);
        noise$V <- V
        e <- mu * (-1 + V) + sigma * sqrt(V) * rnorm(n)
    } else if (noise$noise_type == "normal") {
        sigma <- drop(exp(noise$B_sigma %*% noise$theta_sigma))
        e <- rnorm(n, sd = sigma)
    }

    attr(e, "noise") <- noise
    e
}
