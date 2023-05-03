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
    n <- model$W_size
    sim_noise <- simulate.ngme_noise(model$noise, nsim = n, seed = seed)

    # create operator structure
    if (model$model == "ar1") {
        # assume W_0 = 0
        # W_1 = a*W_0 + noise
        # do the loop
        alpha <- ar1_th2a(model$theta_K)
        # for loop
        W <- Reduce(function(x, y) {y + alpha * x}, sim_noise, accumulate = T)
    } else if (model$model == "matern") {
        # K_a %*% W = noise
        W <- with(model$operator, {
            C.inv <- as(Matrix::diag(1 / Matrix::diag(C)), "sparseMatrix")
            if (ncol(B_K) > 1) {
                kappas <- as.numeric(exp(B_K %*% theta_K))
                Kappa <- diag(kappas)
                # build K_a
                if (alpha == 2) {
                    K_a <- (Kappa %*% C %*% Kappa + G)
                } else if (alpha == 4) {
                    K_a <- (Kappa %*% C %*% Kappa + G) %*% C.inv %*% (Kappa %*% C %*% Kappa + G)
                }
            } else {
                kappa <- exp(theta_K)
                K_a <- if (alpha == 2) kappa^2 * C + G else kappa^4 * C + G
            }
            as.numeric(solve(K_a, sim_noise))
        })
    } else if (model$model == "rw" && model$rw_order == 1) {
        # assume W_0 = 0
        # W_n = W_0 + cumsum(simnoise)
        W <- cumsum(sim_noise)

        W <- c(0, W) # watch out! compensate
        # h <- sim_noise
    } else if (model$model == "rw" && model$rw_order == 2) {
        # assume W_0 = 0, W_1 = 0
        W <- cumsum(cumsum(sim_noise))
        # starting value W_0, W_1
        # W_{n+1} <- W_0 + n * (W_1 - W_0) + cumsum(cumsum(sim_noise))
    } else {
        # general case
        W <- as.numeric(solve(model$operator$K, sim_noise))
    }

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

    if (is.null(seed)) seed <- as.numeric(Sys.time())

    if (!is.null(list(...)$n_noise)) {
        # only for simulate measurement noise
        # means noise h = 1
        n <- list(...)$n_noise
        new_mn <- ngme_noise(
            noise_type = noise$noise_type,
            theta_mu = noise$theta_mu,
            theta_sigma = noise$theta_sigma,
            nu = noise$nu,
            n = n
        )
        return(simulate.ngme_noise(new_mn))
    }

    # create nig noise
    n <- noise$n_noise
    if (noise$noise_type == "nig") {
        mu <- drop(noise$B_mu %*% noise$theta_mu)
        sigma <- drop(exp(noise$B_sigma %*% noise$theta_sigma))
        nu <- noise$nu
        V <- ngme2::rig(n, a=nu, b=nu*(noise$h)^2, seed = seed)
        noise$V <- V
        e <- mu * (V - noise$h) + sigma * sqrt(V) * rnorm(n)
    } else if (noise$noise_type == "normal") {
        sigma <- drop(exp(noise$B_sigma %*% noise$theta_sigma))
        e <- rnorm(n, sd = sigma * sqrt(noise$h))
    } else if (noise$noise_type == "gal") {
        mu <- drop(noise$B_mu %*% noise$theta_mu)
        sigma <- drop(exp(noise$B_sigma %*% noise$theta_sigma))
        nu <- noise$nu
        V <- rgamma(n, shape = noise$h * nu, rate = nu)
        noise$V <- V
        e <- mu * (V - noise$h) + sigma * sqrt(V) * rnorm(n)
    }
    attr(e, "noise") <- noise
    e
}
