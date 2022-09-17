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
ngme_simulate <- function(x, ...) {
    UseMethod("ngme_simulate")
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
ngme_simulate.ngme_model <- function(
    model,
    seed   = NULL
) {
    n <- model$n_mesh
    noise <- ngme_simulate.ngme_noise(model$noise, n = n, seed = seed)

    # create operator structure
    if (model$model == "ar1") {
        alpha <- model$theta_K
        W <- Reduce(function(x, y) {y + alpha * x}, noise, accumulate = T)
    } else if (model$model == "matern") {
        # K_a %*% W = noise
        W <- with(model, {
            C.inv <- as(Matrix::diag(1 / Matrix::diag(C)), "sparseMatrix")
            kappas <- drop(exp(B_kappa %*% theta_kappa))
            Kappa <- diag(kappas)
            # build K_a
            if (alpha == 2) {
                K_a <- (Kappa %*% C %*% Kappa + G)
            } else if (alpha == 4) {
                K_a <- (Kappa %*% C %*% Kappa + G) %*% C.inv %*% (Kappa %*% C %*% Kappa + G)
            }
            drop(solve(K_a, noise))
        })
    } else {
        stop("not implement yet")
    }

    attr(W, "noise") <- attr(noise, "noise")
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
#' ngme.simulate(
#'   noise = ngme.noise(
#'     theta_mu = 2,
#'     theta_sigma = 0,
#'     theta_V = 2
#'   ),
#'   n = 10,
#'   seed=NULL
#' )$realization
ngme_simulate.ngme_noise <- function(noise, n, seed = NULL) {
    if (is.null(seed)) seed <- as.numeric(Sys.time())

    noise <- update.ngme.noise(noise, n = n)
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
