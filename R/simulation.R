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
#' simulate(f(1:10, model="ar1", rho = 0.4, noise = noise_nig()))
#' simulate(f(rnorm(10), model="rw1", noise = noise_normal()))
simulate.ngme_model <- function(
  object,
  nsim   = NULL,
  seed   = NULL,
  ...
) {
  if (is.null(seed)) seed <- Sys.time()
  model <- object
  noise <- model$noise

  # simulate noise
  h <- model$operator$h
  mu    <- as.numeric(noise$B_mu %*% noise$theta_mu)
  sigma <- as.numeric(exp(noise$B_sigma %*% noise$theta_sigma))
  nu    <- as.numeric(exp(noise$B_nu %*% noise$theta_nu))
  n <- length(mu)

  if (length(noise$noise_type) == 2) {
    # bivariate noise
    e1 <- simulate_noise(noise$noise_type[[1]],
      head(h, n/2), head(mu, n/2), head(sigma, n/2), head(nu, n/2),
      seed+1, noise$single_V)
    if (noise$share_V)
      e2 <- e1
    else e2 <- simulate_noise(noise$noise_type[[2]],
      tail(h, n/2), tail(mu, n/2), tail(sigma, n/2), tail(nu, n/2),
      seed+2, noise$single_V)
    e <- c(e1, e2);
    attr(e, "V") <- c(attr(e1, "V"), attr(e2, "V"))
  } else {
    e <- simulate_noise(noise$noise_type, h, mu, sigma, nu, seed+3, noise$single_V)
  }
  W <- as.numeric(solve(model$operator$K, e))
  W <- W - mean(W)

  # attach noise attributes
  attr(W, "noise") <- noise
  attr(W, "V") <- attr(e, "V")

  return (W)
}

simulate_noise <- function(
  noise_type, h_vec, mu_vec, sigma_vec, nu_vec, seed, single_V
) {
  set.seed(seed)
  stopifnot(
    length(mu_vec) == length(sigma_vec),
    length(mu_vec) == length(nu_vec),
    length(mu_vec) == length(h_vec)
  )
  n <- length(mu_vec)

  if (noise_type == "normal") {
    V <- h_vec
  } else if (noise_type == "nig" || noise_type == "normal_nig") {
    V <- if (single_V) h_vec * ngme2::rig(1, a=nu_vec[1], b=nu_vec[1], seed=seed)
      else ngme2::rig(n, a=nu_vec, b=nu_vec * (h_vec)^2, seed = seed)
  } else if (noise_type == "gal") {
    V <- if (single_V) h_vec * rgamma(1, nu_vec[1], nu_vec[1])
      else rgamma(n, shape = h_vec * nu_vec, rate = nu_vec)
  } else {
    stop("unknown noise to simulate")
  }

  e <- mu_vec * (V - h_vec) + sigma_vec * sqrt(V) * rnorm(n)
  attr(e, "V") <- V
  e
}

#' Simulate ngme noise object
#'
#' @param object  ngme noise object
#' @param h should be of same length as nsim
#' @param seed seed
#' @param nsim length of the noise
#' @param ... ignored
#'
#' @return a realization of noise
#' @export
simulate.ngme_noise <- function(
  object,
  nsim = NULL,
  seed = NULL,
  h = NULL,
  ...
) {
  n_noise <- max(nrow(object$B_mu), nrow(object$B_sigma), nrow(object$B_nu), nsim)
  if (is.null(seed)) seed <- Sys.time()
  if (is.null(h)) h <- rep(1, n_noise)
  if (length(h) > n_noise) n_noise <- length(h)
  stopifnot(length(h) == n_noise)

  # return this
  with(object, {
    res <- numeric(n_noise)
    mu_vec    <- as.numeric(B_mu %*% theta_mu)
    sigma_vec <- as.numeric(exp(B_sigma %*% theta_sigma))
    nu_vec <- as.numeric(exp(B_nu %*% theta_nu))
    if (length(mu_vec) == 1) mu_vec <- rep(mu_vec, n_noise)
    if (length(sigma_vec) == 1) sigma_vec <- rep(sigma_vec, n_noise)
    if (length(nu_vec) == 1) nu_vec <- rep(nu_vec, n_noise)

    if (!corr_measurement) {
      res <- simulate_noise(noise_type, h, mu_vec, sigma_vec, nu_vec, seed, single_V)
    } else {
      # simulate correlated noise
      i = 1;
      while (i <= n_noise) {
        if (i==n_noise || index_corr[[i]] != index_corr[[i+1]]) {
          res[i] <- simulate_noise(noise_type, h[i], mu_vec[i], sigma_vec[i], nu_vec[i], seed, single_V)
          i = i + 1
        } else {
          # simulate a pair correlated noise
          if (noise_type == "normal") {
            cov_mat <- matrix(c(sigma_vec[i]^2, sigma_vec[i]*sigma_vec[i+1]*rho,
                                sigma_vec[i]*sigma_vec[i+1]*rho, sigma_vec[i+1]^2),
                              nrow=2)
            res[i:(i+1)] <- mvtnorm::rmvnorm(1, rep(0, 2), cov_mat)
            i = i + 2
          } else {
            # to-do
            warning("Simulation of correlated NIG and GAL is not implemented yet, return only 0 for now.")
            return (res)
          }
        }
      }
    }
    res
  })
}


#' Simulate from a ngme object (possibly with replicates)
#'
#' @param object  ngme object
#' @param posterior whether to simulate from posterior sampling of latent fields
#' @param m_noise whether to add the measurement noise
#' @param seed seed
#' @param nsim ignored
#' @param ... ignored
#'
#' @return a realization of ngme object
#' @export
simulate.ngme <- function(
  object,
  posterior = TRUE,
  m_noise = TRUE,
  seed = NULL,
  nsim = NULL,
  ...
) {
  attr <- attributes(object)
  Y <- numeric(attr$fit$n_data)

  # simulate from different replicates
  replicate <- attr$fit$replicate
  for (repl in levels(replicate)) {
    repl_idx <- replicate == repl
    this_repl <- object$replicate[[repl]]
    Y[repl_idx] <- simulate_1rep(this_repl, posterior, seed)

    # add measurement noise
    if (m_noise) Y[repl_idx] <-
      Y[repl_idx] + simulate(this_repl$noise, nsim=nsim, seed=seed)
  }

  Y
}

# simulate from one replicate
simulate_1rep <- function(ngme_1rep, posterior=TRUE, seed=NULL) {
  # extract A and cbind!
  As <- list();
  for (i in seq_along(ngme_1rep$models)) {
    As[[i]] <- ngme_1rep$models[[i]]$A
  }
  A_block <- Reduce(cbind, x = As)

  if (is.null(seed)) seed <- Sys.time()
  Ws <- sampling_cpp(ngme_1rep, n=1, posterior=posterior, seed=seed)[["W"]] [[1]]

  # return A W + X beta
  as.numeric(A_block %*% Ws + ngme_1rep$X %*% ngme_1rep$feff)
}