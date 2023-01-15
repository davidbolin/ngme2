# This file is for testing the prediction function.
test_that("predict AR(1)", {
  n_obs <<- 500
  x1 <- rexp(n_obs); x2 <- rnorm(n_obs)
  beta <- c(-2, 4, 1)
  mu <- 1.5; sigma <- 2.3; nu <- 2; sigma_eps <- 0.8

  ar1 <- model_ar1(1:n_obs, noise=noise_nig(mu=mu, sigma=sigma, nu=nu))
  W <- simulate(ar1)

  Y <- beta[1] + x1 * beta[2] + x2 * beta[3] + W + rnorm(n_obs, sd = sigma_eps)

  # first we test estimation
  out <- ngme(
    Y ~ x1 + x2 + f(1:n_obs, model="ar1", noise = noise_nig()),
    family = noise_normal(),
    control = ngme_control(
      iteration = 500
   ),
   data = data.frame(Y = Y, x1 = x1, x2 = x2)
  )

  # compare results
  expect_true(
    sum(abs(out$beta - beta)) < 1 &&
    abs(out$noise$theta_sigma - log(sigma_eps)) < 1 &&
    abs(out$latents[[1]]$noise$theta_mu - mu) < 1 &&
    abs(out$latents[[1]]$noise$theta_sigma - log(sigma)) < 1 &&
    abs(out$latents[[1]]$noise$nu - nu) < 1
  )

  # str(out)
  # Next prediction

})
