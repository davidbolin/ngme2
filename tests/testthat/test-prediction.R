# This file is for testing the prediction function.

test_that("test predict general", {
  set.seed(10)
  n_obs <<- 500
  x1 <- rexp(n_obs); x2 <- rnorm(n_obs)
  beta <- c(-2, 4, 1)
  alpha <- 0.75
  mu <- 1.5; sigma <- 2.3; nu <- 2; sigma_eps <- 0.8

  ar1 <- model_ar1(1:n_obs, alpha=alpha, noise=noise_nig(mu=mu, sigma=sigma, nu=nu))
  W <- simulate(ar1)
  Y <- beta[1] + x1 * beta[2] + x2 * beta[3] + W + rnorm(n_obs, sd = sigma_eps)

  # make 1/10 observation NA
  idx_NA <- sample(1:n_obs, size=n_obs/10)
  Y2 <- Y; Y2[idx_NA] <- NA

  # first we test estimation
  out <- ngme(
    Y2 ~ x1 + x2 + f(1:n_obs, model="ar1", noise = noise_nig(
      # fix_V = TRUE, V = attr(W, "noise")$V
    )),
    family = noise_normal(),
    control = control_opt(
      n_parallel_chain = 4,
      estimation = T,
      iteration = 1000,
      print_check_info = FALSE
   ),
   seed = 100,
   data = data.frame(Y2 = Y2, x1 = x1, x2 = x2)
  )
  out

  # evaluation
  traceplot(out, "field1")

  # compare noise
  plot(out$latents[[1]]$noise, attr(W, "noise"))

  # see the prediction
  pd <- attr(out, "prediction")
  Y[pd$index_NA] - pd$lp[pd$index_NA]

  # compare results
  expect_true(sum(abs(out$beta - beta)) < 2)
  expect_true(abs(out$noise$theta_sigma - log(sigma_eps)) < 1)
  expect_true(abs(out$latents[[1]]$noise$theta_mu - mu) < 1)
  expect_true(abs(out$latents[[1]]$noise$theta_sigma - log(sigma)) < 1)
  expect_true(abs(out$latents[[1]]$noise$nu - nu) < 5)

  # str(out)
  # Next is prediction
  # predict(out, locs)
})

test_that("test predict with NA", {
  set.seed(10)
  n_obs <<- 50
  x1 <- rexp(n_obs); x2 <- rnorm(n_obs)
  beta <- c(-2, 4, 1)
  alpha <- 0.75
  mu <- 1.5; sigma <- 2.3; nu <- 2; sigma_eps <- 0.8

  ar1 <- model_ar1(1:n_obs, alpha=alpha, noise=noise_nig(mu=mu, sigma=sigma, nu=nu))
  W <- simulate(ar1)
  Y <- beta[1] + x1 * beta[2] + x2 * beta[3] + W + rnorm(n_obs, sd = sigma_eps)

  # make 1/10 observation NA
  idx_NA <- sample(1:n_obs, size=n_obs/10)
  Y2 <- Y; Y2[idx_NA] <- NA

  # first we test estimation
  out <- ngme(
    Y2 ~ x1 + x2 + f(1:n_obs, model="rw1", noise = noise_nig(
      # fix_V = TRUE, V = attr(W, "noise")$V
    )),
    family = noise_normal(),
    control = control_opt(
      n_parallel_chain = 4,
      estimation = T,
      iteration = 500,
      print_check_info = FALSE
   ),
   seed = 100,
   data = data.frame(Y2 = Y2, x1 = x1, x2 = x2)
  )
  out$X_pred
  out$latents[[1]]$A_pred
  # at NA location
  predict(out)

  # evaluation
  traceplot(out, "field1")

  # compare noise
  plot(out$latents[[1]]$noise, attr(W, "noise"))

  # see the prediction
  pd <- attr(out, "prediction")
  Y[pd$index_NA] - pd$lp[pd$index_NA]

  # compare results
  expect_true(sum(abs(out$beta - beta)) < 3)
  expect_true(abs(out$noise$theta_sigma - log(sigma_eps)) < 2)
  expect_true(abs(out$latents[[1]]$noise$theta_mu - mu) < 2)
  expect_true(abs(out$latents[[1]]$noise$theta_sigma - log(sigma)) < 2)
  expect_true(abs(out$latents[[1]]$noise$nu - nu) < 5)
})
