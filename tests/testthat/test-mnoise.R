test_that("test basic mn", {
  set.seed(100)
  n_obs <- 500
  x <- rexp(n_obs)
  beta <- c(-2, 3)
  y <- beta[[1]] + x * beta[[2]] + rnorm(n_obs, sd = 3)

  # summary(lm(y~x))
  out <- ngme(
    y ~ x,
    data=data.frame(x=x,y=y),
    control_opt=control_opt(
      iterations = 100,
      # verbose = T,
      print_check_info = FALSE,
      n_parallel_chain = 1,
      preconditioner = "fast"
    )
  )
  expect_true(out$replicates[[1]]$noise$theta_sigma - log(1.5) < 1)
  out
  traceplot(out)

  # 2. nig case
  y <- beta[[1]] + x * beta[[2]] +
    rnig(n_obs, delta = -3, mu = 3, nu = 2, sigma = 2, seed=10)

  out <- ngme(y ~ x, data=data.frame(x=x,y=y), family = noise_nig(),
    control_opt = control_opt(
      iterations = 500,
      n_slope_check = 10,
      print_check_info = FALSE,
      preconditioner = "fast"
  ))
  out
  traceplot(out)

  expect_true(out$replicates[[1]]$noise$theta_sigma - log(2) < 1)
  expect_true(out$replicates[[1]]$noise$theta_mu - (-3) < 1)
  expect_true(out$replicates[[1]]$noise$nu - 2 < 2)

  plot(noise_nig(mu=3, nu=2, sigma=2), out$replicates[[1]]$noise)
})
