# unit test measurement noise

test_that("test basic mn", {
  load_all()
  set.seed(100)
  n_obs <<- 500
  x <- rexp(n_obs)
  beta <- c(-2, 3)
  y <- beta[[1]] + x * beta[[2]] + rnorm(n_obs, sd = 1.5)

  # summary(lm(y~x))
  out <- ngme(
    y ~ x,
    data=data.frame(x=x,y=y),
    control_opt=control_opt(
      print_check_info = FALSE
    )
  )
  out
  expect_true(out$replicates[[1]]$noise$theta_sigma - log(1.5) < 1)
  out

  # 2. nig case
  y <- beta[[1]] + x * beta[[2]] +
    rnig(n_obs, delta = 3, mu = -3, nu = 2, sigma = 2, seed=10)

  out <- ngme(y ~ x, data=data.frame(x=x,y=y), family = noise_nig(),
    control_opt = control_opt(
      iterations = 1000,
      n_slope_check = 100,
      print_check_info = FALSE
  ))
  out
  traceplot(out)

  expect_true(out$replicates[[1]]$noise$theta_sigma - log(2) < 1)
  expect_true(out$replicates[[1]]$noise$theta_mu - (-3) < 1)
  expect_true(out$replicates[[1]]$noise$nu - 2 < 2)

  # very close!!
  # plot(noise_nig(mu=-3, nu=2, sigma=2), out$noise)
})
