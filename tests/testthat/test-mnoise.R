# test measurement noise

test_that("test basic mn", {
  load_all()
  n_obs <<- 500
  x <- rexp(n_obs)
  beta <- c(-2, 3)
  y <- beta[[1]] + x * beta[[2]] + rnorm(n_obs, sd = 1.5)

  # summary(lm(y~x))
  out <- ngme(y ~ x, data=list(x=x,y=y))
  expect_true(out$noise$theta_sigma - log(1.5) < 1)

  # 2. nig case
  y <- beta[[1]] + x * beta[[2]] +
    rnig(n_obs, delta = 3, mu = -3, nu = 2, sigma = 2, seed=10)

  out <- ngme(y ~ x, data=list(x=x,y=y), family = noise_nig(),
    control = ngme_control(
      iterations = 1000,
      n_slope_check = 100
    ))
  out

  # traceplot(out, "mn")
  expect_true(out$noise$theta_sigma - log(2) < 1)
  expect_true(out$noise$theta_mu - (-3) < 1)
  expect_true(out$noise$nu - nu < 1.5)

  # very close!!
  # plot(noise_nig(mu=-3, nu=2, sigma=2), out$noise)
})
