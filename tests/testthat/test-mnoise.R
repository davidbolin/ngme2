test_that("test basic mn", {
  set.seed(100)
  n_obs <- 1000
  x <- rexp(n_obs)
  beta <- c(-1, 2)
  y <- beta[[1]] + x * beta[[2]] + rnorm(n_obs, sd = 0.5)

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
  n_obs <- 3000
  x <- rexp(n_obs)
  beta <- c(-2, 3)
  # y <- beta[[1]] + x * beta[[2]] +
  y =  rnig(n_obs, delta = -5, mu = 5, nu = 0.05, sigma = 0.8, seed=10)

  out <- ngme(y ~ 0, data=data.frame(x=x,y=y), family = noise_nig(),
    control_opt = control_opt(
      iterations = 500,
      n_parallel_chain = 4,
      n_slope_check = 10,
      # rao_blackwellization = TRUE,
      print_check_info = FALSE,
      preconditioner = "fast"
  ))
  out
  traceplot(out, hline=c(5, 0.8, 0.05, -2, 3))

  expect_true(out$replicates[[1]]$noise$theta_sigma - log(2) < 1)
  expect_true(out$replicates[[1]]$noise$theta_mu - 5 < 1)
  expect_true(out$replicates[[1]]$noise$theta_nu - log(2) < 2)

  plot(noise_nig(mu=5, nu=0.05, sigma=0.8), out$replicates[[1]]$noise)
})

test_that("test ar1 + NIG noise", {
  set.seed(100)

  n_obs <- 2000
  W = simulate(f(1:n_obs, rho=0.6, model="ar1", noise=noise_nig(mu=-2,sigma=1,nu=1)))[[1]]

  x <- rexp(n_obs)
  beta <- c(-2, 3)
  # y <- beta[[1]] + x * beta[[2]] +
  y = W +
    # rnorm(n_obs, sd=0.1)
    rnig(n_obs, delta = -3, mu = 3, nu = 0.1, sigma = 0.8, seed=10)

  out <- ngme(y ~ 0 + f(1:n_obs, model="ar1", noise=noise_nig()),
    data=data.frame(x=x,y=y),
    family = noise_nig(),
    control_opt = control_opt(
      iterations = 2000,
      n_parallel_chain = 4,
      n_slope_check = 10,
      # rao_blackwellization = TRUE,
      print_check_info = FALSE,
      preconditioner = "fast"
  ))
  out

  traceplot(out, "field1")
  traceplot(out)
  # plot(out$replicates[[1]]$noise, noise_normal(sigma=0.1))

  traceplot(out, hline=c(3, 0.8, 0.1, -2, 3))

  expect_true(out$replicates[[1]]$noise$theta_sigma - log(2) < 1)
  expect_true(out$replicates[[1]]$noise$theta_mu - 5 < 1)
  expect_true(out$replicates[[1]]$noise$theta_nu - log(2) < 2)

  plot(noise_nig(mu=3, nu=0.1, sigma=0.8), out$replicates[[1]]$noise)
})
