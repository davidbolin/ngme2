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
      # optimizer = precond_sgd(),
      print_check_info = FALSE,
      n_parallel_chain = 4
    )
  )
  expect_true(out$replicates[[1]]$noise$theta_sigma - log(1.5) < 1)
  out
  traceplot(out)

  # 2. nig case
  n_obs <- 2000
  x <- rnorm(n_obs); z <- rexp(n_obs)
  # beta <- c(-2, 3)
  mu = 5; sigma = 0.8; nu = 0.1
  sim_noise <- rnig(n_obs, delta = -mu, mu = mu, nu = nu, sigma = sigma, seed=10)
  mean(sim_noise)
  # y <- z*beta[[1]] + x * beta[[2]] + sim_noise - mean(sim_noise)
  y <- sim_noise - mean(sim_noise)

  out <- ngme(
    y ~ 0, 
    data=data.frame(x=x,y=y), family = noise_nig(),
    control_opt = control_opt(
      iterations = 100,
      n_parallel_chain = 4,
      n_slope_check = 10,
      optimizer = precond_sgd(),
      print_check_info = FALSE
  ))
  out
  traceplot(out, hline=c(mu, sigma, nu, beta))
  expect_true(out$replicates[[1]]$noise$theta_sigma - log(2) < 1)
  expect_true(out$replicates[[1]]$noise$theta_mu - mu < 1)
  expect_true(out$replicates[[1]]$noise$theta_nu - log(nu) < 2)
  plot(noise_nig(mu=mu, nu=nu, sigma=sigma), out$replicates[[1]]$noise)
})

test_that("test ar1 + normal noise", {
  set.seed(100)
  n_obs <- 1000
  rho <- 0.6; mu <- -2; sigma <- 1; nu <- 1
  W = simulate(f(1:n_obs, rho=rho, model="ar1", noise=noise_nig(mu=mu,sigma=sigma,nu=nu)))[[1]]
  sd(W)

  x <- rexp(n_obs); sigma_e <- 0.5
  y = W + rnorm(n_obs, sd=sigma_e)

  fit <- ngme(y ~ 0 + f(1:n_obs, model="ar1", noise=noise_nig()),
    data=data.frame(x=x,y=y),
    control_opt = control_opt(
      rao_blackwellization = TRUE,
      iterations = 500,
      optimizer = sgd(),
      verbose=TRUE,
      n_parallel_chain = 4,
      n_slope_check = 10
  ))
  fit 
  traceplot(fit, "field1", hline=c(0.6, -2, 1, 1))
  traceplot(fit, hline=sigma_e)

  ar_result <- as.numeric(ngme_result(fit)$field1)
  score <- abs((ar_result -  c(rho, mu, sigma, nu)) / c(rho, mu, sigma, nu))
  score
  expect_true(all(score < c(0.1, 0.2, 0.3, 0.3)))

  est_sigma_e <- ngme_result(fit)$data$sigma
  expect_true(all(abs((est_sigma_e - sigma_e) / sigma_e) < 0.3))
})


test_that("test ar1 (normal) + NIG noise", {
  set.seed(100)
  n_obs <- 1000
  rho = 0.6; sigma_ar = 5
  W = simulate(f(1:n_obs, rho=rho, model="ar1", noise=noise_normal(sigma=sigma_ar)))[[1]]
  sd(W)

  x <- rexp(n_obs)
  mu = 2; sigma=0.8; nu=0.1
  sim_noise <- rnig(n_obs, delta = -mu, mu = mu, nu = nu, sigma = sigma, seed=10)
  sd(sim_noise)
  y = W + sim_noise
  mean(sim_noise)

  out <- ngme(y ~ 0 + f(1:n_obs, model="ar1", noise=noise_normal()),
    data=data.frame(x=x,y=y),
    family = noise_nig(),
    control_opt = control_opt(
      iterations = 1000,
      n_parallel_chain = 4,
      # optimizer = precond_sgd(),
      n_slope_check = 10,
      verbose=TRUE
  ))
  out
  traceplot(out, "field1", hline=c(rho, sigma_ar))
  traceplot(out, hline=c(mu, sigma, nu))

  expect_true(out$replicates[[1]]$noise$theta_sigma - log(sigma) < 2)
  expect_true(out$replicates[[1]]$noise$theta_mu - mu < 1)
  expect_true(out$replicates[[1]]$noise$theta_nu - log(nu) < 2)

  plot(noise_nig(mu=mu, nu=nu, sigma=sigma), out$replicates[[1]]$noise)
})
