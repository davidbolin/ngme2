# test rw related model

test_that("simulate and estimate of rw with NIG", {
  n_obs <- 500
  x <- rexp(n_obs, rate = 2)
  mu <- -3; sigma <- 5; nu <- 2; sigma_eps <- 0.8
  my_rw <- f(x, model="rw1", noise=noise_nig(mu=mu, sigma=sigma, nu=nu))

  W <- simulate(my_rw, seed = 3)[[1]]
  Y <- as.numeric(my_rw$A %*% W) + rnorm(n=length(W), sd=sigma_eps)

  out <- ngme(
    Y ~ 0 + f(x,
      model="rw1",
      name="rw",
      noise=noise_nig(
        # fix_nu = TRUE, nu = 2,
        # fix_theta_sigma = TRUE, sigma = sigma,
        # fix_V = TRUE, V = V
      ),
      # fix_W = TRUE, W = W,
      debug = FALSE
    ),
    data = data.frame(Y = Y),
    control_opt = control_opt(
      seed = 12,
      estimation = T,
      iterations = 500,
      n_parallel_chain = 4,
      print_check_info = F,
      preconditioner = "fast",
      # preconditioner = "none",
      verbose = F
    )
    # debug = TRUE
  )

  traceplot(out)
  traceplot(out, "rw")
  out$replicates[[1]]$models[[1]]$n_theta_K
  plot(out$replicates[[1]]$models[[1]]$noise,
    noise_nig(mu=mu, sigma=sigma, nu=nu))
})

############################## AR1 case
test_that("test estimation of basic ar with normal measurement noise", {
  n_obs <- 800
  mu <- 3; sigma <- 2; nu <- 1; sigma_eps <- 0.8
  ar1 <- f(1:n_obs, model="ar1", rho=0.7, noise=noise_nig(mu=mu, sigma=sigma, nu=nu))

  W <- simulate(ar1)[[1]]
  x1 = rnorm(n_obs);
  # Y <-  -3 + x1 * 2  + rnorm(n_obs, sd = sigma_eps)
  Y <- W + rnorm(n_obs, sd = sigma_eps)

# plot(Y, type="l")
  out <- ngme(
    Y ~  0 +
    + f(1:n_obs,
      model="ar1",
      name="ar",
      # replicate = rep(1:6, each=100),
      noise=noise_nig(
        # fix_nu = T, nu = 2,
        # h = ar1$noise$h,
        # fix_V = TRUE, V = attr(W, "noise")$V
      ),
      control=control_f(numer_grad = F),
      debug = FALSE,
      # fix_W = T, W = W
    ),
    data = data.frame(Y = Y, x1=x1),
    control_ngme = control_ngme(
      # n_gibbs_samples = 5
    ),
    control_opt = control_opt(
      estimation = T,
      burnin = 100,
      iterations = 500,
      n_parallel_chain = 4,
      print_check_info = F,
      verbose = F
    ),
    # start = out,
    debug = F
  )
  out
  traceplot(out, "ar")
  traceplot(out)
plot(attr(W, "noise"), out$replicates[[1]]$models[[1]]$noise)

  out$replicates[[1]]$models[[1]]$theta_K

  predict(out, list(ar=801:900))
  # cross_validation(out, type="loo")
  expect_true(TRUE)
})

