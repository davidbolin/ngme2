test_that("nig-normal noise", {
  n_obs <- 10
  sigma_eps <- 0.5
  alpha <- 0.5
  mu = 2; delta = -mu
  sigma <- 3
  nu <- 1
  sigma_normal <- 0.8

  # First we generate V. V_i follows inverse Gaussian distribution
  trueV <- ngme2::rig(n_obs, nu, nu, seed = 10)

  # Then generate the nig noise
  mynoise <- delta + mu*trueV + sigma * sqrt(trueV) * rnorm(n_obs)
  # Add normal noise
  mynoise <- mynoise + rnorm(n_obs, mean=0, sd=sigma_normal)

  trueW <- Reduce(function(x,y){y + alpha*x}, mynoise, accumulate = T)
  Y = trueW + rnorm(n_obs, mean=0, sd=sigma_eps)

  # Add some fixed effects
  x1 = runif(n_obs)
  x2 = rexp(n_obs)
  beta <- c(-3, -1, 2)
  X <- (model.matrix(Y ~ x1 + x2))  # design matrix
  Y = as.numeric(Y + X %*% beta)

  load_all()
  ngme_out <- ngme(
    Y ~ x1 + x2 + f(
      1:n_obs,
      name = "my_ar",
      model = "ar1",
      noise = noise_normal_nig(),
      # noise = noise_nig()
    ),
    data=data.frame(x1=x1, x2=x2, Y=Y),
    control_opt = control_opt(
      burnin = 100,
      iterations = 50,
      std_lim = 0.4,
      n_parallel_chain = 4,
      stop_points = 1,
      print_check_info = FALSE,
      seed = 3
      # verbose = T
    )
    # ,start = ngme_out
  )
  ngme_out
  traceplot(ngme_out)
  load_all()
  traceplot(ngme_out, "my_ar")
})
