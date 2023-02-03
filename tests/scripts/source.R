  load_all()
  n_obs <<- 500
  mu <- -3; sigma <- 5; nu <- 2; sigma_eps <- 0.8
  h <- rexp(n_obs)
  # h <- rep(1, n_obs)
  loc <<- c(0, cumsum(h))

  V <- rig(n_obs, a=nu, b=nu*h^2)
  # V <- h
  dW <- -mu*h + mu * V + sigma * sqrt(V) * rnorm(n_obs) # type-G noise

  W <- c(0, cumsum(dW))
  Y <- W + rnorm(n=length(W), sd=sigma_eps)

  # check model specification
  my_rw <- model_rw(loc, order=1)
  expect_true(all(my_rw$K == my_rw$C + my_rw$G))
  expect_true(all(as.numeric(my_rw$K %*% W) - dW < 1e-5))
  expect_true(all(my_rw$nosie$h - h < 1e-5))

  # first we test the gradient of mu
  out <- ngme(
    Y ~ 0 + f(loc,
      model="rw1",
      name="rw1",
      noise=noise_nig(
        # fix_nu = TRUE, nu = 2,
        # fix_theta_sigma = TRUE, sigma = sigma,
        # fix_V = TRUE, V = V
      ),
      # fix_W = TRUE, W = W,
      debug = FALSE
    ),
    data = list(Y = Y),
    contro = ngme_control(
      estimation = T,
      iterations = 500,
      n_parallel_chain = 4,
      print_check_info = TRUE
    ),
    debug = TRUE
  )
  out