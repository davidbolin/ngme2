test_that("general model", {
  mesh2d <- fmesher::fm_mesh_2d(
    loc.domain = cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5),
    max.edge = c(0.3, 10),
    cutoff = 0.1
  )
  mesh2d$n
  # generate random loc for each year
  n_obs <- c(102, 85, 120, 105, 109, 100, 30) # observation for each year
  year <- rep(2001:2007, times = n_obs)
  # 2d coordinate
  x <- runif(sum(n_obs)) * 10;
  y <- runif(sum(n_obs)) * 5
  loc <- cbind(x, y)
  mesh <- mesh2d

  # set the model for simulation
  true_model <- ngme2::f(
    map = list(year, ~x+y),
    model = "tp",
    first = list(model="ar1", rho = 0.5),
    second = list(model="matern", mesh = mesh2d),
    noise = noise_nig(mu=-2, sigma=1, nu=2)
  )
  dim(true_model$A)
  n_Y <- dim(true_model$A)[1]
  true_model$operator$first
  true_model$operator$second
  # W <- simulate(true_model)[[1]]

  W <- rnorm(n_Y, mean = 0, sd = 0.1)
  Y_obs <- W + rnorm(n_Y, sd = 0.5)
  df <- data.frame(year, x, y, Y_obs)

  out3 <- ngme(
    Y ~ 0 + f(loc,
        model="general",
        name="spde",
        # mesh = mesh,
        theta_K = c(rho = 0.5),
        trans = c(rho = "tanh"),
        theta_K2 = c(kappa = 2),
        trans2 = c(kappa = "exp2"),
        h = true_model$operator$first$h,
        h2 = true_model$operator$second$h,
        matrices = list(
          true_model$operator$first$C,
          true_model$operator$first$G
        ),
        matrices2 = list(
          true_model$operator$second$C,
          true_model$operator$second$G
        ),
        interact = "kronecker",
        A = true_model$A,
        noise=noise_normal(),
        control = control_f(),
        debug = F
      ),
      data = data.frame(Y = Y_obs),
      control_opt = control_opt(
        estimation = T,
        iterations = 20,
        verbose = TRUE,
        optimizer = adam(),
        rao_blackwellization = TRUE,
        n_parallel_chain = 1,
        max_num_threads = 1,
        print_check_info = F,
        std_lim = 0.01
      ),
    # start=out,
    debug = F
  )
  out3

  traceplot(out3)
  traceplot(out3, "spde")
})
