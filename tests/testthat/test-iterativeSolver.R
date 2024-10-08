test_that("iterativeSolver works", {
  set.seed(123)
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh_2k <- fmesher::fm_mesh_2d(
    loc.domain = pl01, cutoff = 0.1,
    max.edge = c(0.3, 10)
  )
  mesh_2k$n
  mesh_8k <- fmesher::fm_mesh_2d(
    loc.domain = pl01, cutoff = 0.08,
    max.edge = c(0.1, 10)
  )
  mesh_8k$n
  mesh_17k <- fmesher::fm_mesh_2d(
    loc.domain = pl01, cutoff = 0.01,
    max.edge = c(0.1, 10)
  )
  mesh_17k$n

  mesh <- mesh_17k
  # plot(mesh)
  mesh$n

  n_obs <- 1000
  loc <- cbind(runif(n_obs, 0, 10), runif(n_obs, 0, 5))
  # plot(mesh); points(loc)
  true_noise = noise_nig(mu=-2, sigma=1, nu=0.5)
  true_noise = noise_normal(sigma=4)
  # plot(true_noise)

  true_model <- f(
  map = loc,
  model="matern",
  theta_K = log(4),
  mesh = mesh,
  noise = true_noise
  )

  W <- rnorm(n_obs)
  # W <- simulate(true_model)[[1]]
  attr(W, "noise")
  Y <- W + rnorm(n_obs, sd=0.5)

  # make bubble plot
  # sp_obj <- as.data.frame(mesh$loc); sp_obj[, 3] <- W
  # names(sp_obj) <- c("s1", "s2", "y")
  # coordinates(sp_obj) <- ~ s1 + s2
  # bubble(sp_obj, zcol=3)
  # range(mesh$loc[, 1]); range(mesh$loc[, 2])

  # Matern case
  # load_all()
  out <- ngme(
  Y ~ 0 + f(loc,
    model="matern",
    name="spde",
    mesh = mesh,
    # noise=noise_nig(),
    noise=noise_normal(),
    control = control_f(
      # iterative_solver = TRUE,
      numer_grad = FALSE
    ),
    debug = F
  ),
  data = data.frame(Y = Y),
  control_ngme = control_ngme(
    # iterative_solver = TRUE
  ),
  control_opt = control_opt(
    estimation = T,
    iterations = 20,
    optimizer = adam(),
    rao_blackwellization = TRUE,
    n_parallel_chain = 4,
    print_check_info = F,
    verbose = T,
    std_lim = 0.01
  ),
  # start=out,
  debug = F
)

  out
  traceplot(out, "spde")

})
