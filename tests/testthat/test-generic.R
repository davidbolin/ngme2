test_that("generic model operator K works", {
  A = matrix(1, 2, 2); B = matrix(2, 2, 2); C = matrix(3, 2, 2)
  generic_model <- generic(
    theta_K = c(alpha=1, beta=2, gamma=3),
    trans = c(alpha="identity", beta="identity", gamma="identity"),
    matrices = list(A, B, C),
    h = c(1, 1)
  )
  K = generic_model$K
  expect_equal(K, 1 * A + 2 * B + 3 * C)
})

test_that("generic model == AR1 model", {
  # Test AR1 representation
  seed <- 10; n_obs <- 10; Y <- rnorm(n_obs)
  ar1 <- ar1(1:n_obs, rho=0.5)
  g = name2fun("tanh", inv=TRUE)
  ar1$param_name
  ar1$param_trans

  generic_ar1 <- generic(
    theta_K = c(rho=g(0.5)),
    trans = c(rho="tanh"),
    matrices = list(ar1$C, ar1$G),
    h = ar1$h, 
    mesh = 1:n_obs
  )

  expect_equal(ar1$K, ar1$C * 0.5 + ar1$G)
  expect_equal(generic_ar1$K, ar1$K)
  expect_equal(generic_ar1$matrices, list(ar1$C, ar1$G))

  control <- control_opt(
    seed = seed,
    iterations = 100,
    n_parallel_chain = 4,
    stop_points = 1
  )

  fit_ar1 <- ngme(
    Y ~ 0 + f(
      1:n_obs,
      name = "my_ar",
      model = "ar1",
      rho = 0.5,
    ),
    data = data.frame(Y=Y),
    control_opt = control
  )
  fit_ar1
  est_rho_ar1 <- ar1_th2a(ngme_result(fit_ar1, "my_ar")$operator$theta_K)
  print(est_rho_ar1)
  traceplot(fit_ar1, "my_ar")

  fit_generic <- ngme(
    Y ~ 0 + f(
      1:n_obs,
      model = "generic",
      name = "generic",
      theta_K = c(rho=g(0.5)),
      trans = c(rho="tanh"),
      matrices = list(ar1$C, ar1$G),
      h = ar1$h
    ),
    data = data.frame(Y=Y),
    control_opt = control
  )
  fit_generic
  est_rho_generic <- ar1_th2a(ngme_result(fit_generic, "generic")$operator$theta_K)
  print(est_rho_generic)
  traceplot(fit_generic, "generic")

  expect_equal(est_rho_generic[[1]], est_rho_ar1[[1]])
})

test_that("generic model == Matern model", {
  mesh = fmesher::fm_mesh_1d(seq(0, 1, length.out = 10))
  matern <- matern(mesh, theta_K = 0.7)
  K = matern$K
  expect_equal(K, matern$C * exp(2*0.7) + matern$G)
  
  generic_model_matern <- generic(
    theta_K = c(x=0.7),
    trans = c(x="exp2"),
    matrices = list(matern$C, matern$G),
    h = matern$h
  )

  expect_equal(generic_model_matern$K, matern$K)
})


test_that("generic model == Matern model", {
  seed <- 10; n_obs <- 5; Y <- rnorm(n_obs)
  mesh = fmesher::fm_mesh_1d(seq(0, 1, length.out = 6))

  matern <- matern(
    mesh, theta_K = 0.7
  )
  g = name2fun("exp2", inv=TRUE)
  matern$symmetric
  matern$zero_trace

  generic_matern <- generic(
    theta_K = c(x=0.7),
    trans = c(x="exp2"),
    matrices = list(matern$C, matern$G),
    h = matern$h, 
    mesh = mesh
  )
  expect_equal(generic_matern$symmetric, matern$symmetric)
  expect_equal(generic_matern$zero_trace, matern$zero_trace)
  expect_equal(generic_matern$mesh, matern$mesh)
  expect_equal(generic_matern$h, matern$h)

  expect_equal(matern$K, matern$C * exp(2*0.7) + matern$G)
  expect_equal(generic_matern$K, matern$K)
  expect_equal(generic_matern$matrices, list(matern$C, matern$G))

  control <- control_opt(
    seed = seed,
    iterations = 100,
    n_parallel_chain = 4,
    stop_points = 1,
    # verbose = TRUE
  )

  fit_matern <- ngme(
    Y ~ 0 + f(
      1:n_obs,
      name = "my_ar",
      model = "matern",
      theta_K = 0.7,
      mesh = mesh
    ),
    data = data.frame(Y=Y),
    control_opt = control
  )
  fit_matern
  est_x_matern <- ngme_result(fit_matern, "my_ar")$operator$theta_K
  exp(est_x_matern[[1]])
  traceplot(fit_matern, "my_ar")

  fit_generic <- ngme(
    Y ~ 0 + f(
      1:n_obs,
      model = "generic",
      name = "generic",
      theta_K = c(x=0.7),
      trans = c(x="exp2"),
      matrices = list(matern$C, matern$G),
      h = matern$h,
      mesh = mesh
    ),
    data = data.frame(Y=Y),
    control_opt = control
  )
  fit_generic
  est_x_generic <- ngme_result(fit_generic, "generic")$operator$theta_K
  traceplot(fit_generic, "generic")
  exp(est_x_generic[[1]])
  expect_equal(est_x_generic[[1]], est_x_matern[[1]])
})


test_that("generic model == RW1 model", {
  # Test AR1 representation
  seed <- 10; n_obs <- 10; Y <- rnorm(n_obs)
  rw1 <- rw1(1:n_obs)

  generic_rw1 <- generic(
    matrices = list(rw1$K),
    h = rw1$h, 
    mesh = 1:n_obs
  )
  generic_rw1$K

  expect_equal(generic_rw1$K, rw1$K)
  # expect_equal(generic_rw1$matrices, list(rw1$K))

  control <- control_opt(
    seed = seed,
    iterations = 100,
    n_parallel_chain = 4,
    stop_points = 1
  )

  fit_rw1 <- ngme(
    Y ~ 0 + f(
      1:n_obs,
      name = "my_rw1",
      model = "rw1",
    ),
    data = data.frame(Y=Y),
    control_opt = control
  )
  fit_rw1
  est_sigma_rw1 <- ngme_result(fit_rw1, "my_rw1")$noise$theta_sigma

  fit_generic <- ngme(
    Y ~ 0 + f(
      1:n_obs,
      model = "generic",
      name = "generic",
      matrices = list(rw1$K),
      h = rw1$h
    ),
    data = data.frame(Y=Y),
    control_opt = control
  )
  fit_generic
  est_sigma_generic <- ngme_result(fit_generic, "generic")$noise$theta_sigma

  expect_equal(est_sigma_generic[[1]], est_sigma_rw1[[1]])
})
