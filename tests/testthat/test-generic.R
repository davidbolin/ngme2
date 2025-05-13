test_that("generic model operator K works", {
  A = matrix(1, 2, 2); B = matrix(2, 2, 2); C = matrix(3, 2, 2); D = matrix(4, 2, 2)
  generic_model <- generic(
    theta_K = c(alpha=1, beta=2, gamma=3),
    trans = c(alpha="identity", beta="exp", gamma="exp2"),
    matrices = list(A, B, C, D),
    h = c(1, 1)
  )
  K = generic_model$K
  expect_equal(as.matrix(K), 1 * A + exp(2) * B + exp(2*3) * C + D)
})

test_that("generic model == AR1 model", {
  # Test AR1 representation
  seed <- 10; n_obs <- 10; Y <- rnorm(n_obs)
  ar1 <- ar1(1:n_obs, rho=0.5)
  g = name2fun("tanh", inv=TRUE)
  ar1$param_name
  ar1$param_trans

  generic_ar1 <- generic(
    theta_K = c(x=g(0.5)), # trans(X) = rho
    trans = c(x="tanh"),
    matrices = list(ar1$C, ar1$G),
    h = ar1$h, 
    mesh = 1:n_obs
  )
  generic_ar1
  generic_ar1$param_name
  generic_ar1$param_trans

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

  expect_equal(est_rho_generic[[1]][1], est_rho_ar1[[1]][1])
})

test_that("generic model == Matern model (alpha == 2)", {
  mesh = fmesher::fm_mesh_1d(seq(0, 1, length.out = 10))

  matern <- matern(mesh, theta_K = 0.7)
  matern_2 <- matern(mesh, kappa = exp(0.7))

  expect_equal(matern$K, matern$C * exp(2*0.7) + matern$G)
  expect_equal(matern_2$K, matern$K)
  
  generic_model_matern <- generic(
    theta_K = c(x=0.7),
    trans = c(x="exp2"),
    matrices = list(matern$C, matern$G),
    h = matern$h
  )

  expect_equal(generic_model_matern$K, matern$K)

  # Fitting the matern model
  control <- control_opt(
    seed = 10,
    iterations = 10,
    n_parallel_chain = 4,
    stop_points = 1,
  )

  Y <- rnorm(n_obs)
  fit_matern <- ngme(
    Y ~ 0 + f(
      1:n_obs,
      mesh = mesh,
      model = "matern",
      alpha = 2,
    ),
    data = data.frame(Y=Y),
    control_opt = control
  )
  fit_matern

  est_theta_matern <- ngme_result(fit_matern, "field1")$operator$theta_K  
  est_theta_matern[[1]]

  fit_generic <- ngme(
    Y ~ 0 + f(
      1:n_obs,
      mesh = mesh,
      model = "generic",
      name = "generic",
      theta_K = c(theta=0),
      trans = c(theta="exp2"),
      matrices = list(matern$C, matern$G),
      h = matern$h
    ),
    data = data.frame(Y=Y),
    control_opt = control
  )

  fit_generic
  est_theta_generic <- ngme_result(fit_generic, "generic")$operator$theta_K
  expect_equal(est_theta_generic[[1]], est_theta_matern[[1]])
})

test_that("generic model == Matern model (alpha == 4)", {
  x <- seq(0, 1, length.out = 10); y <- seq(0, 1, length.out = 10)
  mesh = fmesher::fm_mesh_2d(cbind(x, y))

  kappa <- exp(0.7)
  matern <- matern(mesh, alpha = 4, theta_K = log(kappa))
  matern_2 <- matern(mesh, alpha = 4, kappa = kappa)

  # Build Cinv
  C <- matern$C
  G <- matern$G
  Cinv <- C; diag(Cinv) <- 1 / Matrix::diag(C)

  expect_equal(matern$K, kappa^4 * C + 2 * kappa^2 * G + G %*% Cinv %*% G)
  
  generic_model_matern <- generic(
    theta_K = c(theta=log(kappa)),
    trans = list(theta=c("exp4", "exp2", "null")),
    matrices = list(C, 2*G, G %*% Cinv %*% G),
    h = matern$h
  )
  expect_equal(generic_model_matern$symmetric, matern$symmetric)
  expect_equal(generic_model_matern$zero_trace, matern$zero_trace)
  expect_equal(generic_model_matern$h, matern$h)
  expect_equal(generic_model_matern$K, matern$K)

  # Fitting the matern model
  control <- control_opt(
    seed = 10,
    iterations = 100,
    n_parallel_chain = 4,
    stop_points = 1,
  )

  Y <- rnorm(n_obs)
  fit_matern <- ngme(
    Y ~ 0 + f(
      cbind(x, y),
      mesh = mesh,
      model = "matern",
      alpha = 4,
    ),
    data = data.frame(Y=Y),
    control_opt = control
  )
  fit_matern
  est_theta_matern <- ngme_result(fit_matern, "field1")$operator$theta_K
  est_theta_matern[[1]]

  fit_generic <- ngme(
    Y ~ 0 + f(
      cbind(x, y),
      mesh = mesh,
      model = "generic",
      name = "generic",
      theta_K = c(theta=0),
      trans = list(theta=c("exp4", "exp2", "null")),
      matrices = list(C, 2*G, G %*% Cinv %*% G),
      h = matern$h
    ),
    data = data.frame(Y=Y),
    control_opt = control
  )
  fit_generic
  est_theta_generic <- ngme_result(fit_generic, "generic")$operator$theta_K
  expect_equal(est_theta_generic[[1]], est_theta_matern[[1]], tolerance = 1e-4)
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

test_that("test compute coef", {
  # Create 3 random matrices for testing
  set.seed(123)
  n <- 5
  
  # Create random matrices with different patterns
  mat1 <- matrix(rnorm(n*n), nrow = n)
  mat2 <- matrix(runif(n*n), nrow = n)
  mat3 <- matrix(rpois(n*n, lambda = 2), nrow = n)
  
  # Test with the new compute_coef_new function
  theta_K <- c(param1 = 2, param2 = 3, param3 = 4)
  
  # Create trans list in the new format
  trans <- list(
    param1 = c("identity", "null", "null"),
    param2 = c("null", "identity", "null"),
    param3 = c("null", "null", "identity")
  )
  
  coef <- compute_coef_new(theta_K, trans, 3)
  expect_equal(coef, c(2, 3, 4))
  
  # Test with mixed transformations
  trans2 <- list(
    param1 = c("exp", "null", "null"),
    param2 = c("null", "identity", "null"),
    param3 = c("null", "null", "exp2")
  )
  
  coef2 <- compute_coef_new(theta_K, trans2, 3)
  expect_equal(coef2, c(exp(2), 3, exp(2*4)))
  
  # Test with a more complex case where parameters affect multiple matrices
  trans3 <- list(
    param1 = c("identity", "identity", "null"),
    param2 = c("null", "identity", "identity")
  )
  
  theta_K2 <- c(param1 = 2, param2 = 3)
  
  coef3 <- compute_coef_new(theta_K2, trans3, 3)
  expect_equal(coef3, c(2, 2*3, 3))
})

test_that("generic model with new trans format works", {
  # Create test matrices
  A <- matrix(1, 3, 3)
  B <- matrix(2, 3, 3)
  C <- matrix(3, 3, 3)
  D <- matrix(4, 3, 3)
  
  # Set parameters
  theta_K <- c(a=1, b=2, c=3)
  
  # Define transformations for each parameter across matrices
  trans <- list(
    a = c("identity", "exp", "null", "null"),
    b = c("null", "identity", "null", "null"),
    c = c("null", "null", "identity", "null")
  )
  
  # Create the generic model with the new format
  generic_model <- generic(
    theta_K = theta_K,
    trans = trans,
    matrices = list(A, B, C, D),
    h = c(1, 1, 1)
  )
  
  # Calculate expected K manually
  # K = a*A + b*exp(a)*B + c*C + D
  expected_K <- 1*A + 2*exp(1)*B + 3*C + D
  
  expect_equal(as.matrix(generic_model$K), expected_K)
  
  # Test with different transformations
  trans2 <- list(
    a = c("identity", "exp2", "null", "null"),
    b = c("null", "null", "tanh", "null"),
    c = c("null", "null", "null", "exp")
  )
  
  generic_model2 <- generic(
    theta_K = theta_K,
    trans = trans2,
    matrices = list(A, B, C, D),
    h = c(1, 1, 1)
  )
  
  # Calculate expected K manually with new transformations
  # K = a*A + exp(2*a)*B + tanh(b)*C + exp(c)*D
  expected_K2 <- 1*A + exp(2*1)*B + ar1_th2a(2)*C + exp(3)*D
  
  expect_equal(as.matrix(generic_model2$K), expected_K2)
  
  # Test default behavior with no trans provided
  generic_model3 <- generic(
    matrices = list(A, B, C, D),
    h = c(1, 1, 1)
  )
  
  # By default, should be the sum of all matrices with coefficient 1
  expected_K3 <- A + B + C + D
  
  expect_equal(as.matrix(generic_model3$K), expected_K3)
})

test_that("get_param_trans correctly extracts transformation functions", {
  # Case 1: Each parameter has exactly one transformation
  trans1 <- list(
    param1 = c("exp", "null", "null"),
    param2 = c("null", "identity", "null"),
    param3 = c("null", "null", "exp2")
  )
  
  result1 <- get_param_trans(trans1)
  expect_equal(result1$param1(2), exp(2))
  expect_equal(result1$param2(3), 3)
  expect_equal(result1$param3(4), exp(2*4))
  
  # Case 2: Some parameters have multiple transformations
  trans2 <- list(
    a = c("identity", "exp", "null", "null"),
    b = c("null", "identity", "null", "null"),
    c = c("null", "null", "identity", "null")
  )
  
  result2 <- get_param_trans(trans2)
  expect_equal(result2$a(5), 5)  # Should use identity
  expect_equal(result2$b(6), 6)
  expect_equal(result2$c(7), 7)
  
  # Case 3: Mixed cases including no transformations
  trans3 <- list(
    x = c("null", "null", "null", "null"),
    y = c("exp2", "null", "null", "null"),
    z = c("identity", "tanh", "null", "null")
  )
  
  result3 <- get_param_trans(trans3)
  expect_equal(result3$x(8), 8)  # Should use identity
  expect_equal(result3$y(2), exp(2*2))
  expect_equal(result3$z(9), 9)  # Should use identity
})
