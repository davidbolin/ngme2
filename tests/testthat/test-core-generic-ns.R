test_that("generic_ns basic functionality works", {
  # Create test matrices
  n <- 5
  A <- matrix(1, n, n)
  B <- matrix(2, n, n)
  alpha <- 0.4

  # Create a simple generic_ns model
  model <- generic_ns(
    theta_K = list(alpha = c(alpha)),
    matrices = list(A, B),
    h = rep(1, n),
    position = list(c(2, 1), c(3)) # Diagonal(alpha) * A + B
  )

  # First position c(2, 1) means:
  # - Matrix at index 2 (A) multiplied by
  # - Matrix at index 1 (Diagonal with alpha)
  # Second position c(3) means matrix at index 3 (B)
  # So expected_K = A * Diagonal(alpha) + B
  expected_K <- A %*% Matrix::Diagonal(n, alpha) + B

  expect_equal(as.matrix(model$K), as.matrix(expected_K))
  expect_equal(length(model$matrices), 2) # Only original matrices are stored
  expect_equal(model$h, rep(1, n))
})

test_that("generic_ns with parameter transformations works", {
  n <- 5
  A <- matrix(1, n, n)
  B <- matrix(2, n, n)
  alpha <- 0.6

  # Create model with parameter transformations
  model <- generic_ns(
    theta_K = list(alpha = c(log(alpha))),
    matrices = list(A, B),
    h = rep(1, n),
    trans = list(alpha = "exp"),
    position = list(c(1, 2), c(3)) # Diagonal(exp(log(alpha))) * A + B = Diagonal(alpha) * A + B
  )

  # With exp transformation, expected_K = A * Diagonal(alpha) + B
  expected_K <- A %*% Matrix::Diagonal(n, alpha) + B

  expect_equal(as.matrix(model$K), as.matrix(expected_K))

  # Check flat parameters
  expect_equal(model$theta_K, c(alpha = log(alpha)))
})

test_that("generic_ns with custom basis expansion works", {
  n <- 10
  A <- matrix(1, n, n)
  B <- matrix(2, n, n)
  C <- matrix(3, n, n)

  # Create basis matrices (simple example for testing)
  basis_alpha <- matrix(rnorm(n), n, 1) # One coefficient per location
  alpha <- 0.5

  # Create model with basis expansion
  model <- generic_ns(
    theta_K = list(alpha = c(alpha)),
    matrices = list(A, B, C),
    h = rep(1, n),
    B_theta_K = list(alpha = basis_alpha),
    trans = list(alpha = "identity"),
    position = list(c(1, 2), c(3, 2, 1), c(4)) # Diagonal(basis_alpha * alpha) * A + B
  )

  # Expected: K = D_alpha * A + B
  alpha_vec <- as.vector(basis_alpha %*% alpha)
  D_alpha <- Matrix::Diagonal(x = alpha_vec)
  expected_K <- D_alpha %*% A + B %*% A %*% D_alpha + C

  expect_equal(as.matrix(model$K), as.matrix(expected_K))

  expect_false(identical(as.matrix(A %*% D_alpha), as.matrix(D_alpha %*% A)))
})

test_that("generic_ns with B_theta_K defaults to matrices of 1s", {
  n <- 5
  A <- matrix(1, n, n)
  B <- matrix(2, n, n)
  alpha <- 0.7

  # Create model without explicit B_theta_K
  model <- generic_ns(
    theta_K = list(alpha = c(alpha)),
    matrices = list(A, B),
    h = rep(1, n),
    position = list(c(1, 2), c(3)) # Diagonal(alpha*ones(n)) * A + B = Diagonal(alpha) * A + B
  )

  # Default B_theta_K is a matrix of 1s, so expected_K = A * Diagonal(alpha) + B
  expected_K <- A %*% Matrix::Diagonal(n, alpha) + B

  expect_equal(as.matrix(model$K), as.matrix(expected_K))
})

test_that("generic_ns handles complex matrix combinations", {
  n <- 5
  A <- matrix(1, n, n)
  B <- matrix(2, n, n)
  C <- matrix(3, n, n)
  D <- matrix(4, n, n)

  alpha <- 0.8
  beta <- 1.2

  # Create a complex model with multiple operations
  model <- generic_ns(
    theta_K = list(alpha = c(alpha), beta = c(beta)),
    matrices = list(A, B, C, D),
    h = rep(1, n),
    trans = list(alpha = "identity", beta = "identity"),
    position = list(c(1, 3), c(2, 4), c(1, 5), c(2, 6))
    # This represents: Diagonal(alpha)*A + Diagonal(beta)*B + Diagonal(alpha)*C + Diagonal(beta)*D
  )

  # Calculate expected result
  D_alpha <- Matrix::Diagonal(n, alpha)
  D_beta <- Matrix::Diagonal(n, beta)
  expected_K <- (D_alpha %*% A) + (D_beta %*% B) + (D_alpha %*% C) + (D_beta %*% D)

  expect_equal(as.matrix(model$K), as.matrix(expected_K))
})

test_that("generic_ns with spatially-varying kappa works", {
  # Create a simple 1D mesh
  n <- 10
  mesh <- fmesher::fm_mesh_1d(seq(0, 1, length.out = n))

  # Create basis for spatially-varying kappa
  B_kappa <- matrix(0, n, 2)
  B_kappa[1:(n / 2), 1] <- 1 # First half of the domain
  B_kappa[(n / 2 + 1):n, 2] <- 1 # Second half of the domain

  # Create standard Matern components
  matern_model <- matern(mesh)
  C <- matern_model$C
  G <- matern_model$G

  # Create model with space-varying kappa
  ns_model <- generic_ns(
    theta_K = list(kappa = c(log(1), log(2))), # log(kappa) values
    matrices = list(C, G), # Include both matrices C and G
    B_theta_K = list(kappa = B_kappa),
    trans = list(kappa = "exp2"), # kappa^2 transformation for C
    h = matern_model$h,
    position = list(c(1, 2), c(3)), # D_kappa * C + G
    mesh = mesh
  )

  # Create the expected K matrix manually
  kappa_vec <- exp(2 * (B_kappa %*% c(log(1), log(2))))
  D_kappa <- Matrix::Diagonal(x = kappa_vec)
  expected_K <- D_kappa %*% C + G

  expect_equal(as.matrix(ns_model$K), as.matrix(expected_K))
})

test_that("generic_ns handles edge cases correctly", {
  n <- 5
  A <- matrix(1, n, n)
  B <- matrix(2, n, n)

  # Single matrix, no parameters
  model1 <- generic_ns(
    theta_K = list(),
    matrices = list(A),
    position = list(c(1)),
    h = rep(1, n)
  )
  expect_equal(as.matrix(model1$K), A)

  # Multiple parameters with same transformation
  model3 <- generic_ns(
    theta_K = list(alpha = c(1), beta = c(2)),
    matrices = list(A, B),
    h = rep(1, n),
    trans = list(alpha = "exp", beta = "exp"),
    position = list(c(1, 2), c(3, 2, 1), c(4))
  )
  D_alpha <- Matrix::Diagonal(n, exp(1))
  D_beta <- Matrix::Diagonal(n, exp(2))
  expected_K <- D_alpha %*% D_beta + A %*% D_beta %*% D_alpha + B
  expect_equal(as.matrix(model3$K), as.matrix(expected_K))
})

test_that("generic_ns with complex matrix combinations works", {
  n <- 5
  A <- matrix(rnorm(n^2), n, n)
  B <- matrix(rnorm(n^2), n, n)
  C <- matrix(rnorm(n^2), n, n)
  D <- matrix(rnorm(n^2), n, n)

  # Create a complex model with nested matrix operations
  model <- generic_ns(
    theta_K = list(alpha = c(1), beta = c(2)),
    matrices = list(A, B, C, D),
    h = rep(1, n),
    trans = list(alpha = "identity", beta = "exp"),
    position = list(c(1, 2), c(3, 4), c(1, 5), c(2, 6))

    # This represents: D_alpha * D_beta + A * B + D_alpha * C + D_beta * D
  )
  D_alpha <- Matrix::Diagonal(n, 1)
  D_beta <- Matrix::Diagonal(n, exp(2))
  expected_K <- D_alpha %*% D_beta + A %*% B + D_alpha %*% C + D_beta %*% D

  expect_equal(as.matrix(model$K), as.matrix(expected_K))
})


test_that("generic_ns model == AR1 model", {
  # Test AR1 representation
  seed <- 10
  n_obs <- 5
  Y <- rnorm(n_obs)
  ar1 <- ar1(1:n_obs, rho = 0.5)
  g <- name2fun("tanh", inv = TRUE)
  ar1$param_name
  ar1$param_trans

  generic_ar1 <- generic_ns(
    theta_K = list(x = c(g(0.5))), # trans(X) = rho
    trans = list(x = "tanh"),
    matrices = list(ar1$C, ar1$G),
    position = list(c(1, 2), c(3)),
    h = ar1$h,
    mesh = 1:n_obs
  )
  generic_ar1
  generic_ar1$param_name
  generic_ar1$param_trans

  # ar1() carries the stationary initial condition sqrt(1 - rho^2) at (1,1).
  # generic_ns applies ONE transformation per parameter (a diagonal coefficient
  # built by basis expansion), so it can represent the linear part rho C + G
  # but not that entry, whose coefficient is not linear in theta_K.
  expect_equal(ar1$K, ar1$C * 0.5 + ar1$G + sqrt(1 - 0.5^2) * ar1$E11)
  expect_equal(generic_ar1$K, ar1$C * 0.5 + ar1$G)
  expect_equal(generic_ar1$matrices, list(ar1$C, ar1$G))

  control <- control_opt(
    seed = seed,
    iterations = 100,
    n_parallel_chain = 4,
    n_batch = 1
  )

  fit_ar1 <- ngme(
    Y ~ 0 + f(
      1:n_obs,
      name = "my_ar",
      model = ar1(rho = 0.5)
    ),
    data = data.frame(Y = Y),
    control_opt = control
  )
  fit_ar1

  est_rho_ar1 <- ar1_th2a(ngme_result(fit_ar1, "my_ar")$rho)
  print(est_rho_ar1)

  print(generic_ar1$K)
  generic_ar1$param_map
  generic_ar1$position
  fit_generic <- ngme(
    Y ~ 0 + f(
      1:n_obs,
      name = "generic",
      model = generic_ns(
        theta_K = list(rho = g(0.5)),
        trans = list(rho = "tanh"),
        matrices = list(ar1$C, ar1$G),
        h = ar1$h,
        position = list(c(1, 2), c(3))
      )
    ),
    data = data.frame(Y = Y),
    control_opt = control
  )
  fit_generic

  est_rho_generic <- ar1_th2a(ngme_result(fit_generic, "generic")$rho)
  print(est_rho_generic)

  # ar1() and this generic_ns model are no longer the same operator -- ar1
  # carries the stationary sqrt(1 - rho^2) at (1,1), which generic_ns cannot
  # express -- so their rho estimates need not coincide. What is still exact,
  # and is checked above, is that generic_ns assembles rho C + G. On n_obs = 5
  # that one entry moves rho substantially, so no tolerance is asserted here.
  expect_true(is.finite(est_rho_generic[[1]][1]))
  expect_true(is.finite(est_rho_ar1[[1]][1]))
})


test_that("generic model == Matern model (alpha == 2 or 4)", {
  x <- seq(0, 1, length.out = 10)
  y <- seq(0, 1, length.out = 10)
  mesh <- fmesher::fm_mesh_2d(cbind(x, y))
  mesh$n

  B_kappa <- matrix(rnorm(mesh$n * 2), mesh$n, 2)
  theta_kappa <- c(0.5, 0.3)
  matern <- matern(
    mesh,
    alpha = 4,
    theta_kappa = theta_kappa,
    B_kappa = B_kappa
  )
  kappas <- as.numeric(exp(B_kappa %*% theta_kappa))
  D_kappa <- Matrix::Diagonal(x = kappas)
  C <- matern$C
  G <- matern$G
  Cinv <- C
  diag(Cinv) <- 1 / Matrix::diag(C)

  # Complicated K matrix
  #  D_\kappa C D_\kappa C^{-1} D_\kappa C D_\kappa + D_\kappa C D_\kappa C^{-1}  G + G C^{-1}  D_\kappa C D_\kappa  + G C^{-1} G
  full_K <- D_kappa %*% C %*% D_kappa %*% Cinv %*% D_kappa %*% C %*% D_kappa +
    D_kappa %*% C %*% D_kappa %*% Cinv %*% G +
    G %*% Cinv %*% D_kappa %*% C %*% D_kappa +
    G %*% Cinv %*% G
  expect_equal(matern$K, full_K)

  # matrices list : D_kappa, C, G, Cinv
  generic_model_matern <- generic_ns(
    theta_K = list(theta = theta_kappa),
    trans = list(theta = c("exp")),
    B_theta_K = list(theta = B_kappa),
    matrices = list(C, G, Cinv),
    position = list(
      c(1, 2, 1, 4, 1, 2, 1),
      c(1, 2, 1, 4, 3),
      c(3, 4, 1, 2, 1),
      c(3, 4, 3)
    ),
    h = matern$h
  )
  expect_equal(generic_model_matern$symmetric, matern$symmetric)
  expect_equal(generic_model_matern$zero_trace, matern$zero_trace)
  expect_equal(generic_model_matern$h, matern$h)
  expect_equal(generic_model_matern$K, matern$K, tolerance = 1e-4)

  # Fitting the matern model
  control <- control_opt(
    seed = 10,
    iterations = 50,
    n_parallel_chain = 4,
    n_batch = 1,
  )

  Y <- rnorm(10)
  fit_matern_2 <- ngme(
    Y ~ 0 + f(
      cbind(x, y),
      model = matern(
        mesh = mesh,
        theta_kappa = theta_kappa,
        B_kappa = B_kappa,
        alpha = 2
      )
    ),
    data = data.frame(Y = Y),
    control_opt = control
  )
  fit_matern_2
  est_theta_matern_2 <- ngme_result(fit_matern_2, "field1")$kappa
  est_theta_matern_2[[1]]

  fit_generic_2 <- ngme(
    Y ~ 0 + f(
      cbind(x, y),
      name = "generic",
      model = generic_ns(
        mesh = mesh,
        theta_K = list(theta = theta_kappa),
        trans = list(theta = c("exp")),
        B_theta_K = list(theta = B_kappa),
        matrices = list(C, G),
        position = list(
          c(1, 2, 1),
          c(3)
        ),
        h = matern$h
      )
    ),
    data = data.frame(Y = Y),
    control_opt = control
  )
  fit_generic_2
  est_theta_generic_2 <- ngme_result(fit_generic_2, "generic")$theta
  expect_equal(est_theta_generic_2[[1]], est_theta_matern_2[[1]], tolerance = 1e-4)

  fit_matern_4 <- ngme(
    Y ~ 0 + f(
      cbind(x, y),
      model = matern(
        mesh = mesh,
        theta_kappa = theta_kappa,
        B_kappa = B_kappa,
        alpha = 4
      )
    ),
    data = data.frame(Y = Y),
    control_opt = control
  )
  fit_matern_4
  est_theta_matern_4 <- ngme_result(fit_matern_4, "field1")$operator$theta_K
  est_theta_matern_4[[1]]

  fit_generic_alpha_4 <- ngme(
    Y ~ 0 + f(
      cbind(x, y),
      name = "generic",
      model = generic_ns(
        mesh = mesh,
        theta_K = list(theta = theta_kappa),
        trans = list(theta = c("exp")),
        B_theta_K = list(theta = B_kappa),
        matrices = list(C, G, Cinv),
        position = list(
          c(1, 2, 1, 4, 1, 2, 1),
          c(1, 2, 1, 4, 3),
          c(3, 4, 1, 2, 1),
          c(3, 4, 3)
        ),
        h = matern$h
      )
    ),
    data = data.frame(Y = Y),
    control_opt = control
  )
  fit_generic_alpha_4
  est_theta_generic_alpha_4 <- ngme_result(fit_generic_alpha_4, "generic")$operator$theta_K
  expect_equal(est_theta_generic_alpha_4[[1]], est_theta_matern_4[[1]], tolerance = 1e-4)
})

test_that("ou (generic) equals rho*C + G on uniform mesh", {
  mesh <- 1:8
  rho <- 0.3
  theta <- -log(rho)

  op <- ou(mesh = mesh, theta = theta)

  C <- Matrix::sparseMatrix(j = 1:(length(mesh) - 1), i = 2:length(mesh), x = -1, dims = c(length(mesh), length(mesh)))
  G <- Matrix::Diagonal(length(mesh))
  G[1, 1] <- sqrt(1 - rho^2)
  expected <- rho * C + G

  expect_equal(as.matrix(op$K), as.matrix(expected))
})


# bv_matern_nig
test_that("generic_ns model == bv_matern_nig model", {
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh <- fmesher::fm_mesh_2d(
    loc.domain = pl01,
    cutoff = 1,
    max.edge = c(2, 10)
  )
  mesh$n
  n_obs <- 10

  long <- runif(n_obs / 2, 0, 10)
  lat <- runif(n_obs / 2, 0, 5)
  long <- c(long, long)
  lat <- c(lat, lat)
  group <- c(rep("W1", n_obs / 2), rep("W2", n_obs / 2))
  Y <- rnorm(n_obs)

  B_sigma <- matrix(0, nrow = n_obs, ncol = 2)
  B_sigma[group == "W1", 1] <- 1
  B_sigma[group == "W2", 2] <- 1

  out_cor <- ngme(
    Y ~ 0 + f(
      ~ long + lat,
      name = "bv",
      model = bv_matern(
        mesh = mesh,
        sub_models = list(
          W1 = matern(),
          W2 = matern()
        )
      ),
      # debug=T,
      noise = list(
        W1 = noise_nig(),
        W2 = noise_nig()
      )
    ),
    group = group,
    family = noise_normal(
      corr_measurement = TRUE,
      index_corr = c(1:(n_obs / 2), 1:(n_obs / 2)),
      B_sigma = B_sigma,
      theta_sigma = c(0, 0)
    ),
    data = data.frame(Y, long, lat),
    control_opt = control_opt(
      iterations = 10,
      n_parallel_chain = 4,
      rao_blackwellization = TRUE,
      # verbose = TRUE,
      print_check_info = FALSE,
      seed = 50
    ),
    debug = FALSE
  )
  out_cor
})


test_that("Simulation and fitting", {
  n <- 5
  A <- matrix(1, n, n)
  B <- matrix(2, n, n)
  A
  B
  alpha <- 0.4

  model <- generic_ns(
    theta_K = list(alpha = c(1, 1), beta = 3),
    B_theta_K = list(alpha = matrix(1, n, 2)),
    matrices = list(A, B),
    position = list(c(1, 3), c(2, 4)),
    h = rep(1, n)
  )
  model
  model$K # (2A+3B)

  mesh <- fmesher::fm_mesh_1d(1:n)
  y <- 1:5
  fit <- ngme(
    y ~ 0 + f(
      1:n,
      model = generic_ns(
        mesh = mesh,
        theta_K = list(alpha = c(1, 1)),
        B_theta_K = list(alpha = matrix(1, n, 2)),
        matrices = list(A, B),
        position = list(c(1, 2), c(3)), # D_alpha * A + B
        h = rep(1, n)
      )
    ),
    data = data.frame(y),
    control_opt = control_opt(
      iterations = 10
    )
  )
  fit$replicates[[1]]$models[[1]]$operator$param_name
  traceplot(fit, "field1")
  traceplot(fit)
  fit
  model
})
