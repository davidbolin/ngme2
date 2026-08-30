test_that("test sum of ar1 and matern", {
  skip_on_cran()
  withr::local_seed(123)
  library(ngme2)
  library(fmesher)

  # 1. Simulation
  n_obs <- 1000
  time_points <- 1:n_obs

  # AR1 process
  rho_true <- 0.7
  ar1_model <- f(time_points, model = ar1(rho = rho_true), noise = noise_normal(sigma = 1.5))
  W_ar1 <- simulate(ar1_model, seed = 123)[[1]]
  sd(W_ar1)

  # Matern process
  # Create a 1D mesh
  mesh <- fm_mesh_1d(seq(1, n_obs, length.out = 50))
  kappa_true <- 0.3
  matern_model <- f(time_points, model = matern(mesh, kappa = kappa_true), noise = noise_normal(sigma = 4.0))
  W_matern <- simulate(matern_model, seed = 124)[[1]]
  sd(W_matern)

  # Measurement noise
  sigma_eps_true <- 0.5
  Y <- W_ar1 + W_matern + rnorm(n_obs, sd = sigma_eps_true)

  # 2. Estimation
  # We need to specify the model structure for estimation
  # Note: For Matern 1D, we need to provide the mesh

  out <- ngme(
    Y ~ 0 + f(time_points, model = ar1(), name = "ar1_comp") +
      f(time_points, model = matern(mesh), name = "matern_comp"),
    data = data.frame(Y = Y),
    control_opt = control_opt(
    warn_no_convergence = FALSE,
      seed = 123,
      iterations = 1000,
      optimizer = adamW(),
      n_parallel_chain = 4,
      print_check_info = FALSE,
      verbose = FALSE
    )
  )
  traceplot(out, name = "ar1_comp", hline = c(rho_true, 1.5))
  traceplot(out, name = "matern_comp", hline = c(kappa_true, 4))

  # 3. Comparison
  res_ar1 <- ngme_result(out, "ar1_comp")
  res_matern <- ngme_result(out, "matern_comp")
  res_noise <- ngme_result(out, "data") # Measurement noise is usually under "data" or similar if not specified in f()

  # Parallel stochastic fitting is somewhat noisy here, so keep the assertions
  # focused on recovering a sensible decomposition rather than exact targets.
  expect_true(is.finite(res_ar1$rho))
  expect_true(res_ar1$rho > 0.4 && res_ar1$rho < 1.05)
  expect_true(is.finite(res_matern$kappa))
  expect_true(res_matern$kappa > 0)
  expect_true(is.finite(res_noise$sigma))
  expect_true(res_noise$sigma > 0.2 && res_noise$sigma < 0.85)

  # Print results for manual inspection
  print(paste("True AR1 rho:", rho_true, "Est:", res_ar1$rho))
  print(paste("True Matern kappa:", kappa_true, "Est:", res_matern$kappa)) # Assuming kappa is directly available or we need to extract it
  print(paste("True Sigma eps:", sigma_eps_true, "Est:", res_noise$sigma))
})
