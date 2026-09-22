test_that("Matern uses the symmetry of the mesh stiffness matrix", {
  for (alpha in c(2, 4)) {
    free <- matern(fmesher::fm_mesh_1d(1:20, boundary = "free"), alpha = alpha)
    expect_false(Matrix::isSymmetric(free$K))
    expect_false(free$symmetric)
    neumann <- matern(fmesher::fm_mesh_1d(1:20, boundary = "neumann"), alpha = alpha)
    expect_true(Matrix::isSymmetric(neumann$K))
    expect_true(neumann$symmetric)
  }
})

test_that("the CV vignette's free-boundary Matern fit has finite traces", {
  withr::local_seed(123)
  n <- 200
  spec <- f(map = seq_len(n), model = ou(theta = 0.15),
            noise = noise_nig(mu = 0.5, sigma = 2, nu = 1))
  y <- simulate(spec, seed = 123)[[1]] + rnorm(n, sd = 0.5)
  mesh <- fmesher::fm_mesh_1d(seq_len(n), boundary = "free")
  fit <- ngme(
    y ~ 0 + f(idx, model = matern(mesh = mesh, alpha = 2), noise = noise_nig()),
    data = data.frame(y = y, idx = seq_len(n)),
    control_opt = control_opt(
      seed = 456, iterations = 30, n_parallel_chain = 4,
      max_num_threads = 2, solver_backend = "cholmod",
      polish_iterations = 0, warn_no_convergence = FALSE, verbose = FALSE
    )
  )
  expect_true(all(is.finite(unlist(ngme_result(fit)))))
  cv <- cross_validation(fit, type = "custom",
    train_idx = list(1:180), test_idx = list(181:190),
    N_sim = 1, n_gibbs_samples = 20, n_burnin = 10,
    max_num_threads = 2, seed = 123, print = FALSE)
  expect_true(all(is.finite(as.matrix(cv$mean.scores))))
})
