test_that("a non-finite optimizer step is reported as an R error", {
  withr::local_seed(123)
  n <- 200
  spec <- f(map = seq_len(n), model = ou(theta = 0.15),
            noise = noise_nig(mu = 0.5, sigma = 2, nu = 1))
  y <- simulate(spec, seed = 123)[[1]] + rnorm(n, sd = 0.5)
  mesh <- fmesher::fm_mesh_1d(seq_len(n), boundary = "free")
  # Deliberately overflow a finite SGD step to exercise the error path on
  # worker threads without depending on a numerically unstable model.
  expect_error(
    ngme(y ~ 0 + f(idx, model = matern(mesh = mesh, alpha = 2, kappa = 1),
                   noise = noise_nig()),
      data = data.frame(y = y, idx = seq_len(n)),
      control_opt = control_opt(
        seed = 456, iterations = 100, n_parallel_chain = 4,
        optimizer = sgd(stepsize = .Machine$double.xmax),
        max_num_threads = 2, solver_backend = "cholmod",
        polish_iterations = 0, warn_no_convergence = FALSE,
        verbose = FALSE
      )
    ),
    "Non-finite optimizer step"
  )
})
