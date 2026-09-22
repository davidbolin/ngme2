test_that("fractional Matern estimation and prediction work without INLA", {
  skip_if_not_installed("rSPDE")
  # Exercise the rational operator and a free smoothness parameter. The
  # independent INLAbru comparison lives in tests/demo: it does not test ngme2
  # and requires a working external INLA binary as well as R packages.
  x <- seq(0, 1, length.out = 40)
  mesh <- fmesher::fm_mesh_1d(x)
  withr::local_seed(2024)
  y <- sin(2 * pi * x) + rnorm(length(x), sd = 0.1)
  fit <- ngme(
    y ~ 0 + f(x,
      model = matern(mesh = mesh, kappa = 6, alpha = 2.3,
                     rational_order = 1, fix_alpha = FALSE),
      noise = noise_normal(sigma = 2)
    ),
    data = data.frame(x = x, y = y),
    family = noise_normal(sigma = 0.1, fix_sigma = TRUE),
    control_opt = control_opt(
      seed = 2024, iterations = 30, burnin = 10,
      n_parallel_chain = 1, max_num_threads = 1,
      polish_iterations = 0, warn_no_convergence = FALSE,
      verbose = FALSE
    )
  )
  expect_s3_class(fit, "ngme")
  pars <- unlist(ngme_result(fit))
  expect_true(length(pars) > 0)
  expect_true(all(is.finite(pars)))
  pred <- predict(fit, map = list(field1 = x), seed = 2024,
                  sampling_size = 20, burnin_size = 10)
  expect_length(pred$mean, length(y))
  expect_true(all(is.finite(pred$mean)))
})

test_that("fractional Matern supports a small two-dimensional mesh", {
  skip_if_not_installed("rSPDE")
  withr::local_seed(2024)
  loc <- as.matrix(expand.grid(x = seq(0, 1, length.out = 5),
                              y = seq(0, 1, length.out = 5)))
  mesh <- fmesher::fm_mesh_2d(loc = loc, max.edge = c(0.3, 0.6))
  y <- sin(2 * pi * loc[, 1]) + rnorm(nrow(loc), sd = 0.1)
  fit <- ngme(
    y ~ 0 + f(loc,
      model = matern(mesh = mesh, kappa = 6, alpha = 2.3,
                     rational_order = 1, fix_alpha = FALSE),
      noise = noise_normal(sigma = 2)
    ),
    data = data.frame(y = y),
    family = noise_normal(sigma = 0.1, fix_sigma = TRUE),
    control_opt = control_opt(
      seed = 2024, iterations = 20, burnin = 10,
      n_parallel_chain = 1, max_num_threads = 1,
      polish_iterations = 0, warn_no_convergence = FALSE,
      verbose = FALSE
    )
  )
  expect_s3_class(fit, "ngme")
  expect_true(all(is.finite(unlist(ngme_result(fit)))))
  pred <- predict(fit, map = list(field1 = loc), seed = 2024,
                  sampling_size = 20, burnin_size = 10)
  expect_length(pred$mean, length(y))
  expect_true(all(is.finite(pred$mean)))
})
