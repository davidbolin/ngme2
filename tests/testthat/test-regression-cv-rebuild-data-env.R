test_that("cross_validation(data=...) rebuilds non-stationary model without formula globals", {
  skip_if_not_installed("fmesher")
  withr::local_seed(1)

  x <- seq(0, 1, length.out = 80)
  mesh <- fmesher::fm_mesh_1d(seq(0, 1, length.out = 30))
  mesh_loc <- as.numeric(mesh$loc)
  B <- cbind(1, scale(mesh_loc, scale = FALSE))
  n_basis <- ncol(B)

  y <- sin(2 * pi * x) + rnorm(length(x), sd = 0.1)
  train_data <- data.frame(Y = y[1:60], x = x[1:60])
  test_data <- data.frame(Y = y[61:80], x = x[61:80])

  fit <- ngme(
    Y ~ 0 + f(
      x,
      model = matern(
        mesh = mesh,
        B_kappa = B,
        theta_kappa = rep(0, n_basis)
      ),
      noise = noise_normal(
        B_sigma = B,
        theta_sigma = rep(-2, n_basis)
      )
    ),
    data = train_data,
    control_opt = control_opt(estimation = FALSE)
  )

  rm(mesh, mesh_loc, B, n_basis)

  expect_no_error({
    cv <- cross_validation(
      fit,
      data = test_data,
      type = "k-fold",
      k = 2,
      N_sim = 1,
      n_gibbs_samples = 5,
      n_burnin = 2,
      print = FALSE,
      parallel = FALSE
    )
  })

  expect_equal(rownames(cv$mean.scores), "model_1")
  expect_true(all(is.finite(as.matrix(cv$mean.scores))))
})
