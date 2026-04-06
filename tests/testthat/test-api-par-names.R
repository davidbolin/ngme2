test_that("ngme_model par_names drop fixed sigma entries", {
  noise <- noise_normal(
    theta_sigma = c(0, 1),
    B_sigma = diag(2),
    fix_theta_sigma = c(TRUE, FALSE)
  )

  op <- ngme2:::ngme_operator(mesh = NULL, model = "iid", K = diag(1), h = 1)

  model <- ngme2:::ngme_model(
    model = "iid",
    operator = op,
    noise = noise,
    theta_K = op$theta_K
  )

  expect_equal(model$par_names, "sigma_2")
  expect_equal(model$n_params, length(model$par_names))
})

test_that("ngme_replicate par_names stay consistent with fixed sigma", {
  noise <- noise_normal(
    theta_sigma = c(0, 1),
    B_sigma = diag(2),
    fix_theta_sigma = c(TRUE, FALSE)
  )
  ctrl <- control_ngme(beta_init = numeric(0))

  ng <- ngme_replicate(
    Y = NULL,
    X = matrix(numeric(0), nrow = 0, ncol = 0),
    noise = noise,
    models = list(),
    control_ngme = ctrl,
    standardize = FALSE
  )

  expect_equal(ng$par_names, "sigma_2")
  expect_equal(ng$n_params, length(ng$par_names))
})
