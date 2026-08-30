test_that("sgld optimizer object validates inputs", {
  opt <- sgld(stepsize = 0.01, temperature = 0.5)
  expect_s3_class(opt, "ngme_optimizer")
  expect_equal(opt$method, "sgld")
  expect_equal(opt$stepsize, 0.01)
  expect_equal(opt$sgd_parameters, 0.5)

  expect_error(sgld(stepsize = 0, temperature = 1))
  expect_error(sgld(stepsize = 0.01, temperature = -1))
})

test_that("ngme runs with sgld and stores trajectories", {
  skip_on_cran()
  withr::local_seed(1)

  n <- 20
  dat <- data.frame(y = rnorm(n))

  fit <- ngme(
    y ~ 0 + f(1:n, model = ar1(), name = "ar1_field"),
    data = dat,
    family = noise_normal(),
    control_opt = control_opt(
      burnin = 2,
      iterations = 20,
      n_batch = 2,
      n_parallel_chain = 1,
      optimizer = sgld(stepsize = 0.001, temperature = 1),
      stepsize_control = poly_decay(alpha = 0.6, t0 = 10),
      R_hat_conv_check = FALSE,
      trend_std_conv_check = FALSE,
      store_traj = TRUE
    ),
    control_ngme = control_ngme(n_gibbs_samples = 1, n_post_samples = 2)
  )

  expect_s3_class(fit, "ngme")

  tr <- get_trace_trajectories(fit, name = "general", apply_transform = FALSE)
  expect_s3_class(tr, "ngme_trajectories")
  expect_true(tr$n_iterations >= 20)
  expect_true(length(tr$trajectories) >= 1)
})

test_that("sgld warns when stepsize is constant without decay", {
  expect_warning(
    control_opt(
      optimizer = sgld(stepsize = 0.01, temperature = 1)
    ),
    "constant stepsize \\(no decay\\)"
  )

  expect_no_warning(
    control_opt(
      optimizer = sgld(stepsize = 0.01, temperature = 1),
      stepsize_control = poly_decay(alpha = 0.6, t0 = 10)
    )
  )
})
