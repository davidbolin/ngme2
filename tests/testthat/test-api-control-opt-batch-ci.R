test_that("control_opt_batch_ci sets recommended CI defaults", {
  ctl <- control_opt_batch_ci(
    iterations = 400,
    n_batch = 10,
    n_parallel_chain = 3,
    alpha = 0.6,
    t0 = 5,
    schedule_burnin_iter = 12
  )

  expect_s3_class(ctl, "control_opt")
  expect_equal(ctl$iterations, 400)
  expect_equal(ctl$n_batch, 10)
  expect_equal(ctl$n_parallel_chain, 3)

  expect_true(isTRUE(ctl$store_traj))
  expect_false(isTRUE(ctl$trend_std_conv_check))
  expect_false(isTRUE(ctl$R_hat_conv_check))
  expect_equal(ctl$stepsize_decay, "none")

  expect_equal(ctl$stepsize_schedule, "poly")
  expect_equal(ctl$stepsize_schedule_alpha, 0.6)
  expect_equal(ctl$stepsize_schedule_t0, 5)
  expect_equal(ctl$stepsize_schedule_burnin_iter, 12)
})


test_that("control_opt_batch_ci allows overriding defaults via dots", {
  ctl <- control_opt_batch_ci(
    store_traj = FALSE,
    R_hat_conv_check = TRUE,
    stepsize_control = batch_decay(),
    optimizer = adam(stepsize = 0.01)
  )

  expect_false(isTRUE(ctl$store_traj))
  expect_true(isTRUE(ctl$R_hat_conv_check))
  expect_equal(ctl$stepsize_schedule, "constant")
  expect_equal(ctl$stepsize_decay, "grad_norm_plateau")
  expect_equal(ctl$sgd_method, "adam")
})
