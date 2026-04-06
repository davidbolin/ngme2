test_that("stepsize_schedule validates inputs", {
  sch <- stepsize_schedule(method = "constant")
  expect_s3_class(sch, "ngme_stepsize_schedule")
  expect_equal(sch$method, "constant")

  sch_poly <- stepsize_schedule(method = "poly", alpha = 0.6, t0 = 2, burnin_iter = 11)
  expect_equal(sch_poly$method, "poly")
  expect_equal(sch_poly$alpha, 0.6)
  expect_equal(sch_poly$t0, 2)
  expect_equal(sch_poly$burnin_iter, 11)

  expect_error(stepsize_schedule(method = "poly", alpha = 0.4, t0 = 1))
  expect_error(stepsize_schedule(method = "poly", alpha = 1.0, t0 = 1))
  expect_error(stepsize_schedule(method = "poly", alpha = 0.6, t0 = -1))
  expect_error(stepsize_schedule(method = "poly", alpha = 0.6, t0 = 1, burnin_iter = -1))
  expect_error(stepsize_schedule(method = "poly", alpha = 0.6, t0 = 1, burnin_iter = 2.5))
})

test_that("stepsize_control helpers produce expected components", {
  sc <- stepsize_control(
    schedule = stepsize_schedule(method = "poly", alpha = 0.6, t0 = 2),
    decay = stepsize_decay(method = "grad_norm_plateau", patience = 4, gamma = 0.3)
  )
  expect_s3_class(sc, "ngme_stepsize_control")
  expect_equal(sc$schedule$method, "poly")
  expect_equal(sc$schedule$alpha, 0.6)
  expect_equal(sc$schedule$t0, 2)
  expect_equal(sc$decay$method, "grad_norm_plateau")
  expect_equal(sc$decay$patience, 4)
  expect_equal(sc$decay$gamma, 0.3)

  pd <- poly_decay(alpha = 0.61, t0 = 9, burnin_iter = 13)
  expect_equal(pd$schedule$method, "poly")
  expect_equal(pd$schedule$burnin_iter, 13)
  expect_equal(pd$decay$method, "none")

  bd <- batch_decay(patience = 7, gamma = 0.4)
  expect_equal(bd$schedule$method, "constant")
  expect_equal(bd$decay$method, "grad_norm_plateau")
  expect_equal(bd$decay$patience, 7)
  expect_equal(bd$decay$gamma, 0.4)
})


test_that("control_opt accepts stepsize schedule settings", {
  ctl <- control_opt(
    stepsize_control = poly_decay(alpha = 0.6, t0 = 5, burnin_iter = 9)
  )
  expect_equal(ctl$stepsize_schedule, "poly")
  expect_equal(ctl$stepsize_schedule_alpha, 0.6)
  expect_equal(ctl$stepsize_schedule_t0, 5)
  expect_equal(ctl$stepsize_schedule_burnin_iter, 9)

  sch <- stepsize_schedule(method = "poly", alpha = 0.7, t0 = 3)
  ctl2 <- control_opt(stepsize_control = stepsize_control(schedule = sch))
  expect_equal(ctl2$stepsize_schedule, "poly")
  expect_equal(ctl2$stepsize_schedule_alpha, 0.7)
  expect_equal(ctl2$stepsize_schedule_t0, 3)

  ctl3 <- control_opt(
    stepsize_control = stepsize_control(
      schedule = stepsize_schedule(method = "poly", alpha = 0.8, t0 = 10)
    )
  )
  expect_equal(ctl3$stepsize_schedule_alpha, 0.8)
  expect_equal(ctl3$stepsize_schedule_t0, 10)

  expect_error(poly_decay(alpha = 0.4, t0 = 1))
  expect_error(poly_decay(alpha = 1.0, t0 = 1))
  expect_error(poly_decay(alpha = 0.6, t0 = -1))
  expect_error(poly_decay(alpha = 0.6, t0 = 1, burnin_iter = -2))
})


test_that("ngme runs with polynomial stepsize schedule", {
  skip_on_cran()

  set.seed(1)
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
      optimizer = sgd(stepsize = 0.01),
      stepsize_control = poly_decay(alpha = 0.6, t0 = 1),
      pflug_conv_check = FALSE,
      R_hat_conv_check = FALSE,
      trend_std_conv_check = FALSE,
      store_traj = FALSE
    ),
    control_ngme = control_ngme(n_gibbs_samples = 1, n_post_samples = 2)
  )

  expect_s3_class(fit, "ngme")
})
