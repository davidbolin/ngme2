test_that("calibrate_inv_exp_lambda_driven_nig returns coherent calibration", {
  out <- calibrate_inv_exp_lambda_driven_nig(
    r_target = 1.3,
    alpha = 0.1,
    c = 2.5,
    mu = 0,
    sigma = 1,
    h = 1,
    n_samples = 20000,
    nu_lower = 0.05,
    nu_upper = 50,
    tol = 0.25,
    max_iter = 12,
    seed = 123
  )

  expect_true(is.list(out))
  expect_true(inherits(out, "ngme_inv_exp_calibration"))
  expect_gt(out$lambda, 0)
  expect_gt(out$nu_r, 0)
  expect_equal(out$lambda, -out$nu_r * log(out$alpha), tolerance = 1e-12)
  expect_equal(out$alpha_implied, out$alpha, tolerance = 1e-12)
  expect_lt(abs(out$rc_nu_r - out$r_target), 0.6)
})

test_that("calibrate_inv_exp_lambda is a direct alias", {
  out1 <- calibrate_inv_exp_lambda_driven_nig(
    r_target = 1.3,
    alpha = 0.2,
    c = 2.5,
    n_samples = 20000,
    nu_lower = 0.05,
    nu_upper = 50,
    tol = 0.25,
    max_iter = 12,
    seed = 99
  )

  out2 <- calibrate_inv_exp_lambda(
    r_target = 1.3,
    alpha = 0.2,
    c = 2.5,
    n_samples = 20000,
    nu_lower = 0.05,
    nu_upper = 50,
    tol = 0.25,
    max_iter = 12,
    seed = 99
  )

  expect_equal(out2$lambda, out1$lambda, tolerance = 1e-12)
  expect_equal(out2$nu_r, out1$nu_r, tolerance = 1e-12)
})

test_that("calibration reports unattainable tail target clearly", {
  expect_error(
    calibrate_inv_exp_lambda_driven_nig(
      r_target = 2,
      alpha = 0.1,
      c = 2,
      mu = 0,
      sigma = 1,
      h = 2,
      n_samples = 20000,
      seed = 1
    ),
    "Target may be unattainable"
  )
})
