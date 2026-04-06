test_that("compute_ngme_ci refits from existing object and returns CI", {
  skip_on_cran()
  withr::local_seed(42)

  n <- 20
  dat <- data.frame(y = rnorm(n))

  fit0 <- ngme(
    y ~ 0 + f(1:n, model = ar1(), name = "ar1_field"),
    data = dat,
    family = noise_normal(),
    control_opt = control_opt(
      burnin = 2,
      iterations = 20,
      n_batch = 2,
      n_parallel_chain = 1,
      optimizer = adam(stepsize = 0.01),
      pflug_conv_check = FALSE,
      R_hat_conv_check = FALSE,
      trend_std_conv_check = FALSE,
      store_traj = FALSE
    ),
    control_ngme = control_ngme(n_gibbs_samples = 1, n_post_samples = 2)
  )

  ci <- compute_ngme_ci(
    fit = fit0,
    iterations = 120,
    alpha = 0.6,
    optimizer = sgd(stepsize = 0.01),
    burnin = 2,
    n_batch = 4,
    n_parallel_chain = 1,
    t0 = 1,
    seed = 123,
    verbose = FALSE,
    burnin_iter = 7,
    M = 5
  )

  expect_s3_class(ci, "ngme_batch_ci")
  expect_true(inherits(ci$refit, "ngme"))
  expect_true(inherits(ci$refit_source, "ngme"))
  expect_true(inherits(ci$refit_control_opt, "control_opt"))
  expect_equal(ci$alpha, 0.6)
  expect_equal(ci$burnin_iter, 7)
  expect_equal(ci$refit_control_opt$stepsize_schedule_burnin_iter, 7)
  expect_true(is.matrix(ci$conf_int))
  expect_true(all(c("lower", "upper") %in% colnames(ci$conf_int)))
})


test_that("compute_ngme_CI is an alias of compute_ngme_ci", {
  expect_identical(compute_ngme_CI, compute_ngme_ci)
})


test_that("compute_ngme_ci disallows separate schedule_burnin_iter in dots", {
  withr::local_seed(42)
  n <- 10
  dat <- data.frame(y = rnorm(n))
  fit0 <- ngme(
    y ~ 0 + f(1:n, model = ar1(), name = "ar1_field"),
    data = dat,
    family = noise_normal(),
    control_opt = control_opt(
      burnin = 2,
      iterations = 20,
      n_batch = 2,
      n_parallel_chain = 1,
      optimizer = adam(stepsize = 0.01),
      pflug_conv_check = FALSE,
      R_hat_conv_check = FALSE,
      trend_std_conv_check = FALSE,
      store_traj = FALSE
    ),
    control_ngme = control_ngme(n_gibbs_samples = 1, n_post_samples = 2)
  )

  expect_error(
    compute_ngme_ci(
      fit = fit0,
      iterations = 40,
      alpha = 0.6,
      burnin_iter = 5,
      schedule_burnin_iter = 3
    ),
    "unified with `burnin_iter`"
  )
})
