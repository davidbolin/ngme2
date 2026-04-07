test_that("predict and cross_validation support predictive_average over chains", {
  withr::local_seed(1)

  n <- 40
  y <- as.numeric(arima.sim(list(ar = 0.6), n = n)) + rnorm(n, sd = 0.2)

  fit <- ngme(
    y ~ 0 + f(1:n, model = ar1(), name = "ar1_field"),
    data = data.frame(y = y),
    control_opt = control_opt(
      iterations = 30,
      n_batch = 3,
      n_parallel_chain = 3,
      burnin = 1,
      n_min_batch = 1,
      verbose = FALSE,
      print_check_info = FALSE,
      # store_traj = FALSE,
      seed = 1
    )
  )

  chain_outputs <- attr(fit, "chain_outputs")
  expect_true(is.list(chain_outputs))
  expect_equal(length(chain_outputs), 3)

  pred_param_mean <- predict(
    fit,
    map = list(ar1_field = seq_len(n)),
    estimator = "mean",
    sampling_size = 20,
    burnin_size = 5,
    seed = 1,
    chain_combine = "param_mean"
  )
  pred_predictive_avg <- predict(
    fit,
    map = list(ar1_field = seq_len(n)),
    estimator = "mean",
    sampling_size = 20,
    burnin_size = 5,
    seed = 1,
    chain_combine = "predictive_average"
  )

  expect_equal(length(pred_predictive_avg$mean), n)
  expect_true(all(is.finite(pred_predictive_avg$mean)))
  expect_true(ncol(attr(pred_predictive_avg, "samples")) >= ncol(attr(pred_param_mean, "samples")))

  cv_predictive_avg <- cross_validation(
    fit,
    type = "k-fold",
    k = 2,
    N_sim = 1,
    n_gibbs_samples = 5,
    n_burnin = 2,
    seed = 1,
    print = FALSE,
    parallel = FALSE,
    chain_combine = "predictive_average"
  )

  expect_true(all(is.finite(as.matrix(cv_predictive_avg$mean.scores))))
})
