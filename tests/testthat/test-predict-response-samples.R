test_that("predict supports response draws and output-level samples", {
  withr::local_seed(42)

  n <- 30
  w <- as.numeric(arima.sim(list(ar = 0.5), n = n))
  y <- w + rnorm(n, sd = 0.3)

  fit <- ngme(
    y ~ 0 + f(1:n, model = ar1(), name = "ar1_field"),
    data = data.frame(y = y),
    family = noise_normal(),
    control_opt = control_opt(
      iterations = 30,
      burnin = 1,
      n_batch = 3,
      n_parallel_chain = 2,
      n_min_batch = 1,
      verbose = FALSE,
      print_check_info = FALSE,
      store_traj = FALSE,
      seed = 42
    )
  )

  pred <- predict(
    fit,
    map = list(ar1_field = seq_len(n)),
    type = "response",
    estimator = c("mean", "0.25q", "0.75q"),
    sampling_size = 25,
    burnin_size = 5,
    seed = 7,
    return_samples = TRUE
  )

  response_samples <- attr(pred, "samples")
  latent_samples <- attr(pred, "latent_samples")

  expect_true(is.matrix(response_samples))
  expect_equal(dim(response_samples), c(n, 25))
  expect_true(is.matrix(latent_samples))
  expect_equal(ncol(latent_samples), 25)

  expect_equal(pred$mean, rowMeans(response_samples), tolerance = 1e-8)
  expect_equal(
    pred$`0.25q`,
    apply(response_samples, 1, quantile, 0.25),
    tolerance = 1e-8
  )
  expect_equal(
    pred$`0.75q`,
    apply(response_samples, 1, quantile, 0.75),
    tolerance = 1e-8
  )
})


test_that("predict only requires maps for requested latent fields", {
  withr::local_seed(123)

  n <- 20
  w1 <- as.numeric(arima.sim(list(ar = 0.4), n = n))
  w2 <- cumsum(rnorm(n, sd = 0.15))
  y <- w1 + w2 + rnorm(n, sd = 0.2)

  fit <- ngme(
    y ~ 0 +
      f(1:n, model = ar1(), name = "ar1_field") +
      f(1:n, model = rw1(), name = "rw1_field"),
    data = data.frame(y = y),
    family = noise_normal(),
    control_opt = control_opt(
      iterations = 20,
      burnin = 1,
      n_batch = 2,
      n_parallel_chain = 2,
      n_min_batch = 1,
      verbose = FALSE,
      print_check_info = FALSE,
      store_traj = FALSE,
      seed = 123
    )
  )

  pred_mean <- predict(
    fit,
    map = list(ar1_field = seq_len(n)),
    type = "ar1_field",
    estimator = "mean",
    sampling_size = 12,
    burnin_size = 3,
    seed = 11
  )

  pred_samples <- predict(
    fit,
    map = list(ar1_field = seq_len(n)),
    type = "ar1_field",
    estimator = "mean",
    sampling_size = 12,
    burnin_size = 3,
    seed = 11,
    return_samples = TRUE
  )

  expect_length(pred_mean$mean, n)
  expect_length(pred_samples$mean, n)
  expect_equal(dim(attr(pred_samples, "samples")), c(n, 12))

  expect_error(
    predict(
      fit,
      map = list(ar1_field = seq_len(n)),
      type = "lp",
      estimator = "mean",
      sampling_size = 12,
      burnin_size = 3,
      seed = 11
    ),
    "rw1_field"
  )
})
