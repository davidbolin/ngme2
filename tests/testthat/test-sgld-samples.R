test_that("ngme_sgld_samples returns thinned samples as data.frame", {
  skip_on_cran()
  withr::local_seed(2)

  n <- 20
  dat <- data.frame(y = rnorm(n))

  fit <- ngme(
    y ~ 0 + f(1:n, model = ar1(), name = "ar1_field"),
    data = dat,
    family = noise_normal(),
    control_opt = control_opt(
      burnin = 2,
      iterations = 30,
      n_batch = 3,
      n_parallel_chain = 2,
      optimizer = sgld(stepsize = 0.001, temperature = 1),
      stepsize_control = poly_decay(alpha = 0.6, t0 = 10),
      pflug_conv_check = FALSE,
      R_hat_conv_check = FALSE,
      trend_std_conv_check = FALSE,
      store_traj = TRUE
    ),
    control_ngme = control_ngme(n_gibbs_samples = 1, n_post_samples = 2)
  )

  burnin_iter <- 5
  thinning <- 2
  tr <- get_trace_trajectories(fit, name = "general", apply_transform = FALSE)
  kept <- seq.int(burnin_iter + 1L, tr$n_iterations, by = thinning)

  samp <- ngme_sgld_samples(
    ngme = fit,
    name = "general",
    burnin_iter = burnin_iter,
    thinning = thinning,
    apply_transform = FALSE,
    combine_chains = TRUE
  )

  expect_true(is.data.frame(samp))
  expect_true(all(c(".chain", ".draw", ".iter") %in% colnames(samp)))
  expect_true(all(tr$parameter_names %in% colnames(samp)))
  expect_equal(nrow(samp), length(kept) * tr$n_chains)
  expect_equal(sort(unique(samp$.chain)), seq_len(tr$n_chains))

  samp_by_chain <- ngme_sgld_samples(
    ngme = fit,
    name = "general",
    burnin_iter = burnin_iter,
    thinning = thinning,
    apply_transform = FALSE,
    combine_chains = FALSE
  )
  expect_true(is.list(samp_by_chain))
  expect_equal(length(samp_by_chain), tr$n_chains)
  expect_true(all(vapply(samp_by_chain, nrow, integer(1)) == length(kept)))
})

test_that("ngme_sgld_samples validates burnin_iter", {
  skip_on_cran()
  withr::local_seed(3)

  n <- 12
  dat <- data.frame(y = rnorm(n))
  fit <- ngme(
    y ~ 0 + f(1:n, model = ar1(), name = "ar1_field"),
    data = dat,
    family = noise_normal(),
    control_opt = control_opt(
      burnin = 2,
      iterations = 12,
      n_batch = 2,
      n_parallel_chain = 1,
      optimizer = sgld(stepsize = 0.001, temperature = 1),
      stepsize_control = poly_decay(alpha = 0.6, t0 = 10),
      pflug_conv_check = FALSE,
      R_hat_conv_check = FALSE,
      trend_std_conv_check = FALSE,
      store_traj = TRUE
    ),
    control_ngme = control_ngme(n_gibbs_samples = 1, n_post_samples = 2)
  )

  tr <- get_trace_trajectories(fit, name = "general", apply_transform = FALSE)
  expect_error(
    ngme_sgld_samples(
      ngme = fit,
      name = "general",
      burnin_iter = tr$n_iterations,
      thinning = 1
    ),
    "burnin_iter must be smaller"
  )
})

test_that("ngme_sgld_ci returns quantile intervals and covariance", {
  skip_on_cran()
  withr::local_seed(4)

  n <- 20
  dat <- data.frame(y = rnorm(n))
  fit <- ngme(
    y ~ 0 + f(1:n, model = ar1(), name = "ar1_field"),
    data = dat,
    family = noise_normal(),
    control_opt = control_opt(
      burnin = 2,
      iterations = 30,
      n_batch = 3,
      n_parallel_chain = 2,
      optimizer = sgld(stepsize = 0.001, temperature = 1),
      stepsize_control = poly_decay(alpha = 0.6, t0 = 10),
      pflug_conv_check = FALSE,
      R_hat_conv_check = FALSE,
      trend_std_conv_check = FALSE,
      store_traj = TRUE
    ),
    control_ngme = control_ngme(n_gibbs_samples = 1, n_post_samples = 2)
  )

  samp <- ngme_sgld_samples(
    ngme = fit,
    name = "general",
    burnin_iter = 5,
    thinning = 2,
    apply_transform = FALSE
  )
  ci <- ngme_sgld_ci(
    samples = samp,
    lower = 0.1,
    upper = 0.9
  )
  ci_from_list <- ngme_sgld_ci(
    samples = ngme_sgld_samples(
      ngme = fit,
      name = "general",
      burnin_iter = 5,
      thinning = 2,
      apply_transform = FALSE,
      combine_chains = FALSE
    ),
    lower = 0.1,
    upper = 0.9
  )

  expect_s3_class(ci, "ngme_sgld_ci")
  expect_true(is.numeric(ci$estimates))
  expect_true(is.numeric(ci$std_error))
  expect_true(is.matrix(ci$ci))
  expect_true(is.matrix(ci$covariance))
  expect_true(all(ci$ci[, "lower"] <= ci$ci[, "upper"]))
  expect_equal(nrow(ci$covariance), length(ci$estimates))
  expect_equal(ncol(ci$covariance), length(ci$estimates))
  expect_true(is.data.frame(ci$samples))
  expect_equal(length(ci_from_list$estimates), length(ci$estimates))
})

test_that("ngme_sgld_ci validates quantile bounds", {
  expect_error(
    ngme_sgld_ci(
      samples = data.frame(x = rnorm(10)),
      lower = 0.95,
      upper = 0.05
    )
  )
})
