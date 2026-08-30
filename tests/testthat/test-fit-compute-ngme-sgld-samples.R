test_that("compute_ngme_sgld_samples refits and returns samples", {
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
      R_hat_conv_check = FALSE,
      trend_std_conv_check = FALSE,
      store_traj = FALSE
    ),
    control_ngme = control_ngme(n_gibbs_samples = 1, n_post_samples = 2)
  )

  samp <- compute_ngme_sgld_samples(
    fit = fit0,
    iterations = 120,
    optimizer = sgld(stepsize = 0.001, temperature = 1),
    burnin = 2,
    n_batch = 4,
    n_parallel_chain = 1,
    alpha = 0.6,
    t0 = 10,
    seed = 123,
    verbose = FALSE,
    name = "general",
    burnin_iter = 7,
    thinning = 2,
    apply_transform = FALSE
  )

  expect_true(is.data.frame(samp))
  expect_true(all(c(".chain", ".draw", ".iter") %in% colnames(samp)))
  expect_true(nrow(samp) > 0)
  expect_true(inherits(attr(samp, "refit"), "ngme"))
  expect_true(inherits(attr(samp, "refit_source"), "ngme"))
  expect_true(inherits(attr(samp, "refit_control_opt"), "control_opt"))
  expect_equal(attr(samp, "burnin_iter"), 7)
  expect_equal(attr(samp, "thinning"), 2)
})

test_that("compute_ngme_sgld_samples disallows separate schedule_burnin_iter in dots", {
  skip_on_cran()
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
      R_hat_conv_check = FALSE,
      trend_std_conv_check = FALSE,
      store_traj = FALSE
    ),
    control_ngme = control_ngme(n_gibbs_samples = 1, n_post_samples = 2)
  )

  expect_error(
    compute_ngme_sgld_samples(
      fit = fit0,
      iterations = 40,
      burnin_iter = 5,
      schedule_burnin_iter = 3
    ),
    "unified with `burnin_iter`"
  )
})

test_that("compute_ngme_sgld_samples can override n_gibbs_samples for refit", {
  skip_on_cran()
  withr::local_seed(7)

  n <- 12
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
      R_hat_conv_check = FALSE,
      trend_std_conv_check = FALSE,
      store_traj = FALSE
    ),
    control_ngme = control_ngme(n_gibbs_samples = 1, n_post_samples = 2)
  )

  samp_inherit <- compute_ngme_sgld_samples(
    fit = fit0,
    iterations = 40,
    optimizer = sgld(stepsize = 0.001, temperature = 1),
    burnin = 2,
    n_batch = 2,
    n_parallel_chain = 1,
    seed = 321,
    verbose = FALSE,
    name = "general",
    burnin_iter = 5,
    thinning = 2,
    apply_transform = FALSE
  )

  expect_equal(
    attr(samp_inherit, "refit")$replicates[[1]]$control_ngme$n_gibbs_samples,
    1
  )

  samp_override <- compute_ngme_sgld_samples(
    fit = fit0,
    iterations = 40,
    optimizer = sgld(stepsize = 0.001, temperature = 1),
    burnin = 2,
    n_batch = 2,
    n_parallel_chain = 1,
    seed = 321,
    verbose = FALSE,
    name = "general",
    burnin_iter = 5,
    thinning = 2,
    n_gibbs_samples = 3,
    apply_transform = FALSE
  )

  expect_equal(
    attr(samp_override, "refit")$replicates[[1]]$control_ngme$n_gibbs_samples,
    3
  )
})
