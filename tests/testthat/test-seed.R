test_that("seed is deterministic", {
  seed <- 10; n_obs <- 10; Y <- rnorm(n_obs)
  control <- control_opt(
    seed = seed,
    iterations = 5,
    n_parallel_chain = 1,
    stop_points = 1
  )
  fit1 <- ngme(
    Y ~ 0 + f(
      1:n_obs,
      name = "my_ar",
      model = "ar1",
      rho = 0.5,
    ),
    data = data.frame(Y=Y),
    control_opt = control
  )
  fit2 <- ngme(
    Y ~ 0 + f(
      1:n_obs,
      name = "my_ar",
      model = "ar1",
      rho = 0.5,
    ),
    data = data.frame(Y=Y),
    control_opt = control
  )
  expect_equal(
    ngme_result(fit1, "my_ar")$theta_K, 
    ngme_result(fit2, "my_ar")$theta_K
  )
})
