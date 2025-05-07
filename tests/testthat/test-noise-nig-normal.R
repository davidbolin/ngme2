test_that("test-noise-nig-normal", {
  load_all()
  n_obs <- 10
  Y <- rnorm(n_obs)
  fit <- ngme(
    Y ~ 0 + f(
      1:n_obs,
      name = "my_ar",
      model = "ar1",
      noise = noise_normal_nig(),
    ),
    data = data.frame(Y=Y),
    control_opt = control_opt(
      iterations = 100,
      n_parallel_chain = 1,
      verbose = FALSE,
      seed = 12
    )
  )

})