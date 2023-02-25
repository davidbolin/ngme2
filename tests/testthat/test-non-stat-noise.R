test_that("test non-stat mu", {
  n_obs <<- 500
  mu_cov <<- rexp(n_obs)
  ar <- model_ar1(1:n_obs, noise=noise_nig(
    B_mu = cbind(1, mu_cov),
    theta_mu = c(-4, -2),
    sigma = 5,
    nu = 3
  ))

  W <- simulate(ar)
  Y <- W + rnorm(n_obs, sd = 0.8)

  res <- ngme(
    Y ~ 0 + f(1:n_obs, name="ar", noise=noise_nig(
      theta_mu = c(0, 0),
      B_mu = cbind(1, mu_cov)
    )),
    family = "normal",
    data = list(Y=Y),
    control = ngme_control(
      iterations = 1000,
      print_check_info = FALSE
    )
  )
  res
  traceplot(res, "ar")

  # mu result close to theta_mu
  expect_true(abs(res$latents[[1]]$noise$theta_mu[1] - (-4)) < 5)
  expect_true(abs(res$latents[[1]]$noise$theta_mu[2] - (-2)) < 5)
})
