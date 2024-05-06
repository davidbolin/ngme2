test_that("test non-stat mu", {
  n_obs <- 500
  mu_cov <- rexp(n_obs)
  ar <- f(1:n_obs, rho=0.5,model="ar1", noise=noise_nig(
    B_mu = cbind(1, mu_cov),
    theta_mu = c(-4, -2),
    sigma = 5,
    nu = 3
  ))

  W <- simulate(ar)[[1]]
  Y <- W + rnorm(n_obs, sd = 0.8)

  res <- ngme(
    Y ~ 0 + f(1:n_obs, model="ar1", name="ar", noise=noise_nig(
      theta_mu = c(0, 0),
      B_mu = cbind(1, mu_cov)
    )),
    family = "normal",
    data = data.frame(Y=Y),
    control_opt = control_opt(
      iterations = 20,
      print_check_info = FALSE
    )
  )
  res
  traceplot(res, "ar")
})
