# Setup code - accessible to all test_that blocks

Y <- 1:10; n_obs <- length(Y)
x1 <- runif(n_obs)
x2 <- rexp(n_obs)

ngme_out <- ngme(
  Y ~ x1 + x2 + f(
    1:n_obs,
    rho = 0.5,
    name = "my_ar",
    model = "ar1",
    noise = noise_nig(
      sigma=3, 
      nu=2,
      theta_mu=c(1,2),
      B_mu=matrix(c(1,0,0,1), 2, 2)
    )
  ) + f(
    1:n_obs,
    name = "my_ou",
    model = "ou",
    noise = noise_normal(theta_sigma=c(1,2), B_sigma=matrix(c(1,0,0,1), 2, 2))
  ),
  data = data.frame(x1=x1, x2=x2, Y=Y),
  control_opt = control_opt(estimation = FALSE)
)

withr::defer(rm(ngme_out))

test_that("get_noise_param works", {
  noise <- noise_nig(mu=2, sigma=3, nu=2)
  expect_equal(get_noise_param(noise, "mu"), list(mu = 2))
  expect_equal(get_noise_param(noise, "sigma"), list(sigma = 3))
  expect_equal(get_noise_param(noise, "nu"), list(nu = 2))

  noise2 <- noise_nig(
    theta_sigma=c(1,2), 
    B_sigma=matrix(c(1,0,0,1), 2, 2)
  )
  expect_equal(
    get_noise_param(noise2, "sigma"), 
    list(theta_sigma = c(1,2))
  )
})


test_that("extract_parameters (transformed) works", {
  # 1st latent model
  expect_equal(
    extract_parameters(ngme_out)[["transformed"]][["my_ar"]],
    list(
      rho = 0.5,
      theta_mu = c(1,2),
      sigma = 3,
      nu = 2
    )
  )

  # 2nd latent model
  expect_equal(
    extract_parameters(ngme_out)[["transformed"]][["my_ou"]],
    list(
      theta_K1 = 0,
      theta_sigma = c(1,2)
    )
  )

  # fixed effects and measurement noise
  fe <- extract_parameters(ngme_out)[["transformed"]][["data"]][["fixed_effects"]]
  expect_true(all(c("(Intercept)", "x1", "x2") %in% names(fe)))
})

test_that("extract_parameters (raw) works", {
  # 1st latent model
  expect_equal(
    extract_parameters(ngme_out)[["raw"]][["my_ar"]],
    list(
      theta_rho = 1.09861229,
      theta_mu = c(1,2),
      sigma = 3,
      nu = 2
    )
  )

  # 2nd latent model
  expect_equal(
    extract_parameters(ngme_out)[["raw"]][["my_ou"]],
    list(
      theta_K1 = 0,
      theta_sigma = c(1,2)
    )
  )

  # fixed effects and measurement noise
  fe <- extract_parameters(ngme_out)[["raw"]][["data"]][["fixed_effects"]]
  expect_true(all(c("(Intercept)", "x1", "x2") %in% names(fe)))
})


test_that("ngme_result works", {
  expect_equal(
    ngme_result(ngme_out, model = "my_ar", transformed = TRUE),
    list(rho = 0.5, theta_mu = c(1,2), sigma = 3, nu = 2)
  )

  expect_equal(
    ngme_result(ngme_out, model = "my_ou", transformed = TRUE),
    list(theta_K1 = 0, theta_sigma = c(1,2))
  )

  expect_equal(
    ngme_result(ngme_out, model = "my_ar", transformed = FALSE),
    list(theta_rho = 1.09861229, theta_mu = c(1,2), sigma = 3, nu = 2)
  )
})
