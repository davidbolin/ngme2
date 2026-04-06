library(ngme2)
library(testthat)

test_that("control_opt robust parameter works", {
  # Create a simple model
  n_obs <- 100
  sigma_eps <- 0.5
  x <- runif(n_obs)
  beta <- 2
  y <- beta * x + rnorm(n_obs, sd = sigma_eps)

  # Fit model with robust = TRUE
  control <- control_opt(
    iterations = 10,
    robust = TRUE,
    estimation = TRUE,
    verbose = TRUE
  )

  # Just checking if it runs without error
  expect_error(
    ngme(
      formula = y ~ x,
      data = data.frame(y = y, x = x),
      control_opt = control,
      family = "normal"
    ),
    NA
  )
})
