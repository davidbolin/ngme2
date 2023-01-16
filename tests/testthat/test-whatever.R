
test_that("test names", {
  expect_equal(f(1:10, name="ar")$name, "ar")

  bm <- ngme(
    Y ~ f(x=1:10, name="f1", model="ar1") +
      f(x=1:10, name="f2", model="ar1"),
    data = data.frame(Y = rnorm(10)),
    control = ngme_control(estimation = FALSE)
  )

  expect_equal(names(bm$latents), c("f1", "f2"))
})

# run manually
test_that("predict traceplot", {
  n_obs <<- 500
  x1 <- rexp(n_obs); x2 <- rnorm(n_obs)
  beta <- c(-2, 4, 1)
  alpha <- 0.75
  mu <- 1.5; sigma <- 2.3; nu <- 2; sigma_eps <- 0.8

  ar1 <- model_ar1(1:n_obs, alpha=alpha, noise=noise_nig(mu=mu, sigma=sigma, nu=nu))
  W <- simulate(ar1)
  Y <- beta[1] + x1 * beta[2] + x2 * beta[3] + W + rnorm(n_obs, sd = sigma_eps)

  # make 1/10 observation NA
  idx_NA <- sample(1:n_obs, size=n_obs/10)
  Y2 <- Y; Y2[idx_NA] <- NA

  # first we test estimation
  out <- ngme(
    Y2 ~ x1 + x2
      + f(1:n_obs, model="ar1", noise = noise_nig())
      + f(1:n_obs, model="ar1", noise = noise_normal(), name="f2"),
    family = noise_normal(),
    control = ngme_control(
      n_parallel_chain = 4,
      iteration = 20
   ),
   seed = 100,
   data = data.frame(Y2 = Y2, x1 = x1, x2 = x2)
  )
  out

  expect_no_error(traceplot2(out, name = "field1"))
  expect_no_error(traceplot2(out, name = "f2"))
})
