
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


# temp tests
# KW is sim_noise, and fix V gives perfect estimation of mu
test_that("delete that", {
# load_all()
  n_obs <<- 100
  alpha <- 0.3; mu = -3; sigma=2; nu=2; sigma_eps <- 0.5
  my_ar <- model_ar1(1:n_obs, alpha=alpha, noise=noise_nig(mu=mu, sigma=sigma, nu=nu))

  W <- simulate(my_ar)
  Y <- W + rnorm(n=length(W), sd=sigma_eps)

  expect_true(all(my_ar$K == my_ar$C * alpha + my_ar$G))
  # expect_true(all(as.numeric(my_ar$K %*% W) - dW < 1e-5))

  # first we test the gradient of mu
  out <- ngme(
    Y ~ 0 + f(1:n_obs,
      model="ar1",
      fix_theta_K = TRUE,
      alpha=0.3,
      noise=noise_nig(
        fix_nu = TRUE, nu = nu,
        fix_theta_sigma = TRUE, sigma = sigma,
        fix_V = TRUE, V = attr(W, "noise")$V
      ),
      # fix_W = TRUE, W = W,
      debug = TRUE
    ),
    data = list(Y = Y),
    contro = ngme_control(
      estimation = T,
      iterations = 1000,
      n_parallel_chain = 1
    ),
    debug = TRUE
  )
  out
  traceplot2(out, 1)
# out$latents[[1]]$noise
  plot(out$latents[[1]]$noise,
    noise_nig(mu=mu, sigma=sigma, nu=nu))

with(out, {
  expect_true(all(latents[[1]]$noise$V - V < 1e-5))
  expect_true(all(latents[[1]]$noise$h - h < 1e-5))
  expect_true(all(latents[[1]]$W - W < 1e-5))
})

expect_true(all(as.numeric(out$latents[[1]]$K %*% W) - dW < 1e-5))
expect_true(all(diff(loc) - out$latents[[1]]$h < 1e-5))
expect_true(out$latents[[1]]$noise$nu > 0.1)
})

# tests on posterior sampling
test_that("test posterior sampling and model_validation()", {
# load_all()
  n_obs <<- 20
  alpha <- 0.3; mu = -3; sigma=2; nu=2; sigma_eps <- 0.5
  my_ar <- model_ar1(1:n_obs, alpha=alpha, noise=noise_nig(mu=mu, sigma=sigma, nu=nu))
  W <- simulate(my_ar)
  Y <- W + rnorm(n=length(W), sd=sigma_eps)

  out <- ngme(Y ~ 0 + f(model=my_ar), data=list(Y=Y))
  # traceplot2(out, "field1")

  expect_no_error(samples <- sampling_cpp(out, 10, posterior = TRUE))
  # samples[["AW"]]

  model_validation(out, N=100)
})
