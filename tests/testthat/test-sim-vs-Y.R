test_that("test cor model single replicate", {
  library(ngme2)
  library(mvtnorm)
  sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2)
  mu = c(0, 0)
  set.seed(123)
  n = 50
  W = 1:(2*n)
  noise = rmvnorm(n = n, mean = mu, sigma = sigma)

  # collapse Y into a vector
  Y = W + as.vector(noise)
  Y

  index_corr = rep(1:n, 2)
  index_corr

  # fit cor model
  mod = ngme(
    Y ~ f(c(1:(2*n)), model="ar1", noise=noise_nig()),
    data = data.frame(Y = Y),
    family = noise_normal(
      corr_measurement = TRUE,
      index_corr = index_corr
    )
  )

  plot(Y)
  sim1 = simulate(mod)[[1]]
  lines(sim1, col = "red")
  expect_true(cor(sim1, Y) > 0.98)
})


test_that("test cor model single replicate", {
  library(ngme2)
  library(mvtnorm)
  sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2)
  mu = c(0, 0)
  set.seed(123)
  n = 50
  W = 1:(2*n)
  noise = rmvnorm(n = n, mean = mu, sigma = sigma)
  
  repl = c(rep(1, n/2), rep(2, n/2))
  repl = c(repl, repl)
  plot(repl)

  # collapse Y into a vector
  Y = W + as.vector(noise)
  Y

  index_corr = rep(1:n, 2)
  index_corr

  # fit cor model
  mod = ngme(
    Y ~ f(c(1:(2*n)), model="ar1", noise=noise_nig()),
    data = data.frame(Y = Y),
    replicate = repl,
    family = noise_normal(
      corr_measurement = TRUE,
      index_corr = index_corr
    )
  )
  sim1 = simulate(mod)[[1]]
  plot(sim1)
  lines(Y, col = "red")

  expect_true(cor(sim1, Y) > 0.98)
})



test_that("test independent model two replicates", {
  library(ngme2)
  set.seed(123)
  n = 50

  repl = c(rep(1, n/2), rep(2, n/2))
  repl = c(repl, repl)
  plot(repl)

  # collapse Y into a vector
  Y = W + rnorm(2*n)
  Y

  # fit cor model
  mod = ngme(
    Y ~ f(c(1:(2*n)), model="ar1", noise=noise_nig()),
    data = data.frame(Y = Y),
    replicate = repl,
    family = noise_normal()
  )
  sim1 = simulate(mod)[[1]]
  plot(sim1)
  lines(Y, col = "red")
  expect_true(cor(sim1, Y) > 0.98)
})