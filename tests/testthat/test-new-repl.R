
# test the interface in R
test_that("R interface of replicate", {
  library(INLA)
  load_all()
  # test
  y <- 11:17
  x <- 21:27
  repl <- c(1,1,1,2,2,2,2)
  loc = cbind(rnorm(7), rnorm(7))
  mesh2d = inla.mesh.2d(loc=loc, max.n=20)
####### case 1. ar + rw + matern2d
  load_all()
  out <- ngme(
    y ~ 0 + x + f(1:7, model="ar1", replicate=c(1,1,1,2,2,2,2)) +
      f(2:8, model="rw1", replicate=repl) +
      f(loc, model="matern", mesh=mesh2d),
    data = data.frame(
      y=y, x=x, repl=repl, loc = loc
    ),
    control_opt = control_opt(
      estimation = T,
      iterations = 10
    )
  )
  out

####### case 2. ar + rw + matern2d with f_replicate
rep1 <- c(1,1,2,2,3,3,3)
rep2 <- c(1,1,1,1,2,2,2)
rep3 <- c(1,1,3,4,2,2,2)
merge_repls(list(rep1, rep2, rep3))
  out2 <- ngme(
    y ~ 0 + x + f(1:7, model="ar1",       replicate=rep1) +
      f(2:8, model="rw1",                 replicate=rep2) +
      f(loc, model="matern", mesh=mesh2d, replicate=rep3),
    data = data.frame(
      y=y, x=x, repl=repl, loc = loc
    ),
    control_opt = control_opt(
      estimation = T,
      iterations = 10
    )
  )
  out2
  expect_true(TRUE)
})

test_that("basic ar1 case with different length", {
  # load_all()
  n_obs1  <- 500
  n_obs2 <- 100
  alpha <- 0.75
  mu <- -3; sigma <- 2.3; nu <- 2; sigma_eps <- 0.8
  ar_1 <- model_ar1(1:n_obs1, alpha=alpha, noise=noise_nig(mu=mu, sigma=sigma, nu=nu))
  ar_2 <- model_ar1(1:n_obs2, alpha=alpha, noise=noise_nig(mu=mu, sigma=sigma, nu=nu))

  W1 <- simulate(ar_1)
  W2 <- simulate(ar_2)
  Y <- c(W1, W2) + rnorm(n_obs1 + n_obs2, sd = sigma_eps)
  beta <- c(3,1,-1)
  x1 <- rnorm(n_obs1 + n_obs2)
  x2 <- rexp(n_obs1 + n_obs2)
  Y <- beta[1] + beta[2] * x1 + beta[3] * x2 + Y
  # Y <- beta[1] + Y

  load_all()
  out <- ngme(
    Y ~ x1 + x2 + f(time,
      model="ar1",
      name="ar",
      alpha = -0.5,
      replicate = c(rep(1, n_obs1), rep(2, n_obs2)),
      noise=noise_nig(),
      control=control_f(
        numer_grad = T
      ),
      debug = FALSE
    ),
    data = list(Y = Y, time=c(1:n_obs1, 1:n_obs2), x1=x1, x2=x2),
    control_opt = control_opt(
      estimation = T,
      iterations = 500,
      n_parallel_chain = 2,
      print_check_info = FALSE,
      verbose = F,
      exchange_VW = T
    ),
    debug = FALSE
  )
  out
  cross_validation(out)
  traceplot(out)
  traceplot(out, "ar")
  plot(simulate(out[[1]]$latents[["ar"]]), type="l")
  plot(Y, type="l")
  plot(attr(W1, "noise"), out[[1]]$latents[[1]]$noise)

  # test on predict function
  predict(out,
    loc = list(c(1,2,3,4)),
    data = list(cbind(1, rnorm(4), rexp(4)))
  )
  expect_true(TRUE)
})
