
# test the interface in R
test_that("R interface of replicate", {
  library(INLA)
  # load_all()
  # test
  y <- 11:17
  x <- 21:27
  repl <- c(1,1,1,2,2,2,2)
  loc = cbind(rnorm(7), rnorm(7))
  mesh2d = inla.mesh.2d(loc=loc, max.n=20)
####### case 1. ar + rw + matern2d
  out <- ngme(
    y ~ 0 + x + f(1:7, model="ar1", replicate=c(1,1,1,2,2,2,2)) +
      f(2:8, model="rw", order=1, replicate=repl) +
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
      f(2:8, model="rw", order=1,                 replicate=rep2) +
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
    data = data.frame(Y = Y, time=c(1:n_obs1, 1:n_obs2), x1=x1, x2=x2),
    control_opt = control_opt(
      estimation = T,
      iterations = 100,
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
    loc = list(list(c(1,2,3,4)), list(c(1,2,3,4))),
    data = list(cbind(1, rnorm(4), rexp(4)), cbind(1, rnorm(4), rexp(4)))
  )
  expect_true(TRUE)
})

# compare subsampling
test_that("subsampling vs sample all", {
  n_reps <- 10; n_each <- 50; n_obs <- n_reps * n_each
  myar1 <- ar1(1:n_each, alpha=-0.5, noise = noise_nig(mu = -3, sigma=2, nu=1))
  W <- double(); reps <- double()
  for (i in seq_len(n_reps)) {
    W    <- c(W, simulate(myar1))
    reps <- c(reps, rep(i, n_each))
  }
  Y <- W + rnorm(n_obs, sd=0.5)

  # load_all()
  out <- ngme(
    Y ~ 0 + f(1:n_obs, model="ar1", replicate=reps, noise=noise_nig()),
    family = "normal",
    data = data.frame(Y = Y),
    control_opt = control_opt(
      seed = 1,
      iterations = 1000,
      n_parallel_chain = 4,
      sampling_strategy = "ws"
    )
  )
  out
# is : time = 2
# all : time = 9

  traceplot(out, "field1")
})


# test_that("test_f(~1,...)", {
#   n <- 500
#   R1 <- rnorm(n, mean=3, sd=1)
#   R2 <- rnorm(n, mean=-1, sd=2)

#   n_obs <- 2*n
#   Y <- c(R1, R2) + rnorm(n_obs, sd=0.5)

#   fm <- Y ~ 1 + f(~1, family="normal")
#   load_all()
#   out = ngme(fm,
#     data=data.frame(Y=Y),
#     control_opt = control_opt(
#       iterations=10,
#       estimation = F
#     )
#   )
#   traceplot(out)
# })
