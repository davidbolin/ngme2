# test rw related model

test_that("simulate and estimate of rw with NIG", {
  n_obs <<- 500
  mu <- -3; sigma <- 5; nu <- 2; sigma_eps <- 0.8
  # h <- rpois(n_obs, lambda=3) + 1
  # h <- rep(1, n_obs)
  h <- rexp(n_obs, rate = 2)
  loc <<- c(0, cumsum(h))

  # V <- rig(n_obs, a=nu, b=nu*h^2, seed = 3)
  # dW <- -mu*h + mu * V + sigma * sqrt(V) * rnorm(n_obs) # type-G noise
  # W <- c(0, cumsum(dW))

  # plot(W)
  my_rw <- model_rw(loc, order=1, noise=noise_nig(mu=-3, sigma=5, nu=2))
  W <- simulate(my_rw, seed = 3)

# compare W and V
# plot(V, type="l")
# lines(attr(W2, "noise")$V, col="red")
  # points(W2, col="red")

  # check model specification
  Y <- W + rnorm(n=length(W), sd=sigma_eps)
  expect_true(all(my_rw$K == my_rw$C + my_rw$G))
  # expect_true(all(as.numeric(my_rw$K %*% W) - dW < 1e-5))

# ???
# expect_true(all(my_rw$noise$h - h < 1e-5))
# all(c(diff(loc[-1]), mean(diff(loc[-1]))) - my_rw$noise$h < 1e-5)

  # first we test the gradient of mu
  out <- ngme(
    Y ~ 0 + f(loc,
      model="rw",
      name="rw",
      noise=noise_nig(
        # fix_nu = TRUE, nu = 2,
        # fix_theta_sigma = TRUE, sigma = sigma,
        # fix_V = TRUE, V = V
      ),
      # fix_W = TRUE, W = W,
      debug = FALSE
    ),
    data = data.frame(Y = Y),
    control_opt = control_opt(
      estimation = T,
      iterations = 1,
      n_parallel_chain = 4,
      print_check_info = TRUE,
      verbose = F
    ),
    debug = TRUE
  )
  out
  traceplot(out, "rw")
  plot(out$replicates[[1]]$models[[1]]$noise,
    noise_nig(mu=mu, sigma=sigma, nu=nu))

expect_true(all(as.numeric(out$replicates[[1]]$models[[1]]$K %*% W) - dW < 1e-5))
})

######################################################################
test_that("the order of W same as order of index?", {
  library(INLA)
  XX    <- c(1.1, 3.1, 2.2, 2.2, 4.5, 5)
  index_NA <- c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
  myloc <- c(3.2, 1.2)

  mesh1 <- INLA::inla.mesh.1d(XX); mesh1

  rw <- model_rw(XX, order=1, index_NA = index_NA)
  rw$C + rw$G
  rw$A; rw$A_pred

  # use INLA make A? the order is bad
  A <- INLA::inla.spde.make.A(mesh=mesh1, loc = myloc); A

# 2 x 4 sparse Matrix of class "dgCMatrix"
# [1,] .         .         0.7142857 0.2857143
# [2,] 0.4545455 0.5454545 .         .

  # the order is wrong, recover by permutation matrix
  # A %*% as(mesh1$idx$loc, "pMatrix")

# 3 x 4 sparse Matrix of class "dgCMatrix"
# [1,] 1 -1  .  .
# [2,] .  1 -1  .
# [3,] .  .  1 -1

  expect_true(TRUE)
})

############################## AR1 case
test_that("test estimation of basic ar with normal measurement noise", {
  n_obs <- 600
  alpha <- 0.75
  mu <- 4; sigma <- 2; nu <- 1; sigma_eps <- 0.8
  ar1 <- f(1:n_obs, model="ar1", theta_K = ar1_a2th(0.8), noise=noise_nig(mu=mu, sigma=sigma, nu=nu), eval = T)

  W <- simulate(ar1)
  mean(W)
  x1 = rnorm(n_obs);
  # Y <-  -3 + x1 * 2  + rnorm(n_obs, sd = sigma_eps)
  Y <- W + rnorm(n_obs, sd = sigma_eps) - 3 + x1 * 2

# plot(Y, type="l")
load_all()
  out <- ngme(
    Y ~  1 + x1
    + f(1:n_obs,
      model="ar1",
      name="ar",
      noise=noise_nig(
        # fix_nu = T, nu = 2,
        # h = ar1$noise$h,
        # fix_V = TRUE, V = attr(W, "noise")$V
      ),
      control=control_f(numer_grad = T),
      debug = FALSE,
      # fix_W = T, W = W
    ),
    data = data.frame(Y = Y, x1=x1),
    control_ngme = control_ngme(
      # n_gibbs_samples = 5
    ),
    control_opt = control_opt(
      estimation = T,
      iterations = 5000,
      n_parallel_chain = 4,
      print_check_info = FALSE,
      verbose = F,
    ),
    debug = F
  )
  traceplot(out, "ar")
  traceplot(out)
plot(attr(W, "noise"), out$replicates[[1]]$models[[1]]$noise)

  out$replicates[[1]]$models[[1]]$theta_K
  # out$replicates[[1]]$models[[1]]$control$numer_grad
  # load_all()

  # plot(simulate(out$replicates[[1]]$models[["ar"]]), type="l")
  # plot(Y, type="l")
  # prds <- predict(out$replicates[[1]], loc=list(ar=501:600))$mean

  with(out$replicates[[1]]$models[[1]], {
    expect_true(abs(noise$theta_mu - mu) < 4)
    expect_true(abs(noise$theta_sigma - log(sigma)) < 2)
    expect_true(abs(noise$nu - nu) < 4)
    expect_true(abs(ar1_th2a(out$replicates[[1]]$models[[1]]$theta_K) - alpha) < 0.1)
  })
})

# very tricky!!!
test_that("test rw definition", {
  # load_all()
  m1 <- model_rw(rexp(5), order = 1, noise = noise_normal())

  # KW = noise => nrow(K) = length(h)
  expect_equal(nrow(m1$K), length(m1$h))
  # AW = Y
  expect_equal(m1$mesh$n, 5)

  m2 <- model_rw(rnorm(10), order = 2, circular = TRUE)
  expect_equal(nrow(m2$K), length(m2$h))
})


######################################################################
test_that("test ou process", {
  # load_all()
  # ?Diagonal
  # model_ou(1:5)$K
# M1 <- cbind(c(1,2), c(2,3))
# M1 %*% diag(c(2,3))

  simulate(model_ou(1:5))
  n_obs <- 5
  Y <- rnorm(n_obs)
  B_theta_K <- cbind(1, 1:5)
  B_theta_K
  out <- ngme(
    Y ~ 0 + f(Y, model="ou", theta=1, theta_K = c(0.5, 0.5), B_theta_K = B_theta_K),
    data = data.frame(Y = Y),
    control_opt = control_opt(
      estimation = T,
      iterations = n_obs,
      n_parallel_chain = 4,
      print_check_info = TRUE,
      verbose = T
    ),
  )
  out$replicates[[1]]$models[[1]]$theta_K
  # traceplot(out, "field")

  model_ou(c(1,3,5), order=1)$h

  m1 <- model_ou(c(1,3,5), order=2)
  simulate(m1)
  expect_true(TRUE)
})