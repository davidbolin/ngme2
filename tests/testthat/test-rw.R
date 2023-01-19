# test rw related model

test_that("simulation of rw", {
# load_all()
  n_obs <<- 500
  mu <- -3; sigma <- 2.3; nu <- 2; sigma_eps <- 0.8
  # h <- rexp(n_obs)
  h <- rep(1, n_obs)
  loc <- c(0, cumsum(h))

  V <- rig(n_obs, a=nu, b=nu*h^2)
  dW <- -mu + mu * V + sigma * sqrt(V) * rnorm(n_obs) # type-G noise
  W <- c(0, cumsum(dW))

  Y <- W + rnorm(n=length(W), sd=sigma_eps)

  # check model specification
  K <- model_rw(loc, order=1)$K
  expect_true(all(as.numeric(K %*% W) - dW < 1e-5))

  # first we test the gradient of mu
  out <- ngme(
    Y ~ 0 + f(loc,
      model="rw1",
      name="rw1",
      noise=noise_nig(
        fix_nu = TRUE, nu = 2,
        fix_theta_sigma = TRUE, sigma = sigma,
        fix_V = TRUE, V = V
      ),
      fix_W = TRUE, W = W,
      debug = TRUE
    ),
    data = list(Y = Y),
    contro = ngme_control(
      estimation = T,
      iterations = 100,
      n_parallel_chain = 4
    ),
    debug = TRUE
  )
  out
# out$latents[[1]]$noise
  traceplot2(out, 1)

with(out, {
  expect_true(all(latents[[1]]$noise$V - V < 1e-5))
  expect_true(all(latents[[1]]$noise$h - h < 1e-5))
  expect_true(all(latents[[1]]$W - W < 1e-5))
})

  plot(out$latents[[1]]$noise,
    noise_nig(mu=mu, sigma=sigma, nu=nu))

expect_true(all(as.numeric(out$latents[[1]]$K %*% W) - dW < 1e-5))
expect_true(all(diff(loc) - out$latents[[1]]$h < 1e-5))
expect_true(out$latents[[1]]$noise$nu > 0.1)
})

test_that("test basic estimation of rw", {
  # test noise with h
  # load_all()
  n_obs <<- 500
  mu <- -3; sigma <- 2.3; nu <- 2; sigma_eps <- 0.8

  loc <- rnorm(n_obs, sd=1)
  rw1 <- model_rw(loc, noise=noise_nig(mu=mu, sigma=sigma, nu=nu))
  expect_true(mean(rw1$noise$h) != 1)

  W <- simulate(rw1)
  Y <- W + rnorm(n_obs, sd = sigma_eps)

  out <- ngme(
    Y ~ 0 + f(loc,
      model="rw1",
      name="rw1",
      noise=noise_nig(
        fix_V = TRUE,
        V = attr(W, "noise")$V
      ),
      # fix_W = TRUE, W = W,
      debug = TRUE
    ),
    data = list(Y = Y),
    contro = ngme_control(
      estimation = T,
      iterations = 100,
      n_parallel_chain = 1
    ),
    debug = TRUE
  )
  out
  out$latents[[1]]$noise$h

  traceplot2(out, "rw1")
  traceplot2(out, "mn")
  plot(attr(W, "noise"), out$latents[[1]]$noise)

})

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
  n_obs <<- 500
  alpha <- 0.75
  mu <- -3; sigma <- 2.3; nu <- 2; sigma_eps <- 0.8
  ar1 <- model_ar1(1:n_obs, alpha=alpha, noise=noise_nig(mu=mu, sigma=sigma, nu=nu))
  # ar1$K
# ar1$noise$h <- rexp(n_obs)

  W <- simulate(ar1)
  Y <- W + rnorm(n_obs, sd = sigma_eps)

  out <- ngme(
    Y ~ 0 + f(1:n_obs,
      model="ar1",
      name="ar",
      noise=noise_nig(
        # fix_nu = T, nu = 2,
        # h = ar1$noise$h,
        # fix_V = TRUE, V = attr(W, "noise")$V
      ),
      # fix_W = TRUE, W = W,
      debug = TRUE
    ),
    data = list(Y = Y),
    contro = ngme_control(
      estimation = T,
      iterations = 500,
      n_parallel_chain = 4
    ),
    debug = TRUE
  )
  out

  traceplot2(out, "ar")
  traceplot2(out, "mn")
  plot(attr(W, "noise"), out$latents[[1]]$noise)
  # out$latents[[1]]$noise$h

  with(out$latents[[1]], {
    expect_true(abs(noise$theta_mu - mu) < 1)
    expect_true(abs(noise$theta_sigma - log(sigma)) < 1)
    expect_true(abs(noise$nu - nu) < 2)
    expect_true(abs(ar1_th2a(out$latents[[1]]$theta_K) - alpha) < 0.1)
  })
})



# n_obs <- 1000
# rw_mu <- 4
# rw_sigma <- 1
# rw_nu <- 2.3
# rw1_process <- simulate(
#   f(1:n_obs,
#     model = "rw1",
#     noise = noise_gal(
#       mu = rw_mu,
#       sigma = rw_sigma,
#       nu = rw_nu
#     )
#   ),
#   seed = 1
# )
# Y <- ar1_process + rnorm(n_obs)
# attr(rw1_process, "noise")$h

# devtools::load_all()
# ngme_out <- ngme(
#   Y ~ 0 +
#   f(1:n_obs,
#     model = "ar1",
#     theta_K = 0.7,
#     # fix_theta_K = TRUE,
#     # W = as.numeric(ar1_process),
#     # fix_W = TRUE,
#     noise = noise_gal(
#       mu = 1.1,
#       sigma = 1,
#       # V = attr(ar1_process, "noise")$V,
#       # fix_V = TRUE
#     ),
#     control = ngme_control_f(
#       numer_grad       = F,
#       use_precond      = T
#     ),
#     debug = T
#   ),
#   data = data.frame(Y = Y),
#   family = "normal",
#   control = ngme_control(
#     estimation = T,
#     exchange_VW = TRUE,
#     n_parallel_chain = 4,
#     stop_points = 50,
#     burnin = 10,
#     iterations = 2000,
#     gibbs_sample = 5,
#     stepsize = 1,
#     threshold = 1e-4,

#     std_lim = 0.001,
#     trend_lim = 0.001
#   ),
#   seed = 10,
#   debug = TRUE
# )

# ngme_out
# traceplot(ngme_out, f_index = 1, param="alpha")
# traceplot(ngme_out, f_index = 1, param="mu")
# traceplot(ngme_out, f_index = 1, param="sigma")
# traceplot(ngme_out, f_index = 1, param="nu")