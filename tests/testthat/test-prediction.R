# This file is for testing the prediction function.

test_that("test predict general", {
  set.seed(10)
  n_obs <- 500
  x1 <- rexp(n_obs); x2 <- rnorm(n_obs)
  feff <- c(-2, 4, 1)
  alpha <- 0.75
  mu <- 1.5; sigma <- 2.3; nu <- 2; sigma_eps <- 0.8

  ar1 <- f(1:n_obs, model="ar1", rho=alpha, noise=noise_nig(mu=mu, sigma=sigma, nu=nu))
  W <- simulate(ar1)[[1]]
  Y <- feff[1] + x1 * feff[2] + x2 * feff[3] + W + rnorm(n_obs, sd = sigma_eps)

  # first we test estimation
  out <- ngme(
    Y ~ x1 + x2 + f(1:n_obs, model="ar1", noise = noise_nig(
    )),
    family = noise_normal(),
    control_opt = control_opt(
      n_parallel_chain = 4,
      estimation = T,
      iteration = 1000,
      print_check_info = FALSE
   ),
   data = data.frame(Y = Y, x1 = x1, x2 = x2)
  )
  out
  traceplot(out)
  traceplot(out, "field1")
  plot(out$replicates[[1]]$models[[1]]$noise, attr(W, "noise"))

  # compare results
  expect_true(sum(abs(out$replicates[[1]]$feff - feff)) < 2)
  expect_true(abs(out$replicates[[1]]$noise$theta_sigma - log(sigma_eps)) < 1)
  expect_true(abs(out$replicates[[1]]$models[[1]]$noise$theta_mu - mu) < 1)
  expect_true(abs(out$replicates[[1]]$models[[1]]$noise$theta_sigma - log(sigma)) < 1)
  expect_true(abs(out$replicates[[1]]$models[[1]]$noise$theta_nu - log(nu)) < 5)
})


test_that("test posterior sampling and model_validation()", {
  n_obs <- 20
  alpha <- 0.3; mu = -3; sigma=2; nu=2; sigma_eps <- 0.5
  my_ar <- f(1:n_obs, model="ar1", noise=noise_nig(mu=mu, sigma=sigma, nu=nu), eval=TRUE)
  W <- simulate(my_ar)[[1]]
  Y <- W + rnorm(n=length(W), sd=sigma_eps)

  out <- ngme(
    Y ~ 1 + f(1:n_obs, name="ar", model="ar1") + f(1:n_obs, name="rw", model="rw1") + f(1:n_obs, model="matern", mesh=fmesher::fm_mesh_1d(1:10)),
    data=data.frame(Y=Y),
    control_opt=control_opt(print_check_info = FALSE, iteration=10)
  )
  out

  # 1. estimator
  predict(out, map = list(ar = c(2,3,5), rw=c(5,2,1),field1 = c(2,3,4) ))

  # 2. compute idx
  cross_validation(out)

  out$replicates[[1]]$Y
  expect_no_error(samples <- sampling_cpp(out$replicates[[1]], 10, posterior = TRUE, seed=10))
})


test_that("test lpo CV", {
  n_obs <- 100
  ar_mod <- f(1:n_obs, model="ar1", noise=noise_nig())
  yy <- simulate(ar_mod)[[1]] + rnorm(n_obs, sd=0.5)
  ng_100 <- ngme(
    yy~0+f(1:n_obs, model="ar1", noise=noise_nig()),
    data = data.frame(yy=yy),
    control_opt = control_opt(iterations = 100)
  )
  cross_validation(ng_100, type="lpo", times=10)

  ng_1000 <- ngme(
    yy~0+f(1:n_obs, model="ar1", noise=noise_nig()),
    data = data.frame(yy=yy),
    control_opt = control_opt(iterations = 100)
  )
  cross_validation(ng_1000, type="lpo", times=10)

  expect_true(TRUE)
})

test_that("test prediction with repls", {
  # load_all()
  # out <- test_ngme("ar1", n_obs_per_rep = 50, n_replicate = 5)
  # out
  # predict.ngme(out$out, map=list(field1=c(1,2,3)))

  # out2 <- ngme(
  #   formula,
  #   replicate = repl,
  #   group = group,
  #   data = data.frame(Y=Y),
  #   control_ngme = control_ngme(
  #     n_gibbs_samples = n_gibbs_samples
  #   ),
  #   control_opt = control_opt(
  #     estimation = F
  #   ),
  #   start = out,
  #   family = fm_mn_noise,
  #   debug = debug
  # )
  # predict.ngme(out, map=list(field1=c(1,2,3)))
  # predict.ngme(out2, map=list(field1=c(1,2,3)))
})