# test tensor product R structure
test_that("R structure", {
library(INLA)
  load_all()
  tp(first=ar1(1:3), second=ar1(1:3))

  out <- ngme(
    y ~ f(1:6, name="tp", model="tp", first=ar1(1:2), second=ar1(1:3)),
    data = data.frame(y=1:6),
    control_opt = control_opt(
      estimation = T,
      iterations = 1
    )
  )
out
  str(out$replicates[[1]]$latents[[1]]$operator)
  str(out$replicates[[1]]$latents[[1]]$operator$first)


  # f(i_space, model = spde, group = i_time, control.group = list(model="ar1"))
})

test_that("ar x 2d case", {
  set.seed(16)
  library(INLA)

##############################  simulation
  mesh2d <- inla.mesh.2d(
    loc.domain = cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5),
    max.edge = c(0.5, 10)
  )
  nr <- 300
  loc <- cbind(runif(nr, 0, 10), runif(nr, 0, 5))
  # mesh2d$n
  # plot(mesh2d)
  # points(loc[,1], loc[,2], pch=19, col="red")

  # build left model (ar)
  n_ar <- 3
  arr <- model_ar1(1:n_ar, alpha=0.7, debug=F)

  # build right model (matern)
  matern <- model_matern(map=loc, mesh=mesh2d, kappa = 2, debug=T)
  K <- arr$K %x% matern$K
  n <- nrow(K)

  eps <- simulate(noise_nig(mu=-2, sigma=1, nu=1, n = n))
  eps <- simulate(noise_normal(sigma=1, n = n))
  W <- solve(K, eps)
  A <- arr$A %x% matern$A
  dim(A)

  AW <- as.numeric(A %*% W)
  n_obs <- length(AW)
  Y <- AW + rnorm(n_obs, sd=0.5)

  # f(model=matern(mesh), group=ar1(1:3))
  tp <- f(model="tp", right=matern, left=arr, eval=T); tp
  expect_true(all(tp$A == A))
##############################  estimation
# str(f(model=matern, group=ar, noise=noise_nig())$model_right$noise)
  out <- ngme(
    Y ~ 0 + f(model="tp",
      left = ar1(1:n_ar, alpha=0.6),
      right = matern(loc, mesh=mesh2d, kappa=2.5),
      control = control_f(numer_grad = T)
    ),
    data = data.frame(Y=Y),
    family = "normal",
    control_opt = control_opt(
      iterations = 300,
      n_parallel_chain = 4,
      estimation = T,
      verbose = T,
      stepsize = 3
    ),
    debug = TRUE
  )
  out
  traceplot(out)
  traceplot(out, "field1")
})

######################################################################

test_that("iid x ar case", {
set.seed(16); library(INLA);
  n_obs <- 200
  Y1 <- simulate(ar1(1:n_obs, alpha=0.7, noise=noise_normal(n=n_obs)))
  Y2 <- simulate(ar1(1:n_obs, alpha=0.7, noise=noise_normal(n=n_obs)))
  Y3 <- simulate(ar1(1:n_obs, alpha=0.7, noise=noise_normal(n=n_obs)))
  Y <- c(Y1, Y2, Y3) + rnorm(n = 3*n_obs)
  out <- ngme(
  Y ~ 0 + f(model="tp",
    left=iid(1:3),
    right=ar1(1:n_obs, alpha=0.2),
    control = control_f(numer_grad = T),
    # fix_W = T, W = c(Y1, Y2),
    debug=F
  ),
  data = data.frame(Y=Y),
  control_opt = control_opt(
    estimation = T,
    iterations = 200,
    n_parallel_chain = 1,
    verbose = T
  ),
  debug= T
  )
  out; traceplot(out, "field1")
})

##########################################################################

test_that("ar x ar case", {
set.seed(16); library(INLA);

  nl <- 5
  Kl <- ar1(1:nl, alpha=0.3)$K
  Al <- ar1(1:nl, alpha=0.3)$A

  nr <- 300
  Kr <- ar1(1:nr, alpha=0.7)$K
  Ar <- ar1(1:nr, alpha=0.7)$A

  n <- nl * nr
  K <- Kl %x% Kr
  A <- Al %x% Ar

  eps <- simulate(noise_nig(
    mu = -2, sigma = 2, nu=1, n = n
  ))
  # eps <- simulate(noise_normal(sigma=2, n=n))
  W <- solve(K, eps)
  Y <- as.numeric(A %*% W) + rnorm(n, sd=0.5)

tp <- f(model="tp", right=ar1(1:nr, alpha=0.7), left=ar1(1:nl, alpha=0.3), eval=T)
  expect_true(all(tp$A == A))

  ##############################  estimation
  # str(f(model=matern, group=ar, noise=noise_nig())$model_right$noise)

  out <- ngme(
    Y ~ 0 + f(model="tp",
      left=ar1(1:nl, alpha=0.9), # true 0.3
      right=ar1(1:nr, alpha=0.2), # true 0.7
      control_f=control_f(numer_grad = T),
      noise=noise_nig()
    ),
    data = data.frame(Y=Y),
    family = "normal",
    control_opt = control_opt(
      iterations = 500,
      n_parallel_chain = 4,
      estimation = T,
      verbose = T
    ),
    debug = F
  )
  out$replicates[[1]]

  traceplot(out, "field1")
})

