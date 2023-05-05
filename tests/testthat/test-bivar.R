# Y1 Y2
# group = factor()
# f(~x+y, model="matern")
# coordinates()


test_that("test bv(ar1, ar1) with 1 noise", {
  load_all()
  n = 500
  true_model <- f(model="bv",
    first  = ar1(1:n, theta_K = ar1_a2th(0.8)),
    second = ar1(1:n, theta_K = ar1_a2th(0.2)),
    zeta = 1, rho = 0.5,
    noise = noise_nig(mu=-3, sigma=2, nu=1),
    eval=T
  )

  W <- simulate(true_model)
  AW <- as.numeric(true_model$A %*% W)
  n_obs <- length(AW)
  Y <- AW + rnorm(n_obs, sd=0.5)

  out <- ngme(
    Y ~
    # 0 + f(1:5, model="ar1", group="temp") +
    0 + f (model="bv", name="bv",
      zeta = 1,
      first=ar1(1:n),
      second=ar1(1:n),
      noise = noise_nig(),
      control = control_f(numer_grad = F)
    ),
    # noise = noise_nig(correlated = TRUE),
    data = data.frame(Y = Y
      # group=group, loc1=.., loc2=..
    ),
    control_ngme = control_ngme(
      n_gibbs_samples = 5
    ),
    control_opt = control_opt(
      iterations = 100,
      n_parallel_chain = 4,
      estimation = T,
      verbose = T
    )
  )
  out
  traceplot(out, "bv")
  out$replicates[[1]]$latents[[1]]$theta_K
  out$replicates[[1]]$latents[[1]]$operator$theta_K
})

test_that("test on bv(matern, matern)", {
  library(INLA)
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh <- inla.mesh.2d(loc.domain = pl01, cutoff = 0.2, max.edge = c(1,10))
  mesh$n

  n_obs1 <- 300
  loc1 <- cbind(runif(n_obs1, 0, 10), runif(n_obs1, 0, 5))
  n_obs2 <- 400
  loc2 <- cbind(runif(n_obs2, 0, 10), runif(n_obs2, 0, 5))

load_all()
  true_model <- f(model="bv",
    first  = matern(loc1, mesh=mesh, theta_K = log(1)),
    second = matern(loc2, mesh=mesh, theta_K = log(3)),
    # replicate = ...,
    zeta = 1, rho = 0.5,
    noise = noise_nig(mu=-3, sigma=2, nu=1),
    eval=T
  )

  W <- simulate(true_model)
  AW <- as.numeric(true_model$A %*% W)
  n_obs <- length(AW)
  Y <- AW + rnorm(n_obs, sd=0.5)
  length(Y)

  load_all()
  out <- ngme(
    Y ~ 0 +
    # f(x, model="rw1", group=1) +
    f(model="bv", name="bv",
      first =  matern(loc1, mesh=mesh),
      second = matern(loc2, mesh=mesh),
      zeta = 0.8,
      noise = noise_nig(),
      control = control_f(numer_grad = F)
    ),
    # noise = noise_nig(correlated = TRUE),
    data = data.frame(
      Y = Y
      # x =
    ),
    # group=c(111,222),
    control_ngme = control_ngme(
      n_gibbs_samples = 5
    ),
    control_opt = control_opt(
      iterations = 100,
      n_parallel_chain = 4,
      estimation = T,
      verbose = T
    )
  )
  out
  traceplot(out, "bv")
  out$replicates[[1]]$latents[[1]]$theta_K
  out$replicates[[1]]$latents[[1]]$operator$theta_K

})

test_that("Comparing different structures", {
  Y1 <- rnorm(3); Y2 <- rnorm(3)
  x1 <- rnorm(3); x2 <- rnorm(3)
  group_num = factor(rep(c(1,2), each=3))
  group_str = factor(rep(c("sal", "temp"), each=3))

  1 %in% levels(group_num)
  "sal" %in% levels(group_str)

# 1. sal, temp ~ f(x, "rw1", group = sal) + bv(ar1, ar1)
  load_all()
  m1 <- ngme(
    Y ~ f(x, model = "rw1", group = "sal") +
    f(model="bv", first=ar1(1:3), second=ar1(1:3)),
    data = data.frame(Y = c(Y1, Y2), x = c(x1, x2)),
    group = rep(c("sal", "temp"), each = 3),
    control_opt = control_opt(
      estimation = T,
      iterations = 1
    )
  )

  m1$replicates[[1]]$latents[[1]]$A
})
# f(x, "rw", group=1) on Y1
# f(x, "rw") on Y1 and Y2
# f(x, "rw", replicate=group) on Y1 and Y2 (share parameter)
# f(x, "rw", group=1) + f(x, "rw", group=2)
# on Y1 and Y2, indpendent process

# f(x, subset=c(T,T,F)) of length Y
# f(x, subset=c(1,2)) if group provide
