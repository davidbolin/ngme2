test_that("test tp and bv operator", {
  # 1. modify bv
  n = 10
  true_model <- f(
    1:n,
    model="bv",
    theta = 0.5, rho = 0.5,
    sub_models = list(
      first = "ar1", second="ar1"
    ),
    group = c(rep("first", n/2), rep("second", n/2)),
    noise = list(
      first=noise_nig(mu=-3, sigma=2, nu=1),
      second=noise_nig(mu=-3, sigma=2, nu=1)
    )
  )
  true_model
})

test_that("test bv(ar1, ar1) with 2 noise", {
  n <- 500
  true_model <- f(
    1:n,
    model="bv",
    theta = 0, rho = 0.6,
    sub_models = list(
      first = "ar1", second="ar1"
    ),
    group = c(rep("first", n/2), rep("second", n/2)),
    noise = list(
      first=noise_nig(mu=-3, sigma=2, nu=1),
      second=noise_nig(mu=-3, sigma=2, nu=1)
    )
  )

  W <- simulate(true_model)
  AW <- as.numeric(true_model$A %*% W)
  n_obs <- length(AW)
  Y <- AW + rnorm(n_obs, sd=0.5)

  out <- ngme(
    Y ~ 0 + f(
      1:n,
      model="bv",
      sub_models = list(
        first = "ar1", second="ar1"
      ),
      control = control_f(numer_grad = F),
      noise = list(first=noise_nig(), second=noise_nig())
    ),
    group = c(rep("first", n/2), rep("second", n/2)),
    data = data.frame(Y = Y),
    control_ngme = control_ngme(
      n_gibbs_samples = 5
    ),
    control_opt = control_opt(
      iterations = 100,
      n_parallel_chain = 4,
      estimation = T,
      verbose = T,
      print_check_info = F,
      preconditioner = "none"
    )
    # ,start = out
  )
  out
  traceplot(out,"field1")
  out$replicates[[1]]$models[[1]]$theta_K
  out$replicates[[1]]$models[[1]]$operator$theta_K
})

test_that("test on bv(matern, matern)", {
  # library(INLA)
  # pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  # mesh <- inla.mesh.2d(loc.domain = pl01, cutoff = 0.2, max.edge = c(1,10))
  # mesh$n

  # n_obs1 <- 300
  # loc1 <- cbind(runif(n_obs1, 0, 10), runif(n_obs1, 0, 5))
  # n_obs2 <- 400
  # loc2 <- cbind(runif(n_obs2, 0, 10), runif(n_obs2, 0, 5))

  # true_model <- f(model="bv",
  #   first  = matern(loc1, mesh=mesh, theta_K = log(1)),
  #   second = matern(loc2, mesh=mesh, theta_K = log(3)),
  #   # replicate = ...,
  #   zeta = 1, rho = 0.5,
  #   noise = noise_nig(mu=-3, sigma=2, nu=1),
  #   eval=T
  # )

  # W <- simulate(true_model)
  # AW <- as.numeric(true_model$A %*% W)
  # n_obs <- length(AW)
  # Y <- AW + rnorm(n_obs, sd=0.5)
  # length(Y)

  # out <- ngme(
  #   Y ~ 0 +
  #   # f(x, model="rw1", group=1) +
  #   f(model="bv", name="bv",
  #     first =  matern(loc1, mesh=mesh),
  #     second = matern(loc2, mesh=mesh),
  #     zeta = 0.8,
  #     noise = noise_nig(),
  #     control = control_f(numer_grad = F)
  #   ),
  #   # noise = noise_nig(correlated = TRUE),
  #   data = data.frame(
  #     Y = Y
  #     # x =
  #   ),
  #   # group=c(111,222),
  #   control_ngme = control_ngme(
  #     n_gibbs_samples = 5
  #   ),
  #   control_opt = control_opt(
  #     iterations = 100,
  #     n_parallel_chain = 4,
  #     estimation = T,
  #     verbose = T
  #   )
  # )
  # out
  # traceplot(out, "bv")
  # out$replicates[[1]]$models[[1]]$theta_K
  # out$replicates[[1]]$models[[1]]$operator$theta_K
})
