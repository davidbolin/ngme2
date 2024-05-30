test_that("test tp and bv operator", {
  # 1. modify bv
  n = 10
  true_model <- f(
    c(1:(n/2), 1:(n/2)),
    model = "bv",
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
  n <- 1000
  true_model <- f(
    c(1:(n/2), 1:(n/2)),
    model="bv",
    theta = 0.5, rho = 0.8,
    sub_models = list(
      A = list(model="ar1", rho=0.5),
      B = list(model="ar1", rho=-0.5)
    ),
    group = c(rep("A", n/2), rep("B", n/2)),
    noise = list(
      A=noise_nig(mu=3, sigma=2, nu=1),
      B=noise_nig(mu=-3, sigma=2, nu=1)
    )
  )
true_model$operator$K
  W <- simulate(true_model, seed=1)[[1]]
  n_obs <- length(W)
  Y <- W + rnorm(n_obs, sd=0.5)
  # Y = c(Y_A, Y_B)
  # label = c(A, B)

  out <- ngme(
    Y ~ 0 + f(
      c(1:(n/2), 1:(n/2)),
      model="bv",
      # rho = 0.1,
      theta = 0.45,
      sub_models = list(
        A = "ar1",
        B = "ar1"
      ),
      control = control_f(numer_grad = F),
      noise = list(A=noise_nig(), B=noise_nig())
    ),
    group = c(rep("A", n/2), rep("B", n/2)),
    data = data.frame(Y = Y),
    control_ngme = control_ngme(
      n_gibbs_samples = 5
    ),
    control_opt = control_opt(
      optimizer = precond_sgd(
        preconditioner = "full",
        precond_by_diff_chain = TRUE,
        numerical_eps = 1e-5
      ),
      seed = 3,
      iterations = 500,
      n_parallel_chain = 4,
      print_check_info = F
    )
  )
  out
  traceplot(out,"field1")
  out$replicates[[1]]$models[[1]]$theta_K
  out$replicates[[1]]$models[[1]]$operator$theta_K
  predict(out, map=list(field1=c(1,2,3)))
})

test_that("test on bv(matern, matern)", {
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh <- fmesher::fm_mesh_2d(loc.domain = pl01, cutoff = 0.2, max.edge = c(1,10))
  mesh$n

  n_obs1 <- 300; loc1 <- cbind(runif(n_obs1, 0, 10), runif(n_obs1, 0, 5))
  n_obs2 <- 400; loc2 <- cbind(runif(n_obs2, 0, 10), runif(n_obs2, 0, 5))

  true_model <- f(
    map = rbind(loc1, loc2),
    model="bv",
    sub_models = list(
      A = list(model="matern", theta_K = log(1)),
      B = list(model="matern", theta_K = log(3))
    ),
    group = c(rep("A", 300), rep("B", 400)),
    mesh = mesh,
    noise = list(A=noise_nig(), B=noise_nig())
  )

  W <- simulate(true_model)[[1]]
  n_obs <- length(W)
  Y <- W + rnorm(n_obs, sd=0.5)
  length(Y)

  out <- ngme(
    Y ~ 0 + f(map = rbind(loc1, loc2),
      model="bv",
      sub_models = list(
        A = list(model="matern"),
        B = list(model="matern")
      ),
      mesh = mesh,
      noise = list(A=noise_nig(), B=noise_nig())
    ),
    group = c(rep("A", 300), rep("B", 400)),
    data = data.frame(
      Y = Y
    ),
    control_ngme = control_ngme(
      n_gibbs_samples = 5
    ),
    control_opt = control_opt(
      iterations = 10,
      n_parallel_chain = 4,
      estimation = T,
      verbose = T
    )
  )

  out
  traceplot(out, "field1")
  predict(out, map=list(field1=cbind(c(1,2,3), c(2,3,4))))
  # out$replicates[[1]]$models[[1]]$theta_K
  # out$replicates[[1]]$models[[1]]$operator$theta_K
})
