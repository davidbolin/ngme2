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
  n <- 100
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
      iterations = 5,
      stop_points = 5,
      n_parallel_chain = 4,
      print_check_info = F
    )
  )

  out
  traceplot(out,"field1")
  out$replicates[[1]]$models[[1]]$theta_K
  out$replicates[[1]]$models[[1]]$operator$theta_K

  predict(out, 
    map = list(field1=c(1,2,3)),
    group = c("B", "B", "A")
  )
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
  predict(out, 
    map=list(field1=cbind(c(1,2,3), c(2,3,4))),
    group = rep("A", 6)
  )
  # out$replicates[[1]]$models[[1]]$theta_K
  # out$replicates[[1]]$models[[1]]$operator$theta_K
})

test_that("test on bv matern NIG", {
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh <- fmesher::fm_mesh_2d(
    loc.domain = pl01,
    cutoff = 0.2,
    max.edge = c(2, 10)
  )
  mesh$n
  n_obs <- 100

  # generate random locations (same for two fields)
  long <- runif(n_obs/2, 0, 10); lat <- runif(n_obs/2, 0, 5)
  long <- c(long, long); lat <- c(lat, lat)
  group <- c(rep("W1", n_obs/2), rep("W2", n_obs/2))

  # plot the mesh
  # plot(mesh); points(long, lat)
 
  # parameters of the bivariate model
  theta = pi/4; rho = 0.2
  theta_K_1 = log(2); theta_K_2 = log(5);
  mu_1 = -2; sigma_1 = 0.5; nu_1 = 1
  mu_2 = 2; sigma_2 = 0.3; nu_2 = 0.5

  true_model <- f(
    ~ long + lat,
    theta = theta,
    rho = rho,
    mesh = mesh,
    sd1=sigma_1, sd2=sigma_2,
    model = "bv_matern_nig",
    sub_models = list(
      W1 = list(model = "matern", theta_K = theta_K_1),
      W2 = list(model = "matern", theta_K = theta_K_2)
    ),
    group = group,
    noise = list(
      W1 = noise_nig(),
      W2 = noise_nig()
    )
  )
  true_model

  sim_fields <- simulate(true_model)[[1]]

  sd_1 = 0.6; sd_2 = 0.9; rho_e = 0.9
  Cov_same_idx <- matrix(c(sd_1^2, rho_e*sd_1*sd_2, rho_e*sd_1*sd_2, sd_2^2), nrow=2)

  Cov_measurement <- Cov_same_idx %x% diag(n_obs / 2)

  # e ~ N(0, Cov_measurement)
  L <- t(chol(Cov_measurement))
  e <- L %*% rnorm(n_obs)

  # fixed effects
  x1 <- rexp(n_obs)
  x2 <- rnorm(n_obs)
  feff <- c(-3, 1.5)

  Y <- sim_fields + x1 * feff[1] + x2 * feff[2] + as.numeric(e)

  B_sigma = matrix(0, nrow=n_obs, ncol=2)
  B_sigma[group=="W1", 1] = 1
  B_sigma[group=="W2", 2] = 1

# 
  out_cor_gauss <- ngme(
    Y ~ 0 + x1 + x2 + f(
      ~ long + lat,
      mesh = mesh,
      model = "bv_matern_normal",
      name = "bv",
      sub_models = list(
        W1 = list(model = "matern"),
        W2 = list(model = "matern")
      ),
      # debug=T,
      noise = list(
        W1 = noise_normal(),
        W2 = noise_normal()
      )
    ),
    group = group,
    family = noise_normal(
      corr_measurement = TRUE,
      index_corr = c(1:(n_obs/2), 1:(n_obs/2)),
      B_sigma = B_sigma,
      theta_sigma = c(0, 0)
    ),
    data = data.frame(Y, long, lat),
    control_opt = control_opt(
      iterations = 20,
      n_parallel_chain = 4,
      seed = 50
    ),
    debug = FALSE
  )
  out_cor_gauss

  out_cor_nig <- ngme(
    Y ~ 0 + x1 + x2 + f(
      ~ long + lat,
      mesh = mesh,
      model = "bv_matern_nig",
      name = "bv",
      sub_models = list(
        W1 = list(model = "matern"),
        W2 = list(model = "matern")
      ),
      # debug=T,
      noise = list(
        W1 = noise_nig(),
        W2 = noise_nig()
      )
    ),
    group = group,
    family = noise_normal(
      corr_measurement = TRUE,
      index_corr = c(1:(n_obs/2), 1:(n_obs/2)),
      B_sigma = B_sigma,
      theta_sigma = c(0, 0)
    ),
    data = data.frame(Y, long, lat),
    control_opt = control_opt(
      estimation = TRUE,
      iterations = 20,
      print_check_info = FALSE,
      seed = 50
    ),
    debug = FALSE,
    start = out_cor_gauss
  )
  out_cor_nig
})
