# test effiency
test_ngme <- function(
  model,
  n_obs_per_rep,
  n_replicate,
  n_iter,
  n_parallel_chain = 4,
  stop_points = 20,
  numer_grad = FALSE,
  preconditioner = "fast",
  sampling_strategy = "all",
  precond_by_diff_chain = TRUE,
  compute_precond_each_iter = FALSE,
  max.n = 1000,
  print = FALSE
) {
  # create 2d mesh
  if (model %in% c("matern", "bvmatern")) {
    pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
    mesh <- fmesher::fm_mesh_2d(
      loc.domain = pl01, cutoff = 0.2,
      max.edge = c(0.5, 10),
      max.n = max.n
    )
print(paste("nodes of mesh = ", mesh$n))
  }

  # ------- Simulate data for each model --------
  sim_data <- switch(model,
    "ar1" = {
      idx <- 1:n_obs_per_rep
      ar1_model <- f(idx, model="ar1", rho = 0.5,
      noise = noise_nig(mu = -3, sigma = 2, nu=0.4))
      W <- simulate(ar1_model, seed = 16)
      Y <- W + rnorm(n_obs_per_rep, sd = 2)
      list(Y=Y, idx=idx, group=rep(1, n_obs_per_rep))
    },
    "matern" = {
      loc <- cbind(runif(n_obs_per_rep, 0, 10), runif(n_obs_per_rep, 0, 5))
      true_model <- f(
        map = loc,
        model="matern",
        theta_K = log(2),
        mesh = mesh,
        noise = noise_nig(mu=-2, sigma=1, nu=0.5)
      )
      W <- simulate(true_model)
      Y <- as.numeric(true_model$A %*% W) + rnorm(n_obs_per_rep, sd=0.5)
      list(Y=Y, idx=loc, group=rep(1, n_obs_per_rep))
    },
    "bvar1" = {
      group_per_rep <- c(rep("first", n_obs_per_rep/2), rep("second", n_obs_per_rep/2))
      idx_per_rep <- c(1:(n_obs_per_rep/2), 1:(n_obs_per_rep/2))
      true_model <- f(
        idx_per_rep,
        model="bv",
        theta = 0.5, rho = 0.8,
        sub_models = list(
          first = list(model="ar1", rho=0.5),
          second= list(model="ar1", rho=-0.5)
        ),
        group = group_per_rep,
        noise = list(
          first=noise_nig(mu=-3, sigma=2, nu=1),
          second=noise_nig(mu=-3, sigma=2, nu=1)
        )
      )
      W <- simulate(true_model)
      AW <- as.numeric(true_model$A %*% W)
      Y <- AW + rnorm(length(AW), sd=0.5)
      list(Y=Y, idx=idx_per_rep, group=group_per_rep)
    },
    "bvmatern" = {
      loc <- cbind(runif(n_obs_per_rep, 0, 10), runif(n_obs_per_rep, 0, 5))
      group_per_rep <- c(rep("first", n_obs_per_rep/2), rep("second", n_obs_per_rep/2))
      idx_per_rep <- c(1:(n_obs_per_rep/2), 1:(n_obs_per_rep/2))
      true_model <- f(
        idx_per_rep,
        model="bv",
        theta = 0.5, rho = 0.8,
        sub_models = list(
          first = list(model="matern", kappa=0.5),
          second= list(model="matern", kappa=0.5)
        ),
        group = group_per_rep,
        noise = list(
          first=noise_nig(mu=-3, sigma=2, nu=1),
          second=noise_nig(mu=-3, sigma=2, nu=1)
        )
      )
      W <- simulate(true_model)
      AW <- as.numeric(true_model$A %*% W)
      Y <- AW + rnorm(length(AW), sd=0.5)
      list(Y=Y, idx=idx_per_rep, group=group_per_rep)
    },
    stop("Unknown test model")
  )

  # ------- Specify formula for each model -------
  formula <- switch(model,
    "ar1" = Y ~ f(idx,
      model = "ar1",
      noise = noise_nig(),
      control = control_f(numer_grad = numer_grad)
    ),
    "matern" = Y ~ 0 + f(idx,
      model="matern",
      mesh = mesh,
      noise=noise_nig(),
      control = control_f(numer_grad = numer_grad)
    ),
    "bvar1" = Y ~ f(idx,
      model="bv",
      sub_models = list(first = "ar1", second="ar1"),
      control = control_f(numer_grad = numer_grad),
      noise = list(first=noise_nig(), second=noise_nig())
    ),
    "bvmatern" = Y ~ f(idx,
      model="bv",
      sub_models = list(first = "matern", second="matern"),
      control = control_f(numer_grad = numer_grad),
      noise = list(first=noise_nig(), second=noise_nig())
    ),
  )

  # make replicate
  idx <- rep_map(sim_data$idx, n_replicate)
  group <- rep(sim_data$group, n_replicate)
  Y <- rep(sim_data$Y, n_replicate)
  repl <- rep(1:n_replicate, each=n_obs_per_rep)

  # fit
  start <- proc.time()
  out <- ngme(
    formula,
    replicate = repl,
    group = group,
    data = data.frame(Y=Y),
    control_opt = control_opt(
      burnin = 100,
      std_lim = 0.001,
      print_check_info = FALSE,
      verbose = FALSE,
      seed = 3,
      iterations = n_iter,
      precond_by_diff_chain = precond_by_diff_chain,
      compute_precond_each_iter = compute_precond_each_iter,
      n_parallel_chain = n_parallel_chain,
      stop_points = stop_points,
      preconditioner = preconditioner,
      sampling_strategy = sampling_strategy
    )
  )

  print(proc.time() - start)
  list(
    out = out,
    time = proc.time() - start
  )
}
