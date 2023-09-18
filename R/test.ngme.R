#' @title Test ngme function
#'
#' @description
#' Test ngme function for different models
#' @param model model name
#' @param n_obs_per_rep number of observation per replicate
#' @param n_replicate number of replicate
#' @param n_iter number of iteration
#' @param n_parallel_chain number of parallel chains
#' @param stop_points     number of stop points for convergence check
#' @param numer_grad numerical gradient
#' @param preconditioner  preconditioner, can be c("none", "fast", "full")
#' "none" means no preconditioner, "fast" means precondition everything except for the parameter of K matrix (for speed reason), "full" means precondition everything
#' @param sampling_strategy subsampling method of replicates of model, c("all", "is")
#' @param precond_by_diff_chain logical, if TRUE, use different chains to estimate preconditioner (only computed at check points), if FALSE, use the same chain to estimate preconditioner (computed at each iteration)
#' @param compute_precond_each_iter logical, if TRUE, compute preconditioner at each iteration, if FALSE, only compute preconditioner at check points (if has only 1 chain running, it will be set TRUE)
#' @param max.n maximum number for building mesh
#' @param print print the process
#' @param debug debug mode
#' @param debug_f debug mode for latent process
#' @export
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
  num_threads = c(n_parallel_chain, 4),
  n_gibbs_samples = 5,
  f_noise = "nig",
  family = "nig",
  max.n = 1000,
  print = FALSE,
  debug = FALSE,
  debug_f = FALSE,
  verbose = FALSE,
  rao_blackwellization = FALSE,
  start = NULL,
  seed = Sys.time()
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

f_simu_noise <- if (f_noise=="nig") noise_nig(mu = -3, sigma = 2, nu=0.4)
  else noise_normal()
f_form_noise <- if (f_noise=="nig") noise_nig() else noise_normal()

  # ------- Simulate data for each model --------
  sim_data <- switch(model,
    "ar1" = {
      idx <- 1:n_obs_per_rep
      ar1_model <- f(idx, model="ar1", rho = 0.5,
      noise = f_simu_noise,
      debug = debug_f
      )
      W <- simulate(ar1_model, seed = seed)
      Y <- W + rnorm(n_obs_per_rep, sd = 2)
      list(Y=Y, idx=idx, group=rep(1, n_obs_per_rep))
    },
    "matern" = {
      loc <- cbind(runif(n_obs_per_rep, 0, 10), runif(n_obs_per_rep, 0, 5))
      true_model <- f(
        map = loc,
        model="matern",
        theta_K = log(5),
        mesh = mesh,
        noise = f_simu_noise,
        debug = debug_f
      )
      W <- simulate(true_model, seed = seed)
      Y <- as.numeric(true_model$A %*% W) + rnorm(n_obs_per_rep, sd=0.5)
      # Y <- as.numeric(true_model$A %*% W) + rnig(n_obs_per_rep, mu=1,delta=-1,nu=1,sigma=1)
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
          first = f_simu_noise,
          second = f_simu_noise
        ),
        debug = debug_f
      )
      W <- simulate(true_model, seed=seed)
      AW <- as.numeric(true_model$A %*% W)
      Y <- AW + rnorm(length(AW), sd=0.5)
      list(Y=Y, idx=idx_per_rep, group=group_per_rep)
    },
    "bvmatern" = {
      loc <- cbind(runif(n_obs_per_rep, 0, 10), runif(n_obs_per_rep, 0, 5))
      group_per_rep <- c(rep("first", n_obs_per_rep/2), rep("second", n_obs_per_rep/2))
      idx_per_rep <- c(1:(n_obs_per_rep/2), 1:(n_obs_per_rep/2))
      true_model <- f(
        loc,
        model="bv",
        theta = 0, rho = 0,
        sub_models = list(
          first = list(model="matern", kappa=2),
          second= list(model="matern", kappa=4)
        ),
        mesh = mesh,
        group = group_per_rep,
        noise = list(
          first=noise_nig(mu=-3, sigma=2, nu=1),
          second=noise_nig(mu=-3, sigma=2, nu=1)
        ),
        debug = debug_f
      )
      W <- simulate(true_model, seed=seed)
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
      noise = f_form_noise,
      control = control_f(numer_grad = numer_grad)
    ),
    "matern" = Y ~ 0 + f(idx,
      model="matern",
      mesh = mesh,
      noise=f_form_noise,
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
  start_time <- proc.time()
  out <- ngme(
    formula,
    replicate = repl,
    group = group,
    data = data.frame(Y=Y),
    control_ngme = control_ngme(
      n_gibbs_samples = n_gibbs_samples,
      rao_blackwellization = rao_blackwellization
    ),
    control_opt = control_opt(
      burnin = 100,
      std_lim = 0.001,
      print_check_info = FALSE,
      seed = seed,
      num_threads = num_threads,
      iterations = n_iter,
      precond_by_diff_chain = precond_by_diff_chain,
      compute_precond_each_iter = compute_precond_each_iter,
      n_parallel_chain = n_parallel_chain,
      stop_points = stop_points,
      verbose = verbose,
      preconditioner = preconditioner,
      sampling_strategy = sampling_strategy
    ),
    start = start,
    family = family,
    debug = debug
  )

  print(proc.time() - start_time)
  list(
    out = out,
    time = proc.time() - start_time
  )
}
