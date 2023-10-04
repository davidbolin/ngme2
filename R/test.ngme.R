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
  n_iter,
  n_replicate=1,
  n_parallel_chain = 4,
  stop_points = 10,
  stepsize = 0.5,
  numer_grad = TRUE,
  preconditioner = "fast",
  sampling_strategy = "all",
  precond_by_diff_chain = TRUE,
  compute_precond_each_iter = FALSE,
  num_threads = c(n_parallel_chain, 4),
  f_noise = noise_nig(mu = -2, sigma = 1.5, nu=0.5),
  n_gibbs_samples = 5,
  family = "normal",
  max.n = 1000,
  print = FALSE,
  debug = FALSE,
  debug_f = FALSE,
  verbose = FALSE,
  rao_blackwellization = FALSE,
  precond_eps = 1e-5,
  n_trace_iter = 10,
  start = NULL,
  estimation = TRUE,
  fix_theta_K = FALSE,
  fix_theta_mu= FALSE,
  fix_theta_sigma = FALSE,
  fix_nu = FALSE,
  seed = Sys.time()
) {
  set.seed(seed)
  # create 2d mesh
  if (model %in% c("matern", "bvmatern", "ar1+matern")) {
    pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
    mesh <- fmesher::fm_mesh_2d(
      loc.domain = pl01, cutoff = 0.2,
      max.edge = c(0.4, 10),
      max.n = max.n
    )
print(paste("nodes of mesh = ", mesh$n))
  }

  stopifnot(length(f_noise$noise_type) == 1)
  f_fm_noise <- f_noise
  f_fm_noise$fix_theta_mu = fix_theta_mu
  f_fm_noise$fix_theta_sigma = fix_theta_sigma
  f_fm_noise$fix_nu = fix_nu
  f_fm_noise$theta_mu = rep(0, ncol(f_fm_noise$B_mu))
  f_fm_noise$theta_sigma = rep(0, ncol(f_fm_noise$B_sigma))
  f_fm_noise$nu = 1

  mn_noise <- if (family == "nig")
    # rnig(n_obs_per_rep, delta=-3, mu=3, nu=0.8, sigma=0.5, seed=seed)
    rnig(n_obs_per_rep, delta=0, mu=0, nu=0.8, sigma=0.5, seed=seed)
  else if (family == "cor_normal") {
    family = noise_normal(
      corr_measurement=TRUE,
      index_corr=rep(1:(n_obs_per_rep/2), 2)
    )
    stopifnot(n_obs_per_rep %% 2 == 0)
    rho = -0.5
    Cov_kron <- matrix(c(.5, rho*.5, rho*.5, .5), nrow=2) %x% diag(n_obs_per_rep / 2)
    L <- t(chol(Cov_kron))
    as.numeric(L %*% rnorm(n_obs_per_rep))
  } else
    rnorm(n_obs_per_rep, sd=0.2)

  # ------- Simulate data for each model --------
  sim_data <- switch(model,
    "none" = {
      list(Y = mn_noise, group=rep(1, n_obs_per_rep))
    },
    "ar1" = {
      idx <- 1:n_obs_per_rep
      ar1_model <- f(idx, model="ar1", rho = 0.5,
        noise = f_noise
      )
      W <- simulate(ar1_model, seed = seed)
      Y <- W + mn_noise
      list(Y=Y, idx=idx, group=rep(1, n_obs_per_rep))
    },
    "matern" = {
      loc <- cbind(runif(n_obs_per_rep, 0, 10), runif(n_obs_per_rep, 0, 5))
      matern_model <- f(
        map = loc,
        model="matern",
        theta_K = log(4),
        mesh = mesh,
        noise = f_noise
      )
      W <- simulate(matern_model, seed=seed)
      Y <- as.numeric(matern_model$A %*% W) + mn_noise
      list(Y=Y, idx=loc, group=rep(1, n_obs_per_rep))
    },
    "ar1+ar1" = {
      idx <- 1:n_obs_per_rep
      ar1_model_1 <- f(idx, model="ar1", rho = 0.5, noise = f_noise)
      ar1_model_2 <- f(idx, model="ar1", rho = -0.6, noise = f_noise)
      W1 <- simulate(ar1_model_1, seed = seed)
      W2 <- simulate(ar1_model_2, seed = seed)
      Y <- W1 + W2 + mn_noise
      list(Y=Y, idx=idx, group=rep(1, n_obs_per_rep))
    },
    "ar1+matern" = {
      idx <- 1:n_obs_per_rep
      ar1_model <- f(idx, model="ar1", rho = 0.5,noise = f_noise)
      loc <- cbind(runif(n_obs_per_rep, 0, 10), runif(n_obs_per_rep, 0, 5))
      matern_model <- f(
        map = loc, model="matern",
        theta_K = log(4), mesh = mesh,
        noise = f_noise
      )
      W1 <- simulate(ar1_model, seed = seed)
      W2 <- simulate(matern_model, seed = seed)
      Y <- W1 + as.numeric(matern_model$A %*% W2) + mn_noise
      list(Y=Y, idx=list(idx, loc), group=rep(1, n_obs_per_rep))
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
          first = f_noise,
          second = f_noise
        )
      )
      W <- simulate(true_model, seed=seed)
      Y <- as.numeric(true_model$A %*% W) + mn_noise
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
        )
      )
      W <- simulate(true_model, seed=seed)
      Y <- as.numeric(true_model$A %*% W) + mn_noise
      list(Y=Y, idx=idx_per_rep, group=group_per_rep)
    },
    stop("Unknown test model")
  )

  # ------- Specify formula for each model -------
  formula <- switch(model,
    "none" = Y ~ 0,
    "ar1" = Y ~ 1 + f(idx,
      fix_theta_K = fix_theta_K,
      model = "ar1",
      noise = f_fm_noise,
      debug = debug_f,
      control = control_f(numer_grad = numer_grad)
    ),
    "matern" = Y ~ 0 + f(idx,
      fix_theta_K = fix_theta_K,
      model="matern",
      mesh = mesh,
      noise=f_fm_noise,
      debug = debug_f,
      control = control_f(numer_grad = numer_grad)
    ),
    "ar1+ar1" =  Y ~ 0 + f(idx,
      fix_theta_K = fix_theta_K,
      model = "ar1",
      noise = f_fm_noise,
      debug = debug_f,
      name = "field1",
      control = control_f(numer_grad = numer_grad)
      ) + f(idx,
      fix_theta_K = fix_theta_K,
      model = "ar1",
      noise = f_fm_noise,
      debug = debug_f,
      name = "field2",
      control = control_f(numer_grad = numer_grad)
    ),
    "ar1+matern" = Y ~ 0 + f(idx[[1]],
      fix_theta_K = fix_theta_K,
      model = "ar1",
      noise = f_fm_noise,
      debug = debug_f,
      control = control_f(numer_grad = numer_grad)
      ) + f(idx[[2]],
      fix_theta_K = fix_theta_K,
      model="matern",
      mesh = mesh,
      noise=f_fm_noise,
      debug = debug_f,
      control = control_f(numer_grad = numer_grad)
    ),
    "bvar1" = Y ~ f(idx,
      model="bv",
      fix_theta_K = fix_theta_K,
      sub_models = list(first = "ar1", second="ar1"),
      control = control_f(numer_grad = numer_grad),
      debug = debug_f,
      noise = list(first=f_fm_noise, second=f_fm_noise)
    ),
    "bvmatern" = Y ~ f(idx,
      model="bv",
      fix_theta_K = fix_theta_K,
      sub_models = list(first = "matern", second="matern"),
      debug = debug_f,
      control = control_f(numer_grad = numer_grad),
      noise = list(first=f_fm_noise, second=f_fm_noise)
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
      n_gibbs_samples = n_gibbs_samples
    ),
    control_opt = control_opt(
      seed = seed,
      burnin = 0,
      std_lim = 0.001,
      print_check_info = FALSE,
      precond_eps = precond_eps,
      rao_blackwellization = rao_blackwellization,
      n_trace_iter = n_trace_iter,
      num_threads = num_threads,
      iterations = n_iter,
      precond_by_diff_chain = precond_by_diff_chain,
      compute_precond_each_iter = compute_precond_each_iter,
      n_parallel_chain = n_parallel_chain,
      stop_points = stop_points,
      verbose = verbose,
      preconditioner = preconditioner,
      estimation = estimation,
      stepsize = stepsize,
      sampling_strategy = sampling_strategy
    ),
    start = start,
    family = family,
    debug = debug
  )

  print(proc.time() - start_time)
  list(
    out = out,
    time = proc.time() - start_time,
    f_noise = f_noise,
    m_noise = if (inherits(family, "ngme_noise")) family
      else if (family=="nig") noise_nig(mu=0, sigma=0.5, nu=0.8)
      else noise_normal(sigma=0.2)
  )
}
