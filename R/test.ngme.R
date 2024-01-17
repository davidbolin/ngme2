#' @title Test ngme function
#'
#' @description
#' Test ngme function for different models
#' @param model model name
#' @param n_obs_per_rep number of observation per replicate
#' @param n_replicate number of replicate
#' @param numer_grad numerical gradient
#' @param max.n maximum number for building mesh
#' @param debug debug mode
#' @param debug_f debug mode for latent process
#' @param f_noise noise function
#' @param n_gibbs_samples number of gibbs samples
#' @param family family of noise
#' @param seed seed
#' @param start start value for optimization
#' @param fix_theta_K fix theta_K
#' @param fix_theta_mu fix theta_mu
#' @param fix_theta_sigma fix theta_sigma
#' @param fix_nu fix nu
#' @param control_opt control options for optimization, see \code{\link{control_opt}}
#' @export
test_ngme <- function(
  model,
  n_obs_per_rep,
  n_replicate=1,
  control_opt = NULL,
  numer_grad = TRUE,
  f_noise = noise_nig(mu = -2, sigma = 1.5, nu=0.5),
  n_gibbs_samples = 5,
  family = "normal",
  max.n = 1000,
  debug = FALSE,
  debug_f = FALSE,
  start = NULL,
  fix_theta_K = FALSE,
  fix_theta_mu= FALSE,
  fix_theta_sigma = FALSE,
  fix_nu = FALSE,
  seed = Sys.time()
) {
  if (is.null(control_opt)) control_opt <- control_opt()

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

  if (family == "nig") {
    mn_noise <- rnig(n_obs_per_rep, delta=2, mu=-2, nu=1, sigma=0.5, seed=seed)
    real_mn_noise <- noise_nig(mu=-2,sigma=0.5,nu=1)
    fm_mn_noise <- noise_nig(
      # nu=2, fix_nu=TRUE
      # V = attr(mn_noise, "V"), fix_V = TRUE
    )
  } else if (family == "normal") {
    mn_noise <- rnorm(n_obs_per_rep, sd=0.1)
    real_mn_noise <- noise_normal(sigma=0.1)
    fm_mn_noise <- noise_normal()
  } else if (family == "cor_normal") {
    stopifnot(n_obs_per_rep %% 2 == 0)
    rho = -0.5
    Cov_kron <- matrix(c(.5, rho*.5, rho*.5, .5), nrow=2) %x% diag(n_obs_per_rep / 2)
    L <- t(chol(Cov_kron))
    mn_noise <- as.numeric(L %*% rnorm(n_obs_per_rep))
    real_mn_noise <- noise_normal(
      corr_measurement=TRUE,
      index_corr=rep(1:(n_obs_per_rep/2), 2),
      rho = -0.5
    )
    fm_mn_noise = noise_normal(
      corr_measurement=TRUE,
      index_corr=rep(1:(n_obs_per_rep/2), 2)
    )
  } else if (family == "cor_nig") {
    # simulation of cor_nig?
    stopifnot(n_obs_per_rep %% 2 == 0)
    rho = -0.5
    Cov_kron <- matrix(c(.5, rho*.5, rho*.5, .5), nrow=2) %x% diag(n_obs_per_rep / 2)
    L <- t(chol(Cov_kron))
    mn_noise <- as.numeric(L %*% rnorm(n_obs_per_rep))
    real_mn_noise <- noise_nig(
      corr_measurement=TRUE,
      index_corr=rep(1:(n_obs_per_rep/2), 2),
      rho = -0.5
    )
    fm_mn_noise = noise_nig(
      corr_measurement=TRUE,
      index_corr=rep(1:(n_obs_per_rep/2), 2)
    )
  }

  # ------- Simulate data for each model --------
  sim_data <- switch(model,
    "none" = {
      list(Y = mn_noise, group=rep(1, n_obs_per_rep))
    },
    "iid" = {
      idx <- 1:n_obs_per_rep
      iid_model <- f(idx, model="iid",
        noise = f_noise
      )
      W <- simulate(iid_model, seed = seed)
      Y <- W + mn_noise
      list(Y=Y, idx=idx, group=rep(1, n_obs_per_rep))
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
      loc <- cbind(stats::runif(n_obs_per_rep, 0, 10), stats::runif(n_obs_per_rep, 0, 5))
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
      loc <- cbind(stats::runif(n_obs_per_rep, 0, 10), stats::runif(n_obs_per_rep, 0, 5))
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
      group_per_rep <- c(rep("first", n_obs_per_rep/2), rep("second", n_obs_per_rep/2))
      idx_per_rep <- c(1:(n_obs_per_rep/2), 1:(n_obs_per_rep/2))
      loc <- cbind(stats::runif(n_obs_per_rep/2, 0, 10), stats::runif(n_obs_per_rep/2, 0, 5))
      loc <- rbind(loc, loc)
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
      list(Y=Y, idx=loc, group=group_per_rep)
    },
    "graph" = {
      # library(MetricGraph)
      edge1 <- rbind(c(0,0),c(1,0))
      edge2 <- rbind(c(0,0),c(0,1))
      edge3 <- rbind(c(0,1),c(-1,1))
      theta <- seq(from=pi,to=3*pi/2,length.out = 20)
      edge4 <- cbind(sin(theta),1+ cos(theta))
      graph <- MetricGraph::metric_graph$new(edges = list(edge1, edge2, edge3, edge4))
      # graph <- MetricGraph::metric_graph$new(edges = list(edge1))
      graph$build_mesh(h = 0.02)

      matern_graph <- f(
        model="matern",
        theta_K = log(8),
        mesh=graph,
        noise=f_noise
      )

      W <- simulate(matern_graph, seed=seed)

      # build observation and A matrices
      obs.per.edge <- n_obs_per_rep / graph$nE
      obs.loc <- NULL
      for(i in 1:graph$nE) {
        obs.loc <- rbind(obs.loc,
                        cbind(rep(i,obs.per.edge), stats::runif(obs.per.edge)))
      }
      A <- graph$fem_basis(obs.loc)
      Y <- as.numeric(A %*% W) + mn_noise

      df_data <- data.frame(
        Y = Y, edge_number = obs.loc[,1],
        distance_on_edge = obs.loc[,2])

      graph$clear_observations()
      graph$add_observations(data = df_data, normalized = TRUE)
  # browser()

      list(Y=graph$get_data()$Y, idx=NULL, group=rep(1, n_obs_per_rep))
    },
    stop("Unknown test model")
  )

  # ------- Specify formula for each model -------
  formula <- switch(model,
    "none" = Y ~ 0,
    "iid" = Y ~ 0 + f(idx,
      fix_theta_K = fix_theta_K,
      model = "iid",
      noise = f_fm_noise,
      debug = debug_f,
      control = control_f(numer_grad = numer_grad)
    ),
    "ar1" = Y ~ 0 + f(idx,
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
      mesh = mesh,
      debug = debug_f,
      control = control_f(numer_grad = numer_grad),
      noise = list(first=f_fm_noise, second=f_fm_noise)
    ),
    "graph" = Y ~ 0 + f(idx,
      fix_theta_K = fix_theta_K,
      model="matern",
      mesh = graph,
      noise=f_fm_noise,
      debug = debug_f,
      control = control_f(numer_grad = numer_grad)
    )
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
    control_opt = control_opt,
    start = start,
    family = fm_mn_noise,
    debug = debug
  )

  print(proc.time() - start_time)
  list(
    out = out,
    time = proc.time() - start_time,
    f_noise = f_noise,
    m_noise = real_mn_noise
  )
}
