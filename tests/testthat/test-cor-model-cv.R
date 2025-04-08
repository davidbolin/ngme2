test_that("cross validation for AR1 correlated noise model", {
  n_obs <- 2000
  n_each <- n_obs / 2

  sim_fields <- simulate(
    f(1:n_obs, model="ar1", rho=0.6),
    noise = noise_normal()
  )[[1]]

  sd_1 = 1; sd_2 = 1; rho_e = 0.9
  Cov_same_idx <- matrix(c(sd_1^2, rho_e*sd_1*sd_2, rho_e*sd_1*sd_2, sd_2^2), nrow=2)

  print("The covariance matrix for 2 correlated fields: ")
  print(Cov_same_idx)

  tmp <- replicate(n_each, Cov_same_idx, simplify = FALSE)
  Cov_measurement <- Matrix::bdiag(tmp)

  # e ~ N(0, Cov_measurement)
  L <- t(chol(Cov_measurement))
  e <- L %*% rnorm(n_obs)

  y <- sim_fields + as.numeric(e)
  y

  control_opt <- control_opt(
    iterations = 1000
  )

  m_cor <- ngme(
    y ~ 0 + f(1:n_obs, model="ar1", rho=0.6),
    family = noise_normal(
      corr_measurement = TRUE,
      index_corr = rep(1:n_each, each=2)
    ),
    data = data.frame(y=y),
    control_opt = control_opt
  )
  m_cor
  traceplot(m_cor, hline = c(1, 0.9))
  traceplot(m_cor, "field1", hline = c(0.6, 1))

  m_neg_cor <- ngme(
    y ~ 0 + f(1:n_obs, model="ar1", rho=0.6),
    family = noise_normal(
      corr_measurement = TRUE,
      index_corr = rep(1:n_each, each=2),
      rho = -0.9
    ),
    data = data.frame(y=y),
    control_opt = control_opt(estimation = FALSE)
  )
  m_neg_cor


  m_ind <- ngme(
    y ~ 0 + f(1:n_obs, model="ar1", rho=0.6),
    family = noise_normal(),
    data = data.frame(y=y),
    control_opt = control_opt(estimation = FALSE)
  )
  traceplot(m_ind)
  traceplot(m_ind, "field1", hline = c(0.6, 1))

  cv_res <- ngme2::cross_validation(
    list(
      m_cor = m_cor,
      m_neg_cor = m_neg_cor,
      m_ind = m_ind
    ),
    type = "k-fold",
    k = 5,
    seed = 30
  )

  cv_res
})


test_that("cross validation for AR1 correlated noise model", {
  set.seed(42)
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh <- fmesher::fm_mesh_2d(
    loc.domain = pl01, cutoff = 0.1,
    max.edge = c(0.3, 10)
  )

  n_obs <- 1000
  loc <- cbind(runif(n_obs, 0, 10), runif(n_obs, 0, 5))
  true_noise = noise_nig(mu=-2, sigma=1, nu=0.5)

  true_model <- f(
    map = loc,
    model="matern",
    theta_K = log(0.3),
    mesh = mesh,
    noise = true_noise
  )

  W <- simulate(true_model)[[1]]

  n_each <- n_obs / 2

  sd_1 = 1; sd_2 = 1; rho_e = 0.9
  Cov_same_idx <- matrix(c(sd_1^2, rho_e*sd_1*sd_2, rho_e*sd_1*sd_2, sd_2^2), nrow=2)

  print("The covariance matrix for 2 correlated fields: ")
  print(Cov_same_idx)

  tmp <- replicate(n_each, Cov_same_idx, simplify = FALSE)
  Cov_measurement <- Matrix::bdiag(tmp)

  # e ~ N(0, Cov_measurement)
  L <- t(chol(Cov_measurement))
  e <- L %*% rnorm(n_obs)

  y <- W + as.numeric(e)

  control_opt <- control_opt(
    estimation = FALSE
  )

  m_cor <- ngme(
    y ~ 0 + f(loc,
      model="matern",
      name="spde",
      theta_K = log(0.3),
      mesh = mesh,
      noise=true_noise
    ),
    data = data.frame(y = y),
    family = noise_normal(
      corr_measurement = TRUE,
      rho = 0.9,
      index_corr = rep(1:n_each, each=2)
    ),
    control_opt = control_opt
  )
  m_cor

  m_ind <- ngme(
    y ~ 0 + f(loc,
      model="matern",
      name="spde",
      mesh = mesh,
      theta_K = log(0.3),
      noise=true_noise
    ),
    data = data.frame(y=y),
    family = noise_normal(),
    control_opt = control_opt
  )
  m_ind

  cv_res <- cross_validation(
    list(
      m_cor = m_cor,
      m_ind = m_ind
    ),
    type = "k-fold",
    k = 10,
    n_gibbs_samples = 500,
    seed = 42,
    print = TRUE
  )

  cv_res
})