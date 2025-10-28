setup({
  set.seed(42)
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh <- fmesher::fm_mesh_2d(
    loc.domain = pl01, cutoff = 0.1,
    max.edge = c(0.3, 10)
  )
  mesh$n #2163

  n_obs <- 1000
  loc <- cbind(runif(n_obs, 0, 10), runif(n_obs, 0, 5))
  true_noise = noise_nig(mu=-2, sigma=1, nu=0.5)

  # NIG model
  true_model <- f(
    map = loc,
    model="matern",
    theta_K = log(0.3),
    mesh = mesh,
    noise = true_noise
  )

  W <- simulate(true_model)[[1]]
  Y <- W + rnorm(n_obs, sd=0.5)

  control_opt <- control_opt(
    iterations = 2000,
    n_parallel_chain = 4,
    print_check_info = F,
    rao_blackwellization = TRUE,
    optimizer = adam(),
    std_lim = 0.01
  )

  m_nig_gauss <- ngme(
    Y ~ 0 + f(loc,
      model="matern",
      name="spde",
      mesh = mesh,
      noise=noise_nig(),
      control = control_f(),
    ),
    data = data.frame(Y = Y),
    control_ngme = control_ngme(n_gibbs_samples = 3),
    control_opt = control_opt
  )
  # Total time of the estimation is (s): 125 (gibbs = 5)

  m_nig_gauss
  ngme_result(m_nig_gauss)
  traceplot(m_nig_gauss, "spde", hline = c(0.3, -2, 1, 0.5))
  traceplot(m_nig_gauss, hline = 0.5)

  m_gauss_gauss <- ngme(
    Y ~ 0 + f(loc,
      model="matern",
      name="spde",
      mesh = mesh,
      noise=noise_normal(),
    ),
    data = data.frame(Y = Y),
    control_opt = control_opt
  )
  m_gauss_gauss
  traceplot(m_gauss_gauss, "spde")
})

test_that("test CV between NIG and Gauss", {
  cv = ngme2::cross_validation(
    list(
      m_nig_gauss = m_nig_gauss,
      m_gauss_gauss = m_gauss_gauss
    ),
    N_sim=1,
    print=TRUE,
    type = "k-fold",
    k = 5,
    n_gibbs_samples = 500,
    seed = 42
  )
  cv

  expect_true(cv$mean.scores["m_nig_gauss", "MAE"] < cv$mean.scores["m_gauss_gauss", "MAE"])
  expect_true(cv$mean.scores["m_nig_gauss", "MSE"] < cv$mean.scores["m_gauss_gauss", "MSE"])
  expect_true(cv$mean.scores["m_nig_gauss", "neg.CRPS"] < cv$mean.scores["m_gauss_gauss", "neg.CRPS"])
  expect_true(cv$mean.scores["m_nig_gauss", "neg.sCRPS"] < cv$mean.scores["m_gauss_gauss", "neg.sCRPS"])
})
