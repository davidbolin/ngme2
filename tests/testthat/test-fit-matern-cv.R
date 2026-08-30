test_that("test fractional matern specification", {
  mesh <- fmesher::fm_mesh_1d(1:10)
  matern(mesh, alpha = 2.01)
  matern(mesh, alpha = 2.00, fix_alpha = FALSE, kappa = 1)
})


test_that("test fit matern", {
  skip_on_cran()
  set.seed(42)
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh <- fmesher::fm_mesh_2d(
    loc.domain = pl01, cutoff = 0.2,
    max.edge = c(0.5, 10)
  )
  plot(mesh)
  mesh$n

  n_obs <- 1000
  loc <- cbind(runif(n_obs, 0, 10), runif(n_obs, 0, 5))
  true_noise <- noise_nig(mu = -2, sigma = 1, nu = 0.5)

  kappa <- 0.3
  true_model <- f(
    map = loc,
    model = matern(mesh, alpha = 2, kappa = kappa, fix_alpha = TRUE),
    noise = true_noise
  )

  W <- simulate(true_model, seed = 42)[[1]]
  Y <- W + rnorm(n_obs, sd = 0.5)

  control_opt <- control_opt(
    warn_no_convergence = FALSE,
    n_min_batch = 1,
    optimizer = precond_sgd(),
    iterations = 1000,
    n_parallel_chain = 4,
    print_check_info = TRUE,
    n_batch = 10,
    robust = TRUE,
    verbose = TRUE
  )

  m_nig_gauss <- ngme(
    Y ~ 0 + f(loc,
      model = matern(
        mesh,
        alpha = 2.00
        # fix_alpha = FALSE
      ),
      noise = noise_nig()
    ),
    data = data.frame(Y = Y),
    control_opt = control_opt
  )

  m_nig_gauss
  traceplot(m_nig_gauss, "field1", hline = c(kappa, -2, 1, 0.5))
  # traceplot(m_nig_gauss, "field1", hline = c(2, kappa, -2, 1, 0.5))
  traceplot(m_nig_gauss, hline = 0.5)
})

# m_nig_gauss_2 <- ngme(
#   Y ~ 0 + f(loc,
#     model="matern",
#     name="spde",
#     mesh = mesh,
#     alpha=2,
#     noise=noise_nig(),
#     control = control_f(),
#   ),
#   data = data.frame(Y = Y),
#   control_opt = control_opt
# )
# # Total time of the estimation is (s): 25
# m_nig_gauss_2
# traceplot(m_nig_gauss_2, "spde", hline = c(0.3, -2, 1, 0.5))
# traceplot(m_nig_gauss_2, hline = 0.5)


# m_gauss_gauss <- ngme(
#   Y ~ 0 + f(loc,
#     model="matern",
#     name="spde",
#     mesh = mesh,
#     noise=noise_normal(),
#   ),
#   data = data.frame(Y = Y),
#   control_opt = control_opt
# )
# m_gauss_gauss
#   traceplot(m_gauss_gauss, "spde")

# test_that("test CV between NIG and Gauss", {
#   cv = ngme2::cross_validation(
#     list(
#       m_nig_gauss = m_nig_gauss,
#       m_gauss_gauss = m_gauss_gauss
#     ),
#     type = "k-fold",
#     k = 10,
#     n_gibbs_samples = 500,
#     seed = 42
#   )

#   expect_true(cv$mean.scores["m_nig_gauss", "MAE"] < cv$mean.scores["m_gauss_gauss", "MAE"])
#   expect_true(cv$mean.scores["m_nig_gauss", "MSE"] < cv$mean.scores["m_gauss_gauss", "MSE"])
#   expect_true(cv$mean.scores["m_nig_gauss", "neg.CRPS"] < cv$mean.scores["m_gauss_gauss", "neg.CRPS"])
#   expect_true(cv$mean.scores["m_nig_gauss", "neg.sCRPS"] < cv$mean.scores["m_gauss_gauss", "neg.sCRPS"])
