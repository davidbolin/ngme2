# ignored for now, the github action run out of time

# setup({
#   set.seed(42)
#   n_obs <- 1000
#   day <- 1:n_obs
#   ar1_model <- ngme2::f(day, model="ar1", rho = 0.5,
#     noise = noise_nig(mu = -3, sigma = 4, nu=0.4))
#   W <- simulate(ar1_model, seed = 16, nsim=1)[[1]]
#   Y <-  W + rnorm(n_obs, sd = 2)

#   # Fit the model with the AR1 model
#   control_opt = control_opt(
#     optimizer = adam(),
#     burnin = 100,
#     iterations = 2000,
#     std_lim = 0.01,
#     n_parallel_chain = 4,
#     stop_points = 10,
#     verbose = FALSE,
#     print_check_info = FALSE,
#     seed = 3
#   )

#   m_nig_gauss <- ngme(
#     Y ~ 0 + f(
#       1:n_obs,
#       name = "my_ar",
#       model = "ar1",
#       noise = noise_nig(),
#       control = control_f(numer_grad = TRUE, improve_hessian = FALSE)
#     ),
#     data = data.frame(Y=Y),
#     control_opt = control_opt
#   )
#   m_nig_gauss
#   traceplot(m_nig_gauss, hline = c(2))
#   traceplot(m_nig_gauss, "my_ar", hline = c(0.5, -3, 4, 0.4))

#   m_gauss_gauss <- ngme(
#     Y ~ 0 + f(
#       1:n_obs,
#       model = "ar1",
#       name = "my_ar",
#       noise = noise_normal()
#     ),
#     data = data.frame(Y=Y),
#     control_opt = control_opt
#   )
#   m_gauss_gauss
#   traceplot(m_gauss_gauss, hline = c(2))
#   traceplot(m_gauss_gauss, "my_ar")
# })

# # Seem it's not obvious in MAE, and MSE
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
#   cv
#   expect_true(cv$mean.scores["m_nig_gauss", "neg.CRPS"] < cv$mean.scores["m_gauss_gauss", "neg.CRPS"])
#   expect_true(cv$mean.scores["m_nig_gauss", "neg.sCRPS"] < cv$mean.scores["m_gauss_gauss", "neg.sCRPS"])
# })