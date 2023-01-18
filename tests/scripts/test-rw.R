test_that("the order of W same as order of index?", {
  XX    <- c(1.1, 3.1, 2.2, 2.2, 4.5, 5)
  index_NA <- c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
  myloc <- c(3.2, 1.2)

  mesh1 <- INLA::inla.mesh.1d(XX); mesh1

  rw <- model_rw(XX, order=1, index_NA = index_NA)
  rw$C + rw$G
  rw$A; rw$A_pred

  str(rw)

  # use INLA make A? the order is bad
  A <- INLA::inla.spde.make.A(mesh=mesh1, loc = myloc); A
# 2 x 4 sparse Matrix of class "dgCMatrix"
# [1,] .         .         0.7142857 0.2857143
# [2,] 0.4545455 0.5454545 .         .

  # the order is wrong, recover by permutation matrix
  # A %*% as(mesh1$idx$loc, "pMatrix")

# 3 x 4 sparse Matrix of class "dgCMatrix"
# [1,] 1 -1  .  .
# [2,] .  1 -1  .
# [3,] .  .  1 -1

})


# n_obs <- 1000
# rw_mu <- 4
# rw_sigma <- 1
# rw_nu <- 2.3
# rw1_process <- simulate(
#   f(1:n_obs,
#     model = "rw1",
#     noise = noise_gal(
#       mu = rw_mu,
#       sigma = rw_sigma,
#       nu = rw_nu
#     )
#   ),
#   seed = 1
# )
# Y <- ar1_process + rnorm(n_obs)
# attr(rw1_process, "noise")$h

# devtools::load_all()
# ngme_out <- ngme(
#   Y ~ 0 +
#   f(1:n_obs,
#     model = "ar1",
#     theta_K = 0.7,
#     # fix_theta_K = TRUE,
#     # W = as.numeric(ar1_process),
#     # fix_W = TRUE,
#     noise = noise_gal(
#       mu = 1.1,
#       sigma = 1,
#       # V = attr(ar1_process, "noise")$V,
#       # fix_V = TRUE
#     ),
#     control = ngme_control_f(
#       numer_grad       = F,
#       use_precond      = T
#     ),
#     debug = T
#   ),
#   data = data.frame(Y = Y),
#   family = "normal",
#   control = ngme_control(
#     estimation = T,
#     exchange_VW = TRUE,
#     n_parallel_chain = 4,
#     stop_points = 50,
#     burnin = 10,
#     iterations = 2000,
#     gibbs_sample = 5,
#     stepsize = 1,
#     threshold = 1e-4,

#     std_lim = 0.001,
#     trend_lim = 0.001
#   ),
#   seed = 10,
#   debug = TRUE
# )

# ngme_out
# traceplot(ngme_out, f_index = 1, param="alpha")
# traceplot(ngme_out, f_index = 1, param="mu")
# traceplot(ngme_out, f_index = 1, param="sigma")
# traceplot(ngme_out, f_index = 1, param="nu")