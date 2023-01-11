devtools::load_all()

n_obs <- 100
ar_mu <- 4
ar_sigma <- 1
ar_nu <- 2.3

ar1_process <- simulate(
  f(1:n_obs,
    model = "ar1",
    theta_K = 0.8,
    noise = noise_gal(
      mu = ar_mu,
      sigma = ar_sigma,
      nu = ar_nu
    )
  ),
  seed = 1
)

noise_theta_mu    <- -6
noise_theta_sigma <- 2
noise_nu     <- 1.7

nig_noise <- simulate(
  noise_nig(
    theta_mu = noise_theta_mu,
    theta_sigma = noise_theta_sigma,
    nu = noise_nu,
    n = n_obs
  )
)

beta <- c(-3, 1, 2)
x1 <- rnorm(n_obs)
x2 <- rexp(n_obs)

# Y <- ar1_process + nig_noise
Y <- ar1_process + rnorm(n_obs) + beta[1] + x1 * beta[2] + x2 * beta[3]

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
#     iterations = 20,
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


# str(ngme_out$latents[[1]]$noise)
# str(ngme_out$latents[[1]]$W)

# #

# load_all()
#   ngme_out2 <- ngme(
#     Y ~ 0 +
#     f(1:n_obs,
#       model = "ar1",
#       noise = noise_nig(),
#       control = ngme_control_f(
#         numer_grad       = F,
#         use_precond      = T
#       ),
#       debug = T
#     ),
#     data = data.frame(Y = Y),
#     family = noise_nig(),
#     start = ngme_out,
#     control = ngme_control(
#       estimation = FALSE
#     )
#   )
# str(ngme_out$noise)
#   str(ngme_out$latents[[1]]$noise)

#   str(ngme_out2$latents[[1]]$W)
#   ngme_out2$noise$V <- NULL

#   lll <- list(a = NULL , b = 2); lll
#   lll$a <- NULL  ; lll
#   # plot
#   # noise
#   traceplot(ngme_out, parameter = "theta_mu", f_index = 0)
#   traceplot(ngme_out, parameter = "theta_sigma", f_index = 0)
#   traceplot(ngme_out, parameter = "nu", f_index = 0)

#   # ar1 model
#   traceplot(ngme_out, parameter = "theta_K",     f_index = 1)
#   traceplot(ngme_out, parameter = "theta_mu",    f_index = 1)
#   traceplot(ngme_out, parameter = "theta_sigma", f_index = 1)
#   traceplot(ngme_out, parameter = "nu",     f_index = 1)

#   # compare nig noise
#   plot(noise_nig(
#         mu = ar_mu,
#         sigma = ar_sigma,
#         nu = ar_nu
#       ),
#       ngme_out$latents[[1]]$noise
#   )

# devtools::check()


# test_ar1 with normal_nig noise
x3 <- rexp(100)
x4 <- rexp(100)
ngme_out <- ngme(
  Y ~ x1 + x2 + x3 + x4 +
  f(1:n_obs,
    model = "ar1",
    theta_K = 0.7,
    # fix_theta_K = TRUE,
    # W = as.numeric(ar1_process),
    # fix_W = TRUE,
    # noise = noise_nig(
    #   mu=1,
    #   sigma = 0.7,
    #   nu = 0.6
    # ),
    noise = noise_normal_nig(
      mu=1,
      sigma_normal = 1.5,
      sigma_nig = 0.7,
      nu = 0.6
    ),
    control = ngme_control_f(
      numer_grad       = F,
      use_precond      = T
    ),
    debug = T
  ),
  data = data.frame(Y = Y),
  family = "nig",
  control = ngme_control(
    estimation = T,
    exchange_VW = TRUE,
    n_parallel_chain = 1,
    stop_points = 50,
    burnin = 10,
    iterations = 3,
    gibbs_sample = 5,
    stepsize = 1,
    threshold = 1e-4,

    std_lim = 0.001,
    trend_lim = 0.001
  ),
  seed = 10,
  debug = TRUE
)

p1 <- traceplot(ngme_out, f_index = 0, param = "beta", param_index = 1)
p2 <- traceplot(ngme_out, f_index = 0, param = "beta", param_index = 2)
p3 <- traceplot(ngme_out, f_index = 0, param = "beta", param_index = 3)
grid.arrange(p1, p2, p3, ncol=2)

load_all()
traceplot2(ngme_out)
traceplot2(ngme_out, param=1)

ngme_out$latents[[1]]$noise$theta_sigma_normal
ngme_out$beta

lapply(1:3, function(x) traceplot(ngme_out, f_index=0, param="beta", param_index = x))


ngme_out$latents[[1]]

ngme_out$latents[[1]]$theta_K
