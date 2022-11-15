devtools::load_all()

n_obs <- 500
ar_mu <- 4
ar_sigma <- 1
ar_nu <- 1
ar1_process <- simulate(
  f(1:n_obs,
    model = "ar1",
    theta_K = 0.9,
    noise = noise_nig(
      mu = ar_mu,
      sigma = ar_sigma,
      nu = ar_nu
    )
  ),
  seed = 1
)

noise_theta_mu    <- -6
noise_theta_sigma <- 2
noise_theta_V     <- 1.7

nig_noise <- simulate(
  noise_nig(
    theta_mu = noise_theta_mu,
    theta_sigma = noise_theta_sigma,
    theta_V = noise_theta_V,
    n = n_obs
  )
)
# Y <- ar1_process + nig_noise
Y <- ar1_process + rnorm(n_obs)

load_all()
  ngme_out <- ngme(
    Y ~ 0 +
    f(1:n_obs,
      model = "ar1",
      theta_K = 0.2,
      # fix_theta_K = TRUE,
      # W = as.numeric(ar1_process),
      # fix_W = TRUE,
      noise = noise_nig(
        mu = 1.1,
        sigma = 1,
        fix_theta_mu      = F,
        fix_theta_sigma   = F,
        fix_theta_V       = F
      ),
      control = ngme_control_f(
        numer_grad       = F,
        use_precond      = T
      ),
      debug = T
    ),
    data = data.frame(Y = Y),
    family = "normal",
    control = ngme_control(
      estimation = T,
      exchange_VW = TRUE,
      n_parallel_chain = 8,
      stop_points = 50,
      burnin = 200,
      iterations = 100,
      gibbs_sample = 5,
      stepsize = 1,
      threshold = 1e-4,

      std_lim = 0.001,
      trend_lim = 0.001
    ),
    seed = 10,
    debug = TRUE
  )
  str(ngme_out$latents[[1]]$noise)
  str(ngme_out$latents[[1]]$W)

load_all()
  ngme_out2 <- ngme(
    Y ~ 0 +
    f(1:n_obs,
      model = "ar1",
      noise = noise_nig(),
      control = ngme_control_f(
        numer_grad       = F,
        use_precond      = T
      ),
      debug = T
    ),
    data = data.frame(Y = Y),
    family = noise_nig(),
    start = ngme_out,
    control = ngme_control(
      estimation = FALSE
    )
  )
str(ngme_out$noise)
  str(ngme_out$latents[[1]]$noise)

  str(ngme_out2$latents[[1]]$W)
  ngme_out2$noise$V <- NULL

  lll <- list(a = NULL , b = 2); lll
  lll$a <- NULL  ; lll
  # plot
  # noise
  traceplot(ngme_out, parameter = "theta_mu", f_index = 0)
  traceplot(ngme_out, parameter = "theta_sigma", f_index = 0)
  traceplot(ngme_out, parameter = "theta_V", f_index = 0)

  # ar1 model
  traceplot(ngme_out, parameter = "theta_K",     f_index = 1)
  traceplot(ngme_out, parameter = "theta_mu",    f_index = 1)
  traceplot(ngme_out, parameter = "theta_sigma", f_index = 1)
  traceplot(ngme_out, parameter = "theta_V",     f_index = 1)

  # compare nig noise
  plot(noise_nig(
        mu = ar_mu,
        sigma = ar_sigma,
        nu = ar_nu
      ),
      ngme_out$latents[[1]]$noise
  )

devtools::check()
