# Simple scripts for test ngme function
{
  library(devtools);load_all()
  a2th <- function(k) {log((-1-k)/(-1+k))}
  th2a <- function(th) {-1 + (2*exp(th)) / (1+exp(th))}

  seed <- 10
  set.seed(seed)
}

{ ############  1. simulate AR with nig noise
n_obs <- 500
ar_mu <- 2
ar_sigma <- 2
ar_eta <- 2

ar1_process <- simulate(
  f(1:n_obs,
    model = "ar1",
    theta_K = 0.8,
    noise = ngme.noise.nig(
      theta_mu = ar_mu,
      theta_sigma = ar_sigma,
      theta_V = ar_eta
    )
  ),
  seed = 1
)

noise_theta_mu    <- -6
noise_theta_sigma <- 2
noise_theta_V     <- 1.7

nig_noise <- simulate(
  ngme.noise.nig(
    theta_mu = noise_theta_mu,
    theta_sigma = noise_theta_sigma,
    theta_V = noise_theta_V,
  ),
  nsim = n_obs
)

Y <- ar1_process + nig_noise

# use normal noise
  Y <- ar1_process + rnorm(n_obs)
}
load_all()
ngme_out <- ngme(
  Y ~ 0 +
  f(model = "ar1",
    theta_K = 0.7,
    # W = as.numeric(ar1_process),
    # fix_W = TRUE,
    noise = ngme.noise.nig(
      theta_mu = ar_mu,
      theta_sigma = ar_sigma,
      theta_V = ar_eta,
      V = attr(ar1_process, "noise")$V,
      fix_V             = F,
      fix_theta_mu      = FALSE,
      fix_theta_sigma   = F,
      fix_theta_V       = F
    ),
    control = ngme.control.f(
      numer_grad       = FALSE,
      use_precond      = TRUE
    ),
    debug = F
  ),
  data = data.frame(Y = Y),
  control = ngme.control(
    estimation = TRUE,
    n_parallel_chain = 4,
    stop_points = 5,
    burnin = 10,
    iterations = 100,
    gibbs_sample = 5,
    stepsize = 1,
    kill_var = FALSE,
    threshold = 1e-4,

    std_lim = 0.1,
    trend_lim = 0.1
  ),
  noise = ngme.noise.normal(),
  # noise = attr(nig_noise, "noise"),
  seed = 2,
  # , last_fit = ngme_out
  debug = F
)

ngme_out
str(ngme_out)

# noise
traceplot(ngme_out, parameter = "theta_mu", f_index = 0)
traceplot(ngme_out, parameter = "theta_V", f_index = 0)
traceplot(ngme_out, parameter = "theta_sigma", f_index = 0, transform = exp)

# ar1 model
traceplot(ngme_out, parameter = "theta_K",     f_index = 1)
traceplot(ngme_out, parameter = "theta_mu",    f_index = 1)
traceplot(ngme_out, parameter = "theta_sigma", f_index = 1)
traceplot(ngme_out, parameter = "theta_V",     f_index = 1)

# compare ar model
# plot(ngme_out$latents[[1]]$noise, col = "red")
plot(ngme.noise.nig(
      theta_mu = ar_mu,
      theta_sigma = ar_sigma,
      theta_V = ar_eta
    ),
    ngme_out$latents[[1]]$noise
)