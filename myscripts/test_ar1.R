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
  ar_mu <- 4
  ar_sigma <- 3
  ar_eta <- 2

  ar1_process <- simulate(
    f(1:n_obs,
      model = "ar1",
      theta_K = 0.9,
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

  # Y <- ar1_process + nig_noise

  # use normal noise
    Y <- ar1_process + rnorm(n_obs)
}

{ #simulate using old fashion
n_obs <- 500
alpha1 <- 0.9
mu1 = 4; delta = -mu1
sigma1 = 3
nu1 = 1

trueV1 <- ngme2::rig(n_obs, nu1, nu1)
noise1 <- delta + mu1*trueV1 + sigma1 * sqrt(trueV1) * rnorm(n_obs)

trueW1 <- Reduce(function(x,y){y + alpha1*x}, noise1, accumulate = T)
Y = trueW1 + rnorm(n_obs, mean=0, sd = 0.5)
}

ngme_out <- ngme(
  Y ~ 0 +
  f(model = "ar1",
    theta_K = a2th(0.8),
    # fix_theta_K = TRUE,
    # W = as.numeric(ar1_process),
    # fix_W = TRUE,
    noise = ngme.noise.nig(
      theta_mu = 0,
      theta_sigma = ar_sigma,
      theta_V = ar_eta,
      # V = attr(ar1_process, "noise")$V,
      # fix_V = TRUE,
      fix_theta_mu      = FALSE,
      fix_theta_sigma   = F,
      fix_theta_V       = F
    ),
    control = ngme.control.f(
      numer_grad       = F,
      use_precond      = T
    ),
    debug = F
  ),
  data = data.frame(Y = Y),
  control = ngme.control(
    estimation = TRUE,
    n_parallel_chain = 4,
    stop_points = 100,
    burnin = 200,
    iterations = 1000,
    gibbs_sample = 5,
    stepsize = 1,
    kill_var = FALSE,
    threshold = 1e-4,

    std_lim = 0.01,
    trend_lim = 0.01
  ),
  noise = ngme.noise.normal(),
  # noise = attr(nig_noise, "noise"),
  seed = 10,
  # , last_fit = ngme_out
  debug = T
)

ngme_out
str(ngme_out)

# noise
traceplot(ngme_out, parameter = "theta_mu", f_index = 0)
traceplot(ngme_out, parameter = "theta_sigma", f_index = 0, transform = exp)
traceplot(ngme_out, parameter = "theta_V", f_index = 0)

# ar1 model
traceplot(ngme_out, parameter = "theta_K",     f_index = 1, transform = th2a)
traceplot(ngme_out, parameter = "theta_mu",    f_index = 1)
traceplot(ngme_out, parameter = "theta_sigma", f_index = 1, transform = exp)
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