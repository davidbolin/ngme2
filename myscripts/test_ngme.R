# Simple scripts for test ngme function
# sth. wrong with nig measurement noise

{ 
a2th <- function(k) {log((-1-k)/(-1+k))}
th2a <- function(th) {-1 + (2*exp(th)) / (1+exp(th))}
library(devtools); load_all()

seed <- 51
set.seed(seed)
}

{ ############  1. simulate AR with nig noise
n_obs <- 500
ar_mu <- 4
ar_sigma <- 1.3
ar_eta <- 0.8

ar1_process <- ngme.simulate(
  f(1:n_obs,
    model = "ar1",
    noise = ngme.noise.nig(
      theta_mu = ar_mu,
      theta_sigma = ar_sigma,
      theta_V = ar_eta
    )
  ),
  seed = 123
)

noise_theta_mu    <- -6
noise_theta_sigma <- 2
noise_theta_V     <- 1.7

nig_noise <- ngme.simulate.noise(
  ngme.noise.nig(
    theta_mu = noise_theta_mu,
    theta_sigma = noise_theta_sigma,
    theta_V = noise_theta_V,

    fix_sigma = TRUE,
    fix_var = TRUE,
    fix_V = FALSE
  ),
  n = n_obs
)

Y <- ar1_process + nig_noise
}

ngme_control <- ngme.control(
  estimation = TRUE,

  burnin = 200,
  iterations = 500,
  gibbs_sample = 5,
  stepsize = 1,
  kill_var = FALSE,
  threshold = 1e-4,
  opt_beta = T
)

ngme_out <- ngme(
  Y ~ 0 +
  f(model = "ar1",
    theta_K = 0.4,
    W = ar1_process,
    noise = ngme.noise(
      theta_mu = ar_mu,
      theta_sigma = ar_sigma,
      theta_V = ar_eta,
      V = attr(ar1_process, "noise")$V,
      fix_mu           = TRUE,
      fix_sigma        = TRUE,
      fix_var          = TRUE,
      fix_V            = FALSE
    ),
    control = ngme.control.f(
      numer_grad       = FALSE,
      use_precond      = TRUE,
      fix_operator     = TRUE,
      fix_W            = FALSE
    ),
    debug = TRUE
  ),
  data = data.frame(Y = Y),
  control = ngme_control,
  noise = attr(nig_noise, "noise"),
  debug = ngme.debug(
    debug = TRUE,
    not_run = FALSE
  ),
  seed = 1
  # , last_fit = ngme_out
)

str(ngme_out$est_output)
# c(noise_theta_mu, noise_theta_sigma, noise_theta_V)

# trace plot of mu
ngme.traceplot(ngme_out$opt_trajectory, start = 5, n = 1)

plot(attr(nig_noise, "noise")) # measurement noise
plot(create.ngme.noise(ngme_out$est_output$noise), add = TRUE, col="blue")

