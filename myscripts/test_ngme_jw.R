# Simple scripts for test ngme function
# sth. wrong with nig measurement noise
library(ngme2)
{
a2th <- function(k) {log((-1-k)/(-1+k))}
th2a <- function(th) {-1 + (2*exp(th)) / (1+exp(th))}
library(devtools); load_all()


#set.seed(seed)
}

{ ############  1. simulate AR with nig noise

n_obs <- 50

noise_theta_mu    <- -7
noise_theta_sigma <- 1
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

Y <-   nig_noise
}

ngme_control <- ngme.control(
  estimation = TRUE,
  burnin = 2,
  iterations = 1000,
  gibbs_sample = 10,
  stepsize = 1,
  kill_var = FALSE,
  threshold = 1e-4,
  opt_beta = T
)

ngme_out <- ngme(
  Y ~ 0 ,
  data = data.frame(Y = Y),
  control = ngme_control,
  noise = attr(nig_noise, "noise"),
  debug = ngme.debug(
    debug = TRUE,
    not_run = FALSE
  ),
  #seed = 2
  # , last_fit = ngme_out
)

str(ngme_out$est_output)
# c(noise_theta_mu, noise_theta_sigma, noise_theta_V)


plot(attr(nig_noise, "noise")) # measurement noise
sim = 10000
#V <- ngme2::rig(sim, eta <-attr(nig_noise,"noise")$theta_V, eta <-attr(nig_noise,"noise")$theta_V, seed =2)
#e <- attr(nig_noise,"noise")$theta_mu * (-1 + V) + exp(attr(nig_noise,"noise")$theta_sigma) * sqrt(V) * rnorm(sim)
#hist(e,add=T,probability = T,breaks=2000)
plot(create.ngme.noise(ngme_out$est_output$noise), add = TRUE, col="blue")

