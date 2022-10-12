# Simple scripts for test ngme function
# sth. wrong with nig measurement noise
# library(ngme2)
library(devtools);
load_all()
{
a2th <- function(k) {log((-1-k)/(-1+k))}
th2a <- function(th) {-1 + (2*exp(th)) / (1+exp(th))}

seed <- 5
set.seed(seed)
}

{ ############  1. simulate AR with nig noise
n_obs <- 300
ar_mu <- 4
ar_sigma <- 1.3
ar_eta <- 0.8

ar1_process <- simulate(
  f(1:n_obs,
    model = "ar1",
    theta_K = 0.2,
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

    fix_theta_sigma = FALSE,
    fix_theta_V     = FALSE,
    fix_V           = FALSE
  ),
  nsim = n_obs
)

Y <- ar1_process + nig_noise
Y
# Y <- ar1_process + rnorm(n_obs)
}

ngme_out <- ngme(
  Y ~ 0 +
  f(model = "ar1",
    theta_K = 0.7,
    W = as.numeric(ar1_process),
    fix_W = FALSE,
    noise = ngme.noise.nig(
      theta_mu = ar_mu,
      theta_sigma = ar_sigma,
      theta_V = ar_eta,
      V = attr(ar1_process, "noise")$V,

      fix_theta_mu      = FALSE,
      fix_theta_sigma   = FALSE,
      fix_theta_V       = FALSE,
      fix_V             = FALSE
    ),
    control = ngme.control.f(
      numer_grad       = FALSE,
      use_precond      = TRUE
    ),
    debug = FALSE
  ),
  data = data.frame(Y = Y),
  control = ngme.control(
    estimation = TRUE,
    n_parallel_chain = 2,
    stop_points = 2,
    burnin = 200,
    iterations = 1000,
    gibbs_sample = 5,
    stepsize = 1,
    kill_var = FALSE,
    threshold = 1e-4
  ),
  # noise = ngme.noise.normal(
  #   fix_theta_sigma = FALSE
  # ),
  noise = attr(nig_noise, "noise"),
  seed = 2,
  # , last_fit = ngme_out
  debug = FALSE
)

ngme_out
# str(ngme_out)

plot_chains(ngme_out, parameter = "theta_mu", f_index = 0)
plot_chains(ngme_out, parameter = "theta_V", f_index = 0)
plot_chains(ngme_out, parameter = "theta_mu", f_index = 1)
plot_chains(ngme_out, parameter = "theta_K", f_index = 1)

# str(ngme_out)

# trajs <- ngme_out

# trajs
# outputs <- ngme_out$outputs
# outputs[[1]]$latents[[1]]
# attr(outputs[[1]], "trajectory")
# attr(outputs[[1]]$latent[[1]], "trajectory")

# str(ngme_out[[1]] + ngme_out[[2]])
# ngme_out[[1]]
# c(noise_theta_mu, noise_theta_sigma, noise_theta_V)

# # trace plot of mu
# opt_trajectory <- attr(ngme_out, "opt_trajectory")

# traceplot(ngme_out[[1]], start = 1, n = 1)
# traceplot(ngme_out)
# ngme.traceplot(opt_trajectory, start = 3, n = 1, transform = identity)



# posterior sampling
# str(sampling_cpp(ngme_out, 10))

# ngme.sampling(ngme, posterior)
# f(X|I, model="nig", data=list(X=...,I=...))

# str(ngme_out)

# outputs <- ngme_out
# outputs[[1]]

# str(outputs[[1]])

# unlist(attr(outputs[[1]], "trajectory")$beta)
# unlist(attr(outputs[[1]], "trajectory")$theta_mu)
# unlist(attr(outputs[[1]], "trajectory")$theta_sigma)
# unlist(attr(outputs[[1]], "trajectory")$theta_V_traj)

# outputs[[1]]$beta


# lll <- list()
# lll[[1]] = list(a = 1, b=3)
# lll[[2]] = list(a = 4, b=5)
# length(lll)
# unlist((Map(function(i) {lll[[i]]$a}, 1:2)))


# beta <- mean(unlist((Map(function(i) {outputs[[i]]$beta}, 1:2))))
# theta_mu <- mean(unlist((Map(function(i) {outputs[[i]]$theta_mu}, 1:2))))

# outputs[]
# theta_mu <- mean(unlist((Map(function(i) {outputs[[i]]$noise$theta_mu}, 1:2))))

# (outputs[[1]]$noise)
# (outputs[[1]]$latents)


# for (l in lll) {
#   for (j in )
# }

# outputs
#