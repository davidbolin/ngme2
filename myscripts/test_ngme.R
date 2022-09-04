a2th <- function(k) {log((-1-k)/(-1+k))}
th2a <- function(th) {-1 + (2*exp(th)) / (1+exp(th))}

library(devtools); load_all()

############  0. generating fix effects and control
seed <- 51
set.seed(seed)

# fix effects
# X <- (model.matrix(Y1 ~ x1 + x2))  # design matrix
# Y1 = as.numeric(Y1 + X %*% beta)

############  1. test AR with nig noise
n_obs <- 500
ar_mu <- 4
ar_sigma <- 1.3
ar_eta <- 0.8

ar1 <- ngme.simulate(
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

nig_simulation <- ngme.simulate.noise(
  ngme.noise.nig(
    theta_mu = noise_theta_mu,
    theta_sigma = noise_theta_sigma,
    theta_V = noise_theta_V,

    fix_sigma = TRUE,
    fix_var = TRUE,
    fix_V = TRUE
  ),
  n = n_obs
)
Y <- ar1$realization + nig_simulation$realization

# Y <- ar1$realization + rnorm(n = n_obs, sd = 1.2)

# range(nig_simulation$realization)
# str(ar1$noise)

################################################################
# index <- c(1:n_obs1, 1:n_obs)
# replicates = c(rep(1, n_obs1), rep(2, n_obs)
# print(Y)

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
  f(1:n_obs,
    replicates = NULL,
    model = "ar1",
    W = ar1$realization,
    noise = ngme.noise(
      theta_mu = ar_mu,
      theta_sigma = ar_sigma,
      theta_V = ar_eta,
      V = ar1$noise$V,

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
  data = data.frame(),
  control = ngme_control,
  noise = nig_simulation$noise,
  debug = ngme.debug(
    debug = TRUE,
    not_run = F
  ),
  seed = 1
  # , last_fit = ngme_out
)

str(ngme_out$est_output)
c(noise_theta_mu, noise_theta_sigma, noise_theta_V)

# trace plot
plot_out(ngme_out$opt_trajectory, start = 5, n = 3)


# ngme_out2 <- ngme(
#   Y ~ 0 +
#   f(1:n_obs,
#     replicates = NULL,
#     model = "ar1",
#     theta_K = 0.5,
#     W = ar1$realization,
#     noise = ngme.noise(
#       theta_mu = 2,
#       theta_sigma = 0,
#       theta_V = 2,
#       V = ar1$noise$V
#     ),
#     control = ngme.control.f(
#       numer_grad       = FALSE,
#       use_precond      = TRUE,

#       fix_operator     = FALSE,
#       fix_mu           = FALSE,
#       fix_sigma        = FALSE,
#       fix_var          = FALSE,
#       fix_V            = FALSE,
#       fix_W            = TRUE
#     ),
#     debug = FALSE
#   ),
#   data = data.frame(),
#   control = ngme.control(
#     estimation = FALSE,
#     fix_V = FALSE,

#     burnin = 200,
#     iterations = 5,
#     gibbs_sample = 5,
#     stepsize = 1,
#     kill_var = FALSE,
#     threshold = 1e-4,
#     opt_beta = T
#   ),
#   noise = nig_simulation$noise,
#   debug = ngme.debug(
#     debug = TRUE,
#     not_run = FALSE
#   ),
#   seed = 1,
#   last_fit = ngme_out
# )

# str(ngme_out2$output)


# ?dnorm
