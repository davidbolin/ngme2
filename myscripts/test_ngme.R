a2th <- function(k) {log((-1-k)/(-1+k))}
th2a <- function(th) {-1 + (2*exp(th)) / (1+exp(th))}

library(devtools); load_all()

############  0. generating fix effects and control
seed <- 8
set.seed(seed)

# fix effects
# X <- (model.matrix(Y1 ~ x1 + x2))  # design matrix
# Y1 = as.numeric(Y1 + X %*% beta)

############  1. test AR with nig noise
n_obs <- 500
ar1 <- ngme.simulate(
  f(1:n_obs, model = "ar1", theta_K = 0.4),
  noise = ngme.noise(
    theta_mu = 2,
    theta_sigma = 0,
    theta_V = 2
  ),
  seed = NULL
)

nig_simulation <- ngme.simulate.noise(
  ngme.noise.nig(
    theta_mu = 2,
    theta_sigma = 1.5,
    theta_V = 0.8
  ),
  n = n_obs
)

str(ar1$noise)
Y <- ar1$realization + nig_simulation$realization

################################################################
# index <- c(1:n_obs1, 1:n_obs)
# replicates = c(rep(1, n_obs1), rep(2, n_obs))
# print(Y)

ngme_control <- ngme.control(
  estimation = FALSE,
  fix_V = TRUE,

  burnin = 200,
  iterations = 5,
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
    theta_K = 0.5,
    noise = ngme.noise(
      theta_mu = 2,
      theta_sigma = 0,
      theta_V = 2,
      V = ar1$noise$V
    ),
    control = ngme.control.f(
      numer_grad       = FALSE,
      use_precond      = TRUE,

      fix_operator     = FALSE,
      fix_mu           = FALSE,
      fix_sigma        = FALSE,
      fix_var          = FALSE,
      fix_V            = TRUE,
      fix_W            = FALSE,

      init_W           = NULL
    ),
    debug = FALSE
  ),
  data = data.frame(),
  control = ngme_control,
  noise = nig_simulation$noise,
  debug = ngme.debug(
    debug = TRUE
  ),
  seed = 1
)

head(nig_simulation$noise$V)
str(ngme_out)
# test grad. of V in block model