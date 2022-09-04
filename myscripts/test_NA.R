library(devtools); load_all()

set.seed(7)
set.seed(Sys.time())

n_obs <- 1000
ar_mu <- 0
ar_sigma <- 1.3
ar_eta <- 0.8

ar1 <- ngme.simulate(
  f(1:n_obs,
    model = "ar1",
    theta_K = 0.6,
    noise = ngme.noise.nig(
      theta_mu = ar_mu,
      theta_sigma = ar_sigma,
      theta_V = ar_eta
    )
  ),
  seed = 123
)

# normal noise
Y <- ar1$realization + rnorm(n_obs)

# fix effects
# x1 <- rexp(n_obs); x2 <- rnorm(n_obs)
# beta <- c(-3, 1, 2)
# X <- (model.matrix(Y ~ x1 + x2))  # design matrix
# Y = as.numeric(Y + X %*% beta)


################################################################
# index <- c(1:n_obs1, 1:n_obs)
# replicates = c(rep(1, n_obs1), rep(2, n_obs)
# print(Y)

ngme_control <- ngme.control(
  estimation = TRUE,
  fix_V = FALSE,

  burnin = 200,
  iterations = 200,
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
    theta_K = 0.6,
    W = ar1$realization,
    noise = ngme.noise(
      theta_mu = ar_mu,
      theta_sigma = ar_sigma,
      theta_V = ar_eta,
      V = ar1$noise$V
    ),
    control = ngme.control.f(
      numer_grad       = FALSE,
      use_precond      = TRUE,

      fix_operator     = TRUE,
      fix_mu           = FALSE,
      fix_sigma        = TRUE,
      fix_var          = FALSE,
      fix_V            = FALSE,
      fix_W            = FALSE
    ),
    debug = TRUE
  ),
  data = data.frame(),
  control = ngme_control,
  noise = ngme.noise.normal(),
  debug = ngme.debug(
    debug = TRUE,
    not_run = F
  ),
  seed = 1
)

str(ngme_out$est_output)
# plot_out(ngme_out$opt_trajectory, start = 5)
plot_out(ngme_out$opt_trajectory, start = 2)


ngme.noise.nig(
  theta_mu = 0,
  theta_sigma = 0,
  theta_V = 1,
  V = NULL,
  B_mu = matrix(1),
  B_sigma = matrix(1),
  fix_mu = TRUE
)
