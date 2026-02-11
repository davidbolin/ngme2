if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all()
} else {
  library(ngme2)
}
seed <- 500
set.seed(seed)

n_obs <- 500
day <- 1:n_obs
sigma_eps <- 0.5
mu <- 2
delta <- -mu
sigma <- 3
nu <- 1
rho <- 0.5

ar1_model <- f(day,
  model = ar1(rho = rho),
  noise = noise_nig(mu = mu, sigma = sigma, nu = nu)
)
# noise = noise_normal(sigma = sigma))

# ar1_model <- f(day, model="ar1", rho = 0.5,
#   noise = noise_normal())

process <- simulate(ar1_model, seed = seed, nsim = 1)[[1]]
Y <- process + rnorm(n_obs, sd = sigma_eps)

# fit the non-Gaussian model with grad.norm plateau step decay
control_decay <- control_opt(
  optimizer = precond_sgd(stepsize = 1),
  burnin = 100,
  iterations = 1000,
  n_parallel_chain = 4,
  verbose = TRUE,
  rao_blackwellization = TRUE,
  seed = seed,
  print_check_info = TRUE,
  n_min_batch = 3,
  n_batch = 50,
  n_trace_iter = 30,
  start_sd = 0.0,
  trend_std_conv_check = FALSE,
  std_lim = 0.01,
  trend_lim = 0.01,
  R_hat_conv_check = FALSE,
  max_R_hat = 1.1,
  pflug_conv_check = FALSE,
  pflug_alpha = 1,
  stepsize_decay = stepsize_decay(
    method = "grad_norm_plateau",
    patience = 1,
    gamma = 0.8,
    min_delta = 1e-2,
    warmup = 1,
    min_stepsize = 1e-4
  )
)

system.time({
  ret_nig_decay <- ngme(
    Y ~ 0 + f(
      1:n_obs,
      name = "my_ar",
      model = ar1(),
      noise = noise_nig()
    ),
    family = noise_normal(),
    data = data.frame(Y = Y),
    control_opt = control_decay
  )
})

ret_nig_decay
traceplot(ret_nig_decay, "my_ar", hline = c(rho, mu, sigma, nu))
traceplot(ret_nig_decay)


cross_validation(
  ret_nig_decay,
  k = 5
)
