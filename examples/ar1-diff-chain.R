load_all()

seed <- 500
set.seed(seed)

n_obs <- 800
day <- 1:n_obs
sigma_eps <- 0.5
mu = 2; delta = -mu
sigma <- 3
nu <- 1
rho <- 0.5

ar1_model <- f(day, model="ar1", rho = 0.5,
  noise = noise_nig(mu = mu, sigma = sigma, nu=nu))

# ar1_model <- f(day, model="ar1", rho = 0.5,
#   noise = noise_normal())

W <- simulate(ar1_model, seed = seed, nsim=1)[[1]]
Y <- W + rnorm(n_obs, sd = 1)
Y <- 1:n_obs

# First we generate V. V_i follows inverse Gaussian distribution
trueV <- ngme2::rig(n_obs, nu, nu, seed = 10)

# Then generate the nig noise
mynoise <- delta + mu*trueV + sigma * sqrt(trueV) * rnorm(n_obs)
trueW <- Reduce(function(x,y){y + rho*x}, mynoise, accumulate = T)
Y = trueW + rnorm(n_obs, mean=0, sd=sigma_eps)
plot(Y, type="l")

# fit the non-Gaussian model
control_same <- control_opt(
  optimizer = precond_sgd(),
  burnin = 100,
  # estimation = FALSE,
  iterations = 100,
  n_parallel_chain = 4,
  verbose = FALSE,
  rao_blackwellization = TRUE,
  seed = seed,
  print_check_info = TRUE,
  std_lim = 0.01,
  trend_lim = 0.01,
  n_slope_check = 3
)

ret_nig_same <- ngme(
  Y ~ 0 + f(
    1:n_obs,
    rho = 0.5,
    name = "my_ar",
    model = "ar1",
    noise = noise_nig(),
  ),
  data = data.frame(Y=Y),
  control_opt = control_same
)

ret_nig_same
traceplot(ret_nig_same, "my_ar", hline=c(0.5, mu, sigma, nu))


# ------------------------------------------------------------
control_diff <- control_opt(
  optimizer = precond_sgd(
    precond_by_diff_chain = TRUE
  ),
  burnin = 100,
  # estimation = FALSE,
  iterations = 100,
  n_parallel_chain = 4,
  verbose = FALSE,
  rao_blackwellization = TRUE,
  seed = seed,
  print_check_info = TRUE,
  std_lim = 0.01,
  trend_lim = 0.01,
  n_slope_check = 3
)

ret_nig_diff <- ngme(
  Y ~ 0 + f(
    1:n_obs,
    rho = 0.5,
    name = "my_ar",
    model = "ar1",
    noise = noise_nig(),
  ),
  data = data.frame(Y=Y),
  control_opt = control_diff
)

ret_nig_diff
traceplot(ret_nig_diff, "my_ar", hline=c(0.5, mu, sigma, nu))
