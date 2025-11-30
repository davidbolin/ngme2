# library(ngme2)
load_all()

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

# # First we generate V. V_i follows inverse Gaussian distribution
# trueV <- ngme2::rig(n_obs, nu, nu, seed = 10)

# # Then generate the nig noise
# mynoise <- delta + mu*trueV + sigma * sqrt(trueV) * rnorm(n_obs)
# trueW <- Reduce(function(x,y){y + rho*x}, mynoise, accumulate = T)
# Y = trueW + rnorm(n_obs, mean=0, sd=sigma_eps)
# # plot(Y, type="l")

# fit the non-Gaussian model
control_same <- control_opt(
  optimizer = precond_sgd(),
  burnin = 100,
  iterations = 100,
  n_parallel_chain = 4,
  # verbose = TRUE,
  rao_blackwellization = TRUE,
  seed = seed,
  print_check_info = TRUE,
  std_lim = 0.01,
  trend_lim = 0.01,
  n_slope_check = 3,
  n_trace_iter = 30,
  start_sd = 0.0,
  max_R_hat = 1.3
)

load_all()
system.time({
  ret_nig_same <- ngme(
    Y ~ 0 + f(
      1:n_obs,
      name = "my_ar",
      model = ar1(),
      noise = noise_nig()
    ),
    family = noise_normal(),
    data = data.frame(Y = Y),
    control_opt = control_same
  )
})

ret_nig_same
traceplot(ret_nig_same, "my_ar", hline = c(rho, mu, sigma, nu))
traceplot(ret_nig_same, hline = c(sigma_eps))
