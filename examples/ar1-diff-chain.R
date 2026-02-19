library(ngme2)
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
  iterations = 500,
  n_parallel_chain = 4,
  verbose = TRUE,
  rao_blackwellization = TRUE,
  seed = seed,
  print_check_info = TRUE,
  n_min_batch = 3,
  n_trace_iter = 30,
  start_sd = 0.0,
  trend_std_conv_check = TRUE,
  std_lim = 0.01,
  trend_lim = 0.01,
  R_hat_conv_check = TRUE,
  max_R_hat = 1.1,
  pflug_conv_check = TRUE,
  pflug_alpha = 1
)
t <- 1:n_obs
ret_nig_same <- ngme(
  Y ~ 1 + t + f(
    1:n_obs,
    name = "my_ar",
    model = ar1(),
    noise = noise_nig(
      prior = priors(
        nu = prior_normal(0, 0.5)
      )
    )
  ),
  family = noise_normal(),
  data = data.frame(Y = Y, t = t),
  control_opt = control_same
)
summary(ret_nig_same)
traceplot(ret_nig_same)


pred <- predict(
  ret_nig_same,
  map = list(my_ar = 1:100),
  data = data.frame(t = 1:100),
  estimator = c("mean", "sd", "0.05q", "0.95q", "median", "mode"),
  sampling_size = 500,
  burnin_size = 100,
  seed = seed,
  train_idx = NULL,
  # chain_combine = c("predictive_average")
)

plot(1:100, pred$mean, type = "l")
lines(1:100, pred$`0.05q`, col = "blue")
lines(1:100, pred$`0.95q`, col = "blue")
