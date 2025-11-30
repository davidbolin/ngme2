seed <- 500
set.seed(seed)

n_obs <- 800
day <- 1:n_obs
sigma_eps <- 0.5
mu <- 2
delta <- -mu
sigma <- 3
nu <- 1
rho <- 0.5

ar1_model <- f(day,
  model = "ar1", rho = 0.5,
  noise = noise_nig(mu = mu, sigma = sigma, nu = nu)
)

# ar1_model <- f(day, model="ar1", rho = 0.5,
#   noise = noise_normal())

W <- simulate(ar1_model, seed = seed, nsim = 1)[[1]]
Y <- W + rnorm(n_obs, sd = 1)
Y <- 1:n_obs

# First we generate V. V_i follows inverse Gaussian distribution
trueV <- ngme2::rig(n_obs, nu, nu, seed = 10)

# Then generate the nig noise
mynoise <- delta + mu * trueV + sigma * sqrt(trueV) * rnorm(n_obs)
trueW <- Reduce(function(x, y) {
  y + rho * x
}, mynoise, accumulate = T)
Y <- trueW + rnorm(n_obs, mean = 0, sd = sigma_eps)

# fit the Gaussian model
ret_gauss <- ngme(
  Y ~ 0 + f(
    1:n_obs,
    name = "my_ar",
    model = "ar1",
    rho = rho,
    noise = noise_normal(),
    control = control_f(numer_grad = TRUE, improve_hessian = FALSE)
  ),
  data = data.frame(Y = Y),
  control_opt = control_opt(
    burnin = 100,
    iterations = 1000,
    n_parallel_chain = 4,
    verbose = FALSE,
    seed = seed,
    rao_blackwellization = TRUE,
    print_check_info = TRUE,
    std_lim = 0.01,
    trend_lim = 0.01,
    n_min_batch = 10
  )
)
ret_gauss
traceplot(ret_gauss, "my_ar")
# saveRDS(ret_gauss, "examples/models/ret_gauss.rds")
# ret_gauss_real <- readRDS("examples/models/ret_gauss.rds")
traceplot(ret_gauss, "my_ar")

# fit the non-Gaussian model
ret_nig <- ngme(
  Y ~ 0 + f(
    1:n_obs,
    rho = 0.5,
    name = "my_ar",
    model = "ar1",
    noise = noise_nig(mu = mu, sigma = sigma, nu = nu),
    control = control_f(numer_grad = TRUE, improve_hessian = FALSE)
  ),
  data = data.frame(Y = Y),
  # family = noise_normal(sigma=2),
  control_opt = control_opt(
    burnin = 100,
    # estimation = FALSE,
    iterations = 1000,
    n_parallel_chain = 4,
    verbose = FALSE,
    rao_blackwellization = TRUE,
    seed = seed,
    print_check_info = TRUE,
    std_lim = 0.01,
    trend_lim = 0.01,
    n_min_batch = 3
  )
  # , start = ret_nig
)
ret_nig

traceplot(ret_nig, "my_ar", hline = c(0.5, 2, 3, 1))
saveRDS(ret_nig, "examples/models/ret_nig.rds")
ret_nig <- readRDS("examples/models/ret_nig.rds")
ret_nig_real <- readRDS("examples/models/ret_nig.rds")
traceplot(ret_nig, "my_ar", hline = c(0.5, mu, sigma, nu))

load_all()
seed <- 500

tmp <- make_time_series_cv_index(1:n_obs, train_length = 10, test_length = 1)
train_idx <- tmp$train
test_idx <- tmp$test

cv <- cross_validation(
  list(
    gauss = ret_gauss,
    nig = ret_nig
  ),
  type = "custom",
  # k=10,
  parallel = TRUE,
  n_burnin = 100,
  n_gibbs_samples = 100,
  train_idx = train_idx,
  test_idx = test_idx,
  print = TRUE,
  N_sim = 1,
  seed = seed + 50
)
cv

plot(y_data, type = "l")
lines(AW_N_2[, 100], col = "red")
lines(AW_N_1[, 100], col = "blue")

lines(pred_N_1[, 1], col = "red")
lines(pred_N_2[, 1], col = "blue")
lines(pred, col = "green")


sum(ngme_1rep$models[[i]]$A)
cv
cv_500
cv_1000

ret_nig
ret_nig_2
cv


######## create loo prediction ########
seed <- 1234
idx <- 1:100
# idx <- sample(1:100)
# idx <- sample(1:n_obs, 81)

pred_gauss <- double(length = length(idx))
pred_nig <- double(length = length(idx))
load_all()
for (i in seq_along(idx)) {
  print(i)
  ret_gauss_loo <- ret_gauss$replicates[[1]]
  # ret_gauss_loo$models[[1]]$A[idx[i], idx[i]] <- 0
  ret_nig_loo <- ret_nig$replicates[[1]]
  # ret_nig_loo$models[[1]]$A[idx[i], idx[i]] <- 0

  n_samples <- 100
  W_gauss <- sampling_cpp(
    ret_gauss_loo,
    n = n_samples,
    n_burnin = 100,
    posterior = TRUE,
    seed = seed + i
  )
  W_nig <- sampling_cpp(
    ret_nig_loo,
    n = n_samples,
    n_burnin = 100,
    posterior = TRUE,
    seed = seed + i
  )
  pred_gauss[i] <- mean(sapply(W_gauss$W, function(x) x[idx[i]]))
  pred_nig[i] <- mean(sapply(W_nig$W, function(x) x[idx[i]]))
}

{
  plot(Y[idx], type = "l")
  lines(pred_gauss, col = "red")
  lines(pred_nig, col = "blue")
}

{
  mae_gauss <- mean(abs(Y[idx] - pred_gauss))
  mae_nig <- mean(abs(Y[idx] - pred_nig))
  mse_gauss <- mean((Y[idx] - pred_gauss)^2)
  mse_nig <- mean((Y[idx] - pred_nig)^2)

  table <- data.frame(
    mae = c(mae_gauss, mae_nig),
    mse = c(mse_gauss, mse_nig),
    model = c("gauss", "nig")
  )
  table
}




load_all()
seed <- 50
cv <- cross_validation(
  list(
    gauss = ret_gauss,
    nig = ret_nig
  ),
  k = 20,
  parallel = FALSE,
  n_burnin = 100,
  # n_gibbs_samples=10000,
  print = TRUE,
  N_sim = 1,
  seed = seed + 50
)
cv
dim(pred_N_1)
plot(y_data, type = "l", lty = 2, col = "blue")
for (i in 1:100) {
  print(i)
  lines(pred_N_1[, i], type = "l")
}
lines(y_data, type = "l", col = "blue")
pred <- rowMeans(as.matrix(pred_N_1))
lines(pred, type = "l", col = "red")
