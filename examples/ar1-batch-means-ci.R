if (requireNamespace("devtools", quietly = TRUE)) {
  devtools::load_all(".")
} else {
  library(ngme2)
}

# Example: Xi-style batch-means inference for AR1 with polynomial stepsize
set.seed(2026)
seed <- 2026

n_obs <- 500
idx <- 1:n_obs

# Ground truth used for simulation
rho_true <- 0.55
sigma_latent_true <- 3
sigma_eps_true <- 0.4

# Simulate latent AR1 process and observations
ar1_model <- f(
  idx,
  model = ar1(rho = rho_true),
  noise = noise_normal(sigma = sigma_latent_true)
)

latent <- simulate(ar1_model, seed = seed, nsim = 1)[[1]]
y <- latent + rnorm(n_obs, sd = sigma_eps_true)
dat <- data.frame(y = y)

# CI-focused control helper:
# - fixed number of iterations
# - trajectory storage enabled
# - early-stop diagnostics disabled
# - polynomial stepsize schedule
ctl <- control_opt_batch_ci(
  optimizer = sgd(stepsize = 0.02),
  burnin = 300,
  iterations = 10000,
  n_batch = 10,
  n_parallel_chain = 1,
  alpha = 0.51,
  t0 = 100,
  # schedule_burnin_iter = 500,
  start_sd = 0.2,
  seed = seed,
  verbose = TRUE
)

fit <- ngme(
  y ~ 0 + f(idx, name = "ar1_field", model = ar1(), noise = noise_normal()),
  data = dat,
  family = noise_normal(),
  control_opt = ctl,
  control_ngme = control_ngme(n_gibbs_samples = 5)
)

traceplot(fit, hline = c(rho_true, sigma_latent_true, sigma_eps_true))

# Xi-style batch-means CI from all parameters jointly
ci <- ngme_batch_ci(
  fit,
  name = "all",
  level = 0.95,
  alpha = 0.51,
  # burnin_iter = 1000,
  apply_transform = TRUE
)

res <- data.frame(
  parameter = names(ci$estimates),
  estimate = as.numeric(ci$estimates),
  lower = as.numeric(ci$conf_int[, "lower"]),
  upper = as.numeric(ci$conf_int[, "upper"]),
  row.names = NULL,
  check.names = FALSE
)

truth <- c(
  "rho (ar1_field)" = rho_true,
  "sigma (ar1_field)" = sigma_latent_true,
  "sigma" = sigma_eps_true
)
res$truth <- unname(truth[match(res$parameter, names(truth))])
res$covered <- ifelse(
  is.na(res$truth),
  NA,
  res$lower <= res$truth & res$truth <= res$upper
)

cat("\nBatch-means estimates and 95% CI:\n")
print(res, row.names = FALSE)

cat("\nSelected settings:\n")
cat("  n_chains =", ci$n_chains, "\n")
cat("  n_iterations =", ci$n_iterations, "\n")
cat("  n_eff =", ci$n_eff, "\n")
cat("  alpha =", ci$alpha, "\n")
cat("  M =", ci$M, "\n")






# new approach
fit <- ngme(
  y ~ 0 + f(idx, name = "ar1_field", model = ar1(), noise = noise_normal()),
  data = dat,
  control_opt = control_opt(
    iterations = 1000
  ),
  family = noise_normal(),
  control_ngme = control_ngme(n_gibbs_samples = 5)
)
traceplot(fit)

CIs <- compute_ngme_CI(
  fit = fit,
  iterations = 8000,
  alpha = 0.51,
  optimizer = sgd(stepsize = 0.02),
  burnin = 300,
  n_batch = 10,
  n_parallel_chain = 1,
  t0 = 100,
  burnin_iter = 1000, # 同时用于 schedule warmup + CI 丢弃前段轨迹
  M = 10,
  verbose = TRUE,
  apply_transform = TRUE
)
traceplot(CIs$refit)
summary(CIs)
