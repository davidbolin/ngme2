library(rSPDE)

s <- seq(from = 0, to = 1, length.out = 1001)
mesh <- fmesher::fm_mesh_1d(s)
kappa <- 20
sigma <- 10
tau <- 1 / sigma
nu <- 0.8
alpha <- nu + 1 / 2
r <- sqrt(8 * nu) / kappa

# integer case
kappa <- 2
sigma <- 5
tau <- 1 / sigma
nu <- 1.5
alpha <- nu + 1 / 2
r <- sqrt(8 * nu) / kappa

# alpha=nu + d/2
C <- matern(mesh, kappa = kappa)$C
G <- matern(mesh, kappa = kappa)$G

# op <- matern.operators(
#   C = C,
#   G = G,
#   sigma = sigma,
#   range = r,
#   # nu = nu,
#   nu = 2 - 1/2,
#   # loc_mesh = s,
#   d = 1,
#   m = 1,
#   type = "operator",
#   parameterization = "matern"
# )
alpha
op <- matern.operators(
  C = C,
  G = G,
  tau = tau,
  kappa = kappa,
  alpha = alpha,
  d = 1,
  m = 2,
  type = "operator",
  parameterization = "spde"
)
# op$beta
# op$Pl
# op$Pr

w <- simulate(op, seed = 10)[, 1]
var(w)
y <- w + rnorm(length(w), sd = 0.1)
dat <- data.frame(y = y, s = s)

# 2) 用 rspde_lme 拟合（rspde_order 与 m 一致）
rspde_fit <- rspde_lme(
  formula = y ~ 0,
  loc = "s",
  data = dat,
  model = matern.operators(
    mesh = mesh,
    kappa = kappa, tau = tau, alpha = alpha,
    d = 1, m = 2,
    type = "operator",
    parameterization = "spde"
  ),
  rspde_order = 2,
  model_options = list(
    start_sigma_e = 0.1,
    fix_alpha = alpha
    # 可选：fix_alpha = alpha, fix_kappa = kappa, fix_tau = tau
  )
)
summary(rspde_fit)

# options(error = recover)
fit <- ngme(
  y ~ 0 + f(
    1:length(y),
    model = "matern",
    name = "spde",
    kappa = 20,
    mesh = mesh,
    alpha = alpha,
    fix_alpha = TRUE,
    rational_order = 2,
    noise = noise_normal(sigma = 5)
  ),
  data = data.frame(y = y),
  control_opt = control_opt(
    estimation = TRUE,
    iterations = 1000,
    optimizer = adam(),
    n_batch = 10,
    rao_blackwellization = F,
    n_parallel_chain = 4,
    print_check_info = F,
    verbose = T,
    std_lim = 0.01,
    start_sd = 0.01
  )
)
fit
traceplot(fit, "spde")
traceplot(fit, hline = 0.1)
