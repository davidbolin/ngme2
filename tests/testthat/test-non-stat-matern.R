
test_that("matern offset_K triggers non-stationary operator", {
  mesh <- fmesher::fm_mesh_1d(seq(0, 1, length.out = 5))
  varying_offset <- seq(0, 0.4, length.out = mesh$n)

  model_ns <- matern(
    mesh = mesh,
    theta_K = 0,
    offset_K = varying_offset
  )

  model_stationary <- matern(
    mesh = mesh,
    theta_K = 0,
    offset_K = 0
  )

  expect_identical(model_ns$generic_type, "generic_ns")
  expect_equal(model_ns$offset_K$theta_K, as.numeric(varying_offset))

  diff_norm <- max(abs(as.matrix(model_ns$K) - as.matrix(model_stationary$K)))
  expect_gt(diff_norm, 0)
})

test_that("test Matern", {
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh_2k <- fmesher::fm_mesh_2d(
    loc.domain = pl01, cutoff = 0.1,
    max.edge = c(0.3, 10)
  )
  mesh_2k$n
  mesh_8k <- fmesher::fm_mesh_2d(
    loc.domain = pl01, cutoff = 0.08,
    max.edge = c(0.1, 10)
  )
  mesh_8k$n
  mesh_17k <- fmesher::fm_mesh_2d(
    loc.domain = pl01, cutoff = 0.01,
    max.edge = c(0.1, 10)
  )
  mesh_17k$n

  mesh <- mesh_2k
  # plot(mesh)
  mesh$n

  Basis <- cbind(1, (mesh$loc[,1] - 5) / 10)

  set.seed(123)

  n_obs <- 1000
  loc <- cbind(runif(n_obs, 0, 10), runif(n_obs, 0, 5))
  library(ngme2)
  true_noise = noise_nig(mu=-2, sigma=1, nu=0.5)
  true_noise = noise_normal(
    B_sigma = Basis,
    theta_sigma = c(1.5, 1.0)
    # sigma = 4
  )

  true_model <- ngme2::f(
    map = loc,
    model="matern",
    theta_K = c(0.2, 0.5),
    B_theta_K = Basis,
    mesh = mesh,
    noise = true_noise
  )

  # W <- rnorm(n_obs)
  W <- simulate(true_model)[[1]]
  attr(W, "noise")
  Y <- W + rnorm(n_obs, sd=0.2)
  var(W)

  # Matern case
  out <- ngme(
    Y ~ 1 + f(loc,
      model="matern",
      name="spde",
      mesh = mesh,
      B_theta_K = Basis,
      theta_K = c(0, 0),
      noise=noise_normal(
        B_sigma = Basis,
        theta_sigma = c(0, 0)
      ),
    ),
    data = data.frame(Y = Y),
    control_ngme = control_ngme(
      # use_iterative_solver = TRUE
    ),
    control_opt = control_opt(
      iterations = 5000,
      optimizer = precond_sgd(
        preconditioner = "full"
      ),
      rao_blackwellization = TRUE,
      n_parallel_chain = 4,
      print_check_info = F,
      verbose = T,
    )
  )
  out
  # compute_log_like(out)

  # load_all()
  # traceplot(out, "spde", hline=c(0.5, 0.5, 4))
  traceplot(out, "spde", hline=c(0.2, 0.5, 1.5, 1.0))
  traceplot(out, hline=c(0.2, 0))
  plot(true_noise,
    out$replicates[[1]]$models[[1]]$noise)

  # Now let's do some prediction
  # new_xs <- c(3, 5, 7)
  # new_ys <- c(3, 5, 3)
  # coo <- matrix(c(new_xs, new_ys), ncol=2)
  # predict(out, loc=list(coo))
  # plot(mesh)
  # points(x=new_xs, y = new_ys, type = "p", col="red", pch=16, cex=2)
  expect_true(TRUE)

skip_if_not_installed("INLA")
skip_if_not_installed("rSPDE")

library(INLA)
library(rSPDE)

## 1. 保持与 ngme 例子一致的网格、协变量与观测
mesh    <- mesh_2k                                   # 直接复用已有的 mesh
# B_kappa <- cbind(0, 0, 1, (mesh$loc[, 1] - 5) / 10)
# B_tau <- cbind(rep(0, mesh$n), 1, 0, 0)              
B_kappa <- cbind(0, 1, (mesh$loc[, 1] - 5) / 10, 0, 0)
B_tau <- cbind(0, 0, 0, 1, (mesh$loc[, 1] - 5) / 10)
n_obs   <- 1000
# loc     <- cbind(runif(n_obs, 0, 10), runif(n_obs, 0, 5))

## 2. 构造非平稳 Matern SPDE（log κ 含偏移列）
spde <- inla.spde2.matern(
  mesh  = mesh,
  alpha = 2,
  B.kappa = B_kappa,              # 和 ngme 一致：包含 offset + Basis
  B.tau   = B_tau  
  # prior.range = c(1, 0.5),           # 可按需要调整先验
  # prior.sigma = c(1, 0.5)
)

## 3. 建立观测矩阵 & SPDE 索引
A     <- inla.spde.make.A(mesh = mesh, loc = loc)

## 4. 构建数据 stack
stack <- inla.stack(
  data = list(Y = Y),
  A    = list(A),
  effects = list(
    i = 1:spde$n.spde
  )
)

## 5. 设置测量噪声（正态）以及公式
formula <- Y ~ 1 + f(i, model = spde)

result <- inla(
  formula,
  data    = inla.stack.data(stack),
  family  = "gaussian",
  control.predictor = list(A = inla.stack.A(stack), compute = TRUE),
  verbose = TRUE
)

## 6. 结果查看：非平稳 κ 的后验
summary(result)
out

# 将 κ 的后验平均映射回空间
kappa_post <- spde$fitted.values(result, what = "kappa")
plot(mesh, main = "Posterior mean of kappa (INLA)")
points(loc, pch = 16, cex = 0.3)
inla.mesh.projector(mesh)$plot.field(kappa_post$mean)

rspde_model <- rspde.matern(
  mesh = mesh,
  # alpha = 2,
  nu = 1.5,
  parameterization = "spde",
  B.tau = B_tau,
  B.kappa = B_kappa
)

rspde_stack <- inla.stack(
  data = list(Y = Y),
  A = list(rspde.make.A(mesh = mesh, loc = loc)),
  effects = list(i = 1:rspde_model$n.spde)
)

rspde_result <- inla(
  Y ~ 1 + f(i, model = rspde_model),
  data = inla.stack.data(rspde_stack),
  family = "gaussian",
  control.predictor = list(A = inla.stack.A(rspde_stack), compute = TRUE),
  verbose = TRUE
)

expect_true("(Intercept)" %in% rownames(result$summary.fixed))
expect_true("(Intercept)" %in% rownames(rspde_result$summary.fixed))
expect_equal(
  rspde_result$summary.fixed["(Intercept)", "mean"],
  result$summary.fixed["(Intercept)", "mean"],
  tolerance = 1e-3
)

common_hyper <- intersect(
  rownames(result$summary.hyperpar),
  rownames(rspde_result$summary.hyperpar)
)
expect_true(length(common_hyper) > 0)
for (name in common_hyper) {
  expect_equal(
    rspde_result$summary.hyperpar[name, "mean"],
    result$summary.hyperpar[name, "mean"],
    tolerance = 1e-3
  )
}


spde <- inla.spde2.matern(
  mesh  = mesh,
  alpha = 2,
  B.kappa = B_kappa,              # 和 ngme 一致：包含 offset + Basis
  B.tau   = B_tau  
  # prior.range = c(1, 0.5),           # 可按需要调整先验
  # prior.sigma = c(1, 0.5)
)

})

# test_that("test matern Non-stationary", {
#   # library(devtools); library(INLA); load_all()
#   { # First we create mesh
#     library(INLA)
#     pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
#     mesh <- fmesher::fm_mesh_2d(
#       loc.domain = pl01,
#       max.edge = 0.5
#     )
#     # plot(mesh)

#   n_obs <- 500; index_obs <- sample(1:mesh$n, n_obs)
#   loc_obs <- mesh$loc[index_obs, c(1, 2)]
#   A <- inla.spde.make.A(mesh = mesh, loc = loc_obs)
#   sigma.e <- 0.01

#   B_kappa <- matrix(c(rep(1, mesh$n), sin(1:mesh$n / 20)), ncol = 2)
#   theta_kappa <- c(1.5, 2.5)
#   W <- simulate(model_matern(loc_obs,
#       mesh = mesh,
#       B_kappa = B_kappa,
#       theta_kappa = theta_kappa,
#       # theta_kappa = 2,
#       noise = noise_nig(mu=-3, sigma = 2, nu = 1.5)))
#   # mean(attr(W, "noise")$h)
#   Y_obs <- as.numeric(A %*% W) + sigma.e * rnorm(n_obs)
#   # points(loc_obs[,1], loc_obs[,2], col="red", pch=16, cex=2)
#   }

#   # plot(B_kappa %*% theta_kappa)

#   #  { # 1d matern case
#   #   sigma.e = 0.01
#   #   loc_obs <- 1:500
#   #   B_kappa <- matrix(c(rep(1, 500), sin(1:500 / 20)), ncol = 2)
#   #   mesh = inla.mesh.1d(loc=1:500)
#   #   W <- simulate(
#   #       model_matern(loc_obs,
#   #       mesh = inla.mesh.1d(loc=loc_obs),
#   #       theta_kappa = 2,
#   #       # B_kappa = B_kappa,
#   #       # theta_kappa = c(1.1, 2),
#   #       noise = noise_nig(mu=-5, sigma = 2, nu = 1.5)))
#   #   A <- inla.spde.make.A(loc=loc_obs, mesh = inla.mesh.1d(loc=1:500))
#   #   Y <- as.numeric(A %*% W) + sigma.e * rnorm(n_obs)
#   # }
#   plot(Y, col="red"); abline(h=0)

#   trueV = attr(W,"noise")$V

#   ngme_out <- ngme(
#     Y ~ 0 + f(
#       loc_obs,
#       model="matern",
#       mesh = mesh,
#       theta_kappa = c(0.5, 0.5),
#       B_kappa = B_kappa,
#       # theta_kappa = 2,
#       # fix_theta_K = T,
#       # W = as.numeric(W), fix_W = TRUE,
#       noise = noise_nig(
#         # V = trueV, fix_V = TRUE
#         # fix_theta_nu = T, fix_theta_sigma=T, fix_theta_mu=T
#       ),
#       control = control_f(numer_grad = F),
#       debug = F,
#     ),
#     data = data.frame(Y = Y),
#     family = noise_normal(),
#     control_opt = control_opt(
#       estimation = T,
#       iterations = 100,
#       n_parallel_chain = 4
#       # max_relative_step = 5,
#       # max_absolute_step = 10
#     ),
#     debug = FALSE
#   )
#   ngme_out
#   traceplot(ngme_out, "field1")
#   ngme_out$replicates[[1]]$models[[1]]$theta_K
#   plot(B_kappa %*% theta_kappa, type="l")
#   points(B_kappa %*% ngme_out$replicates[[1]]$models[[1]]$theta_K)

# # compare noise
#   plot(ngme_out$replicates[[1]]$models[[1]]$noise,  noise_nig(mu=-3, sigma = 2, nu = 1.5))

#   expect_true(TRUE)
# })

##################################################
test_that("test 1d matern with numerical g", {
  n_mesh <- 800
  mesh <- fmesher::fm_mesh_1d(runif(n_mesh) * 800)

  n_obs <- 500
  loc <- runif(n_obs, 1, n_mesh)
  real_noise <- noise_nig(mu = -5, sigma=3, nu=2)
  spde1d <- f(map = loc, model="matern", mesh=mesh, theta_K=1,
    noise = real_noise)

# simulate
  W <- simulate(spde1d)[[1]]
  # eps <- simulate.ngme_noise(real_noise)
  # K <- 4 * spde1d$C + spde1d$G
  # W <- as.numeric(solve(K, eps))
  Y <- W + rnorm(n_obs, sd=0.5)
# plot(Y, type="l")

  out <- ngme(
    Y ~ 0 + f(loc, model="matern", mesh=mesh, noise=noise_nig()),
    data = data.frame(Y = Y),
    control_opt = control_opt(
      estimation = T,
      iterations = 500,
      n_parallel_chain = 4
    )
    # ,start = out
  )
  out
traceplot(out, "field1")
plot(out$replicates[[1]]$models[[1]]$noise, real_noise)

  expect_true(TRUE)
})
