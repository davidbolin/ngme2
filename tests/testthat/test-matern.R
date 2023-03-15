# test spatial Matern model

# 1. test estimation of matern model
# 2. test predict(out, loc=new_loc)
test_that("test Matern", {
  load_all()
  library(INLA)
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh <- inla.mesh.2d(
    loc.domain = pl01, cutoff = 0.2,
    max.edge = c(0.5,10)
  )
# mesh$n
# plot(mesh)
  n_obs <- 500
  loc <- cbind(runif(n_obs, 0, 10), runif(n_obs, 0, 5))

  true_model <- model_matern(
    kappa = 3, mesh = mesh, map = loc, noise = noise_nig(
      mu=-2, sigma=1.5, nu=1, n = mesh$n
    )
  )
  W <- simulate(true_model)

  A <- inla.spde.make.A(loc=loc, mesh=mesh)
  Y <- as.numeric(A %*% W) + rnorm(n_obs, sd=0.5)

  # make bubble plot
  # sp_obj <- as.data.frame(mesh$loc); sp_obj[, 3] <- W
  # names(sp_obj) <- c("s1", "s2", "y")
  # coordinates(sp_obj) <- ~ s1 + s2
  # bubble(sp_obj, zcol=3)
  # range(mesh$loc[, 1]); range(mesh$loc[, 2])

# Matern case
load_all()
  out <- ngme(
    Y ~ 0 + f(loc,
      model="matern",
      name="spde",
      mesh = mesh,
      noise=noise_nig(),
      control = control_f(numer_grad = F),
      # fix_theta_K = T, theta_kappa = log(3),
      # fix_W = TRUE,
      # W = W,
      debug = F
    ),
    data = data.frame(Y = Y),
    control_opt = control_opt(
      estimation = T,
      iterations = 500,
      n_parallel_chain = 4,
      print_check_info = T,
      verbose = T
    ),
    debug = F
  )
  # out[[1]]$latents[[1]]$noise$h == out[[1]]$latents[[1]]$h
  out
  traceplot(out)
  traceplot(out, "spde")
  plot(noise_nig(mu=-2,sigma=1.5,nu=1),
    out[[1]]$latents[[1]]$noise)

  # Now let's do some prediction
  new_xs <- c(3, 5, 7)
  new_ys <- c(3, 5, 3)
  coo <- matrix(c(new_xs, new_ys), ncol=2)
  predict(out, loc=coo)
  # plot(mesh)
  # points(x=new_xs, y = new_ys, type = "p", col="red", pch=16, cex=2)
  expect_true(TRUE)
})

# test more
test_that("test matern ns", {
  # library(devtools); library(INLA); load_all()
  { # First we create mesh
    pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
    mesh <- inla.mesh.2d(
      loc.domain = pl01,
      max.edge = c(1)
    )
    B_kappa <- matrix(c(rep(1, mesh$n), rexp(mesh$n)), ncol = 2)
    W <- simulate(
      f(model = model_matern(
        mesh = mesh,
        B_kappa = B_kappa,
        theta_kappa = c(1.1, 0.7)),
        noise = noise_nig()
      )
    )
  }
  # plot(mesh)

  n_obs <- 100; index_obs <- sample(1:mesh$n, n_obs)
  loc_obs <- mesh$loc[index_obs, c(1, 2)]
  A <- inla.spde.make.A(mesh = mesh, loc = loc_obs)
  sigma.e <- 0.7
  Y <- drop(A %*% W + sigma.e * rnorm(n_obs))

  ngme_out <- ngme(
    Y ~ 0 + f(
      model = model_matern(
        loc = loc_obs,
        mesh = mesh,
        theta_kappa = c(0.5, 0.5),
        B_kappa = B_kappa
      ),
      fix_theta_K = FALSE,
      # W = as.numeric(W),
      # fix_W = TRUE,
      noise = noise_nig(),
      debug = F,
    ),
    data = list(Y = Y),
    family = noise_normal(),
    control = control_opt(
      estimation = T,
      iterations = 100,
      n_parallel_chain = 1
    ),
    debug = FALSE
  )
  ngme_out
  expect_true(TRUE)
})


# ####### test for kappa
# library(devtools); library(INLA);load_all()

# pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
# mesh <- inla.mesh.2d(
#   loc.domain = pl01, cutoff = 1,
#   max.edge = c(0.3, 1), offset = c(0.5, 1.5)
# )
# # plot(mesh)

# ############################ Stationary Case ######################################
# sigma <- 1
# alpha <- 2
# mu <- 2;
# delta <- -mu
# nu <- 1
# ngme.matern(mesh = mesh)

# kappa = 3
# Kappa <- diag(rep(kappa, mesh$n))
# sigma.e = 0.25

# fem <- inla.mesh.fem(mesh)
# C = fem$c0 ; G = fem$g1
# # # C = diag(rowSums(C))

# # K = (Kappa %*% C %*% Kappa + G)
# C.sqrt.inv <- as(diag(sqrt(1/diag(C))), "sparseMatrix")
# C.inv <- as(diag(1/diag(C)), "sparseMatrix")

# if (alpha==2) {
#   K_a = (Kappa %*% C %*% Kappa + G)
# } else if (alpha==4) {
#   K_a = (Kappa %*% C %*% Kappa + G) %*% C.inv %*% (Kappa %*% C %*% Kappa + G)
# } else {
#   stop("alpha = 2 or 4")
# }

# # # W|V ~ N(solve(K, delta + mu*V), sigma^2*K^(-1)*diag(V)*K^(-1) )

# trueV1 <- ngme2::rig(n_mesh, nu, nu)
# noise1 = -mu + mu * trueV1 + sigma * sqrt(trueV1) * rnorm(mesh$n)
# trueW1 = drop(solve(K_a, noise1))

# trueV2 <- ngme2::rig(n_mesh, nu, nu)
# noise2 = -mu + mu * trueV2 + sigma * sqrt(trueV2) * rnorm(mesh$n)
# trueW2 = drop(solve(K_a, noise2))

# # get n.samples1 locations with 2 observations
# n.samples1 <- 50
# index.samples1 = sample(1:mesh$n, n.samples1)

# loc1 = mesh$loc[index.samples1, c(1,2)]
# A1 = inla.spde.make.A(mesh=mesh, loc=loc1)
# Y1 = A1%*%trueW1 + sigma.e * rnorm(n.samples1); Y1 = drop(Y1)
# Y2 = A1%*%trueW2 + sigma.e * rnorm(n.samples1); Y2 = drop(Y2)

# # get n.samples2 locations with 1 observations
# left.index <- setdiff(1:n_mesh, index.samples1)
# n.samples2 <- 30
# index.samples2 <- sample(left.index, n.samples2)
# loc2 = mesh$loc[index.samples2, c(1,2)]
# A2 = inla.spde.make.A(mesh=mesh, loc=loc2)
# Y3 = A2%*%trueW1 + sigma.e * rnorm(n.samples2); Y3 = drop(Y3)

# A <- bdiag(rbind(A1, A2), A1)
# Y <- c(Y1, Y3, Y2)

# # ?inla.spde.make.index
# # input

# replicate <- c(rep(1, n.samples1 + n.samples2), rep(2, n.samples1))
# index <- c(1:(n.samples1 + n.samples2), 1:n.samples1)

# # inla.spde.make.index(name = "field", n.spde=n_mesh, mesh = mesh, n.repl=2)

# # A <- bdiag(A1, A2)
# # A = inla.spde.make.A(mesh=mesh, loc=rbind(loc1,loc2), repl = c(rep(1,n.samples1), rep(2,n.samples2)))
# # ?ngme.spde.matern

# # ff <- f(1:mesh$n, model=spde, A=A, noise=ngme.noise(type="nig", theta.noise=1)); str(ff)

# formula1 <- Y ~ 0 + f(
#   1:mesh$n,
#   model=matern <- ngme.matern(
#     alpha=alpha,
#     mesh=mesh,
#     kappa=1
#   ),
#   A=A,
#   debug=TRUE,
#   theta.mu=mu,
#   theta.sigma=log(sigma),
#   noise = "nig",
#   control=control_f(
#     numer_grad       = FALSE,
#     use_precond      = TRUE,

#     fix_operator     = FALSE,
#     fix_mu           = FALSE,
#     fix_sigma        = FALSE,
#     fix_noise        = FALSE,
#     use_iter_solver  = FALSE
#   )
# )


# ngme_out1 <- ngme(
#   formula = formula1,
#   data=data.frame(
#     Y=Y
#   ),
#   family = "normal",
#   control=ngme.control(
#     burnin=100,
#     iterations=100,
#     n_gibbs_samples = 5
#   ),
#   debug=ngme.debug(
#     debug = TRUE
#   ),
#   start=ngme.start(
#     # W=trueW
#   )
# )

##################################  Non-stationary Case ######################################
#
# sigma = 1
# alpha = 2
# mu = 2;
# delta = -mu
# nu = 1
#
# n_mesh <- mesh$n
# trueV <- ngme2::rig(n_mesh, nu, nu)
# noise <- delta + mu*trueV + sigma * sqrt(trueV) * rnorm(n_mesh)
# theta.kappa <- c(-1,1)
# B.kappa = cbind(1, -1 * (mesh$loc[,1] - 5) / 10)
# kappa <- drop(exp(B.kappa %*% theta.kappa)); max(kappa)
# Kappa <- diag(kappa)
#
# fem <- inla.mesh.fem(mesh)
# C = fem$c0 ; G = fem$g1
#
# # Kappa = diag(drop(exp(B.kappa %*% c(1, theta))))
# K = (Kappa %*% C %*% Kappa + G)
# C.sqrt.inv <- as(diag(sqrt(1/diag(C))), "sparseMatrix")
#
# if (alpha==2) {
#   K_a =  C.sqrt.inv %*% (Kappa %*% C %*% Kappa + G)
# } else if (alpha==4) {
#   K_a = C.sqrt.inv %*% (Kappa %*% C %*% Kappa + G) %*% C %*% (Kappa %*% C %*% Kappa + G)
# }
#
# # # W|V ~ N(solve(K, delta + mu*V), sigma^2*K^(-1)*diag(V)*K^(-1) )
# trueW = solve(K_a, sqrt(trueV)) * rnorm(mesh$n) + solve(K_a, delta+mu*trueV)
# trueW = drop(trueW)
#
# n.samples = 1000
# loc = mesh$loc[sample(1:mesh$n, n.samples), c(1,2)]
# A = inla.spde.make.A(mesh=mesh, loc=loc)
# dim(A)
#
# sigma.e = 0.25
# Y = A%*%trueW + sigma.e * rnorm(n.samples); Y = drop(Y)
#
# ################# fitting
#
# spde <- ngme.spde.matern(alpha=2,
#                          mesh=mesh,
#                          theta.kappa=c(-1, -1),
#                          B.kappa = B.kappa)
#
# str(spde)
# ff <- f(1:mesh$n, model=spde, A=A, noise=ngme.noise(type="nig", theta.noise=1)); str(ff)
#
# ngme_out <- ngme(
#   formula = Y ~ 0 + f(
#     1:mesh$n,
#     model=spde,
#     A=A,
#     debug=TRUE,
#     theta.mu=2,
#     theta.sigma=log(1),
#     control=control.f(
#       use_num_hess     = FALSE,
#       opt_operator     = TRUE,
#       opt_mu           = FALSE,
#       opt_sigma        = FALSE,
#       opt_var          = FALSE,
#       use_precond      = FALSE
#     )
#   ),
#   data=data.frame(Y=Y),
#   family = "normal",
#   control=control.ngme(
#     burnin=100,
#     iterations=100,
#     n_gibbs_samples = 5
#   ),
#   debug=debug.ngme(fixW = FALSE)
# )
#
# # plot theta.kappa
# plot_out(ngme_out$trajectory, start=1, n=2)
# ngme_out$estimates
#
