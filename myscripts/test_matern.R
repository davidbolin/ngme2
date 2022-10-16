library(devtools); library(INLA); load_all()
{ # First we create mesh
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh <- inla.mesh.2d(
    loc.domain = pl01, cutoff = 1,
    max.edge = c(0.3, 1), offset = c(0.5, 1.5)
  )

  W <- simulate(
    f(model = ngme.matern(mesh = mesh, theta_kappa = 1.1),
      noise = ngme.noise.nig()
    )
  )
}

# generate A and A_pred
n_obs <- 100; index_obs <- sample(1:mesh$n, n_obs)
loc_obs <- mesh$loc[index_obs, c(1, 2)]
A <- inla.spde.make.A(mesh = mesh, loc = loc_obs)

sigma.e <- 0.7
Y <- drop(A %*% W + sigma.e * rnorm(n_obs))

load_all()
ngme_out <- ngme(
  Y ~ 0 + f(
    model = ngme.matern(mesh = mesh, theta_kappa = 0.4),
    fix_theta_K = FALSE,
    # W = as.numeric(W),
    # fix_W = TRUE,

    noise = ngme.noise.nig(
      fix_theta_mu    = F,
      fix_theta_sigma = F,
      fix_theta_V     = F
    ),
    A = A,
    debug = TRUE,
    control = ngme.control.f(
      numer_grad = F,
      use_precond = F
    )
  ),
  data = list(Y = Y),
  noise = ngme.noise.normal(),
  control = ngme.control(
    iterations = 500,
    estimation = TRUE,
    n_parallel_chain = 2
  ),
  debug = TRUE
)
ngme_out
str(ngme_out)

plot_chains(ngme_out, parameter = "theta_sigma", f_index = 0)

# matern model
plot_chains(ngme_out, parameter = "theta_K",     f_index = 1)
plot_chains(ngme_out, parameter = "theta_mu",    f_index = 1)
plot_chains(ngme_out, parameter = "theta_sigma", f_index = 1)
plot_chains(ngme_out, parameter = "theta_V",     f_index = 1)


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

# replicates <- c(rep(1, n.samples1 + n.samples2), rep(2, n.samples1))
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
#   control=ngme.control.f(
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
#   control = ngme.control(
#     estimation = TRUE,
#     n_parallel_chain = 1,
#     burnin = 200,
#     iterations = 100,
#     gibbs_sample = 5,
#     stepsize = 1,
#     kill_var = FALSE,
#     threshold = 1e-4
#   ),
#   debug = TRUE,
#   noise = ngme.noise.normal(
#     fix_theta_sigma = FALSE
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
#     gibbs_sample = 5
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
#     gibbs_sample = 5
#   ),
#   debug=debug.ngme(fixW = FALSE)
# )
#
# # plot theta.kappa
# plot_out(ngme_out$trajectory, start=1, n=2)
# ngme_out$estimates
#
