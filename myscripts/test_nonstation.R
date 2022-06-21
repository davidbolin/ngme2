######## 1. test AR with non-stationary simga
n_obs <- 100
sigma_eps = 0.5

alpha <- 0.5
mu = 2;
delta = -mu1
nu = 1
B.sigma <- cbind(1, rexp(n_obs))
sigma <- drop(B.sigma %*% c(2, 3))

trueV <- ngme2::rig(n_obs, nu, nu)
noise <- delta + mu*trueV + sigma * sqrt(trueV) * rnorm(n_obs)

trueW <- Reduce(function(x,y){y + alpha*x}, noise, accumulate = T)
Y = trueW + rnorm(n_obs, mean=0, sd=sigma_eps)


# fitting
control = control.ngme(burnin=100, iterations = 100,
                       gibbs_sample = 5, stepsize = 1,
                       kill_var = FALSE, threshold = 1e-4,
                       opt_fix_effect = T)

debug = debug.ngme(fixW = FALSE)

ngme_out = ngme(Y ~ 0 +
                  f(1:length(Y),
                    model = "ar1",
                    var = "nig",
                    control = control.f(numer_grad       = FALSE,
                                      init_operator    = k2th(0.5),
                                      init_var         = 1,
                                      opt_operator     = F,
                                      opt_mu           = F,
                                      opt_sigma        = T,
                                      opt_var          = F),
                    theta.mu = 2,
                    B.sigma = B.sigma,
                    theta.sigma = log(2),
                    debug = T),
                family = "normal",
                data = data.frame(Y=(as.numeric(Y))),
                control = control)

# ###### test for non-stationary case
# library(INLA)
# library(devtools)
# load_all()

# ################### 1. create mesh
# pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
# mesh <- inla.mesh.2d(loc.domain = pl01, cutoff = 0.3,
#                      max.edge = c(0.3, 1), offset = c(0.5, 1.5))
# plot(mesh)
# mesh$n

# ################### 2. simulation nig (alpha = 2)
# # sigma = 1
# alpha = 2
# mu = 2;
# delta = -mu
# nu = 1

# # noise <- delta + mu*trueV + sigma * sqrt(trueV) * rnorm(n_obs)
# trueV <- ngme2::rig(mesh$n, nu, nu)
# theta <- c(-1,0,1)

# logtau0 <- logkappa0 <- 1
# B.kappa = cbind(logkappa0, 0, -1, -1 * (mesh$loc[,1] - 5) / 10)

# fem <- inla.mesh.fem(mesh)
# C = fem$c0 ; G = fem$g1
# # C = diag(rowSums(C))

# Tau = diag(drop(exp(B.tau %*% c(1, theta))))
# Kappa = diag(drop(exp(B.kappa %*% c(1, theta))))
# # K = (Kappa %*% C %*% Kappa + G)

# C.sqrt.inv <- as(diag(sqrt(1/diag(C))), "sparseMatrix")

# if (alpha==2) {
#   K_a = Tau %*% (Kappa %*% C %*% Kappa + G) %*% C.sqrt.inv
# } else if (alpha==4) {
#   K_a = Tau %*% (Kappa %*% C %*% Kappa + G) %*% C %*% (Kappa %*% C %*% Kappa + G) %*% C.sqrt.inv
# }

# # W|V ~ N(solve(K, delta + mu*V), sigma^2*K^(-1)*diag(V)*K^(-1) )
# trueW = diag(sqrt(diag(1/Tau))) %*% solve(K_a, sqrt(trueV)) * rnorm(mesh$n) + solve(K_a, delta+mu*trueV)
# trueW = drop(trueW)

# n.samples = 300
# loc = mesh$loc[sample(1:mesh$n, n.samples), c(1,2)]
# A = inla.spde.make.A(mesh=mesh, loc=loc)
# dim(A)

# sigma.e = 0.25
# Y = A%*%trueW + sigma.e * rnorm(n.samples); Y = drop(Y)

# ###################### 3. NGME
# # ?ngme.spde.matern

# spde <- ngme.spde.matern(alpha=2,
#                          mesh=mesh,
#                          theta.init=c(0,2,3),
#                          B.tau = cbind(logtau0, -1, nu, nu * (mesh$loc[,1] - 5) / 10),
#                          B.kappa = cbind(logkappa0, 0, -1, -1 * (mesh$loc[,1] - 5) / 10))



# str(spde)
# str(ff <- f(1:mesh$n, model=spde, A=A,
#             B.sigma=,

#             B.mu))

# ngme.noise.spec(noise="nig", )

# class(ff$operator_in$B.tau)

# res <- ngme(formula = Y ~ 0 + f(
#   1:mesh$n, model=spde, A=A, debug=TRUE,
#   # noise="nig",
#   # noise=ngme.noise.spec(),
#   control=control.f(use_num_hess = FALSE)),
#             data=data.frame(Y=Y),
#             family = "normal",
#             debug=debug.ngme(fixW = FALSE),
#             control=control.ngme(
#               burnin=10,
#               iterations=10,
#               gibbs_sample = 5),
#   # ngme.start(result from ngme object(lastW, lastV)),
#             )

# # plot operator
# plot_out(res$trajectory, start=1, n=3)
# # plot mu
# plot_out(res$trajectory, start=4, n=1, ylab="mu")
# # plot nu
# plot_out(res$trajectory, start=6, n=1, transform = exp, ylab="nu")
# # plot sigma.e
# plot_out(res$trajectory, start=7, n=1, transform = exp, ylab="sigma_eps")

# # results
# res$estimates


