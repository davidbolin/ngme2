####### test for non-stationary kappa
library(devtools)
library(INLA)
load_all()

# ################### 1. create mesh
pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
mesh <- inla.mesh.2d(loc.domain = pl01, cutoff = 0.2,
                     max.edge = c(0.3, 1), offset = c(0.5, 1.5))
plot(mesh)

# ################### 2. simulation nig (alpha = 2)
sigma = 1
alpha = 2
mu = 2;
delta = -mu
nu = 1

n_mesh <- mesh$n
trueV <- ngme2::rig(n_mesh, nu, nu)
noise <- delta + mu*trueV + sigma * sqrt(trueV) * rnorm(n_mesh)

kappa = 4
Kappa <- diag(rep(kappa, mesh$n))

# kappa <- drop(exp(B.kappa %*% theta.kappa)); max(kappa)

fem <- inla.mesh.fem(mesh)
C = fem$c0 ; G = fem$g1
# # C = diag(rowSums(C))

# Tau = diag(drop(exp(B.tau %*% c(1, theta))))
# Kappa = diag(drop(exp(B.kappa %*% c(1, theta))))
# K = (Kappa %*% C %*% Kappa + G)
C.sqrt.inv <- as(diag(sqrt(1/diag(C))), "sparseMatrix")

if (alpha==2) {
  K_a =  C.sqrt.inv %*% (Kappa %*% C %*% Kappa + G)
} else if (alpha==4) {
  K_a = C.sqrt.inv %*% (Kappa %*% C %*% Kappa + G) %*% C %*% (Kappa %*% C %*% Kappa + G)
}

# # W|V ~ N(solve(K, delta + mu*V), sigma^2*K^(-1)*diag(V)*K^(-1) )
trueW = solve(K_a, sqrt(trueV)) * rnorm(mesh$n) + solve(K_a, delta+mu*trueV)
trueW = drop(trueW)

n.samples = 500
loc = mesh$loc[sample(1:mesh$n, n.samples), c(1,2)]
A = inla.spde.make.A(mesh=mesh, loc=loc)
dim(A)

sigma.e = 0.25
Y = A%*%trueW + sigma.e * rnorm(n.samples); Y = drop(Y)

# ###################### 3. NGME
# ?ngme.spde.matern

# ff <- f(1:mesh$n, model=spde, A=A, noise=ngme.noise(type="nig", theta.noise=1)); str(ff)

ngme_out <- ngme(
  formula = Y ~ 0 + f(
    1:mesh$n,
    model=ngme.matern(
      alpha=2,
      mesh=mesh,
      kappa=1
    ),
    A=A,
    debug=TRUE,
    theta.mu=2,
    theta.sigma=log(1),
    control=ngme.control.f(
      use_num_hess     = FALSE,
      opt_operator     = TRUE,
      opt_mu           = FALSE,
      opt_sigma        = FALSE,
      opt_var          = FALSE,
      use_precond      = FALSE
    )
  ),
  data=data.frame(Y=Y),
  family = "normal",
  control=ngme.control(
    burnin=100,
    iterations=500,
    gibbs_sample = 5
  ),
  debug=ngme.debug(fixW = FALSE)
)
# results

# plot mu
plot_out(ngme_out$trajectory, start=2, n=1)
# plot sigma
plot_out(ngme_out$trajectory, start=3, n=1, transform = exp)
# plot var
plot_out(ngme_out$trajectory, start=4, n=1, transform = exp)
# plot m err
plot_out(ngme_out$trajectory, start=5, n=1, transform = exp)

# plot kappa
plot_out(ngme_out$trajectory, start=1, n=1, transform=exp)
ngme_out$estimates

#   # ngme.start(result from ngme object(lastW, lastV)),
# # plot operator
# # plot mu
# plot_out(res$trajectory, start=4, n=1, ylab="mu")
# # plot nu
# plot_out(res$trajectory, start=6, n=1, transform = exp, ylab="nu")
# # plot sigma.e
# plot_out(res$trajectory, start=7, n=1, transform = exp, ylab="sigma_eps")

# # results
# res$estimates



