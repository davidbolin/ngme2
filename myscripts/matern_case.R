library(devtools)
library(INLA)
load_all(reset = FALSE, recompile = FALSE)

data("mcycle", package="MASS")
# mydata = data.frame(times=c(1,2,3), accel=c(2,3,4))

# Matern 2d case
pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
mesh <- inla.mesh.2d(loc.domain = pl01, cutoff = 0.5,
                     max.edge = c(0.5, 1), offset = c(0.5, 1.5))
plot(mesh)
points(pl01)

################## simulation (alpha = 2)
n_obs = 10
alpha = 2
mu = 2; delta = -mu
nu = 1
sigma = 2

trueV <- ngme2::rig(n_obs, nu, nu)
noise <- delta + mu*trueV + sigma * sqrt(trueV) * rnorm(n_obs)
trueW1 <- Reduce(function(x,y){y + alpha1*x}, noise1, accumulate = T)

# w|V = delta + mu*V + sigma*sqrt(V)*Z ~ N(delta+mu*V, sigma^2 diag(V))


Y1 = trueW1 + rnorm(n_obs, mean=0, sd=sigma_eps)


nu <- 1
alpha <- nu + 2 / 2
# log(kappa)
logkappa0 <- log(8 * nu) / 2
# log(tau); in two lines to keep code width within range
logtau0 <- (lgamma(nu) - lgamma(alpha) -1 * log(4 * pi)) / 2
logtau0 <- logtau0 - logkappa0
# SPDE model
spde <- inla.spde2.matern(mesh,
                          B.tau = cbind(logtau0, -1, nu, nu * (mesh$loc[,1] - 5) / 10),
                          B.kappa = cbind(logkappa0, 0, -1, -1 * (mesh$loc[,1] - 5) / 10),
                          theta.prior.mean = rep(0, 3),
                          theta.prior.prec = rep(1, 3))
theta <- c(-1, 0, 1)
Q <- inla.spde2.precision(spde, theta = theta)
Q
