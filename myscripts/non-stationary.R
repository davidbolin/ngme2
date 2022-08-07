theta <- c(-1, 0, 1)
par(mfrow = c(1, 1), mar = c(3, 3, 1, 1), mgp = 2:0)
plot(function(x) exp(theta[2] + theta[3] * (x - 5) / 10), 0, 10,
     lwd = 2, xlab = 'first coordinate', ylab = expression(rho(s)))

## ----poly----------------------------------------------------------------
pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)

## ----mesh----------------------------------------------------------------
mesh <- inla.mesh.2d(loc.domain = pl01, cutoff = 0.1,
                     max.edge = c(0.3, 1), offset = c(0.5, 1.5))

## ----spde----------------------------------------------------------------
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
cbind(logtau0, -1, nu, nu * (mesh$loc[,1] - 5) / 10)
## ----Q-------------------------------------------------------------------
?inla.spde2.precision
Q <- inla.spde2.precision(spde, theta = theta)

## ----samples-------------------------------------------------------------
sample <- as.vector(inla.qsample(1, Q, seed = 1))

## ----likehy--------------------------------------------------------------
clik <- list(hyper = list(theta = list(initial = 20,
                                       fixed = TRUE)))

## ----res1----------------------------------------------------------------
formula <- y ~ 0 + f(i, model = spde)

res1 <- inla(formula, control.family = clik,
             data = data.frame(y = sample, i = 1:mesh$n))



################### 1. create mesh
pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
mesh <- inla.mesh.2d(loc.domain = pl01, cutoff = 0.1,
                     max.edge = c(0.3, 1), offset = c(0.5, 1.5))
plot(mesh)

################### 2. simulation nig (alpha = 2)
# sigma = 1
n_obs = 10
alpha = 2
mu = 2; delta = -mu
nu = 1

# noise <- delta + mu*trueV + sigma * sqrt(trueV) * rnorm(n_obs)
trueV <- ngme2::rig(mesh$n, nu, nu)
theta <- c(-1,0,1)

logtau0 <- logkappa0 <- 1
B.tau = cbind(logtau0, -1, nu, nu * (mesh$loc[,1] - 5) / 10)
B.kappa = cbind(logkappa0, 0, -1, -1 * (mesh$loc[,1] - 5) / 10)

fem <- inla.mesh.fem(mesh)
C = fem$c0 ; G = fem$g1
# C = diag(rowSums(C))

Tau = diag(drop(exp(B.tau %*% c(1, theta))))
Kappa = diag(drop(exp(B.kappa %*% c(1, theta))))


if (alpha==2) {
  K_a = (Kappa %*% C %*% Kappa + G)
  # K_a = Tau %*% (Kappa %*% C %*% Kappa + G) %*% sqrt(C)
} else if (alpha==4) {
  # K_a = Tau %*% (Kappa %*% C %*% Kappa + G) %*% sqrt(C)
}

# W|V ~ N(solve(K, delta + mu*V), sigma^2*K^(-1)*diag(V)*K^(-1) )
trueW = diag(sqrt(diag(1/Tau))) %*% solve(K_a, sqrt(trueV)) * rnorm(mesh$n) + solve(K_a, delta+mu*trueV)

loc = mesh$loc[sample(1:mesh$n, 1000), c(1,2)]
inla.spde.make.A(mesh=mesh, loc=loc)

A = inla.spde.make.A(mesh=mesh, loc=loc)
dim(A)
sigma.e = 0.1
Y = A%*%trueW + sigma.e * rnorm(1000)


###################### 3. NGME
?ngme.spde.matern
spde <- ngme.spde.matern(alpha=2,
                         mesh=mesh,
                         B.tau = cbind(logtau0, -1, nu, nu * (mesh$loc[,1] - 5) / 10),
                         B.kappa = cbind(logkappa0, 0, -1, -1 * (mesh$loc[,1] - 5) / 10))

str(spde)
f(1:mesh$n, model=spde)
ngme(formula = Y ~ 0 + f(1:mesh$n, model=spde, var="nig"),
     data,
     family = "normal")






# inla.models("rw")


Mat = matrix(rep(1, 9), nrow=3); diag(Mat)=2; Mat
chol(Mat)


