## ----settings,echo=FALSE,results='hide',message=FALSE,warning=FALSE------
source("R/initial_setup.R")
opts_chunk$set(
  fig.path = 'figs/nonstationar-'
)
library(fields)

## ----theta, echo = FALSE, fig.height=4, out.width="69%", fig.cap="Range as function of the first coordinate."----
theta <- c(-1, 0, 1)
par(mfrow = c(1, 1), mar = c(3, 3, 1, 1), mgp = 2:0)
plot(function(x) exp(theta[2] + theta[3] * (x - 5) / 10), 0, 10,
     lwd = 2, xlab = 'first coordinate', ylab = expression(rho(s)))

## ----poly----------------------------------------------------------------
pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)

## ----mesh----------------------------------------------------------------
mesh <- inla.mesh.2d(loc.domain = pl01, cutoff = 0.1,
                     max.edge = c(0.3, 1), offset = c(0.5, 1.5))
?inla.mesh.2d
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
str(Q)
## ----samples-------------------------------------------------------------
sample <- as.vector(inla.qsample(1, Q, seed = 1))

str(sample)
## ----plotsamples, fig.cap = "The simulated random field with increasing range along the horizontal coordinate (top) and the correlation at two location points: $(1,2.5)$ and $(7,2.5)$ (bottom)."----
# Plot parameters
par(mfrow = c(2, 1), mar = c(0, 0, 0, 0))
#Plot Field
proj <- inla.mesh.projector(mesh, xlim = 0:1 * 10, ylim = 0:1 * 5,
                            dims = c(200, 100))
book.plot.field(sample, projector = proj)
# Compute spatial autocorrelation
cx1y2.5 <- book.spatial.correlation(Q, c(1, 2.5), mesh)
cx7y2.5 <- book.spatial.correlation(Q, c(7, 2.5), mesh)
# Plot spatial autocorrelation
book.plot.field(cx1y2.5, projector = proj, zlim = c(0.1, 1))
book.plot.field(cx7y2.5, projector = proj, zlim = c(0.1, 1),
                add = TRUE)

## ----likehy--------------------------------------------------------------
clik <- list(hyper = list(theta = list(initial = 20,
                                       fixed = TRUE)))

## ----res1----------------------------------------------------------------
formula <- y ~ 0 + f(i, model = spde)

res1 <- inla(formula, control.family = clik,
             data = data.frame(y = sample, i = 1:mesh$n))

## ----label = "tabhy1", echo = FALSE--------------------------------------
res1$summary.hyperpar[c(2:4), c(1:3,5)]
tab.res1 <- cbind(true = theta, res1$summary.hyperpar[c(2:4), c(1:3,5)])
tab.res1 <- cbind(
  Parameter = paste("$\\theta_", 1:3, "$", sep = ""), tab.res1)

names(tab.res1) <- c("Parameter", "True", "Mean", "St. Dev.",
                     "2.5\\% quant.", "97.5\\% quant.")

knitr::kable(tab.res1,
             row.names = FALSE,
             caption = "Summary of the posterior distributions of the parameters in the non-stationary example for the data simulated at the mesh nodes.",
             format = "pandoc")

## ----loc-----------------------------------------------------------------
set.seed(2)
n <- 200
loc <- cbind(runif(n) * 10, runif(n) * 5)

## ----projloc-------------------------------------------------------------
projloc <- inla.mesh.projector(mesh, loc)

## ----projectsamples------------------------------------------------------
x <- inla.mesh.project(projloc, sample)

## ----stacks--------------------------------------------------------------
stk <- inla.stack(
  data = list(y = x),
  A = list(projloc$proj$A),
  effects = list(data.frame(i = 1:mesh$n)),
  tag = 'd')

## ----fitt----------------------------------------------------------------
res2 <- inla(formula, data = inla.stack.data(stk),
             control.family = clik,
             control.predictor = list(compute = TRUE, A = inla.stack.A(stk)))

## ----label = "tabres2", echo = FALSE-------------------------------------
tab.res2 <- cbind(true = theta, res2$summary.hyperpar[, c(1:3, 5)])
tab.res2 <- cbind(
  Parameter = paste("$\\theta_", 1:3, "$", sep = ""), tab.res2)

names(tab.res2) <- c("Parameter", "True", "Mean", "St. Dev.",
                     "2.5\\% quant.", "97.5\\% quant.")

knitr::kable(tab.res2,
             row.names = FALSE,
             caption = "Summary of the posterior distributions of the parameters in the non-stationary example for the data simulated at the locations.",
             format = "pandoc")

## ----projpostrf, echo = FALSE--------------------------------------------
x.mean <- inla.mesh.project(proj, field = res2$summary.random$i$mean)
x.sd <- inla.mesh.project(proj, field = res2$summary.random$i$sd)

## ----visprojp, echo = FALSE, fig.height = 9, fig.cap = "Simulated field (top), posterior mean (mid) and posterior standard deviation with the location points added (bottom)."----
par(mfrow=c(3, 1), mar = c(0, 0, 0, 2))
book.plot.field(sample, projector = proj)
book.plot.field(list(x = proj$x, y = proj$y, z = x.mean))
book.plot.field(list(x = proj$x, y = proj$y, z = x.sd),
                col = book.color.c2())
points(loc, cex = 0.5 + 2 * (sample - min(sample)) / diff(range(sample)))

