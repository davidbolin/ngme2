# Doing comparison with INLA and ngme2
# using dataset mcycle
# in this script, mainly show how to fit the model in INLA, and ngme2
set.seed(17)
library(MASS)
library(INLA)

data(mcycle)
str(mcycle)
with(mcycle, {plot(times, accel)})
mesh <- inla.mesh.1d(loc=mcycle$times, max.edge=c(1, 10))
mesh$n
points(mesh$loc, rep(0, mesh$n), col=2)

############################## 1. model estimation inla / ngme2
# fit use INLA
spde <- inla.spde2.matern(mesh, alpha=2)
A <- inla.spde.make.A(mesh, loc=mcycle$times)
ind <- inla.spde.make.index("time", spde$n.spde)
data <- list(accel = mcycle$accel, time = ind$time)

result_inla <- inla(
  accel ~ -1 + f(time, model=spde),
  data = data,
  control.predictor = list(A = A),
  control.compute = list(config=TRUE)
)
summary(result_inla)
spde_res <- inla.spde.result(result_inla, "time", spde)
summary(spde_res)
sqrt(1 / 0.002)
# estimation results
summary(result_inla)
plot(spde_res$marginals.kappa$kappa.1,
  type="l", main="kappa")
mean(spde_res$marginals.kappa$kappa.1[, 1])

plot(spde_res$marginals.tau$tau.1, type="l", main="tau")
mean(spde_res$marginals.tau$tau.1[, 1])
names(spde_res)
exp(spde_res$summary.log.kappa)
exp(spde_res$summary.log.tau)
with(mcycle, {plot(times, accel)})
lines(mesh$loc, result_inla$summary.random$time[, "mean"], col=2)
############################## using ngme2
library(ngme2)
spde_ngme <- model_matern(alpha=2, mesh = mesh, loc = mcycle$times)
result_ngme <- ngme(
  accel ~ -1 + f(model = spde_ngme, name="myspde"),
  data = mcycle,
  family = "normal",
  control = control_opt(
    iterations = 1000
  )
)
result_ngme
traceplot(result_ngme, "myspde")
traceplot(result_ngme, "mn")

with(mcycle, {plot(times, accel)})
lines(mesh$loc, result_ngme$latents[["myspde"]]$W, lwd=2)
lines(mesh$loc, result_inla$summary.random$time[, "mean"], col=2, lwd=2)
pred_W <- predict(result_ngme, loc=list(myspde = mesh$loc))
# by default we compute, a bunch of statistics for the given location
str(pred_W)

with(mcycle, {plot(times, accel)})
lines(mesh$loc, pred_W[["median"]], col=3, lwd=2)
lines(mesh$loc, pred_W[["5quantile"]], col=5, lwd=2)
lines(mesh$loc, pred_W[["95quantile"]], col=4, lwd=2)

# refit the model using nig noise
result_ngme2 <- ngme(
  accel ~ 0 + f(model = spde_ngme, name="myspde", noise=noise_nig()),
  data = mcycle,
  family = "normal",
  control = control_opt(
    n_parallel_chain = 4,
    iterations = 2000
  ),
  start = result_ngme
)
result_ngme2
traceplot(result_ngme2, "myspde")
plot(result_ngme2$latents[["myspde"]]$noise)

############################# model prediction using inla / ngme2
rg <- range(mcycle$times)
rg
locs <- seq(from=rg[1], to=rg[2], length = 100)
locs

# predict with INLA
stk.dat <- inla.stack(
  data = list(accel = mcycle$accel),
  A = list(A),
  effects = list(c(ind)),
  tag = "est"
)
A.prd <- inla.spde.make.A(mesh, loc = locs)
stk.prd <- inla.stack(
  data = list(accel = NA),
  A = list(A.prd),
  effects = list(c(ind)),
  tag = "pred"
)
stk <- inla.stack(stk.dat, stk.prd)
# str(stk)

# doing prediction
result_inla_prd <- inla(
  accel ~ -1 + f(time, model=spde),
  data = inla.stack.data(stk),
  control.predictor = list(A = inla.stack.A(stk)),
  control.compute = list(config=TRUE)
)
# summary(result_inla_prd)

prd_idx <- inla.stack.index(stk, "pred")$data
prd_inla <- result_inla_prd$summary.fitted.values[prd_idx, 1] # mean

# predict with ngme
prd_ngme <- predict(result_ngme2, loc = list(myspde=locs))[["mean"]]

with(mcycle, {plot(times, accel)})
lines(locs, prd_inla)
lines(locs, prd_ngme, col=2)
