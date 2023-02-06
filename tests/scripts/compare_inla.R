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

spde_res <- inla.spde.result(result_inla, "time", spde)

# estimation results
summary(result_inla)
plot(spde_res$marginals.kappa$kappa.1,
  type="l", main="kappa")
mean(spde_res$marginals.kappa$kappa.1[, 1])

plot(spde_res$marginals.tau$tau.1, type="l", main="tau")
1 / mean(spde_res$marginals.tau$tau.1[, 1])

with(mcycle, {plot(times, accel)})
lines(mesh$loc, result_inla$summary.random$time[, "mean"], col=2)

############################## using ngme2
load_all()
spde_ngme <- model_matern(mesh = mesh, loc = mcycle$times)
result_ngme <- ngme(
  accel ~ 0 + f(model = spde_ngme, name="myspde"),
  data = mcycle,
  family = "normal",
  control = ngme_control(
    iterations = 1000
  )
)
result_ngme
traceplot(result_ngme, "myspde")
traceplot(result_ngme, "mn")

with(mcycle, {plot(times, accel)})
lines(mesh$loc, result_ngme$latents[["myspde"]]$W, lwd=2)
lines(mesh$loc, result_inla$summary.random$time[, "mean"], col=2, lwd=2)
postW <- predict(result_ngme, loc=list(myspde = mesh$loc))
postW_mode <- predict(result_ngme, loc=list(myspde = mesh$loc), estimator="mode")
postW_median <- predict(result_ngme, loc=list(myspde = mesh$loc), estimator="median")
postW_q25 <- predict(result_ngme, loc=list(myspde = mesh$loc), estimator="quantile", q=0.25)
postW_q75 <- predict(result_ngme, loc=list(myspde = mesh$loc), estimator="quantile", q=0.75)

lines(mesh$loc, postW_median, col=3, lwd=2)
lines(mesh$loc, postW_q25, col=4, lwd=2)
lines(mesh$loc, postW_q75, col=5, lwd=2)
lines(mesh$loc, postW, col=3, lwd=2)
lines(mesh$loc, postW_mode, col=6, lwd=2)

# refit the model using nig noise
result_ngme2 <- ngme(
  accel ~ 0 + f(model = spde_ngme, name="myspde", noise=noise_nig()),
  data = mcycle,
  family = "normal",
  control = ngme_control(
    n_parallel_chain = 4,
    iterations = 2000
  ),
  start = result_ngme
)
result_ngme2
traceplot(result_ngme2, "myspde")

############################# model prediction using inla / ngme2
rg <- range(mcycle$times)
locs <- seq(from=rg[1], to=rg[2], length = 100)

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

# doing prediction
result_inla_prd <- inla(
  accel ~ -1 + f(time, model=spde),
  data = inla.stack.data(stk),
  control.predictor = list(A = inla.stack.A(stk)),
  control.compute = list(config=TRUE)
)
# summary(result_inla_prd)

str(stk)
prd_idx <- inla.stack.index(stk, "pred")$data
prd_inla <- result_inla_prd$summary.fitted.values[prd_idx, 1] # mean

# predict with ngme
prd_ngme <- predict(result_ngme2, loc = list(spde=locs))
with(mcycle, {plot(times, accel)})
lines(locs, prd_inla)
lines(locs, prd_ngme, col=2)
