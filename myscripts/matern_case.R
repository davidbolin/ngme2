# Matern1d case
data("mcycle", package="MASS")
# mydata = data.frame(times=c(1,2,3), accel=c(2,3,4))

library(devtools)
load_all(reset = FALSE, recompile = FALSE)

########### 3. run ngme
# specify the control for ngme
# args(control.ngme)
# args(control.f)
# args(control.ngme)

control = control.ngme(burnin=100, iterations = 1000,
                       gibbs_sample = 5, stepsize = 1,
                       kill_var = FALSE)
debug = debug.ngme(debug=TRUE)

##### ngme for matern1d
# f(x=mcycle$times, model="matern1d")$operator_in$G # 94*94
# f(x=mcycle$times, model="ar1")$operator_in$G      # 133*133
#
# f(x=mcycle$times, model="matern1d")$operator_in$C

# f(x=mydata$times, model="matern1d")


ngme_out = ngme(accel ~ -1 +
                  f(times, model="matern1d", var="nig",
                    control=control.f(numer_grad = FALSE,
                                      init_kappa    = 0.5,
                                      init_mu       = 0,
                                      init_sigma    = 1,
                                      init_nu       = 1),
                    debug=TRUE
                  ),
                family="normal",
                data=mcycle,
                control=control)

########### 4. results

# nig ar1: a=0.5, mu=2, sigma=2, nu=1
# normal ar1: a=0.7 sigma=3

# plot(ngme_out, param = "fe", type = "traj")
plot(ngme_out, param = "me", type = "traj")
plot(ngme_out, param = "la", type = "traj", which=1)
# plot(ngme_out, param = "la", type = "traj", which=2)
ngme_out$model.types == "matern"
summary(ngme_out)

# mcycle$times
# mesh <- INLA::inla.mesh.1d(mcycle$times)
# mesh
# fem <- INLA::inla.mesh.1d.fem(mesh)
# mesh$loc
# A <- INLA::inla.spde.make.A(mesh, loc=mcycle$times)
# A


# ngme.spde.matern
loc = matrix(runif(10 * 2), ncol=2)
mesh_2d = inla.mesh.2d(loc=loc, max.edge = c(1, 10))
plot(mesh_2d)
points(loc[, 1], loc[, 2])

fem <- inla.mesh.fem(mesh_2d, order=2)
str(fem)


spde.model = ngme.spde.matern(mesh=mesh_2d)
inherits(spde.model, "ngme.spde")
args(ngme)

ngme(formula=Y~0+f(1:length(mesh_2d), model=spde))
load_all()

# simulate data and write test case
