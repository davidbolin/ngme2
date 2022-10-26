library(INLA)
library(devtools)

load_all()

data(PRprec)
data(PRborder)

# dim(PRprec)
# dim(PRborder)
# PRprec[1:6, 1:8]

prdomain <- inla.nonconvex.hull(as.matrix(PRprec[, 1:2]),
  convex = -0.03, concave = -0.05,
  resolution = c(100, 100))

coords <- as.matrix(PRprec[, 1:2])

dim(PRprec)

#
Y1 <- rowMeans(PRprec[, 12 + 1:31]) # 2 + Octobor
Y2 <- apply(PRprec[, 12 + 1:31], 1, max) # 2 + Octobor

ind <- !is.na(Y1)

Y1 <- Y1[ind]
Y2 <- Y2[ind]

coords <- as.matrix(PRprec[ind, 1:2])
alt <- PRprec$Altitude[ind]
seaDist <- apply(spDists(coords, PRborder[1034:1078, ],
  longlat = TRUE
), 1, min)


# 1 + long + lat

prdomain <- inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(100, 100))
prmesh <- inla.mesh.2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.2)


A <- inla.spde.make.A(mesh = prmesh, loc = coords)

mesh.index <- inla.spde.make.index(
  name = "field",
  mesh = prmesh,
  n.spde = prmesh$n
)

data <- data.frame(
  Y_mean  = Y1,
  Y_max   = Y2,
  long    = coords[, 1],
  lat     = coords[, 2]
)

load_all()
out <- ngme(
  formula = Y_mean ~ 1 +
    f(inla.group(seaDist), model="rw1", noise=noise_normal()) +
    f(index = mesh.index$field, model = model_matern(
      A = A,
      mesh = prmesh,
      noise = noise_nig()
    )) +
    f(index = mesh.index$field, model = model_matern(
      A = A,
      mesh = prmesh,
      noise = noise_normal()
    )),
  data = data,
  family = noise_nig(),
  control = ngme_control(
    estimation = T,
    iterations = 10,
    n_slope_check = 4,
    stop_points = 100,
    n_parallel_chain = 8
  ),
  debug = T
)

str(out)


# plots
# beta
# traceplot(out, parameter = "beta",    f_index = 0, param_index = 1)
# traceplot(out, parameter = "beta",    f_index = 0, param_index = 2)



# 1st model
traceplot(out, parameter = "theta_K",     f_index = 1, transform = exp)
traceplot(out, parameter = "theta_mu",    f_index = 1)
traceplot(out, parameter = "theta_sigma", f_index = 1, transform = exp)
traceplot(out, parameter = "theta_V",     f_index = 1, transform = log)

# 2nd model
traceplot(out, parameter = "theta_K",     f_index = 2, transform = exp)
traceplot(out, parameter = "theta_mu",    f_index = 2)
traceplot(out, parameter = "theta_sigma", f_index = 2, transform = exp)
traceplot(out, parameter = "theta_V",     f_index = 2, transform = log)

# 3rd model
traceplot(out, parameter = "theta_K",     f_index = 3, transform = exp)
traceplot(out, parameter = "theta_mu",    f_index = 3)
traceplot(out, parameter = "theta_sigma", f_index = 3, transform = exp)
traceplot(out, parameter = "theta_V",     f_index = 3, transform = log)


# fixed effects
traceplot(out, parameter = "beta", f_index = 0, param_index = 1)

# measurement noise
traceplot(out, parameter = "theta_mu", f_index = 0)
traceplot(out, parameter = "theta_sigma", f_index = 0, transform = exp)
traceplot(out, parameter = "theta_V", f_index = 0)
