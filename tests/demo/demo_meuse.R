# spatial data
library(sp)
data(meuse)
head(meuse)
str(meuse)
coordinates(meuse) <- ~x + y
class(meuse)
summary(meuse)

plot(meuse)
title("points")

#  --- grid ---

data(meuse.grid)
head(meuse.grid)
str(meuse.grid)

coordinates(meuse.grid) <- ~ x + y
meuse.grid <- as(meuse.grid, "SpatialPixels")
image(meuse.grid, col = "grey")
title("grid")

# --- river ---

data(meuse.riv)
meuse.lst <- list(Polygons(list(Polygon(meuse.riv)), "meuse.riv"))
meuse.pol <- SpatialPolygons(meuse.lst)
plot(meuse.pol, col="grey")


image(meuse.grid, col = "lightgrey")
plot(meuse.grid, col = "lightgrey", add = TRUE)
plot(meuse.pol, col="grey", add = TRUE)
plot(meuse, add = TRUE)

##############################
# create mesh
library(INLA)
mesh <- inla.mesh.2d(loc=meuse@coords, max.n=200, max.edge = c(1,2))
plot(mesh)
mesh$n

str(mesh1)
points(meuse, col = "red")

head(meuse)

prdomain <- inla.nonconvex.hull(meuse@coords, -0.03, -0.05, resolution = c(100, 100))
prmesh <- inla.mesh.2d(meuse@coords, boundary = prdomain, max.n = 400, cutoff = 50)
prmesh$n
plot(prmesh)

##############################
# fit using ngme
load_all()
bubble(meuse, "copper")

spde <- model_matern(loc = meuse@coords, mesh=mesh)

out1 <- ngme(
  zinc ~ 1 # + f(dist.m, model="rw1", noise=noise_normal(), name="rw")
    + f(model=spde, name="spde", noise=noise_nig()),
  data = meuse@data,
  control = control_opt(
    iterations = 100,
    n_parallel_chain = 4,
    print_check_info = TRUE
  ),
  debug=TRUE
)

out1$latents[[1]]$h

traceplot(out1, "rw")
traceplot(out1, "rw")
traceplot(out1, "spde")

class(meuse)

library(INLA)

n_meuse <- nrow(meuse)
spde_inla <- inla.spde2.matern(mesh = prmesh, loc = meuse@coords)
res_inla <- inla(
  formula = zinc ~ 0  +
    f(dist.m, model = "rw1") +
    f(i, model = spde_inla),
  data = list(dist.m = meuse$dist.m, i = 1:spde_inla$n.spde, zinc = meuse$zinc),
  family = "gaussian",
  control.predictor = list(1, A = inla.spde.make.A(mesh=prmesh, loc=meuse@coords))
)
list(dist.m = meuse$dist.m, i = 1:spde_inla$n.spde, zinc = meuse$zinc)
summary(res_inla)

meuse@data$copper
