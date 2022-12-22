# NGME2 test on parana data / similar in ngme2.Rmd

library(INLA)
library(splancs)
library(lattice)
library(ggplot2)
library(grid)
library(gridExtra)
library(viridis)
load_all()

# read data
data(PRprec)
data(PRborder)

############################### Create INLA mesh
coords <- as.matrix(PRprec[, 1:2])
prdomain <- inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(100, 100))
prmesh <- inla.mesh.2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.2)

# monthly mean at each location
Y <- rowMeans(PRprec[, 12 + 1:31]) # 2 + Octobor

ind <- !is.na(Y) # non-NA index
Y <- Y_mean <- Y[ind]
coords <- as.matrix(PRprec[ind, 1:2])
seaDist <- apply(spDists(coords, PRborder[1034:1078, ],
  longlat = TRUE
), 1, min)

ggplot() +
  geom_point(aes(
    x = coords[, 1], y = coords[, 2],
    colour = Y
  ), size = 2, alpha = 1) +
  scale_color_gradientn(colours = viridis(100)) +
  geom_path(aes(x = PRborder[, 1], y = PRborder[, 2])) +
  geom_path(aes(x = PRborder[1034:1078, 1], y = PRborder[
    1034:1078,
    2
  ]), colour = "red")

############################### fit the matern model

matern_spde <- model_matern(
  loc = coords,
  mesh = prmesh
)

out <- ngme(
  formula = Y ~ 1 +
    f(inla.group(seaDist), model = "rw1", noise = noise_normal()) +
    f(model = matern_spde, noise = noise_nig()),
  data = list(
    Y = Y
  ),
  family = noise_nig(),
  control = ngme_control(
    estimation = T,
    iterations = 10,
    n_slope_check = 4,
    stop_points = 10,
    std_lim = 0.1,
    n_parallel_chain = 1,
    print_check_info = FALSE
  ),
  debug = TRUE,
  seed = 16
)
