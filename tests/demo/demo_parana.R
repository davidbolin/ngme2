# This script is used for fitting different matern models using Parana data
############################### Create INLA mesh
{
library(INLA)
library(splancs)
library(lattice)
library(ggplot2)
library(grid)
library(gridExtra)
library(viridis)

# read data
data(PRprec)
data(PRborder)

str(PRprec)
str(PRborder)

# maps
plot(PRborder[, 1], PRborder[, 2], type="l")
# plot the seashore
lines(PRborder[1034:1078, 1], PRborder[1034:1078, 2], col=2, lwd=2)

# create mesh
coords <- as.matrix(PRprec[, 1:2])
prdomain <- inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(100, 100))
prmesh <- inla.mesh.2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.2)
prmesh2 <- inla.mesh.2d(coords, max.edge = c(0.45, 1), cutoff = 0.2)
plot(prmesh)

# monthly mean at each location
colnames(PRprec)
Y <- rowMeans(PRprec[, 3 + 1:31]) # January
any(is.na(Y))

ind <- !is.na(Y) # non-NA index
sum(ind)
Y <- Y[ind]

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
}
############################### Matern model and predict location
{
# matern model
library(ngme2)
matern_spde <- model_matern(
  loc = coords,
  mesh = prmesh
)

# make prediction at locations
projgrid <- inla.mesh.projector(
  prmesh,
  xlim = range(PRborder[, 1]),
  ylim = range(PRborder[, 2]),
  dims = c(300, 200)
)

# get the points inside the boundary idx
xy.in <- inout(projgrid$lattice$loc, cbind(PRborder[, 1], PRborder[, 2]))
coord.prd <- projgrid$lattice$loc[xy.in, ]

plot(coord.prd, type = "p", cex = 0.1)
lines(PRborder)
points(coords[, 1], coords[, 2], pch = 19, cex = 0.5, col = "red")

coord.prd.df <- data.frame(
  x1 = coord.prd[,1],
  x2 = coord.prd[,2]
)
coordinates(coord.prd.df) <- c("x1", "x2")

# compute the distance to the shore
seaDist.prd <- apply(spDists(coord.prd,
  PRborder[1034:1078, ],
  longlat = TRUE
), 1, min)
coord.prd.df$seaDist <- seaDist.prd
str(coord.prd.df)
head(coord.prd.df)
}


################ case 1. spde_nig + normal measurement nosie
{
out_nig_norm <- ngme(
  formula = Y ~ 1 +
    f(inla.group(seaDist), model = "rw1", noise = noise_normal(), name="rw") +
    f(model = matern_spde, noise = noise_nig(), name="spde1"),
  data = list(
    Y = Y
  ),
  family = noise_normal(),
  control = ngme_control(
    estimation = T,
    iterations = 1000,
    n_slope_check = 4,
    stop_points = 10,
    std_lim = 0.1,
    n_parallel_chain = 4,
    print_check_info = TRUE
  ),
  # debug = TRUE,
  seed = 16
)
traceplot(out_nig_norm, name="rw")
traceplot(out_nig_norm, name="spde1")

{ # make prediction
str(coord.prd.df@coords)
preds <- predict(out_nig_norm, loc = list(
  rw = seaDist.prd,
  spde1 = coord.prd.df@coords
))

pred_df <- data.frame(
  x1 = coord.prd[, 1],
  x2 = coord.prd[, 2],
  preds = preds
)
ggplot(pred_df, aes(x = x1, y = x2, fill = preds)) +
  geom_raster() +
  scale_fill_viridis()
}
}

################ case 2. spde_nig + spde_normal + normal measurement nosie
{
out_nig_norm_norm <- ngme(
  formula = Y ~ 1 +
    f(inla.group(seaDist), model = "rw1", noise = noise_normal(), name="rw") +
    f(model = matern_spde, noise = noise_nig(), name="spde1") +
    f(model = matern_spde, noise = noise_normal(), name="spde2"),
  data = list(
    Y = Y
  ),
  family = noise_normal(),
  control = ngme_control(
    estimation = T,
    iterations = 1000,
    n_slope_check = 4,
    stop_points = 10,
    std_lim = 0.1,
    n_parallel_chain = 4,
    print_check_info = TRUE
  ),
  # debug = TRUE,
  seed = 16
)
traceplot(out_nig_norm_norm, name="rw")
traceplot(out_nig_norm_norm, name="spde1")

{ # make prediction
str(coord.prd.df@coords)
preds <- predict(out_nig_norm_norm, loc = list(
  rw = seaDist.prd,
  spde1 = coord.prd.df@coords,
  spde2 = coord.prd.df@coords
))

pred_df <- data.frame(
  x1 = coord.prd[, 1],
  x2 = coord.prd[, 2],
  preds = preds
)
ggplot(pred_df, aes(x = x1, y = x2, fill = preds)) +
  geom_raster() +
  scale_fill_viridis()
}
}

################ case 3. spde_nig + spde_normal + normal measurement nosie
{
out_normnig_norm <- ngme(
  formula = Y ~ 1 +
    f(inla.group(seaDist), model = "rw1", noise = noise_normal(), name="rw") +
    f(model = matern_spde, noise = noise_normal_nig(), name="spde1"),
  data = list(
    Y = Y
  ),
  family = noise_normal(),
  control = ngme_control(
    estimation = T,
    iterations = 1000,
    n_slope_check = 4,
    stop_points = 10,
    std_lim = 0.1,
    n_parallel_chain = 4,
    print_check_info = TRUE
  ),
  # debug = TRUE,
  seed = 16
)
traceplot(out_normnig_norm, name="rw")
traceplot(out_normnig_norm, name="spde1")

{ # make prediction
str(coord.prd.df@coords)
preds <- predict(out_normnig_norm, loc = list(
  rw = seaDist.prd,
  spde1 = coord.prd.df@coords
))

pred_df <- data.frame(
  x1 = coord.prd[, 1],
  x2 = coord.prd[, 2],
  preds = preds
)
ggplot(pred_df, aes(x = x1, y = x2, fill = preds)) +
  geom_raster() +
  scale_fill_viridis()
}
}

# cross validation
cross_validation(out_normnig_norm, k=5, N=500)
