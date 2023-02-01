# This script is used for fitting different matern models using Parana data
# full version as in ngme2.Rmd

############################### Create INLA mesh
{
  load_all()
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
}
############################### Matern model and predict location
{
# matern model
matern_spde <- model_matern(
  loc = coords,
  mesh = prmesh
)

# make prediction at locations
nxy <- c(150, 100)
projgrid <- inla.mesh.projector(
  prmesh,
  xlim = range(PRborder[, 1]),
  ylim = range(PRborder[, 2]),
  dims = nxy
)

xy.in <- inout(projgrid$lattice$loc, cbind(PRborder[, 1], PRborder[, 2]))

coord.prd <- projgrid$lattice$loc[xy.in, ]
plot(coord.prd, type = "p", cex = 0.1)
lines(PRborder)
points(coords[, 1], coords[, 2], pch = 19, cex = 0.5, col = "red")

coord.prd.df <- data.frame(x1 = coord.prd[,1],
                          x2 = coord.prd[,2])
coordinates(coord.prd.df) <- c("x1", "x2")

# compute the distance to the shore
seaDist.prd <- apply(spDists(coord.prd,
  PRborder[1034:1078, ],
  longlat = TRUE
), 1, min)
coord.prd.df$seaDist <- seaDist.prd
}

################ case 1. spde_nig + nig measurement nosie
{
  # load_all()
out_nig_nig <- ngme(
  formula = Y ~ 1 +
    f(inla.group(seaDist), model = "rw1", noise = noise_normal(), name="rw") +
    f(model = matern_spde, noise = noise_nig(), name="spde1"),
  data = list(
    Y = Y
  ),
  family = noise_nig(),
  control = ngme_control(
    estimation = T,
    iterations = 1000,
    n_slope_check = 4,
    stop_points = 10,
    std_lim = 0.1,
    n_parallel_chain = 4,
    print_check_info = FALSE
  ),
  debug = TRUE,
  seed = 16
)
traceplot(out_nig_nig, "spde1")
# model CV
cross_validation(out_nig_nig, k = 5, N = 100)
# traceplot(out_nig_nig, name="rw")
# traceplot(out_nig_nig, name="spde1")
{ # make prediction
preds <- predict(out_nig_nig, loc = list(
  rw_loc = seaDist.prd,
  spde_loc = coord.prd.df
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

################ case 2. spde_nig + normal measurement nosie
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
    print_check_info = FALSE
  ),
  # debug = TRUE,
  seed = 16
)
traceplot(out_nig_norm, name="rw")
traceplot(out_nig_norm, name="spde1")
{ # make prediction
preds <- predict(out_nig_norm, loc = list(
  rw_loc = seaDist.prd,
  spde_loc = coord.prd.df
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

################ case 3. spde_normal + spde_nig (share same parameter) + normal measurement nosie
{
  load_all()
out_nignorm_norm <- ngme(
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
    print_check_info = FALSE
  ),
  # debug = TRUE,
  seed = 16
)
# traceplot(out_nignorm_norm, name="rw")
# traceplot(out_nignorm_norm, name="spde1")
{ # make prediction
preds <- predict(out_nignorm_norm, loc = list(
  rw_loc = seaDist.prd,
  spde_loc = coord.prd.df
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


################ case 4. spde_normal + spde_nig (diff kappa parameter) + normal measurement nosie
{
out_nig_norm_norm <- ngme(
  formula = Y ~ 1 +
    f(inla.group(seaDist), model = "rw1", noise = noise_normal(), name="rw") +
    f(model = matern_spde, noise = noise_nig(), name = "spde1") +
    f(model = matern_spde, noise = noise_normal(), name = "spde2"),
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
    print_check_info = FALSE
  ),
  debug = TRUE,
  seed = 16
)
out_nig_norm_norm
# traceplot(out_nig_norm_norm, name="rw")
# traceplot(out_nig_norm_norm, name="spde1")
{ # make prediction
preds <- predict(out_nig_norm_norm, loc = list(
  rw_loc = seaDist.prd,
  spde_loc = coord.prd.df,
  spde_loc = coord.prd.df
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
sampling_cpp(out_nig_nig, n=10, TRUE)
