# this script is try to fit the Argo data at North Atlantic Ocean
# large-scale mapping using global sinusoidal or Mollweide projection.
load("tests/data/argo_2018_10.Rdata")
argo_al$JulDay <- NULL
nrow(argo_df)
argo_df$Long <- (argo_df$Long + 180) %% 360 - 180
names(argo_df)

library(maps)
library(ggplot2)
library(splancs)
library(viridis)

# North Atlantic Ocean
world_map <- ggplot2::map_data("world")

sel_lon = c(-80, -6, -6, -80)
sel_lat = c(50, 50, 0, 0)
sel_poly <- cbind(sel_lon, sel_lat)
argo_idx_sel <- inout(with(argo_df, cbind(Long, Lat)), sel_poly)
idx <- inout(with(world_map, cbind(long, lat)), sel_poly)

# Plot
gg <- ggplot(world_map[idx, ], aes(x = long, y = lat)) +
  geom_polygon(aes(group=group), fill = "grey70", color = "white") +
  geom_point(
    data = as.data.frame(argo_df[argo_idx_sel,]),
    aes(x = Long, y = Lat, col = temp_10), size = 0.7
  ) +
  scale_color_gradientn(colours = viridis(100)) +
  coord_fixed(1.3) +
  theme_void()
print(gg)

# Project longlat -> x and y
library(sp)
coordinates(argo_df) <- ~ Long + Lat
proj4string(argo_df) <- CRS("+proj=longlat +datum=WGS84")
moll_crs <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=km +no_defs"
argo_df_xy <- spTransform(argo_df, CRS(moll_crs))
head(argo_df_xy)

# INLA build mesh
library(INLA)
argo_al <- argo_df_xy[argo_idx_sel, ]
mesh <- inla.mesh.2d(coordinates(argo_al),
  max.n = 1000, cutoff=100, max.edge = c(200, 1000))
mesh$n
plot(mesh)
points(coordinates(argo_al))
head(argo_al)
str(argo_al)
# plot salinity + temperature
gtemp <- ggplot(data = as.data.frame(argo_al)) +
  geom_point(aes(x=Long, y=Lat, col=temp_10)) +
  scale_color_gradientn(colours = viridis(100))

gpsal <- ggplot(data = as.data.frame(argo_al)) +
  geom_point(aes(x=Long, y=Lat, col=psal_10)) +
  scale_color_gradientn(colours = viridis(100))

library(gridExtra)
grid.arrange(gtemp, gpsal, ncol = 2)

##### ngme fit
library(ngme2)
myspde <- model_matern(
  mesh=mesh,
  loc=coordinates(argo_al),
  name="spde"
)

############################## Model 1
m0 <- y ~ 1 + lat +
# f(lat, model="rw1", noise=noise_normal()) +
  f(model=myspde, noise=noise_nig())

res0 <- ngme(
  formula = m0,
  data = list(y = argo_al$temp_10,
    lat = scale(argo_al@coords[, 2])),
  control = control_opt(
    iterations = 1000,
    stepsize = 0.5,
    n_parallel_chain = 4
  ),
  family = "normal"
)

res0
traceplot(res0, "mn")
traceplot(res0, "beta")
traceplot(res0, "spde")

w <- res0$models[["field1"]]$W
length(w)
plot(w)


############################## Model 2
# 1. interpolate salinity first
mod_sal <- ngme(
  formula = sal ~ 1 + f(model=myspde, noise=noise_nig()),
  data = list(
    sal = argo_al$psal_10
  ),
  control = control_opt(
    iterations = 1000,
    n_parallel_chain = 4
  )
)

mod_sal
sal_intp <- predict(mod_sal, data=rep(1, mesh$n), loc=list(myspde=mesh$loc[, 1:2]))

# 2. fit the non-stationary model with intepolated salinity

# non-stationary mu
m1 <- y ~ 1 + f(model=myspde, noise=noise_nig(
  theta_mu = c(0, 0),
  # B_mu = cbind(1, scale(mesh$loc[, 2]))
  B_mu = cbind(1, sal_intp)
))

mod_temp <- ngme(
  formula = m1,
  data = list(y = argo_al$temp_10),
  control = control_opt(
    iterations = 1000,
    n_parallel_chain = 4
  ),
  family = "normal"
)

mod_temp
traceplot(mod_temp, "mn")
traceplot(mod_temp, "spde")

traceplot2(mod_temp, "spde", "mu", 2)

p1 <- traceplot2(mod_temp, "spde", "sigma")
p1 + ylim(0, 0.1)
#
m1 <- y ~ 1 + f(model=myspde, noise=noise_nig(
  theta_mu = c(0, 0),
  # B_mu = cbind(1, scale(mesh$loc[, 2]))
  B_mu = cbind(1, sal_intp-mean(sal_intp))
))

mod_temp <- ngme(
  formula = m1,
  data = list(y = argo_al$temp_10),
  control = control_opt(
    iterations = 1000,
    n_parallel_chain = 4,
    max_relative_step = 20,
    max_absolute_step = 20
  ),
  family = "normal"
)

mod_temp
traceplot(mod_temp, "mn")
traceplot(mod_temp, "spde")
