#
load("tests/data/argo_2018_10.Rdata")
argo_df$Long <- (argo_df$Long + 180) %% 360 - 180

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

# INLA build mesh
library(INLA)
argo_al <- argo_df[argo_idx_sel, ]
coords <- with(argo_al, cbind(Long, Lat))
mesh <- inla.mesh.2d(coords, max.edge=c(4,10), cutoff=1)
mesh$n
plot(mesh)
points(coords)


##### ngme fit
library(ngme2)
myspde <- model_matern(mesh=mesh, loc=coords)
str(argo_al)

load_all()
# B_mu <- cbind(1, mesh$loc[, 2])
m0 <- temp_10 ~ 1 + f(model=myspde, noise=noise_nig(
  # theta_mu = c(0, 0),
  # B_mu = B_mu
))

res0 <- ngme(
  formula = m0,
  data = argo_al,
  control = ngme_control(
    iterations = 1000
  )
)

res0
traceplot(res0, "field1")

