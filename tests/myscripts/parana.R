## Paraná dataset

library(INLA)
library(splancs)
library(lattice)
library(ggplot2)
library(grid)
library(gridExtra)
library(viridis)



data(PRprec)
data(PRborder)

# Create INLA mesh
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
  ]))

library(ngme2)
# leave 0.1 Y as prediction area

# n <- length(Y)
# ind_pred <- sample(1:n, size = 0.1 * n)
# Y_pred <- Y[ind_pred]
# Y[ind_pred] <- NA

matern_spde <- model_matern(
  loc = coords,
  mesh = prmesh
)

out <- ngme(
  formula = Y ~ 1 +
    f(inla.group(seaDist), model = "rw1", noise = noise_normal()) +
    f(model = matern_spde, noise = noise_nig()) +
    f(model = matern_spde, noise = noise_normal()),
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
    print_check_info = TRUE
  ),
  seed = 416
)

out

# traceplots
## merr, beta
traceplot(out, f_index = 0, param = "sigma")
traceplot(out, f_index = 0, param = "beta")

## rw
traceplot(out, f_index = 1, param = "mu")
traceplot(out, f_index = 1, param = "sigma")
traceplot(out, f_index = 1, param = "nu")

## spde
traceplot(out, f_index = 2, param = "kappa")
traceplot(out, f_index = 2, param = "mu")
traceplot(out, f_index = 2, param = "sigma")
traceplot(out, f_index = 2, param = "nu")

## spde
traceplot(out, f_index = 3, param = "kappa")
traceplot(out, f_index = 3, param = "sigma")

### Prediction

nxy <- c(150, 100)
projgrid <- rSPDE::rspde.mesh.projector(prmesh,
  xlim = range(PRborder[, 1]),
  ylim = range(PRborder[, 2]), dims = nxy
)

xy.in <- inout(projgrid$lattice$loc, cbind(PRborder[, 1], PRborder[, 2]))

coord.prd <- projgrid$lattice$loc[xy.in, ]
plot(coord.prd, type = "p", cex = 0.1)
lines(PRborder)
points(coords[, 1], coords[, 2], pch = 19, cex = 0.5, col = "red")

seaDist.prd <- apply(spDists(coord.prd,
  PRborder[1034:1078, ],
  longlat = TRUE
), 1, min)

# construct data
n_prd <- nrow(coord.prd)
Y2 <- c(Y, rep(NA, n_prd))
seaDist2 <- c(seaDist, seaDist.prd)

matern_spde <- model_matern(
  loc = rbind(coords, coord.prd),
  mesh = prmesh,
  index_NA = is.na(Y2)
)

out2 <- ngme(
  formula = Y ~ 1 +
    f(inla.group(seaDist), model = "rw1", noise = noise_normal()) +
    f(model = matern_spde, noise = noise_nig()) +
    f(model = matern_spde, noise = noise_normal()),
  data = list(
    Y = Y2,
    seaDist = seaDist2
  ),
  family = noise_nig(),
  control = ngme_control(
    estimation = F,
    iterations = 1000,
    n_slope_check = 4,
    stop_points = 10,
    std_lim = 0.1,
    n_parallel_chain = 8,
    print_check_info = TRUE
  ),
  seed = 416,
  start = out
)
out2

# Plot prediction
lp <- attr(out2, "prediction")$lp
ggplot() +
  geom_point(aes(
    x = coord.prd[, 1], y = coord.prd[, 2],
    colour = lp[is.na(Y2)]
  ), size = 2, alpha = 1) +
  geom_point(aes(
    x = coords[, 1], y = coords[, 2],
    colour = Y_mean
  ), size = 2, alpha = 1) +
  scale_color_gradientn(colours = viridis(100)) +
  geom_path(aes(x = PRborder[, 1], y = PRborder[, 2])) +
  geom_path(aes(x = PRborder[1034:1078, 1], y = PRborder[
    1034:1078,
    2
  ]), colour = "red")


# Model estimation:
res1 <- data.frame(
  intercept    = format(out$beta, digits=3),
  noise_mu     = format(out$noise$theta_mu, digits=3),
  noise_sigma  = format(exp(out$noise$theta_sigma), digits=3),
  noise_nu     = format(out$noise$theta_V, digits=3),
  rw_sigma     = format(out$latents[[1]]$noise$theta_sigma, digits=3),
  ma_kappa     = format(exp(out$latents[[2]]$theta_K), digits=3),
  ma_mu        = format(out$latents[[2]]$noise$theta_mu, digits=3),
  ma_sigma     = format(exp(out$latents[[2]]$noise$theta_sigma), digits=3),
  ma_nu        = format(out$latents[[2]]$noise$theta_V, digits=3)
)
knitr::kable(res1, caption = "Estimations for the model")

# Result of the optimization trajectory of parameters for the Matern model:

pl <- lapply(c("kappa", "mu", "sigma", "nu"), function(.x)
  traceplot(out, param = .x, f_index = 2));
marrangeGrob(pl, nrow=2, ncol=2)