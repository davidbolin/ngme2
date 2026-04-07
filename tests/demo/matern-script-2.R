library(fmesher)
library(rSPDE)
library(ngme2)

m <- 200
loc_2d_mesh <- matrix(runif(m * 2), m, 2)
mesh_2d <- fm_mesh_2d(
  loc = loc_2d_mesh,
  cutoff = 0.05,
  offset = c(0.1, 0.4),
  max.edge = c(0.05, 0.5)
)

plot(mesh_2d, main = "")
points(loc_2d_mesh[, 1], loc_2d_mesh[, 2])
nu <- 0.8
sigma <- 1.3
range <- 0.15
op <- matern.operators(
  range = range, sigma = sigma,
  nu = nu, m = 2, mesh = mesh_2d,
  parameterization = "matern"
)

n.rep <- 1
u <- simulate(op, nsim = n.rep)
A <- fm_basis(
  x = mesh_2d,
  loc = loc_2d_mesh
)
sigma.e <- 0.1
Y <- A %*% u + matrix(rnorm(m * n.rep), ncol = n.rep) * sigma.e

library(viridis)
library(ggplot2)
proj <- fm_evaluator(mesh_2d, dims = c(70, 70))

df_field <- data.frame(
  x = proj$lattice$loc[, 1],
  y = proj$lattice$loc[, 2],
  field = as.vector(fm_evaluate(proj,
    field = as.vector(u[, 1])
  )),
  type = "field"
)

df_loc <- data.frame(
  x = loc_2d_mesh[, 1],
  y = loc_2d_mesh[, 2],
  field = as.vector(Y[, 1]),
  type = "locations"
)
df_plot <- rbind(df_field, df_loc)

ggplot(df_plot) +
  aes(x = x, y = y, fill = field) +
  facet_wrap(~type) +
  xlim(0, 1) +
  ylim(0, 1) +
  geom_raster(data = df_field) +
  geom_point(
    data = df_loc, aes(colour = field),
    show.legend = FALSE
  ) +
  scale_fill_viridis() +
  scale_colour_viridis()

op_obj <- matern.operators(
  m = 1,
  type = "operator", mesh = mesh_2d
)

y_vec <- as.vector(Y)
repl <- rep(1:n.rep, each = m)
df_data_2d <- data.frame(
  y = y_vec, x_coord = loc_2d_mesh[, 1],
  y_coord = loc_2d_mesh[, 2]
)

# fit_2d <- rspde_lme(y ~ -1, model = op_obj,
#           data = df_data_2d, repl = repl,
#           loc = c("x_coord", "y_coord"),
#           parallel = TRUE)
# glance(fit_2d)
# print(data.frame(
#   sigma = c(sigma, fit_2d$alt_par_coeff$coeff["sigma"]),
#   range = c(range, fit_2d$alt_par_coeff$coeff["range"]),
#   nu = c(nu, fit_2d$alt_par_coeff$coeff["nu"]),
#   row.names = c("Truth", "Estimates")
# ))

# # Total time
# print(fit_2d$fitting_time)

fit <- ngme(
  y ~ 0 + f(
    ~ x_coord + y_coord,
    model = "matern",
    name = "spde",
    kappa = 20,
    mesh = mesh_2d,
    fix_alpha = FALSE,
    rational_order = 2,
    noise = noise_normal(sigma = 5)
  ),
  data = df_data_2d,
  replicate = repl,
  control_opt = control_opt(
    estimation = TRUE,
    iterations = 300,
    optimizer = adam(),
    n_batch = 10,
    rao_blackwellization = F,
    n_parallel_chain = 4,
    print_check_info = F,
    verbose = T,
    std_lim = 0.01,
    start_sd = 0.01
  )
)
fit
traceplot(fit, "spde")
traceplot(fit)
