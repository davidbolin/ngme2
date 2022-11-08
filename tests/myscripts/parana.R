{ # read data
  set.seed(10)
  library(fields)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(INLA)
  library(devtools); load_all()

  data(PRprec)
  data(PRborder)
  coords <- as.matrix(PRprec[, 1:2])

  Y <- rowMeans(PRprec[, 12 + 1:31]) # 2 + Octobor
  # Y2 <- apply(PRprec[, 12 + 1:31], 1, max) # 2 + Octobor
  # remove NA
  ind <- !is.na(Y)
  coords <- as.matrix(PRprec[ind, 1:2])
  Y <- Y_mean <- Y[ind]

  # Covariates
  alt <- PRprec$Altitude[ind]
  seaDist <- apply(spDists(coords, PRborder[1034:1078, ],
    longlat = TRUE
  ), 1, min)

  prdomain <- inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(100, 100))
  prmesh <- inla.mesh.2d(boundary = prdomain, max.edge = c(0.45, 1), cutoff = 0.2)
}

{ # plot data
  ggplot() +
  geom_point(aes(
    x = coords[, 1], y = coords[, 2],
    colour = Y
  ), size = 2, alpha = 1) +
  scale_color_gradientn(colours = tim.colors(100)) +
  geom_path(aes(x = PRborder[, 1], y = PRborder[, 2])) +
  geom_path(aes(x = PRborder[1034:1078, 1], y = PRborder[
    1034:1078,
    2
  ]), colour = "red")
}

# leave the set for prediction
n <- length(Y)
ind_pred <- sample(1:n, size = 0.1 * n)
A <- inla.spde.make.A(mesh = prmesh, loc = coords[-ind_pred, ])
A_pred <- inla.spde.make.A(mesh = prmesh, loc = coords[ind_pred, ])
Y_pred <- Y[ind_pred]; Y[ind_pred] <- NA

mesh.index <- inla.spde.make.index(
  name = "field",
  mesh = prmesh,
  n.spde = prmesh$n
)

# plot(prmesh)
# points(coords)
load_all()
out <- ngme(
  formula = Y ~ 1 +
    # f(inla.group(seaDist), model = "rw1", noise = noise_normal()) +
    # f(seaDist, model = "rw1", noise = noise_normal()) +
    # f(seaDist, model = "ar1", noise = noise_normal()) +
    f(index = mesh.index$field,
      model = model_matern(A = A, A_pred = A_pred, mesh = prmesh, noise = noise_nig())
    ),
  data =  data.frame(
    Y  = Y,
    # Y_max   = Y2,
    long    = coords[, 1],
    lat     = coords[, 2]
  ),
  family = noise_nig(),
  control = ngme_control(
    estimation = TRUE,
    iterations = 1000,
    n_slope_check = 4,
    stop_points = 20,
    n_parallel_chain = 8
  )
)

# Comparing our prediction
mean(Y_mean)
mean(abs(attr(out, "prediction")$lp - Y_mean))

out
str(out)

# plots
# beta
# traceplot(out, parameter = "beta",    f_index = 0, param_index = 1)
# traceplot(out, parameter = "beta",    f_index = 0, param_index = 2)



# 1st model
traceplot(out, parameter = "theta_K",     f_index = 1)
traceplot(out, parameter = "theta_mu",    f_index = 1)
traceplot(out, parameter = "theta_sigma", f_index = 1)
traceplot(out, parameter = "theta_V",     f_index = 1)

# 2nd model
traceplot(out, parameter = "theta_K",     f_index = 2)
traceplot(out, parameter = "theta_mu",    f_index = 2)
traceplot(out, parameter = "theta_sigma", f_index = 2)
traceplot(out, parameter = "theta_V",     f_index = 2)

# 3rd model
traceplot(out, parameter = "theta_K",     f_index = 3)
traceplot(out, parameter = "theta_mu",    f_index = 3)
traceplot(out, parameter = "theta_sigma", f_index = 3)
traceplot(out, parameter = "theta_V",     f_index = 3)


# fixed effects
traceplot(out, parameter = "beta", f_index = 0, param_index = 1)

# measurement noise
traceplot(out, parameter = "theta_mu", f_index = 0)
traceplot(out, parameter = "theta_sigma", f_index = 0)
traceplot(out, parameter = "theta_V", f_index = 0)


# plot(1:3, 1:3)
# library(grid)
# library(gridExtra)
# grid.arrange(rectGrob(), rectGrob())
# ## Not run:
# library(ggplot2)
# pl <- lapply(1:11, function(.x) qplot(1:10, rnorm(10), main=paste("plot", .x)))
# ml <- marrangeGrob(pl, nrow=2, ncol=2)
# ml
# ## non-interactive use, multipage pdf
# ggsave("multipage.pdf", ml)
# ## interactive use; open new devices
# ml

# library(pagedown)
# chrome_print("http://127.0.0.1:3000/myscripts/poster.html",
#   format="pdf")

# #  -------
# load_all()
# out_4k <- ngme(
#   formula = Y_mean ~ 1 +
#     f(inla.group(seaDist), model = "rw1", noise=noise_normal()) +
#     f(index = mesh.index$field,
#       model = model_matern(A = A, mesh = prmesh),
#       noise = noise_nig()
#     ),
#   data = data,
#   family = noise_nig(),
#   control = ngme_control(
#     estimation = T,
#     iterations = 2000,
#     n_slope_check = 4,
#     stop_points = 20,
#     n_parallel_chain = 8
#   ),
#   seed = 5,
#   start = out_2k
# )

# library(grid)
# library(gridExtra)
# pl <- lapply(c("kappa", "mu", "sigma", "nu"), function(.x)
#   traceplot(out_4k, parameter = .x, f_index = 2));
# marrangeGrob(pl, nrow = 2, ncol = 2)

# xxx <- seq(1, 100, length=100)
# plot(xxx, xxx^0.95)


# (1 / (1:10)) ** 0.01


# 1 ** 0.1
# ??knitr


# out <- ngme(
#   formula = Y_mean ~ 1 +
#     f(inla.group(seaDist), model = "rw1", noise=noise_normal()) +
#     f(index = mesh.index$field,
#       model = model_matern(A = A, mesh = prmesh),
#       noise = noise_nig()
#     ),
#   data = data,
#   family = noise_nig(),
#   control = ngme_control(
#     estimation = T,
#     iterations = 2000,
#     n_slope_check = 4,
#     stop_points = 20,
#     n_parallel_chain = 8
#   ),
#   seed = 5
# )

# predict(
#   ngme = out,
#   formula = ~ 1 +
#     # f(inla.group(seaDist_pred), model = "rw1", noise = noise_normal()) +
#     f(seaDist_pred, model = "rw1", noise = noise_normal()) +
#     f(index = mesh.index$field,
#       model = model_matern(A_pred = A_pred, mesh = prmesh),
#       noise = noise_nig()
#     ),
#   data = list(seaDist_pred = seaDist_pred, A_pred = A_pred)
# )

# ?inla.spde.make.A()

# load_all()
# YY <- c(1.2, 2.3, 4.3, NA)

# XX <- c(1.1, 3.1, 2.2, 4.5)
# mesh1 <- inla.mesh.1d(XX); mesh1
# inla.spde.make.A(mesh=mesh1, loc=c(1.7)) %*% mesh1$loc

# ngme(
#   YY ~ f(XX, model = "rw1"),
#   data = data.frame(YY=YY, XX=XX),
#   control = ngme_control(estimation = FALSE)
# )

