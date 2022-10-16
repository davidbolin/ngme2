library(devtools); library(INLA); load_all()
{ # First we create mesh
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh <- inla.mesh.2d(
    loc.domain = pl01, cutoff = 0.3,
    max.edge = c(0.2, 0.7), offset = c(0.5, 1.5)
  )

  W <- simulate(
    f(model = ngme.matern(mesh = mesh, kappa = 1),
      noise = ngme.noise.nig()
    )
  )
}


# generate A and A_pred
n_obs <- 10; index_obs <- sample(1:mesh$n, n_obs)
loc_obs <- mesh$loc[index_obs, c(1, 2)]
A <- inla.spde.make.A(mesh = mesh, loc = loc_obs)

sigma.e <- 0.7
Y1 <- drop(A %*% W + sigma.e * rnorm(n_obs))
Y2 <- drop(A %*% W + sigma.e * rnorm(n_obs))

#bivaraite formula
#formula = (y1 | y2 ~ x1 + f(x1, model="SPDE2D", var="nig") | x2 + f(X2, model="SPDE2D", var="nig"))
fm <- Y1 | Y2 ~ 0 + f(
    model = ngme.matern2D(mesh = mesh, theta_kappa = 0.4, alpha = c(2, 2)),
    fix_theta_K = FALSE,
    # W = as.numeric(W),
    # fix_W = TRUE,

    noise = ngme.noise.nig(
      fix_theta_mu    = F,
      fix_theta_sigma = F,
      fix_theta_V     = F
    ),
    A = A,
    debug = TRUE,
    control = ngme.control.f(
      numer_grad = F,
      use_precond = F
    )
  ) | 0 + f(
    model = ngme.matern2D(mesh = mesh, theta_kappa = 0.4, alpha = c(2, 2)),
    fix_theta_K = FALSE,
    # W = as.numeric(W),
    # fix_W = TRUE,

    noise = ngme.noise.nig(
      fix_theta_mu    = F,
      fix_theta_sigma = F,
      fix_theta_V     = F
    ),
    A = A,
    debug = TRUE,
    control = ngme.control.f(
      numer_grad = F,
      use_precond = F
    )
  )


ngme_out <- ngme(fm,
  data = list(Y1 = Y1, Y2 = Y2),
  noise = ngme.noise.normal(),
  control = ngme.control(
    iterations = 500,
    estimation = TRUE,
    n_parallel_chain = 2
  ),
  debug = TRUE
)
ngme_out
# str(ngme_out)
# plot_chains(ngme_out, parameter = "theta_sigma", f_index = 0)

# # matern model
# plot_chains(ngme_out, parameter = "theta_K",     f_index = 1)
# plot_chains(ngme_out, parameter = "theta_mu",    f_index = 1)
# plot_chains(ngme_out, parameter = "theta_sigma", f_index = 1)
# plot_chains(ngme_out, parameter = "theta_V",     f_index = 1)

# plot(ngme.noise.nig(
#       theta_mu = 0,
#       theta_sigma = 0,
#       theta_V = 1
#     ), add = FALSE)
# plot(ngme_out$latents[[1]]$noise, col = "red", add=TRUE)