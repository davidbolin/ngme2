library(devtools); library(INLA); load_all()
{ # First we create mesh
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh <- inla.mesh.2d(
    loc.domain = pl01, cutoff = 1,
    max.edge = c(0.2, 0.7), offset = c(0.5, 1.5)
  )

  B_kappa <- matrix(c(rep(1, mesh$n), rexp(mesh$n)), ncol = 2)
  W <- simulate(
    f(model = ngme.matern(
      mesh = mesh,
      B_kappa = B_kappa,
      theta_kappa = c(1.1, 0.7)),
      noise = noise_nig()
    )
  )
}
plot(mesh)

# generate A and A_pred
n_obs <- 100; index_obs <- sample(1:mesh$n, n_obs)
loc_obs <- mesh$loc[index_obs, c(1, 2)]
A <- inla.spde.make.A(mesh = mesh, loc = loc_obs)
sigma.e <- 0.7
Y <- drop(A %*% W + sigma.e * rnorm(n_obs))

load_all()
ngme_out <- ngme(
  Y ~ 0 + f(
    model = ngme.matern(
      mesh = mesh,
      theta_kappa = c(0.5, 0.5),
      B_kappa = B_kappa
    ),
    fix_theta_K = FALSE,
    # W = as.numeric(W),
    # fix_W = TRUE,
    noise = noise_nig(
      fix_theta_mu    = F,
      fix_theta_sigma = F,
      fix_nu     = F
    ),
    A = A,
    debug = TRUE,
    control = ngme.control.f(
      numer_grad = T,
      use_precond = F
    )
  ),
  data = list(Y = Y),
  noise = noise_normal(),
  control = ngme.control(
    estimation = T,
    iterations = 100,
    n_parallel_chain = 1
  ),
  debug = TRUE
)
ngme_out
str(ngme_out)

traceplot2(ngme_out, parameter = "theta_sigma", f_index = 0)

# matern model
traceplot2(ngme_out, parameter = "theta_K",     f_index = 1, param_index = 1)
traceplot2(ngme_out, parameter = "theta_K",     f_index = 1, param_index = 2)
traceplot2(ngme_out, parameter = "theta_mu",    f_index = 1)
traceplot2(ngme_out, parameter = "theta_sigma", f_index = 1)
traceplot2(ngme_out, parameter = "nu",     f_index = 1)

plot(noise_nig(
      theta_mu = 0,
      theta_sigma = 0,
      nu = 1
    ), add = FALSE)
plot(ngme_out$latents[[1]]$noise, col = "red", add=TRUE)