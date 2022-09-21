######## test NA with matern ########
library(devtools); library(INLA); load_all()

# test ar1












test_matern <- function() {
  { # First we create mesh
    pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
    mesh <- inla.mesh.2d(
      loc.domain = pl01, cutoff = 1,
      max.edge = c(0.3, 1), offset = c(0.5, 1.5)
    )
    plot(mesh)

    W <- ngme.simulate(
      f(model = ngme.matern(mesh = mesh, theta_kappa = c(0.4, 1)),
        noise = ngme.noise.nig()
      )
    )
  }

  # generate A and A_pred
  n_obs <- 100; index_obs <- sample(1:mesh$n, n_obs)
  loc_obs <- mesh$loc[index_obs, c(1, 2)]
  A <- inla.spde.make.A(mesh = mesh, loc = loc_obs)
  Y_obs <- drop(A %*% W + rnorm(n_obs, sd = 0.4))

  n_NA <- 30; index_NA <- sample(1:mesh$n, n_NA)
  loc_NA <- mesh$loc[index_obs, c(1, 2)]
  A_pred <- inla.spde.make.A(mesh = mesh, loc = loc_NA)

  Y <- c(Y_obs, rep(NA, n_NA))

  # finally
  ngme_out <- ngme(
    formula = Y ~ 0 + f(
      model = "ar1",
      theta_K = 0.4,
      noise = ngme.noise.nig(),
      A = A,
      A_pred = A_pred,
      debug = TRUE
    ),
    data = data.frame(Y = Y),
    noise = ngme.noise.normal(),
    control = ngme.control(
      iteration = 2
    ),
    debug = ngme.debug(not_run = F)
  )
  str(ngme_out)

  index <- c(1,2,3)
  index_NA <- 2

  Filter(index, index %in% index_NA)

  # model.frame(Y~0, na.action=NULL)
}


