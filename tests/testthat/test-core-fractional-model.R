test_that("test fractional model", {
  skip_if_not_installed("fmesher")
  skip_if_not_installed("rSPDE")
  skip_if_not_installed("inlabru")
  skip_if_not_installed("sf")
  library(fmesher)
  library(rSPDE)
  n_loc <- 1000
  loc_2d_mesh <- matrix(runif(n_loc * 2), n_loc, 2)
  mesh_2d <- fm_mesh_2d(
    loc = loc_2d_mesh,
    cutoff = 0.05,
    max.edge = c(0.1, 0.5)
  )
  plot(mesh_2d)
  mesh_2d$n

  sigma <- 2
  range <- 0.25
  nu <- 1.3
  kappa <- sqrt(8 * nu) / range
  op <- matern.operators(
    mesh = mesh_2d, nu = nu,
    range = range, sigma = sigma,
    m = 1,
    parameterization = "matern"
  )
  tau <- op$tau

  u <- simulate(op)
  A <- spde.make.A(
    mesh = mesh_2d,
    loc = loc_2d_mesh
  )
  sigma.e <- 0.1
  y <- A %*% u + rnorm(n_loc) * sigma.e

  ############ Fitting with ngme2 ############
  # library(ngme2)
  toy_df <- data.frame(
    coord1 = loc_2d_mesh[, 1],
    coord2 = loc_2d_mesh[, 2],
    y = as.vector(y)
  )
  ngme_fit <- ngme(
    y ~ -1 + f(
      ~ coord1 + coord2,
      model = matern(
        mesh = mesh_2d,
        kappa = 12,
        alpha = 2.3,
        rational_order = 1,
        fix_alpha = FALSE
      ),
      noise = noise_normal(sigma = 200, fix_theta_sigma = FALSE)
    ),
    data = toy_df,
    family = noise_normal(sigma = 0.1, fix_sigma = FALSE),
    # debug = TRUE,
    control_opt = control_opt(
      n_parallel_chain = 4,
      start_sd = 0.01,
      optimizer = precond_sgd(),
      iterations = 200,
      robust = TRUE,
      verbose = TRUE
    )
    # start = ngme_fit
  )

  summary(ngme_fit)
  traceplot(ngme_fit)
  traceplot(ngme_fit, "field1")

  ############ Fitting with INLAbru ############
  library(inlabru)
  rspde_model <- rspde.matern(
    mesh = mesh_2d,
    nu.upper.bound = 2,
    parameterization = "spde"
  )

  library(sf)
  toy_df <- data.frame(
    coord1 = loc_2d_mesh[, 1],
    coord2 = loc_2d_mesh[, 2],
    y = as.vector(y)
  )
  toy_df <- st_as_sf(toy_df, coords = c("coord1", "coord2"))

  cmp <-
    y ~ -1 + field(geometry,
      model = rspde_model
    )

  rspde_bru_fit <-
    bru(cmp,
      data = toy_df,
      options = list(
        family = "gaussian",
        num.threads = "1:1"
      )
    )

  summary(rspde_bru_fit)

  result_fit <- rspde.result(rspde_bru_fit, "field", rspde_model)
  summary(result_fit)
  tau <- op$tau
  result_df <- data.frame(
    parameter = c("tau", "kappa", "nu"),
    true = c(tau, kappa, nu), mean = c(
      result_fit$summary.tau$mean,
      result_fit$summary.kappa$mean,
      result_fit$summary.nu$mean
    ),
    mode = c(
      result_fit$summary.tau$mode,
      result_fit$summary.kappa$mode,
      result_fit$summary.nu$mode
    )
  )
  print(result_df)
})
