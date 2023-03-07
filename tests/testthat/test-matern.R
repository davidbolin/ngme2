# test spatial Matern model

# 1. test estimation of matern model
# 2. test predict(out, loc=new_loc)

test_that("test Matern", {
  # load_all()
  library(INLA)
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh <- inla.mesh.2d(
    loc.domain = pl01, cutoff = 0.3,
    max.n = 100
  )
  # mesh$n
  # # generate A and A_pred
  # n_obs <- 500; index_obs <- sample(1:mesh$n, n_obs)
  # loc_obs <- mesh$loc[index_obs, c(1, 2)]
  # A <- inla.spde.make.A(mesh = mesh, loc = loc_obs)
  # sigma.e <- 0.7
  # Y <- drop(A %*% W + sigma.e * rnorm(n_obs))

  true_model <- model_matern(
    mesh = mesh,
    kappa = 3,
    noise = noise_nig(
      mu=3, sigma=1.5, nu=2
    )
  )

  # plot(mesh)
  W <<- simulate(true_model)

  n_obs <<- mesh$n
  Y <- W + rnorm(n_obs, sd=0.5)

  # make bubble plot
  sp_obj <- as.data.frame(mesh$loc); sp_obj[, 3] <- W
  names(sp_obj) <- c("s1", "s2", "y")
  coordinates(sp_obj) <- ~ s1 + s2
  bubble(sp_obj, zcol=3)
  range(mesh$loc[, 1]); range(mesh$loc[, 2])

  spde1 <<- model_matern(mesh=mesh, noise=noise_nig())
  Y2 <- Y
  # Y2[1:100] <- NA

# load_all()
  out <- ngme(
    Y ~ 0 + f(
      model=spde1,
      name="spde",
      noise=noise_nig(
        # fix_V = TRUE, V = attr(W, "noise")$V
      ),
      # fix_W = TRUE, W = W,
      debug = FALSE
    ),
    data = list(Y = Y2),
    control = ngme_control(
      estimation = T,
      post_samples_size = 50,
      iterations = 50,
      n_parallel_chain = 4,
      print_check_info = F
    ),
    debug = F
  )
  out
  traceplot(out, "spde")
  plot(attr(W, "noise"), out$latents[[1]]$noise)

  # Now let's do some prediction
  new_xs <- c(3, 5, 7)
  new_ys <- c(3, 5, 3)
  coo <- matrix(c(new_xs, new_ys), ncol=2)
  predict(out, loc=coo)
  # plot(mesh)
  # points(x=new_xs, y = new_ys, type = "p", col="red", pch=16, cex=2)
  expect_true(TRUE)
})

# test more
test_that("test matern ns", {
  # library(devtools); library(INLA); load_all()
  { # First we create mesh
    pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
    mesh <- inla.mesh.2d(
      loc.domain = pl01,
      max.edge = c(1)
    )
    B_kappa <- matrix(c(rep(1, mesh$n), rexp(mesh$n)), ncol = 2)
    W <- simulate(
      f(model = model_matern(
        mesh = mesh,
        B_kappa = B_kappa,
        theta_kappa = c(1.1, 0.7)),
        noise = noise_nig()
      )
    )
  }
  # plot(mesh)

  n_obs <- 100; index_obs <- sample(1:mesh$n, n_obs)
  loc_obs <- mesh$loc[index_obs, c(1, 2)]
  A <- inla.spde.make.A(mesh = mesh, loc = loc_obs)
  sigma.e <- 0.7
  Y <- drop(A %*% W + sigma.e * rnorm(n_obs))

  ngme_out <- ngme(
    Y ~ 0 + f(
      model = model_matern(
        loc = loc_obs,
        mesh = mesh,
        theta_kappa = c(0.5, 0.5),
        B_kappa = B_kappa
      ),
      fix_theta_K = FALSE,
      # W = as.numeric(W),
      # fix_W = TRUE,
      noise = noise_nig(),
      debug = F,
    ),
    data = list(Y = Y),
    family = noise_normal(),
    control = ngme_control(
      estimation = T,
      iterations = 100,
      n_parallel_chain = 1
    ),
    debug = FALSE
  )
  ngme_out
  expect_true(TRUE)
})

##################################################
test_that("test 1d matern with numerical g", {
  library(INLA)
  load_all()
  n_mesh <<- 800
  mesh <- inla.mesh.1d(1:n_mesh)
  n_obs <<- 500
  loc <- runif(n_obs, 1, n_mesh)
  real_noise <- noise_nig(mu = -5, sigma=3, nu=2, n=n_mesh)
  spde1d <- model_matern(loc = loc, mesh=mesh, noise=real_noise, kappa=2)

# simulate
  # W <- simulate(spde1d)
  eps <- simulate.ngme_noise(real_noise)
  K <- 4 * spde1d$C + spde1d$G
  W <- as.numeric(solve(K, eps))
  Y <- drop(spde1d$A %*% W) + rnorm(n_obs, sd=0.5)
plot(Y)

  out <- ngme(
    Y ~ 0 + f(model=spde1d, noise=noise_nig(),
      fix_W = T, W =W
    ),
    data = list(Y = Y),
    control = ngme_control(
      estimation = T,
      iterations = 500,
      n_parallel_chain = 4
    )
  )

  out
  traceplot(out, "field1")
  plot(out$latents[[1]]$noise, real_noise)

  expect_true(T)
})