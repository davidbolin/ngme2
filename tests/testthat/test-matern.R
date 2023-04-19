# test spatial Matern model

# 1. test estimation of matern model
# 2. test predict(out, loc=new_loc)

test_that("test Matern", {
  library(INLA)
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh <- inla.mesh.2d(
    loc.domain = pl01, cutoff = 0.2,
    max.edge = c(0.5,10)
  )
# mesh$n
# plot(mesh)
  n_obs <- 500
  loc <- cbind(runif(n_obs, 0, 10), runif(n_obs, 0, 5))

  true_model <- model_matern(
    theta_kappa = log(3), mesh = mesh, map = loc,
    noise = noise_nig(
      mu=-2, sigma=1.5, nu=1, n = mesh$n
    )
  )
  W <- simulate(true_model)
  mean(W)
  mean(attr(W, "noise")$h)

  A <- inla.spde.make.A(loc=loc, mesh=mesh)
  Y <- as.numeric(A %*% W) + rnorm(n_obs, sd=0.5)

  # make bubble plot
  # sp_obj <- as.data.frame(mesh$loc); sp_obj[, 3] <- W
  # names(sp_obj) <- c("s1", "s2", "y")
  # coordinates(sp_obj) <- ~ s1 + s2
  # bubble(sp_obj, zcol=3)
  # range(mesh$loc[, 1]); range(mesh$loc[, 2])

# Matern case
  load_all()
  out <- ngme(
    Y ~ 0 + f(loc,
      model="matern",
      name="spde",
      mesh = mesh,
      noise=noise_nig(
        # fix_nu = T, nu=1
      ),
      control = control_f(numer_grad = T),
      # fix_theta_K = T, theta_kappa = log(3),
      fix_W = TRUE,
      W = W,
      debug = T
    ),
    data = data.frame(Y = Y),
    control_opt = control_opt(
      estimation = T,
      iterations = 1,
      n_parallel_chain = 4,
      print_check_info = F,
      verbose = F,
      std_lim = 0.001
    ),
    debug = F
  )
  # out$replicates[[1]]$latents[[1]]$noise$h == out$replicates[[1]]$latents[[1]]$h
  out
  traceplot(out, "spde")
  traceplot(out)
  plot(noise_nig(mu=-2,sigma=1.5,nu=1),
    out$replicates[[1]]$latents[[1]]$noise)

  # Now let's do some prediction
  new_xs <- c(3, 5, 7)
  new_ys <- c(3, 5, 3)
  coo <- matrix(c(new_xs, new_ys), ncol=2)
  predict(out, loc=list(coo))
  # plot(mesh)
  # points(x=new_xs, y = new_ys, type = "p", col="red", pch=16, cex=2)
  expect_true(TRUE)
})

# test more
test_that("test matern ns", {
  # library(devtools); library(INLA); load_all()
  { # First we create mesh
    library(INLA)
    pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
    mesh <- inla.mesh.2d(
      loc.domain = pl01,
      max.edge = 0.5
    )
    # plot(mesh)

  n_obs <- 500; index_obs <- sample(1:mesh$n, n_obs)
  loc_obs <- mesh$loc[index_obs, c(1, 2)]
  A <- inla.spde.make.A(mesh = mesh, loc = loc_obs)
  sigma.e <- 0.01

  B_kappa <- matrix(c(rep(1, mesh$n), sin(1:mesh$n / 20)), ncol = 2)
  theta_kappa <- c(1.5, 2.5)
  W <- simulate(model_matern(loc_obs,
      mesh = mesh,
      B_kappa = B_kappa,
      theta_kappa = theta_kappa,
      # theta_kappa = 2,
      noise = noise_nig(mu=-3, sigma = 2, nu = 1.5)))
  # mean(attr(W, "noise")$h)
  Y <- as.numeric(A %*% W) + sigma.e * rnorm(n_obs)
  # points(loc_obs[,1], loc_obs[,2], col="red", pch=16, cex=2)
  }

  # plot(B_kappa %*% theta_kappa)

  #  { # 1d matern case
  #   sigma.e = 0.01
  #   loc_obs <- 1:500
  #   B_kappa <- matrix(c(rep(1, 500), sin(1:500 / 20)), ncol = 2)
  #   mesh = inla.mesh.1d(loc=1:500)
  #   W <- simulate(
  #       model_matern(loc_obs,
  #       mesh = inla.mesh.1d(loc=loc_obs),
  #       theta_kappa = 2,
  #       # B_kappa = B_kappa,
  #       # theta_kappa = c(1.1, 2),
  #       noise = noise_nig(mu=-5, sigma = 2, nu = 1.5)))
  #   A <- inla.spde.make.A(loc=loc_obs, mesh = inla.mesh.1d(loc=1:500))
  #   Y <- as.numeric(A %*% W) + sigma.e * rnorm(n_obs)
  # }
  plot(Y, col="red"); abline(h=0)

  trueV = attr(W,"noise")$V

  ngme_out <- ngme(
    Y ~ 0 + f(loc_obs,
      model="matern",
      mesh = mesh,
      theta_kappa = c(0.5, 0.5),
      B_kappa = B_kappa,
      # theta_kappa = 2,
      # fix_theta_K = T,
      # W = as.numeric(W), fix_W = TRUE,
      noise = noise_nig(
        # V = trueV, fix_V = TRUE
        # fix_nu = T, fix_theta_sigma=T, fix_theta_mu=T
      ),
      control = control_f(numer_grad = F),
      debug = F,
    ),
    data = data.frame(Y = Y),
    family = noise_normal(),
    control_opt = control_opt(
      estimation = T,
      iterations = 100,
      n_parallel_chain = 4
      # max_relative_step = 5,
      # max_absolute_step = 10
    ),
    debug = FALSE
  )
  ngme_out
  traceplot(ngme_out, "field1")
  ngme_out$replicates[[1]]$latents[[1]]$theta_K
  plot(B_kappa %*% theta_kappa, type="l")
  points(B_kappa %*% ngme_out$replicates[[1]]$latents[[1]]$theta_K)

# compare noise
  plot(ngme_out$replicates[[1]]$latents[[1]]$noise,  noise_nig(mu=-3, sigma = 2, nu = 1.5))

  expect_true(TRUE)
})

##################################################
test_that("test 1d matern with numerical g", {
  library(INLA)
  n_mesh <<- 800
  mesh <- inla.mesh.1d(1:n_mesh)
  n_obs <<- 500
  loc <- runif(n_obs, 1, n_mesh)
  real_noise <- noise_nig(mu = -5, sigma=3, nu=2, n=n_mesh)
  spde1d <- model_matern(map = loc, mesh=mesh, kappa=2,
    noise = noise_nig(mu = -5, sigma=3, nu=2, n=n_mesh))
# spde1d$noise$h

# simulate
  # W <- simulate(spde1d)
  eps <- simulate.ngme_noise(real_noise)
  K <- 4 * spde1d$C + spde1d$G
  W <- as.numeric(solve(K, eps))
  Y <- as.numeric(spde1d$A %*% W) + rnorm(n_obs, sd=0.5)
# plot(Y, type="l")

  out <- ngme(
    Y ~ 0 + f(loc, model="matern", mesh=mesh, noise=noise_nig(),
      # fix_W = T, W =W
    ),
    data = data.frame(Y = Y),
    control_opt = control_opt(
      estimation = F,
      iterations = 100,
      n_parallel_chain = 4
    )
  )
  out
# traceplot(out, "field1")
# plot(out$replicates[[1]]$latents[[1]]$noise, real_noise)

  expect_true(TRUE)
})
