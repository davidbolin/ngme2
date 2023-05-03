# Y1 Y2
# group = factor()
# f(~x+y, model="matern")
# coordinates()


test_that("test bv R interface", {
  # test on Matern 2d
  load_all()
  Y1 <- rnorm(5); Y2 <- rnorm(5)
  Y <- cbind(Y1, Y2);
  bv(first = ar1(1:3), second = ar1(1:3))

load_all()

  n = 500
  true_model <- f(model="bv",
    first  = ar1(1:n, theta_K = ar1_a2th(0.8)),
    second = ar1(1:n, theta_K = ar1_a2th(0.2)),
    rho = 0.5, zeta = 1,
    noise = noise_nig(mu=-3, sigma=2, nu=1),
    eval=T
  )

  W <- simulate(true_model)
  AW <- as.numeric(true_model$A %*% W)
  n_obs <- length(AW)
  Y <- AW + rnorm(n_obs, sd=0.5)

out <- ngme(
  Y ~
  # 0 + f(1:5, model="ar1", group="temp") +
   0 + f (model="bv", name="bv",
    first=ar1(1:n),
    second=ar1(1:n),
    noise = noise_nig(),
    control = control_f(numer_grad = F)
  ),
  # noise = noise_nig(correlated = TRUE),
  data = data.frame(Y = Y
    # group=group, loc1=.., loc2=..
  ),
  control_ngme = control_ngme(
    n_gibbs_samples = 5
  ),
  control_opt = control_opt(
    iterations = 100,
    n_parallel_chain = 4,
    estimation = F,
    verbose = T
  )
)
out
traceplot(out, "bv")
out$replicates[[1]]$latents[[1]]$theta_K
out$replicates[[1]]$latents[[1]]$operator$theta_K

ar1_th2a(0.023)


# ngme2
  fm <- Y ~ 0 + f(~1, "fe", group=2) + f(~x1, "re", group=1) + f(
    model="matern2d",
    map=list(1:5, 2:6),
    rho = 0.5, theta = 1,
    kappa1 = 1, kappa2 =1,
    mesh = INLA::inla.mesh.1d(1:10),
    noise = list(noise_nig(), noise_nig()),
    group = c(1,2)
  )

# data =

  out <- ngme(
    fm,
    data = data.frame(
      Y = 1:10,
      x1 = c(1:5, rep(0, 5))
    ),
    group = c(rep(1,5), rep(1,5))
  )

})

test_that("Test on mutlivariate model", {
  library(INLA)
  mesh <- inla.mesh.1d(1:10)
  model_2d(
    a = matern(mesh=mesh),
    b = matern(mesh=mesh),

  )
})