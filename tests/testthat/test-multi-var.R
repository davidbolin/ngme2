test_that("Matern2d formula works", {
  # test on Matern 2d
  Y1 <- rnorm(5); Y2 <- rnorm(5)
  Y <- cbind(Y1, Y2); group <- rep(1:2, each=5)

  load_all()
  mod <- matern_nd(
    map = list(1:5, 1:5),
    names = list("sal", "temp"),
    mesh = INLA::inla.mesh.1d(1:10),
  )

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
