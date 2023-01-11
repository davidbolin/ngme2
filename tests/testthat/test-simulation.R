# This file is test for simulate functions

test_that("simulate noises", {
  out <- simulate(noise_nig(mu=1, sigma = 0.5, nu=1, n=10))
  expect_equal(length(out), 10)
})

test_that("simulate f models", {
  # 1. AR(1)
  out <- simulate(model_ar1(
    1:100,
    noise = noise_nig(mu=1, sigma=2, nu=0.5)
  ))
  expect_equal(length(out), 100)

  # 2. Matern with normal_nig noise (working on)
  # First we create mesh
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh <- INLA::inla.mesh.2d(
    loc.domain = pl01, cutoff = 0.3,
    max.edge = c(0.2, 0.7), offset = c(0.5, 1.5)
  )

  W <- simulate(
    f(model = model_matern(mesh = mesh, kappa = 1),
      noise = noise_normal()
    )
  )

})


