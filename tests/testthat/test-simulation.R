# # This file is test for simulate functions

# test_that("simulate noises", {
#   out <- simulate(noise_nig(mu=1, sigma = 0.5, nu=1, n=10))
#   expect_equal(length(out), 10)
# })

# test_that("simulate f models", {
#   # 1. AR(1)
#   out <- simulate(model_ar1(
#     1:100,
#     noise = noise_nig(mu=1, sigma=2, nu=0.5)
#   ))
#   expect_equal(length(out), 100)

#   # 2. Matern with normal_nig noise (working on)
#   # First we create mesh
#   pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
#   mesh <- fmesher::fm_mesh_2d(
#     loc.domain = pl01, cutoff = 0.3,
#     max.edge = c(0.2, 0.7), offset = c(0.5, 1.5)
#   )

#   # W <- simulate(
#   #   f(model = model_matern(mesh = mesh, kappa = 1),
#   #     noise = noise_normal()
#   #   )
#   # )
# })


test_that("simulate ngme object (with repls)", {
  out = test_ngme("ar1", n_obs_per_rep = 10, n_replicate = 2)
  simulate(out$out)[[1]]

  out2 = ngme(
    y ~ 0+f(c(1:5, 1:5), model="ar1"),
    replicate = repls,
    data = data.frame(y=1:10),
    family = noise_normal()
  )
  out2
  simulate.ngme(out2, posterior = FALSE)
  out2$replicate[["A"]]$Y
})


test_that("simulate correlated noise", {
  out = ngme(
    y ~ 0+f(c(1:5), model="ar1"),
    data = data.frame(y=1:5),
    family = noise_nig(
      corr_measurement = TRUE,
      index_corr = c(1,1,2,3,3)
    ),
    control_opt = control_opt(
      estimation = F
    )
  )
  simulate(out)[[1]]
})