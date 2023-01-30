# test for replicate data

# share the same hyperparamter
test_that("test replicate on 1d data", {
  # load_all()
  rw <- model_rw(c(1:4, 3:4), replicate = c(rep(1,4), rep(2, 2)))
  expect_true(dim(rw$C)[[2]] == dim(rw$A)[[2]])
  expect_true(dim(rw$C)[[1]] == length(rw$h))

  ar1 <- model_ar1(c(1:4, 1:2), replicate = c(rep(1,4), rep(2, 2)))
  expect_true(dim(ar1$C)[[2]] == dim(ar1$A)[[2]])
  expect_true(dim(ar1$C)[[1]] == length(ar1$h))
})

test_that("test replicate on 2d data", {
  library(INLA)
  pl01 <- cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5)
  mesh <- inla.mesh.2d(
    loc.domain = pl01,
    max.n = 10
  )

  locs <- cbind(c(4,5,6,4,5), c(4,5,6,4,5))
  ma <- model_matern(loc=locs, mesh=mesh, replicate = c(1,1,1, 2,2))

  expect_true(dim(ma$C)[[2]] == dim(ma$A)[[2]])
  expect_true(dim(ma$C)[[1]] == length(ma$h))
})
