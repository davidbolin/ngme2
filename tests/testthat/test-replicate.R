# test for replicate data

# share the same hyperparamter
test_that("test replicate on 1d data", {
  # load_all()
  rw <- model_rw(1:5, replicate = c(1,1,1, 2,2))
  expect_true(ncol(rw$C) * 2 == ncol(rw$A))

  ar1 <- model_ar1(1:5, replicate = c(1,1,1, 2,2))
  expect_true(ncol(ar1$C) * 2 == ncol(ar1$A))

  matern <- model_matern(1:5, replicate = c(1,1,1, 2,2), mesh=1:10)
  expect_error(f(model=matern, replicate = c(1,1,1, 2,2)))
  expect_true(ncol(matern$C) * 2 == ncol(matern$A))
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
# dim(ma$A)
# dim(ma$C)
# ma$noise$h

  expect_true(ncol(ma$C) * 2 == ncol(ma$A))
})

test_that("test estimation on AR(1) and RW process", {
  # load_all()
  n <<- 500
  z1 = arima.sim(n, model = list(ar = 0.5), sd = 0.5)
  z2 = arima.sim(n, model = list(ar = 0.5), sd = 0.5)

  out1 <- ngme(
    z ~ f(1:n, model="ar1", noise=noise_normal(),
        control = control_f(numer_grad = TRUE)
      ),
    data = list(z = c(z1)),
    control = control_opt(
      iterations = 500,
      n_parallel_chain = 1
    )
  )
  out1
  traceplot(out1, "field1")

  out <- ngme(
    z ~ f(c(1:n ,1:n), model="ar1", replicate=rep(1:2,each=n), noise=noise_normal(),
        control = control_f(numer_grad = TRUE)),
    data = list(z = c(z1, z2)),
    control = control_opt(
      iterations = 20,
      n_parallel_chain = 1
    )
  )
  out
  # traceplot(out, "field1")

  out2 <- ngme(
    z ~ f(c(1:n ,1:n), model="rw1", replicate=rep(1:2,each=n), noise=noise_normal()),
    data = list(z = c(z1, z2)),
    control = control_opt(
      estimation = TRUE,
      iterations = 20
    )
  )
  prds2 <- predict(out2, loc=list(field1 = c(2,4,5)))
  expect_equal(length(prds2), 3)
})

