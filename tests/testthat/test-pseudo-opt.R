# works with pseudo.opt R class

# structure
# W1, W2
# V1, V2

# test_that("test structure", {
  # load_all()
  # library(INLA)
  # n1 <<- 5
  # n2 <<- 10
  # z1 = arima.sim(n1, model = list(ar = 0.5), sd = 0.5)
  # z2 = arima.sim(n2, model = list(ar = 0.5), sd = 0.5)

  # idx <<- c(1:n1, 1:n2); idx
  # rep <<- c(rep(1, n1), rep(2, n2)); rep

  # # 1. create mesh using largest idx
  # mesh = inla.mesh.1d(loc = 1:10)
  # inla.spde.make.A(loc = idx, mesh=mesh, repl=rep)

  # out <- ngme(
  #   z ~ f(idx, model="ar1", replicate=rep, noise=noise_normal()),
  #   data = list(z = c(z1, z2), rep=rep),
  #   control = control_opt(
  #     iterations = 20
  #   )
  # )

  # out$latents[[1]]$A

  # m1 <- model_ar1(idx, replicate=rep)
  # str(m1)
  # # W_size is 15 or 20??

  # # out
  # # traceplot(out, "field1")
  # prds <- predict(out, loc=list(field1 = c(2,4,5)))
  # expect_true(length(prds$mean) == 3)
# })
