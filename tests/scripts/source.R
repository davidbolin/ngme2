  # library(INLA)
  # load_all()
  # n1 <<- 5
  # n2 <<- 10
  # z1 = arima.sim(n1, model = list(ar = 0.5), sd = 0.5)
  # z2 = arima.sim(n2, model = list(ar = 0.5), sd = 0.5)

  # idx <- c(1:n1, 1:n2); idx
  # rep <- c(rep(1, n1), rep(2, n2)); rep

  # # 1. create mesh using largest idx
  # mesh = inla.mesh.1d(loc = 1:10)
  # inla.spde.make.A(loc = idx, mesh=mesh, repl=rep)

  # out <- ngme(
  #   z ~ f(idx, model="ar1", replicate=rep, noise=noise_normal()),
  #   data = list(z = c(z1, z2)),
  #   control = ngme_control(
  #     iterations = 20
  #   )
  # )

  # m1 <- model_ar1(idx, replicate=rep)
  # str(m1)
  # # W_size is 15 or 20??


  # # out
  # # traceplot(out, "field1")
  # prds <- predict(out, loc=list(field1 = c(2,4,5)))
  # expect_equal(length(prds), 6)

  # out2 <- ngme(
  #   z ~ f(c(1:n ,1:n), model="rw1", replicate=rep(1:2,each=n), noise=noise_normal()),
  #   data = list(z = c(z1, z2)),
  #   control = ngme_control(
  #     estimation = TRUE,
  #     iterations = 20
  #   )
  # )
  # prds2 <- predict(out2, loc=list(field1 = c(2,4,5)))
  # expect_equal(length(prds2), 6)

load_all()
####### test nrep = 1
# n1 <- 20
# idx <<- 1:n1
# z1 <- arima.sim(n1, model = list(ar = 0.5), sd = 0.5)
# ar <- model_ar1(idx)

# ar$V_size
# ar$W_size
# ar$n_rep
# dim(ar$A)
# dim(ar$C)

# load_all()
# ngme(
#   z ~ f(idx, model="ar1", noise=noise_normal(), debug=TRUE),
#   data = list(z = z1),
#   control = ngme_control(
#     iterations = 20,
#     n_parallel_chain = 1
#   ),
#   debug = TRUE
# )

######## test nrep=2
load_all()
n1 <- 5; n2 <- 10
idx <<- c(1:n1, 1:n2)
z1 <- arima.sim(n1, model = list(ar = 0.5), sd = 0.5)
z2 <- arima.sim(n2, model = list(ar = 0.5), sd = 0.5)
ar <- model_ar1(idx, replicate = c(rep(1, n1), rep(2, n2)))
ar$noise$n_noise

# ar$V_size
# ar$n_rep

out <- ngme(
  z ~ f(idx, model="ar1", replicate = c(rep(1, n1), rep(2, n2)), noise=noise_nig(), debug=TRUE),
  data = list(z = c(z1, z2)),
  control = ngme_control(
    estimation = TRUE,
    iterations = 20,
    n_parallel_chain = 1
  ),
  debug = TRUE
)

str(out$latent$field1$noise)
