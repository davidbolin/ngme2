# # test tensor product R structure
# test_that("R structure", {
# library(INLA)
#   load_all()
#   tp(first=ar1(1:3), second=ar1(1:3))
#   true_model <- f(model="tp", first=ar1(1:3), second=ar1(1:3))

#   out <- ngme(
#     y ~ f(1:6, name="tp", model="tp", first=ar1(1:2), second=ar1(1:3)),
#     data = data.frame(y=1:6),
#     control_opt = control_opt(
#       estimation = T,
#       iterations = 1
#     )
#   )
# out
#   str(out$replicates[[1]]$models[[1]]$operator)
#   str(out$replicates[[1]]$models[[1]]$operator$first)


#   # f(i_space, model = spde, group = i_time, control.group = list(model="ar1"))
# })

# test_that("ar x 2d case", {
#   set.seed(16)
#   library(INLA)

# ##############################  simulation
#   mesh2d <- fmesher::fm_mesh_2d(
#     loc.domain = cbind(c(0, 1, 1, 0, 0) * 10, c(0, 0, 1, 1, 0) * 5),
#     max.edge = c(1.2, 10)
#   )
#   plot(mesh2d)
#   mesh2d$n
#   n1 <- 8
#   n2 <- 300
#   loc <- cbind(runif(n2, 0, 10), runif(n2, 0, 5))
# loc


# load_all()
# library(ngme2)
#   tensor_model <- f(model="tp",
#     first = ar1(rep(1:n1, each=n2), theta_K = ar1_a2th(0.7)),
#     second = matern(map=loc, mesh=mesh2d, theta_K = 1),
#     noise = noise_nig(mu=-3, sigma=2, nu=1),
#     eval = T
#   )

#   expect_true(all(tensor_model$operator$K ==
#      ar1(1:n1, theta_K = ar1_a2th(0.7))$K %x% matern(map=loc, mesh=mesh2d, theta_K = 1)$K))
#   W <- simulate(tensor_model)

#   AW <- as.numeric(tensor_model$A %*% W)
#   n_obs <- length(AW)
#   Y <- AW + rnorm(n_obs, sd=0.5)
#   mean(attr(W, "noise")$V)

# ##############################  estimation
#   load_all()
#   # loc is a list
#   out <- ngme(
#     Y ~ 0 +
#     f(
#       model="tp",
#       name="tp",
#       first = iid(1:n1),
#       # first = ar1(1:n1),
#       second = matern(map=loc, mesh=mesh2d),
#       # control = control_f(numer_grad = T),
#       noise = noise_nig(),
#       # fix_W = TRUE, W=W
#     ),
#     # f(
#     #   loc, model="matern", replicate=rep(1:8, each=length_map(loc)),
#     #   control = control_f(numer_grad = T),
#     #   noise = noise_nig()
#     # ),
#     data = data.frame(Y=Y),
#     family = "normal",
#     control_opt = control_opt(
#       iterations = 1,
#       n_parallel_chain = 1,
#       estimation = T,
#       verbose = T,
#       stepsize = 1,
#       max_absolute_step = 10,
#       max_relative_step = 10
#     ),
#     debug = FALSE
#   )
#   out
#   out$replicates[[1]]$models[[1]]$operator
#   traceplot(out, "tp")
#   traceplot(out)
# })

# ######################################################################

test_that("iid x ar case", {
  set.seed(16)
  n_obs <- 200
  Y1 <- simulate(f(1:n_obs, model="ar1", rho=0.7, noise=noise_normal()))
  Y2 <- simulate(f(1:n_obs, model="ar1", rho=0.7, noise=noise_normal()))
  Y3 <- simulate(f(1:n_obs, model="ar1", rho=0.7, noise=noise_normal()))
  Y <- c(Y1, Y2, Y3) + rnorm(n = 3*n_obs)
  time <- rep(1:3, each=n_obs)
  loc <- rep(1:n_obs, 3)

  out <- ngme(
    Y ~ 0 + f(
      map = list(time, loc),
      model = "tp",
      first = iid(n=3),
      second = ar1(loc),
    ),
    data = data.frame(Y=Y),
    control_opt = control_opt(
      iterations = 10
    )
  )

  out <- ngme(
    Y ~ 0 + f(
      map = list(time, loc),
      model = "tp",
      first = list(model="iid", n=3),
      second = list(model="ar1"),
      # control = control_f(numer_grad = T),
      debug=F
    ),
    data = data.frame(Y=Y),
    control_opt = control_opt(
      estimation = T,
      iterations = 500,
      n_parallel_chain = 1,
      verbose = F
    ),
    debug= F
  )
  out
  traceplot(out, "field1")
})

# ##########################################################################

test_that("ar x ar case (simulation)", {
  time_idx <- 1:3; loc_index <- 4:7
  map <- list(
    time = c(1,3,2,3),
    loc = c(4,4,4,4)
  )

  tensor_model <- f(
    map = map,
    model="tp",
    first  = list(model="ar1", rho = 0.7, mesh=time_idx),
    second = list(model="ar1", rho = 0.1, mesh=loc_index),
    noise = noise_nig(mu=-3, sigma=2, nu=1)
  )
  # So the mesh$n is 3 * 4 = 12
  with(tensor_model, {
      expect_true(nrow(A)==4 && ncol(A)==12)
      expect_true(nrow(operator$K)==12)
  })
})

test_that("ar x ar case", {
  set.seed(16)
  time_len <- 8
  n_obs <- 300
  time <- rep(1:time_len, each=n_obs)
  loc <- rep(1:n_obs, time_len)

  tensor_model <- f(
    map = list(time, loc),
    model="tp",
    first  = list(model="ar1", rho = 0.7),
    second = list(model="ar1", rho = 0.1),
    noise = noise_nig(mu=-3, sigma=2, nu=1)
  )
  tensor_model
  W <- simulate(tensor_model)
  AW <- as.numeric(tensor_model$A %*% W)
  n_obs <- length(AW)
  Y <- AW + rnorm(n_obs, sd=0.5)

##############################  estimation
  out <- ngme(
    Y ~ 0 + f(
      map = list(time, loc),
      model="tp", name="tp",
      first = list(model="ar1"),
      second = list(model="ar1"),
      control = control_f(numer_grad = T),
      noise = noise_nig(
        # mu=-3, sigma=2, nu=1,
        # fix_V = T, V = attr(W, "noise")$V
      ),
      # theta_K =  c(ar1_a2th(0.7), ar1_a2th(-0.2)),
      # fix_W = TRUE, W=W
    ),
    data = data.frame(Y=Y),
    family = "normal",
    control_opt = control_opt(
      iterations = 100,
      n_parallel_chain = 4,
      estimation = T,
      verbose = T,
      stepsize = 1,
      max_absolute_step = 10,
      max_relative_step = 10
    ),
    debug = FALSE
  )
  out
  traceplot(out, "tp")
  traceplot(out)
})

test_that("build A for tp model", {
  time = c(1, 1, 3, 2, 3, 2, 1)
  loc = c(3, 2, 1, 1, 2, 2, 1)
  m <- f(
    map = list(time, loc),
    model = "tp",
    first = ar1(mesh=time),
    second = rw1(mesh=loc)
  )
  m$A
})