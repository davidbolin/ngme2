test_that("ou discretisation matches exact AR(1) form", {
  loc <- c(0, 1, 2.5, 3.5)
  theta <- 0.5

  op <- ou(mesh = loc, theta = theta)

  dt <- diff(loc)
  rho <- exp(-theta * dt)
  diag_main <- c(sqrt(1 - rho[1]^2), rep(1, length(loc) - 1))
  sub_diag <- -rho
  expected <- Matrix::bandSparse(
    n = length(loc),
    k = c(0, -1),
    diagonals = list(diag_main, sub_diag)
  )

  expect_equal(as.matrix(op$K), as.matrix(expected))
})

test_that("test tp", {
  n <- 10
  time <- sample(2001:2004, n, replace = TRUE)
  loc <- cbind(runif(n), runif(n)) * 10

  # Display the space-time observation structure
  data.frame(
    observation = 1:n,
    time = time,
    x = round(loc[, 1], 2),
    y = round(loc[, 2], 2)
  )

  # Create spatial mesh
  mesh_spatial <- fmesher::fm_mesh_2d(
    loc.domain = cbind(
      c(0, 10, 10, 0, 0), # x range: 0 to 10
      c(0, 0, 10, 10, 0) # y range: 0 to 10 (to cover all points)
    ),
    max.edge = c(1, 10),
    cutoff = 0.1,
    offset = c(0.5, 2) # Offset to ensure boundary coverage
  )

  tp_model <- f(
    map = list(time, loc),
    model = tp(
      first = ar1(mesh = time, rho = 0.5),
      second = matern(mesh = mesh_spatial, kappa = 2)
    )
  )
  tp_model
})

test_that("test tp-bv", {
  time_len <- 20
  half_bv_len <- 50

  # specify indices for 2 parts
  time_idx <- rep(1:time_len, each = half_bv_len * 2)
  bv_idx <- rep(c(1:half_bv_len, 1:half_bv_len), time_len)

  # group to indicate two different fields for bivariate model
  group <- rep(rep(c("f1", "f2"), each = half_bv_len), time_len)

  tp(
    first = ar1(mesh = time_idx),
    second = bv(
      mesh = bv_idx,
      sub_models = list(
        f1 = ar1(),
        f2 = ar1()
      )
    )
  )
})

test_that("test bivariate (bv)", {
  temp <- c(32, 33, 35.5, 36)
  year_temp <- c(2001, 2002, 2003, 2004)
  precip <- c(0.1, 0.2, 0.5, 1, 0.2)
  year_pre <- c(2001, 2002, 2003, 2004, 2005)

  # bind 2 fields in one vector, and make labels for them
  y <- c(temp, precip)
  year <- c(year_temp, year_pre)
  labels <- c(rep("temp", 4), rep("precip", 5)) # group is label for 2 fields

  x1 <- 1:9
  data <- data.frame(y, year, x1, labels)

  load_all()
  bv_model <- bv(
    year,
    theta = pi / 8,
    rho = 0.5,
    sub_models = list(A = ar1(), B = rw1())
  )

  bv_model
})
