test_that("Check K for all cases", {
  n <- 7
  # rw1
  rw1_constr <- rw1(1:n, constr = TRUE)
  rw1_start_0 <- rw1(1:n, constr = FALSE)
  rw1_cyc <- rw1(1:n, constr = TRUE, cyclic = TRUE)
  print(rw1_constr$K)
  print(rw1_start_0$K)
  print(rw1_cyc$K)

  # rw2
  rw2_constr <- rw2(1:n)
  rw2_cyc <- rw2(1:n, cyclic = TRUE)
  print(rw2_constr$K)
  print(rw2_cyc$K)
})


test_that("test RW1 constraint (non-cyclic)", {
  set.seed(123)
  n <- 500
  h = rw1(1:n)$h
  x <- 1:n

  # default constraint (constr=TRUE)
  w <- simulate(f(x, model="rw1"), nsim = 1)[[1]]
  expect_true(abs(sum(h*w)) < 1e-4)
  y <- w + rnorm(n, sd=0.5)
  plot(x, y, type="l")
  fit <- ngme(
    y ~ 0+f(x, model="rw1", noise=noise_normal()),
    data = data.frame(x=x, y=y),
    control_opt=control_opt(
      iterations=100,
      verbose=TRUE
    )
  )
  fit
  traceplot(fit, "field1")
  traceplot(fit)
  expect_true(abs(ngme_result(fit, "field1")$sigma - 1) < 0.2)
  expect_true(abs(ngme_result(fit, "data")$sigma - 0.5) < 0.3)
})

test_that("test RW1 (start 0) (non-cyclic)", {
  set.seed(123)
  n <- 500
  h = rw1(1:n)$h
  x <- 1:n

  # default constraint (constr=TRUE)
  w <- simulate(f(x, model="rw1", constr=FALSE), nsim = 1)[[1]]
  expect_true(abs(w[1]) < 1e-4)
  
  y <- w + rnorm(n, sd=0.5)
  plot(x, y, type="l")
  fit <- ngme(
    y ~ 0+f(x, model="rw1", constr=FALSE, noise=noise_normal()),
    data = data.frame(x=x, y=y),
    control_opt=control_opt(
      iterations=100
    )
  )
  fit
  expect_true(abs(ngme_result(fit, "field1")$sigma - 1) < 0.2)
  expect_true(abs(ngme_result(fit, "data")$sigma - 0.5) < 0.3)
})

test_that("test RW1 constraint (cyclic)", {
  set.seed(123)
  n <- 500
  x <- 1:n
  m1 <- f(x, model="rw1", cyclic = TRUE)
  h <- m1$operator$h
  w <- simulate(m1, nsim = 1)[[1]]
  # Notice that len(h) == len(y) + 1
  length(h); length(w);
  expect_true(abs(sum(h[-1]*w)) < 1e-4)

  y <- w + rnorm(n, sd=0.5)
  fit <- ngme(
    y ~ 0+f(x, model="rw1", cyclic = TRUE, noise=noise_normal()),
    data = data.frame(x=x, y=y),
    control_opt=control_opt(
      rao_blackwell=TRUE,
      verbose=TRUE,
      iterations=100
    )
  )
  fit
  expect_true(abs(ngme_result(fit, "field1")$sigma - 1) < 0.2)
  expect_true(abs(ngme_result(fit, "data")$sigma - 0.5) < 0.3)
})


test_that("test RW2 (non-cyclic)", {
  set.seed(123)
  n <- 500
  x <- 1:n
  m1 <- f(x, model="rw2")
  w <- simulate(m1, nsim = 1)[[1]]
  plot(x, w, type="l")

  # default constraint (constr=TRUE)
  expect_true(abs(sum(m1$operator$h*w)) < 1e-2)
  expect_true(abs(sum(x*m1$operator$h*w)) < 1e-2)

  y <- w + rnorm(n, sd=0.5)
  fit <- ngme(
    y ~ 0+f(x, model="rw2", noise=noise_normal()),
    data = data.frame(x=x, y=y),
    control_opt=control_opt(
      rao_blackwell=TRUE,
      optimizer = precond_sgd(),
      iterations=100
    )
  )
  fit
  traceplot(fit, "field1")
  expect_true(abs(ngme_result(fit, "field1")$sigma - 1) < 0.2)
  expect_true(abs(ngme_result(fit, "data")$sigma - 0.5) < 0.3)
})


test_that("test RW2 (cyclic)", {
  set.seed(123)
  n <- 500
  x <- 1:n
  m1 <- f(x, model="rw2", cyclic = TRUE)
  m1
  w <- simulate(m1, nsim = 1)[[1]]
  plot(x, w, type="l")

  # default constraint (constr=TRUE)
  expect_true(abs(sum(m1$operator$h[-1]*w)) < 1e-2)

  y <- w + rnorm(n, sd=0.5)
  fit <- ngme(
    y ~ 0+f(x, model="rw2", cyclic = TRUE, noise=noise_normal()),
    data = data.frame(x=x, y=y),
    control_opt=control_opt(
      iterations=100,
      optimizer = precond_sgd()
    )
  )
  fit
  traceplot(fit, "field1")
  expect_true(abs(ngme_result(fit, "field1")$sigma - 1) < 0.2)
  expect_true(abs(ngme_result(fit, "data")$sigma - 0.5) < 0.3)
})


test_that("test RW1 with fix_theta_sigma", {
  set.seed(123)
  n <- 50
  x <- 1:n
  
  # Create B_sigma for RW1 (first parameter fixed, second estimated)
  B_sigma <- matrix(0, nrow = n, ncol = 2)
  B_sigma[1, 1] <- 1      # First observation uses first parameter
  B_sigma[2, 1] <- 1      # First observation uses first parameter
  B_sigma[2:n, 2] <- 1    # Rest use second parameter
  B_sigma
  
  set.seed(42)
  y <- cumsum(rnorm(n, 0, 0.1)) + rnorm(n, 0, 0.3)
  
  fit <- ngme(
    y ~ 0 + f(x, 
      model = "rw1", 
      cyclic = TRUE,
      noise = noise_nig(
        theta_sigma = c(-2, 0),
        B_sigma = B_sigma,
        fix_theta_sigma = c(TRUE, FALSE)  # Fix first, estimate second
      )
    ),
    data = data.frame(x = x, y = y),
    control_opt = control_opt(iterations = 20, verbose = FALSE)
  )
  fit
  
  # Check that first parameter stayed fixed
  result_sigma <- ngme_result(fit, "field1")$theta_sigma
  expect_equal(result_sigma[1], -2, tolerance = 1e-8)
  
  # Check that second parameter was estimated (should not be exactly 0)
  expect_false(abs(result_sigma[2]) < 1e-6)
})

test_that("test RW2 with fix_theta_sigma", {
  n <- 30
  x <- 1:n
  
  # Create B_sigma for RW2 with 3 parameters
  B_sigma <- matrix(0, nrow = n, ncol = 3)
  B_sigma[1:10, 1] <- 1
  B_sigma[11:20, 2] <- 1
  B_sigma[21:n, 3] <- 1
  
  set.seed(123)
  y <- cumsum(cumsum(rnorm(n, 0, 0.05))) + rnorm(n, 0, 0.2)
  
  fit <- ngme(
    y ~ 0 + f(x, 
      model = "rw2",
      noise = noise_normal(
        theta_sigma = c(-1, 0, 0.5),
        B_sigma = B_sigma,
        fix_theta_sigma = c(TRUE, FALSE, TRUE)  # Fix 1st and 3rd
      )
    ),
    data = data.frame(x = x, y = y),
    control_opt = control_opt(iterations = 10, verbose = FALSE)
  )
  
  result_sigma <- ngme_result(fit, "field1")$theta_sigma
  
  # Check fixed parameters stayed fixed
  expect_equal(result_sigma[1], -1, tolerance = 1e-8)
  expect_equal(result_sigma[3], 0.5, tolerance = 1e-8)
  
  # Check estimated parameter changed
  expect_false(abs(result_sigma[2]) < 1e-6)
})
