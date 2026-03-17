## Unit tests for the VAR(1) bivariate operator (var1()) -- Cayley reparameterization
##
## The operator takes four UNCONSTRAINED parameters (p1, p2, p3, p4) and maps
## them to a stationary 2x2 VAR coefficient matrix A via the Cayley transform:
##
##   J = [[0, p1], [-p1, 0]]              (skew-symmetric)
##   L = [[p2, 0], [p3, p4]]              (lower-triangular)
##   R = L L^T + eps*I                    (positive definite, eps = 1e-5)
##   S = J - R                            (stable: S + S^T = -2R < 0)
##   A = (I + S)(I - S)^{-1}             (Cayley: spectral radius < 1 guaranteed)
##
## Precision matrix:  K = M0 + a11*M11 + a22*M22 + a12*M12 + a21*M21
## Sign convention:   K[t, t-1] = -a  (because C[t, t-1] = -1)
##
## Tests cover:
##   1.  Deferred construction (mesh = NULL)
##   2.  Operator meta-fields: model name, param_name, n_theta_K
##   3.  K has dimension 2T x 2T for various T
##   4.  h vector is rep(1, 2*T)
##   5.  Diagonal entries of K are always 1
##   6.  K sub-diagonals match -A entries (generic params)
##   7.  Off-block entries are zero when cross-lags are zero (p1=0, p3=0)
##   8.  Stationarity guarantee: spectral radius < 1 for random unconstrained params
##   9.  Default params give A close to -eps/(2+eps)*I (near-zero A)
##  10.  update_K() reproduces the expected K after a parameter change
##  11.  cayley_to_A() is stored in operator and is callable
##  12.  f() with group builds a 2T x 2T identity A matrix (integer map)
##  13.  f() model can be simulated without error

## ---- helper ----------------------------------------------------------------
make_op <- function(T, p1 = 0, p2 = 1, p3 = 0, p4 = 1) {
  var1(mesh = 1:T, p1 = p1, p2 = p2, p3 = p3, p4 = p4)
}
## ---------------------------------------------------------------------------

test_that("var1: deferred construction returns ngme_operator_def", {
  op_def <- var1()
  expect_s3_class(op_def, "ngme_operator_def")
  expect_equal(op_def$model, "var1")

  # args are stored correctly
  op_def2 <- var1(mesh = NULL, p1 = 0.7, p2 = 1.2)
  expect_s3_class(op_def2, "ngme_operator_def")
  expect_equal(op_def2$args$p1, 0.7)
  expect_equal(op_def2$args$p2, 1.2)
})

test_that("var1: operator meta-fields are correct", {
  T  <- 8
  op <- var1(mesh = 1:T, p1 = 0.5, p2 = 1.0, p3 = 0.2, p4 = 0.8)

  expect_s3_class(op, "ngme_operator")
  expect_equal(op$model, "var1")
  expect_equal(op$n_theta_K, 4L)
  expect_equal(op$param_name, c("p1", "p2", "p3", "p4"))
  # theta_K stores the unconstrained params
  expect_equal(unname(op$theta_K), c(0.5, 1.0, 0.2, 0.8))
})

test_that("var1: h is rep(1, 2*T)", {
  for (T in c(5, 10, 20)) {
    op <- var1(mesh = 1:T)
    expect_equal(op$h, rep(1, 2 * T))
  }
})

test_that("var1: K has dimension 2T x 2T", {
  for (T in c(5, 10, 20)) {
    op <- var1(mesh = 1:T)
    expect_equal(nrow(op$K), 2L * T)
    expect_equal(ncol(op$K), 2L * T)
  }
})

test_that("var1: diagonal entries of K are always 1", {
  # The diagonal is 1 regardless of A, because K = I + A*C and C has zeros on diagonal
  param_sets <- list(
    c(0,   1, 0, 1),       # default  (A near 0)
    c(1.5, 0.5, 0.3, 0.8), # moderate
    c(-2,  2,  -1, 3),     # large
    c(0,   0.1, 0, 0.1)    # small L  (A near -1)
  )
  for (params in param_sets) {
    T  <- 7
    op <- make_op(T, params[1], params[2], params[3], params[4])
    K  <- as.matrix(op$K)
    expect_equal(diag(K), rep(1, 2 * T), tolerance = 1e-12,
                 label = paste("p =", paste(round(params, 3), collapse = ",")))
  }
})

## Sign convention note:
## The VAR(1) model is  y_t = A * y_{t-1} + eps_t.
## Rearranged: K * y = eps  where K[t, t-1] = -a  (NOT +a).
## This is because the sub-diagonal matrix C has C[t, t-1] = -1, so
##   K = I + a * C  =>  K[t, t-1] = a * (-1) = -a.
## All sub-diagonal expectations below use -a accordingly.

test_that("var1: K sub-diagonals match -A entries from cayley_to_A", {
  T  <- 6
  # Use params that give non-trivial A with all four entries non-zero
  op <- make_op(T, p1 = 1.2, p2 = 0.8, p3 = 0.3, p4 = 1.5)
  A  <- op$cayley_to_A(op$theta_K)
  a11 <- A[1, 1]; a12 <- A[1, 2]
  a21 <- A[2, 1]; a22 <- A[2, 2]
  K   <- as.matrix(op$K)

  for (t in 2:T) {
    expect_equal(K[t,         t - 1],     -a11, tolerance = 1e-12,
                 label = sprintf("K11[%d,%d]", t, t - 1))
    expect_equal(K[T + t, T + t - 1],     -a22, tolerance = 1e-12,
                 label = sprintf("K22[%d,%d]", t, t - 1))
    expect_equal(K[t,     T + t - 1],     -a12, tolerance = 1e-12,
                 label = sprintf("K12[%d,%d]", t, t - 1))
    expect_equal(K[T + t,     t - 1],     -a21, tolerance = 1e-12,
                 label = sprintf("K21[%d,%d]", t, t - 1))
  }

  # Also verify upper-triangular entries in each T x T block are zero
  K11 <- K[1:T, 1:T]
  K22 <- K[(T + 1):(2 * T), (T + 1):(2 * T)]
  expect_equal(sum(abs(K11[upper.tri(K11)])), 0, tolerance = 1e-12)
  expect_equal(sum(abs(K22[upper.tri(K22)])), 0, tolerance = 1e-12)
})

test_that("var1: off-block entries are zero when p1 = 0 and p3 = 0 (diagonal A)", {
  ## With J = 0 (p1 = 0) and L diagonal (p3 = 0):
  ##   R = L L^T + eps*I is diagonal
  ##   S = J - R = -R is diagonal
  ##   A = (I+S)(I-S)^{-1} is diagonal  => a12 = a21 = 0
  T  <- 6
  op <- make_op(T, p1 = 0, p2 = 1.5, p3 = 0, p4 = 2.0)
  A  <- op$cayley_to_A(op$theta_K)
  K  <- as.matrix(op$K)

  expect_equal(A[1, 2], 0, tolerance = 1e-12)
  expect_equal(A[2, 1], 0, tolerance = 1e-12)

  # K12 and K21 blocks should be all zero
  expect_equal(sum(abs(K[1:T, (T + 1):(2 * T)])), 0, tolerance = 1e-12)
  expect_equal(sum(abs(K[(T + 1):(2 * T), 1:T])), 0, tolerance = 1e-12)
})

test_that("var1: stationarity guaranteed for arbitrary unconstrained params", {
  set.seed(2025)
  for (i in 1:30) {
    p <- rnorm(4, sd = 3)    # deliberately large / varied values
    # Keep p2 and p4 positive so that R = LL^T + eps*I stays positive definite
    p[2] <- abs(p[2]) + 0.01
    p[4] <- abs(p[4]) + 0.01

    op <- make_op(5, p[1], p[2], p[3], p[4])
    A  <- op$cayley_to_A(op$theta_K)
    sr <- max(abs(eigen(A, only.values = TRUE)$values))

    expect_lt(sr, 1,
              label = sprintf("spectral radius (%.3f) for p = (%.2f,%.2f,%.2f,%.2f)",
                              sr, p[1], p[2], p[3], p[4]))
  }
})

test_that("var1: default params give A close to -eps/(2+eps)*I", {
  ## With (p1=0, p2=1, p3=0, p4=1):
  ##   J = 0, L = I, R = (1+eps)*I
  ##   S = -(1+eps)*I
  ##   I+S = -eps*I,  I-S = (2+eps)*I
  ##   A = solve((2+eps)*I, -eps*I) = -eps/(2+eps) * I
  eps      <- 1e-5
  A_diag   <- -eps / (2 + eps)   # expected diagonal value (≈ -5e-6)
  op <- var1(mesh = 1:5)
  A  <- op$cayley_to_A(op$theta_K)

  expect_equal(A[1, 1], A_diag, tolerance = 1e-10)
  expect_equal(A[2, 2], A_diag, tolerance = 1e-10)
  expect_equal(A[1, 2], 0,      tolerance = 1e-12)
  expect_equal(A[2, 1], 0,      tolerance = 1e-12)
})

test_that("var1: update_K rebuilds K consistently when theta_K changes", {
  T  <- 6
  op <- make_op(T, p1 = 0.3, p2 = 1.0, p3 = 0.2, p4 = 0.9)

  # New unconstrained params (arbitrary)
  new_theta <- c(p1 = -0.8, p2 = 1.5, p3 = -0.4, p4 = 1.2)
  K_new <- as.matrix(op$update_K(new_theta))

  # Diagonal entries still 1
  expect_equal(diag(K_new), rep(1, 2 * T), tolerance = 1e-12)

  # Sub-diagonals must match -A entries for the new params
  A_new <- op$cayley_to_A(new_theta)
  for (t in 2:T) {
    expect_equal(K_new[t,         t - 1], -A_new[1, 1], tolerance = 1e-12)
    expect_equal(K_new[T + t, T + t - 1], -A_new[2, 2], tolerance = 1e-12)
    expect_equal(K_new[t,     T + t - 1], -A_new[1, 2], tolerance = 1e-12)
    expect_equal(K_new[T + t,     t - 1], -A_new[2, 1], tolerance = 1e-12)
  }

  # Must match a fresh operator with the same params
  op2 <- make_op(T, p1 = new_theta["p1"], p2 = new_theta["p2"],
                    p3 = new_theta["p3"], p4 = new_theta["p4"])
  expect_equal(K_new, as.matrix(op2$K), tolerance = 1e-12)
})

test_that("var1: cayley_to_A is stored in operator and is callable", {
  op <- var1(mesh = 1:5, p1 = 1.0, p2 = 0.5, p3 = 0.0, p4 = 0.5)

  expect_true(is.function(op$cayley_to_A))

  A <- op$cayley_to_A(op$theta_K)
  expect_equal(dim(A), c(2L, 2L))

  # stationarity from the stored function
  sr <- max(abs(eigen(A, only.values = TRUE)$values))
  expect_lt(sr, 1)
})

test_that("var1: f() builds a 2T x 2T identity A matrix for integer map + group", {
  T     <- 8
  n_obs <- 2 * T
  group <- factor(c(rep("y1", T), rep("y2", T)), levels = c("y1", "y2"))
  map   <- c(1:T, 1:T)

  f_model <- f(
    map,
    model = var1(mesh = 1:T),
    group = group,
    noise = noise_normal()
  )

  A <- as.matrix(f_model$A)
  expect_equal(dim(A), c(n_obs, n_obs))
  expect_equal(A, diag(n_obs), tolerance = 1e-12)
})

test_that("var1: f() model can be simulated without error", {
  set.seed(42)
  T     <- 20
  group <- factor(c(rep("s1", T), rep("s2", T)), levels = c("s1", "s2"))

  f_model <- f(
    c(1:T, 1:T),
    model = var1(mesh = 1:T, p1 = 0.5, p2 = 1.0, p3 = 0.3, p4 = 0.8),
    group = group,
    noise = noise_nig(mu = 0, sigma = 1, nu = 1)
  )

  sim <- simulate(f_model, seed = 1)
  W   <- sim[[1]]

  expect_length(W, 2 * T)
  expect_true(all(is.finite(W)))
  expect_false(anyNA(W))
})
