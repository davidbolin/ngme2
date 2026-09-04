# The non-separable spacetime operator must start from its STATIONARY
# distribution. Before this was fixed the first time slice had no operator block
# at all -- the rw1 factor's first row is a trapezoid spanning every time slice,
# so slice 1 was pinned only through that aggregate. It left the slice
# effectively unconstrained, degenerated K as cc -> 0, and made K non-triangular
# so the traces needed probes.

make_op <- function(stat, nt = 10, cc = 1, kappa = 2) {
  skip_if_not_installed("fmesher")
  set.seed(9)
  loc <- cbind(stats::runif(200), stats::runif(200))
  mesh_s <- fmesher::fm_mesh_2d(loc, max.edge = c(0.3, 0.7), cutoff = 0.15)
  op <- spacetime(
    mesh = list(fmesher::fm_mesh_1d(1:nt), mesh_s), alpha = 2,
    fix_gamma = TRUE, shared_theta_gamma = FALSE,
    theta_gamma_x = 0, theta_gamma_y = 0,
    stabilization = FALSE, stationary_init = stat)
  list(K = as.matrix(op$update_K(c(log(cc), log(kappa)))), ns = mesh_s$n, nt = nt)
}

test_that("the stationary initial condition leaves K block lower-triangular", {
  skip_on_cran()
  o <- make_op(TRUE)
  occupied <- vapply(seq_len(o$nt), function(j)
    max(abs(o$K[seq_len(o$ns), ((j - 1) * o$ns + 1):(j * o$ns)])) > 1e-10,
    logical(1))
  # only the diagonal block of the first row may be occupied
  expect_true(occupied[1])
  expect_false(any(occupied[-1]))
})

test_that("time slice 1 has approximately the stationary marginal variance", {
  skip_on_cran()
  # The first block is gamma * M with gamma^2 = 1 - tr(N'N)/tr(M'M), which is
  # exactly 1 - rho^2 in the scalar AR(1) case but only an approximation to the
  # matrix correction M'M - N'N. It is used in preference to an exact factor
  # because a factor's ROWS are arbitrary combinations of nodes, and under
  # non-Gaussian noise each row of K carries its own mixing variable -- an exact
  # factor would make the fitted model depend on the factorization ordering.
  # The tolerance reflects that: a few percent, against ratios of 6190 (cc = 1)
  # and 512507 (cc = 0.1) with no first block at all.
  for (cc in c(1, 0.1)) {
    o <- make_op(TRUE, cc = cc)
    S <- solve(o$K); Sig <- S %*% t(S)
    mid <- ((o$nt %/% 2 - 1) * o$ns + 1):((o$nt %/% 2) * o$ns)
    ratio <- mean(diag(Sig[seq_len(o$ns), seq_len(o$ns)])) /
             mean(diag(Sig[mid, mid]))
    expect_equal(ratio, 1, tolerance = 0.05)
  }
})

test_that("the first block has the same sparsity as an interior block", {
  skip_on_cran()
  # This is the property that makes the construction safe for non-Gaussian
  # noise: rows of K_1 are the same kind of object as rows of an interior block,
  # so each carries a mixing variable tied to a mesh node rather than to a
  # combination of them.
  o <- make_op(TRUE)
  b1 <- o$K[seq_len(o$ns), seq_len(o$ns)]
  b2 <- o$K[(o$ns + 1):(2 * o$ns), (o$ns + 1):(2 * o$ns)]
  expect_equal(which(abs(b1) > 1e-12), which(abs(b2) > 1e-12))
})

test_that("K no longer degenerates as cc goes to zero", {
  skip_on_cran()
  rc <- vapply(c(1, 1e-2, 1e-4),
               function(cc) rcond(make_op(TRUE, cc = cc)$K), numeric(1))
  # flat in cc, rather than falling linearly with it as it did before
  expect_gt(min(rc) / max(rc), 0.5)
  rc_old <- vapply(c(1, 1e-4),
                   function(cc) rcond(make_op(FALSE, cc = cc)$K), numeric(1))
  expect_lt(rc_old[2] / rc_old[1], 1e-3)
})

test_that("tr(K^-1 dK) is exact from the diagonal blocks", {
  skip_on_cran()
  cc <- 1.3; e <- 1e-6
  K  <- make_op(TRUE, cc = cc)$K
  o  <- make_op(TRUE, cc = cc)
  dK <- (make_op(TRUE, cc = cc * exp(e))$K -
         make_op(TRUE, cc = cc * exp(-e))$K) / (2 * e)
  truth <- sum(diag(solve(K, dK)))
  blocks <- sum(vapply(seq_len(o$nt), function(t) {
    i <- ((t - 1) * o$ns + 1):(t * o$ns)
    sum(diag(solve(K[i, i], dK[i, i])))
  }, numeric(1)))
  expect_equal(blocks, truth, tolerance = 1e-8)
})

test_that("stationary_init = FALSE restores the previous operator", {
  skip_on_cran()
  o <- make_op(FALSE)
  occupied <- vapply(seq_len(o$nt), function(j)
    max(abs(o$K[seq_len(o$ns), ((j - 1) * o$ns + 1):(j * o$ns)])) > 1e-10,
    logical(1))
  expect_true(all(occupied))   # the trapezoid row spans every time slice
})

test_that("C++ and R assemble the same number of filled time blocks", {
  skip_on_cran()
  skip_if_not_installed("fmesher")
  # The C++ assembly used Ct.rows() as the block count. Ct is the temporal mass
  # matrix and is (nt-1) x (nt-1) for a degree-2 temporal mesh, so the last time
  # block silently lost its spatial operator while R filled it. A degree-2 mesh
  # is what the non-separable space-time models are normally built on, and the
  # disagreement was invisible until the block-trace identity failed.
  set.seed(1)
  loc <- cbind(stats::runif(400), stats::runif(400))
  mesh_s <- fmesher::fm_mesh_2d(loc, max.edge = c(0.2, 0.5), cutoff = 0.1)
  mesh_t <- fmesher::fm_mesh_1d(0:7, boundary = c("neumann", "neumann"),
                                degree = 2)
  op <- spacetime(mesh = list(mesh_t, mesh_s), alpha = 2, fix_gamma = TRUE,
    shared_theta_gamma = FALSE, theta_gamma_x = 0, theta_gamma_y = 0,
    stabilization = FALSE, stationary_init = TRUE)
  expect_lt(nrow(op$Ct), op$nt)   # the configuration that exposed the bug
  ns <- mesh_s$n; nt <- op$nt
  K <- as.matrix(op$update_K(c(0, log(2))))
  b1 <- K[(ns + 1):(2 * ns), (ns + 1):(2 * ns)]
  for (t in 2:(nt - 1)) {
    bt <- K[(t * ns + 1):((t + 1) * ns), (t * ns + 1):((t + 1) * ns)]
    expect_equal(max(abs(bt - b1)), 0, tolerance = 1e-12)
  }
})

test_that("the spacetime advection term is LINEAR in gamma", {
  skip_on_cran()
  skip_if_not_installed("fmesher")
  # It was quadratic -- diag(gamma) B diag(gamma) -- which is not an advection
  # term. Two things follow from that and both are checked here, because either
  # one makes the free-gamma model unusable:
  #   * the derivative vanishes at gamma = 0, so zero advection is an exact
  #     stationary point of the likelihood that a chain can never leave;
  #   * gamma and -gamma give the same operator, so the direction of flow --
  #     the thing the free-gamma model exists to estimate -- is unidentifiable.
  set.seed(9)
  loc <- cbind(stats::runif(300), stats::runif(300))
  mesh_s <- fmesher::fm_mesh_2d(loc, max.edge = c(0.25, 0.6), cutoff = 0.12)
  mesh_t <- fmesher::fm_mesh_1d(0:7, boundary = c("neumann", "neumann"),
                                degree = 1)
  mk <- function(gx, gy) {
    op <- spacetime(mesh = list(mesh_t, mesh_s), alpha = 2, fix_gamma = FALSE,
      shared_theta_gamma = FALSE, theta_gamma_x = gx, theta_gamma_y = gy,
      stabilization = FALSE, stationary_init = TRUE)
    as.matrix(op$update_K(c(0, log(2), gx, gy)))
  }
  # reversing the advection must change the operator
  expect_gt(max(abs(mk(1.3, 0.7) - mk(-1.3, 0.7))), 1e-8)
  # and the gradient must not vanish at gamma = 0
  e <- 1e-5
  d0 <- max(abs((mk(e, 0.7) - mk(-e, 0.7)) / (2 * e)))
  expect_gt(d0, 1e-6)
  # linear in gamma: the derivative is the same at 0 and away from it
  d1 <- max(abs((mk(1 + e, 0.7) - mk(1 - e, 0.7)) / (2 * e)))
  expect_equal(d0, d1, tolerance = 1e-6)
})
