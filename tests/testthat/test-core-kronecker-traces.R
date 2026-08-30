# Tensor-product operators get their tr(K^-1 dK) from the Kronecker identities
#     tr(K^-1 (dK_1 (x) K_2)) = n_2 tr(K_1^-1 dK_1)
#     tr(K^-1 (K_1 (x) dK_2)) = n_1 tr(K_2^-1 dK_2)
# instead of factorizing the full product. These tests pin the algebra the
# implementation relies on, and the observable consequences of using it.

make_st <- function(seed = 11, n = 300, nt = 5, cutoff = 0.30) {
  set.seed(seed)
  loc <- cbind(stats::runif(n), stats::runif(n))
  list(
    d = data.frame(y = stats::rnorm(n), c1 = loc[, 1], c2 = loc[, 2],
                   time = sample(seq_len(nt), n, TRUE)),
    mesh = fmesher::fm_mesh_2d(loc = loc, cutoff = cutoff,
                               max.edge = c(cutoff * 2, cutoff * 4)),
    nt = nt
  )
}

test_that("tp operator K is exactly first (x) second", {
  st <- make_st()
  op <- f(list(st$d$time, cbind(st$d$c1, st$d$c2)),
          model = tp(first = ar1(mesh = seq_len(st$nt)),
                     second = matern(mesh = st$mesh)),
          noise = noise_normal())$operator
  expect_equal(nrow(op$K), nrow(op$first$K) * nrow(op$second$K))
  expect_equal(as.matrix(op$K),
               as.matrix(kronecker(op$first$K, op$second$K)),
               tolerance = 1e-12)
})

test_that("Kronecker trace identities reproduce the full-product traces", {
  skip_on_cran()
  st <- make_st()
  # `ou` rather than `ar1` for the time factor: ar1 has tr(K_1^-1 dK_1) = 0
  # identically, so it would not exercise the n_2 multiplier at all.
  for (first_model in list(ar1(mesh = seq_len(st$nt)), ou(mesh = seq_len(st$nt)))) {
    op <- f(list(st$d$time, cbind(st$d$c1, st$d$c2)),
            model = tp(first = first_model, second = matern(mesh = st$mesh)),
            noise = noise_normal())$operator
    n1 <- nrow(op$first$K); n2 <- nrow(op$second$K)
    p1 <- op$first$n_theta_K
    eps <- 1e-4 # UpdateOptions::eps_dK, forward difference as in update_all()
    trinv <- function(A, dA) sum(diag(solve(as.matrix(A), as.matrix(dA))))
    set.seed(1)
    for (rep in 1:3) {
      th <- op$theta_K + stats::rnorm(length(op$theta_K), 0, 0.5)
      K1 <- op$first$update_K(th[seq_len(p1)])
      K2 <- op$second$update_K(th[-seq_len(p1)])
      K <- kronecker(K1, K2)
      for (j in seq_along(th)) {
        thp <- th; thp[j] <- thp[j] + eps
        K1p <- op$first$update_K(thp[seq_len(p1)])
        K2p <- op$second$update_K(thp[-seq_len(p1)])
        full <- trinv(K, (kronecker(K1p, K2p) - K) / eps)
        ident <- if (j <= p1) n2 * trinv(K1, (K1p - K1) / eps)
                 else         n1 * trinv(K2, (K2p - K2) / eps)
        expect_equal(ident, full, tolerance = 1e-8)
      }
    }
  }
})

test_that("a tp fit no longer depends on the non-symmetric solver choice", {
  skip_on_cran()
  # The identities remove the factorization of the (non-symmetric) product, so
  # nonsym_solver has nothing left to act on for a tensor product.
  st <- make_st(seed = 12)
  d <- st$d; msh <- st$mesh
  ctl <- function(ns) control_opt(
    seed = 4321, burnin = 10, iterations = 20, n_batch = 10,
    n_parallel_chain = 2, max_num_threads = 2, nonsym_solver = ns,
    verbose = FALSE, print_check_info = FALSE, R_hat_conv_check = FALSE,
    trend_std_conv_check = FALSE)
  a <- ngme(y ~ f(list(time, cbind(c1, c2)),
      model = tp(first = ar1(mesh = 1:5), second = matern(mesh = msh)),
      noise = noise_normal()),
    data = d, family = "normal", control_opt = ctl("lu"))
  b <- ngme(y ~ f(list(time, cbind(c1, c2)),
      model = tp(first = ar1(mesh = 1:5), second = matern(mesh = msh)),
      noise = noise_normal()),
    data = d, family = "normal", control_opt = ctl("normal_equations"))
  drop_setting <- function(dg) dg[!grepl("nonsym_solver", names(dg))]
  expect_equal(ngme_digest_max_diff(drop_setting(ngme_numeric_digest(a)),
                                    drop_setting(ngme_numeric_digest(b))),
               0, tolerance = 1e-8)
})

test_that("the triangular time factor of a tp needs no probes at all", {
  skip_on_cran()
  # The 1d factor (ar1) is triangular, so its trace is analytic and the probe
  # budget never reaches it. The spatial factor is not, so a tp as a whole does
  # still depend on n_trace_iter -- only the time half is probe-free. Check the
  # factor directly rather than through the fit.
  st <- make_st()
  op <- f(list(st$d$time, cbind(st$d$c1, st$d$c2)),
          model = tp(first = ar1(mesh = seq_len(st$nt)),
                     second = matern(mesh = st$mesh)),
          noise = noise_normal())$operator
  expect_true(all(Matrix::triu(op$first$K, k = 1) == 0))   # triangular -> analytic
  expect_false(all(Matrix::triu(op$second$K, k = 1) == 0)) # not -> Hutchinson
})

test_that("the ar1 factor of a tp contributes its analytic trace", {
  # ar1 carries the stationary initial condition, so K[1,1] = sqrt(1 - rho^2)
  # and tr(K_1^-1 dK_1) = -rho/2 exactly (see test-core-ar1-stationary.R). The
  # tp trace for that parameter is therefore n_2 * (-rho/2).
  st <- make_st()
  op <- f(list(st$d$time, cbind(st$d$c1, st$d$c2)),
          model = tp(first = ar1(mesh = seq_len(st$nt)),
                     second = matern(mesh = st$mesh)),
          noise = noise_normal())$operator
  eps <- 1e-6
  for (th in c(-0.8, 0, 0.9)) {
    K1 <- as.matrix(op$first$update_K(th))
    dK1 <- (as.matrix(op$first$update_K(th + eps)) - K1) / eps
    expect_equal(sum(diag(solve(K1, dK1))), -ngme2:::ar1_th2a(th) / 2,
                 tolerance = 1e-5)
  }
})

# ---------------------------------------------------------------------------
# Triangular operators: tr(K^-1 dK) = sum_i (dK)_ii / K_ii, exactly.

test_that("triangular operators have analytic traces", {
  eps <- 1e-4 # UpdateOptions::eps_dK
  for (spec in list(list(nm = "ar1", op = ar1(mesh = 1:12, rho = 0.4)),
                    list(nm = "ou",  op = ou(mesh = 1:12, theta = 1)))) {
    op <- spec$op
    expect_true(all(Matrix::triu(op$K, k = 1) == 0),
                info = paste(spec$nm, "is lower triangular"))
    set.seed(5)
    for (r in 1:3) {
      th <- op$theta_K + stats::rnorm(op$n_theta_K, 0, 0.6)
      K <- as.matrix(op$update_K(th))
      dK <- (as.matrix(op$update_K(th + eps)) - K) / eps
      expect_equal(sum(diag(dK) / diag(K)),
                   sum(diag(solve(K, dK))), tolerance = 1e-9,
                   info = spec$nm)
    }
  }
})

test_that("the ar1 trace is exactly -rho/2 for every rho and every n", {
  eps <- 1e-6
  for (n in c(20, 100)) {
    op <- ar1(mesh = seq_len(n), rho = 0.5)
    for (th in c(-1.2, 0, 0.7, 2.0)) {
      K <- as.matrix(op$update_K(th))
      dK <- (as.matrix(op$update_K(th + eps)) - K) / eps
      rho <- ngme2:::ar1_th2a(th)
      expect_equal(sum(diag(solve(K, dK))), -rho / 2, tolerance = 1e-5)
      expect_equal(sum(diag(dK) / diag(K)), -rho / 2, tolerance = 1e-5)
    }
  }
})

test_that("rw1 and rw2 have no K parameters, so need no trace at all", {
  expect_equal(rw1(mesh = 1:12)$n_theta_K, 0)
  expect_equal(rw2(mesh = 1:12)$n_theta_K, 0)
})

# ---------------------------------------------------------------------------
# Hutchinson probe seeding.

test_that("operator traces redraw their probes every iteration", {
  skip_on_cran()
  sp <- withr::with_seed(4, {
    n <- 400
    loc <- matrix(stats::runif(n * 2), n, 2)
    list(mesh = fmesher::fm_mesh_2d(loc = loc, cutoff = 0.09,
                                    max.edge = c(0.22, 0.7)),
         data = data.frame(y = stats::rnorm(n), c1 = loc[, 1], c2 = loc[, 2]))
  })
  msh <- sp$mesh; dat <- sp$data
  ctl <- function(sd) control_opt(
    seed = sd, burnin = 20, iterations = 60, n_batch = 10, n_parallel_chain = 1,
    max_num_threads = 1, n_trace_iter = 10, start_sd = 0, verbose = FALSE,
    print_check_info = FALSE, R_hat_conv_check = FALSE,
    trend_std_conv_check = FALSE)
  k1 <- as.numeric(ngme_result(ngme(
    y ~ 0 + f(~ c1 + c2, model = matern(mesh = msh), noise = noise_normal()),
    data = dat, family = "normal", control_opt = ctl(11)))$field1$kappa)
  k2 <- as.numeric(ngme_result(ngme(
    y ~ 0 + f(~ c1 + c2, model = matern(mesh = msh), noise = noise_normal()),
    data = dat, family = "normal", control_opt = ctl(22)))$field1$kappa)
  expect_true(is.finite(k1) && is.finite(k2))
  expect_gt(abs(k1 - k2), 0)   # trace noise is now sampled, not frozen
})

test_that("a fit is still reproducible from control_opt(seed = )", {
  skip_on_cran()
  # Redrawing probes must not cost reproducibility: the seeds come from the
  # latent's own RNG, which control_opt(seed = ) determines.
  d <- withr::with_seed(8, {
    n <- 200
    w <- as.numeric(stats::arima.sim(n = n, list(ar = 0.5)))
    data.frame(y = w + stats::rnorm(n, sd = 0.3), idx = seq_len(n),
               c1 = stats::runif(n), c2 = stats::runif(n))
  })
  msh <- fmesher::fm_mesh_2d(loc = cbind(d$c1, d$c2), cutoff = 0.12,
                             max.edge = c(0.3, 0.8))
  ctl <- control_opt(seed = 77, burnin = 10, iterations = 30, n_batch = 10,
    n_parallel_chain = 1, max_num_threads = 1, verbose = FALSE,
    print_check_info = FALSE, R_hat_conv_check = FALSE,
    trend_std_conv_check = FALSE)
  a <- ngme(y ~ 0 + f(~ c1 + c2, model = matern(mesh = msh), noise = noise_normal()),
            data = d, family = "normal", control_opt = ctl)
  b <- ngme(y ~ 0 + f(~ c1 + c2, model = matern(mesh = msh), noise = noise_normal()),
            data = d, family = "normal", control_opt = ctl)
  expect_equal(ngme_digest_max_diff(ngme_numeric_digest(a),
                                    ngme_numeric_digest(b)), 0)
})
