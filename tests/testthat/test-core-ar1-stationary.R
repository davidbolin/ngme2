# The AR(1) operator carries its stationary initial condition: the first
# observation has marginal sd 1/sqrt(1 - rho^2), so K[1,1] = sqrt(1 - rho^2)
# must follow the *current* rho rather than the value ar1() was constructed
# with. It lives in its own matrix (E11) because the generic operator builds K
# as a linear combination of fixed matrices, so the only place a non-linear
# function of theta_K can sit is a coefficient -- the "sech" transformation.

test_that("K[1,1] tracks the current rho, not the constructor's", {
  op <- ar1(mesh = 1:8, rho = 0.5)
  for (th in c(-2, -1, 0, 0.5493, 1.5, 3)) {
    rho <- ngme2:::ar1_th2a(th)
    expect_equal(as.numeric(op$update_K(th)[1, 1]), sqrt(1 - rho^2),
                 tolerance = 1e-12, info = paste("theta =", th))
  }
  # and it must not depend on the rho ar1() was built with
  a <- ar1(mesh = 1:8, rho = -0.9)$update_K(0.7309)
  b <- ar1(mesh = 1:8, rho = 0.4)$update_K(0.7309)
  expect_equal(as.matrix(a), as.matrix(b), tolerance = 1e-12)
})

test_that("the sech transformation is sqrt(1 - tanh^2) on both sides", {
  # A silent fall-through to identity in the C++ transform table would make
  # K[1,1] = theta, so pin the R side against the closed form here and let
  # test_that above pin the assembled K.
  f <- ngme2:::name2fun("sech")
  for (th in c(-4, -1.3863, 0, 0.7309, 1.7346, 4))
    expect_equal(f(th), sqrt(1 - ngme2:::ar1_th2a(th)^2), tolerance = 1e-12)
})

test_that("K stays lower triangular and its trace is analytic (-rho/2)", {
  # rho = tanh(theta/2), so d/dtheta log|K| = d/dtheta log sqrt(1 - rho^2)
  #   = -rho rho' / (1 - rho^2) = -rho/2  since rho' = (1 - rho^2)/2.
  op <- ar1(mesh = 1:12, rho = 0.3)
  eps <- 1e-6
  for (th in c(-1.2, -0.4, 0.3, 0.9, 1.6)) {
    K <- op$update_K(th)
    expect_true(all(Matrix::triu(K, k = 1) == 0))
    Km <- as.matrix(K)
    dK <- (as.matrix(op$update_K(th + eps)) - Km) / eps
    rho <- ngme2:::ar1_th2a(th)
    expect_equal(sum(diag(solve(Km, dK))), -rho / 2, tolerance = 1e-5)
    # the O(n) form the operator actually uses
    expect_equal(sum(diag(dK) / diag(Km)), -rho / 2, tolerance = 1e-5)
  }
})

test_that("an ar1 fit is invariant to the Hutchinson probe count", {
  skip_on_cran()
  # The *operator* trace is taken analytically, so n_trace_iter cannot move the
  # fit. It used to: the trace was estimated stochastically and rho moved by
  # ~0.05. This holds with Rao-Blackwellisation on too: for a banded QQ its
  # conditional-variance corrections are taken from the selected inverse rather
  # than estimated, so no probes enter the fit at all.
  d <- withr::with_seed(8, {
    n <- 400
    w <- as.numeric(stats::arima.sim(n = n, list(ar = 0.6)))
    data.frame(y = w + stats::rnorm(n, sd = 0.3), idx = seq_len(n))
  })
  ctl <- function(nti) control_opt(
    seed = 5, burnin = 20, iterations = 60, n_batch = 10, n_parallel_chain = 1,
    max_num_threads = 1, n_trace_iter = nti, start_sd = 0, verbose = FALSE,
    print_check_info = FALSE, R_hat_conv_check = FALSE,
    trend_std_conv_check = FALSE)
  a <- ngme(y ~ 0 + f(idx, model = ar1(), noise = noise_normal()),
            data = d, family = "normal", control_opt = ctl(5))
  b <- ngme(y ~ 0 + f(idx, model = ar1(), noise = noise_normal()),
            data = d, family = "normal", control_opt = ctl(200))
  expect_equal(as.numeric(ngme_result(a)$field1$rho),
               as.numeric(ngme_result(b)$field1$rho), tolerance = 1e-12)
})
