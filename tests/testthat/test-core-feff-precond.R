## Preconditioner blocks for the fixed effects and the measurement noise.
## With the latent field integrated out they must match the exact marginal
## information of a Gaussian AR1 model, not the complete-data curvature (which
## overstates it many times over here).

test_that("par_names follow the C++ order: latent, measurement noise, fixed effects", {
  set.seed(3)
  n <- 100
  dat <- data.frame(y = rnorm(n), x = rnorm(n), t = 1:n)
  fit <- ngme(
    y ~ 1 + x + f(t, model = ar1(), noise = noise_normal(), name = "field"),
    data = dat, family = "normal",
    control_opt = control_opt(estimation = FALSE)
  )
  expect_equal(
    fit$replicates[[1]]$par_names,
    c("rho", "sigma_1", "meas_sigma_1", "feff_1", "feff_2")
  )
})

test_that("control_opt validates the measurement sigma preconditioner options", {
  ctl <- control_opt()
  expect_identical(ctl$precond_meas_sigma, 0L)
  expect_identical(ctl$fisher_refresh_every, 10L)
  expect_identical(control_opt(precond_meas_sigma = "fisher")$precond_meas_sigma, 1L)
  expect_identical(control_opt(precond_meas_sigma = "complete")$precond_meas_sigma, 2L)
  expect_error(control_opt(precond_meas_sigma = "exact"))
  expect_error(control_opt(fisher_refresh_every = 0))
  expect_error(control_opt(fisher_refresh_every = 2.5))
})

test_that("precond_meas_sigma = 'complete' warns for non-Gaussian measurement noise", {
  set.seed(4)
  n <- 100
  dat <- data.frame(y = rnorm(n), t = 1:n)
  expect_warning(
    ngme(y ~ 1 + f(t, model = ar1(), noise = noise_normal()),
         data = dat, family = noise_nig(),
         control_opt = control_opt(estimation = FALSE, precond_meas_sigma = "complete")),
    "only supported for Gaussian"
  )
})

ar1_fixture <- function(n = 400) {
  set.seed(16)
  x <- rnorm(n)
  W <- simulate(
    f(1:n, model = ar1(rho = 0.7), noise = noise_normal(sigma = 1.5)),
    seed = 16
  )[[1]]
  list(n = n, x = x,
       dat = data.frame(y = 2 + 1.5 * x + W + rnorm(n, sd = 0.5), x = x, t = 1:n))
}

fit_ar1 <- function(fx, ...) {
  suppressWarnings(ngme(
    y ~ 1 + x + f(t, model = ar1(), noise = noise_normal(), name = "field"),
    data = fx$dat, family = "normal",
    control_opt = control_opt(
      iterations = 600, n_parallel_chain = 4, seed = 16,
      optimizer = precond_sgd(), ...
    )
  ))
}

# Diagonal of the preconditioner at the last checkpoint, and the chain-mean
# parameters on the internal scale.
last_info <- function(fit) {
  cd <- attr(fit, "conv_diag")
  last <- cd[cd$iteration == max(cd$iteration), ]
  cp <- colMeans(attr(fit, "chain_params"))
  names(cp) <- attr(fit, "par_names")
  list(info = setNames(last$info_kk, last$param), par = cp)
}

test_that("fixed effects and measurement sigma preconditioners match the exact AR1 information", {
  fx <- ar1_fixture()
  n <- fx$n
  li <- last_info(fit_ar1(fx, precond_meas_sigma = "fisher"))
  info <- li$info
  r_hat <- ngme2:::ar1_th2a(li$par[["rho"]])
  s_ar <- exp(li$par[["sigma_1"]])
  s_eps <- exp(li$par[["meas_sigma_1"]])

  K <- Matrix::bandSparse(n, k = c(0, -1),
    diagonals = list(c(sqrt(1 - r_hat^2), rep(1, n - 1)), rep(-r_hat, n - 1)))
  QQ <- Matrix::crossprod(K) / s_ar^2 + Matrix::Diagonal(n, 1 / s_eps^2)
  marg_info <- function(v) {
    Dv <- v / s_eps^2
    sum(v * Dv) - sum(Dv * as.numeric(Matrix::solve(QQ, Dv)))
  }
  # Internally the intercept column is ones and x is centred and SVD-scaled.
  xc <- fx$x - mean(fx$x)
  exact <- c(marg_info(rep(1, n)), marg_info(xc / sqrt(sum(xc^2))))

  expect_equal(unname(info[c("feff_1", "feff_2")]), exact, tolerance = 0.02)
  # And well below the complete-data value for the intercept.
  expect_lt(info[["feff_1"]], 0.1 * n / s_eps^2)

  # Measurement log sigma: the marginal Fisher information
  # 1/2 tr(Sigma^-1 dSigma Sigma^-1 dSigma) = 2 s_eps^4 ||Sigma_Y^-1||_F^2,
  # estimated with Hutchinson probes, so compared loosely.
  Si <- as.matrix(Matrix::Diagonal(n, 1 / s_eps^2)) -
    as.matrix(Matrix::solve(QQ, Matrix::Diagonal(n, 1 / s_eps^2))) / s_eps^2
  fisher_sigma <- 2 * s_eps^4 * sum(Si^2)
  expect_equal(info[["meas_sigma_1"]], fisher_sigma, tolerance = 0.15)
  expect_lt(info[["meas_sigma_1"]], 0.5 * 2 * n)
})

test_that("precond_meas_sigma = 'auto' is the complete-data block for Gaussian noise", {
  fx <- ar1_fixture()
  fit_auto <- fit_ar1(fx)
  fit_complete <- fit_ar1(fx, precond_meas_sigma = "complete")
  fit_fisher <- fit_ar1(fx, precond_meas_sigma = "fisher")
  # Same seed, same blocks: identical optimisation.
  expect_identical(attr(fit_auto, "chain_params"), attr(fit_complete, "chain_params"))
  expect_identical(last_info(fit_auto)$info, last_info(fit_complete)$info)
  # The Fisher block is a different preconditioner, so the path differs.
  expect_false(identical(attr(fit_auto, "chain_params"), attr(fit_fisher, "chain_params")))
})
