# control_opt(nonsym_solver = ) chooses how the operator matrix K of a
# non-symmetric model is factorized when estimating tr(K^-1 dK):
#   "lu"               a sparse LU of K            (default)
#   "normal_equations" a Cholesky of t(K) K        (the behaviour before 0.9.9)
# The two are different unbiased estimators of the same trace, so they give
# different numbers for a non-symmetric operator -- and must give identical
# ones for a symmetric operator, where the option does not apply at all.

library(ngme2)
library(testthat)

# the fitted object records the setting itself, so drop it before comparing
drop_setting <- function(dg) dg[!grepl("nonsym_solver", names(dg))]

nonsym_control <- function(mode, iterations = 15) {
  control_opt(seed = 4321, burnin = 5, iterations = iterations, n_batch = 1,
              n_parallel_chain = 1, max_num_threads = 1,
              solver_backend = "cholmod", nonsym_solver = mode,
              verbose = FALSE, print_check_info = FALSE,
              R_hat_conv_check = FALSE, pflug_conv_check = FALSE,
              trend_std_conv_check = FALSE)
}

test_that("nonsym_solver is validated and defaults to lu", {
  expect_equal(control_opt()$nonsym_solver, 0L)
  expect_equal(control_opt(nonsym_solver = "lu")$nonsym_solver, 0L)
  expect_equal(control_opt(nonsym_solver = "normal_equations")$nonsym_solver, 1L)
  expect_error(control_opt(nonsym_solver = "bogus"))
})

test_that("nonsym_solver changes a non-symmetric fit and not a symmetric one", {
  skip_on_cran()
  skip_if_not_installed("fmesher")
  # The operator must be non-symmetric AND have no exact-trace shortcut for the
  # option to bite: ar1/ou are triangular and tensor products use the Kronecker
  # identities, so neither factorizes K at all any more. spacetime does.
  st <- withr::with_seed(1, {
    n <- 250
    loc <- matrix(runif(n * 2), n, 2)
    list(mesh_s = fmesher::fm_mesh_2d(loc = loc, cutoff = 0.25,
                                      max.edge = c(0.5, 1)),
         mesh_t = fmesher::fm_mesh_1d(1:5),
         data = data.frame(y = rnorm(n), c1 = loc[, 1], c2 = loc[, 2],
                           time = sample(1:5, n, TRUE)))
  })
  mesh_st <- list(st$mesh_t, st$mesh_s)
  st_data <- st$data
  # ngme() resolves its formula in the calling frame, so these are inlined
  # rather than wrapped in a closure.
  a <- drop_setting(ngme_numeric_digest(ngme(
    y ~ f(list(time, cbind(c1, c2)), model = spacetime(mesh = mesh_st, alpha = 2),
          noise = noise_normal()),
    data = st_data, family = "normal", control_opt = nonsym_control("lu"))))
  b <- drop_setting(ngme_numeric_digest(ngme(
    y ~ f(list(time, cbind(c1, c2)), model = spacetime(mesh = mesh_st, alpha = 2),
          noise = noise_normal()),
    data = st_data, family = "normal",
    control_opt = nonsym_control("normal_equations"))))
  expect_identical(names(a), names(b))
  expect_gt(ngme_digest_max_diff(a, b), 0)

  # matern (alpha = 2) IS symmetric, so the option must be inert
  sp <- withr::with_seed(2, {
    n <- 120
    loc <- matrix(runif(n * 2), n, 2)
    z <- rnorm(n)
    list(mesh = fmesher::fm_mesh_2d(loc = loc, cutoff = 0.12,
                                    max.edge = c(0.3, 0.8)),
         data = data.frame(y = 1 + 0.5 * z + rnorm(n, sd = 0.3), z = z,
                           c1 = loc[, 1], c2 = loc[, 2]))
  })
  sp_mesh <- sp$mesh
  sp_data <- sp$data
  s_lu <- drop_setting(ngme_numeric_digest(ngme(
    y ~ z + f(~ c1 + c2, model = matern(mesh = sp_mesh, alpha = 2),
              noise = noise_normal()),
    data = sp_data, family = "normal", control_opt = nonsym_control("lu"))))
  s_ne <- drop_setting(ngme_numeric_digest(ngme(
    y ~ z + f(~ c1 + c2, model = matern(mesh = sp_mesh, alpha = 2),
              noise = noise_normal()),
    data = sp_data, family = "normal",
    control_opt = nonsym_control("normal_equations"))))
  expect_identical(ngme_digest_max_diff(s_lu, s_ne), 0)
})

test_that("triangular operators ignore nonsym_solver: they never factorize K", {
  skip_on_cran()
  # ar1 and ou are lower triangular, so tr(K^-1 dK) = sum_i (dK)_ii / K_ii is
  # taken exactly and no solver is involved.
  d <- withr::with_seed(1, {
    n <- 150
    w <- as.numeric(arima.sim(n = n, list(ar = 0.5)))
    data.frame(y = w + rnorm(n, sd = 0.4), idx = seq_len(n))
  })
  for (mod in c("ar1", "ou")) {
    a <- drop_setting(ngme_numeric_digest(ngme(
      y ~ f(idx, model = if (mod == "ar1") ar1() else ou(mesh = 1:150),
            noise = noise_normal()),
      data = d, family = "normal", control_opt = nonsym_control("lu"))))
    b <- drop_setting(ngme_numeric_digest(ngme(
      y ~ f(idx, model = if (mod == "ar1") ar1() else ou(mesh = 1:150),
            noise = noise_normal()),
      data = d, family = "normal",
      control_opt = nonsym_control("normal_equations"))))
    expect_equal(ngme_digest_max_diff(a, b), 0, tolerance = 1e-10,
                 info = mod)
  }
})
