test_that("tp-bv-matern operator structure and simulation", {
  withr::local_seed(42)

  time_len <- 30
  half_bv_len <- 10
  rho <- 0.6

  # specify indices for 2 parts
  time_idx <- rep(1:time_len, each = half_bv_len * 2)
  bv_idx <- rep(c(1:half_bv_len, 1:half_bv_len), time_len)

  # group to indicate two different fields for bivariate model
  group <- rep(rep(c("f1", "f2"), each = half_bv_len), time_len)

  mesh_t <- fmesher::fm_mesh_1d(1:time_len)
  mesh_s <- fmesher::fm_mesh_1d(1:half_bv_len)

  f_model <- f(
    map = list(time_idx, bv_idx),
    model = tp(
      first = ar1(
        mesh = mesh_t,
        rho = rho
      ),
      second = bv_matern(
        mesh = mesh_s,
        group = group,
        sub_models = list(
          f1 = matern(),
          f2 = matern()
        )
      )
    ),
    group = group,
    noise = noise_nig(
      mu = -2,
      sigma = 1,
      nu = 2
    )
  )

  op <- f_model$operator
  expect_equal(op$model, "tp")
  expect_equal(op$first$model, "ar1")
  expect_equal(op$second$model, "bv_matern")
  expect_equal(op$second$model_names, c("f1", "f2"))

  rho_from_theta <- ar1_th2a(op$first$theta_K)
  expect_equal(as.numeric(rho_from_theta), rho, tolerance = 1e-10)

  K_ar1 <- op$first$K
  expect_equal(as.numeric(K_ar1[1, 1]), sqrt(1 - rho^2), tolerance = 1e-10)
  expect_equal(as.numeric(K_ar1[2, 1]), -rho, tolerance = 1e-10)
  expect_equal(as.numeric(K_ar1[3, 2]), -rho, tolerance = 1e-10)
  expect_equal(as.numeric(K_ar1[2, 2]), 1, tolerance = 1e-10)

  expected_tp_K <- op$first$K %x% op$second$K
  expect_equal(as.matrix(op$K), as.matrix(expected_tp_K))

  expected_bv_size <- 2 * mesh_s$n
  expect_equal(nrow(op$second$K), expected_bv_size)
  expect_equal(ncol(op$second$K), expected_bv_size)
  expect_equal(nrow(op$K), time_len * expected_bv_size)

  sim <- simulate(f_model, seed = 42)
  expect_length(sim, 1)

  W <- sim[[1]]
  expect_length(W, length(time_idx))
  expect_true(is.numeric(W))
  expect_false(anyNA(W))
  expect_true(all(is.finite(W)))

  # Estimate rho from latent W (not A %*% W)
  W_latent <- as.numeric(attr(sim, "W_sim")[[1]])
  expect_length(W_latent, time_len * 2 * mesh_s$n)

  w_mat <- matrix(W_latent, nrow = time_len, byrow = TRUE)
  w_centered <- sweep(w_mat, 2, colMeans(w_mat), FUN = "-")
  x <- w_centered[1:(time_len - 1), , drop = FALSE]
  y <- w_centered[2:time_len, , drop = FALSE]

  rho_by_col <- colSums(x * y) / colSums(x * x)
  rho_by_col <- rho_by_col[is.finite(rho_by_col)]
  rho_hat <- median(rho_by_col, na.rm = TRUE)
  expect_true(is.finite(rho_hat))
  expect_true(rho_hat > 0.35 && rho_hat < 0.9)
})
