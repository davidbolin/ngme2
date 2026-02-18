test_that("batch_means_estimator returns consistent moments", {
  traj <- cbind(theta1 = as.numeric(1:80), theta2 = as.numeric(81:160))

  est <- batch_means_estimator(
    trajectory = traj,
    alpha = 0.6,
    M = 4,
    drop_burnin = TRUE
  )

  expect_s3_class(est, "ngme_batch_means")
  expect_equal(length(est$mean), 2)
  expect_equal(dim(est$covariance), c(2, 2))
  expect_true(all(est$batch_sizes > 0))
  expect_equal(sum(est$batch_sizes), est$n_eff)

  # Weighted mean from retained batches should equal pooled mean
  weighted_mean <- colSums(est$batch_means * est$batch_sizes) / sum(est$batch_sizes)
  expect_equal(as.numeric(est$mean), as.numeric(weighted_mean), tolerance = 1e-10)

  # Covariance must be symmetric positive semidefinite (up to numerical jitter)
  eig <- eigen(est$covariance, symmetric = TRUE, only.values = TRUE)$values
  expect_true(min(eig) > -1e-10)
})


test_that("batch_means_ci pools multiple chains and returns CI", {
  chain1 <- cbind(theta1 = as.numeric(1:120), theta2 = as.numeric(2:121))
  chain2 <- cbind(theta1 = as.numeric(3:122), theta2 = as.numeric(4:123))

  ci <- batch_means_ci(
    chain_trajectories = list(chain1, chain2),
    level = 0.95,
    alpha = 0.6,
    M = 5
  )

  expect_s3_class(ci, "ngme_batch_ci")
  expect_equal(length(ci$estimates), 2)
  expect_equal(length(ci$std_error), 2)
  expect_equal(dim(ci$conf_int), c(2, 2))
  expect_equal(colnames(ci$conf_int), c("lower", "upper"))
  expect_true(all(ci$std_error >= 0))
  expect_true(all(ci$conf_int[, "upper"] >= ci$conf_int[, "lower"]))

  c1 <- batch_means_estimator(chain1, alpha = 0.6, M = 5)
  c2 <- batch_means_estimator(chain2, alpha = 0.6, M = 5)
  expected_mean <- colMeans(rbind(c1$mean, c2$mean))

  expect_equal(as.numeric(ci$estimates), as.numeric(expected_mean), tolerance = 1e-10)
})


test_that("batch_means_ci truncates unequal chain lengths", {
  chain1 <- cbind(theta1 = as.numeric(1:100), theta2 = as.numeric(101:200))
  chain2 <- cbind(theta1 = as.numeric(1:80), theta2 = as.numeric(81:160))

  expect_warning(
    ci <- batch_means_ci(
      chain_trajectories = list(chain1, chain2),
      alpha = 0.6,
      M = 4
    ),
    "different lengths"
  )

  expect_equal(ci$n_iterations, 80)
})


test_that("batch_means_estimator burnin_iter matches manual trajectory trimming", {
  traj <- cbind(theta1 = as.numeric(1:120), theta2 = as.numeric(121:240))
  burnin_iter <- 30

  est_burn <- batch_means_estimator(
    trajectory = traj,
    alpha = 0.6,
    M = 5,
    burnin_iter = burnin_iter
  )

  est_manual <- batch_means_estimator(
    trajectory = traj[(burnin_iter + 1):nrow(traj), , drop = FALSE],
    alpha = 0.6,
    M = 5,
    burnin_iter = 0
  )

  expect_equal(est_burn$mean, est_manual$mean, tolerance = 1e-12)
  expect_equal(est_burn$covariance, est_manual$covariance, tolerance = 1e-12)
  expect_equal(est_burn$n_iterations, nrow(traj) - burnin_iter)
  expect_equal(est_burn$n_iterations_total, nrow(traj))
})


test_that("batch_means_ci burnin_iter propagates across chains", {
  chain1 <- cbind(theta1 = as.numeric(1:140), theta2 = as.numeric(2:141))
  chain2 <- cbind(theta1 = as.numeric(3:142), theta2 = as.numeric(4:143))
  burnin_iter <- 20

  ci_burn <- batch_means_ci(
    chain_trajectories = list(chain1, chain2),
    alpha = 0.6,
    M = 6,
    burnin_iter = burnin_iter
  )

  ci_manual <- batch_means_ci(
    chain_trajectories = list(
      chain1[(burnin_iter + 1):nrow(chain1), , drop = FALSE],
      chain2[(burnin_iter + 1):nrow(chain2), , drop = FALSE]
    ),
    alpha = 0.6,
    M = 6,
    burnin_iter = 0
  )

  expect_equal(ci_burn$estimates, ci_manual$estimates, tolerance = 1e-12)
  expect_equal(ci_burn$covariance, ci_manual$covariance, tolerance = 1e-12)
  expect_equal(ci_burn$n_iterations, nrow(chain1) - burnin_iter)
  expect_equal(ci_burn$n_iterations_total, nrow(chain1))
})


test_that("post-delta transform maps raw CI outputs in transformed scale", {
  chain1 <- cbind(
    theta = as.numeric(seq(-0.2, 0.2, length.out = 120)),
    sigma = as.numeric(seq(log(1.5), log(2.5), length.out = 120))
  )
  chain2 <- cbind(
    theta = as.numeric(seq(-0.1, 0.3, length.out = 120)),
    sigma = as.numeric(seq(log(1.6), log(2.6), length.out = 120))
  )

  ci_raw <- batch_means_ci(
    chain_trajectories = list(chain1, chain2),
    alpha = 0.6,
    M = 5,
    drop_burnin = TRUE
  )

  ci_trans <- .delta_transform_batch_ci(
    ci = ci_raw,
    transforms = list(identity, exp)
  )

  d <- exp(unname(ci_raw$estimates["sigma"]))
  expect_equal(unname(ci_trans$estimates["sigma"]), d, tolerance = 1e-6)

  cov_expected <- d^2 * unname(ci_raw$covariance["sigma", "sigma"])
  expect_equal(unname(ci_trans$covariance["sigma", "sigma"]), cov_expected, tolerance = 1e-6)

  se_expected <- sqrt(cov_expected / (ci_raw$n_eff * ci_raw$n_chains))
  expect_equal(unname(ci_trans$std_error["sigma"]), se_expected, tolerance = 1e-6)
})


test_that("summary.ngme_batch_ci returns parameter table and covariance", {
  chain1 <- cbind(theta1 = as.numeric(1:120), theta2 = as.numeric(2:121))
  chain2 <- cbind(theta1 = as.numeric(3:122), theta2 = as.numeric(4:123))

  ci <- batch_means_ci(
    chain_trajectories = list(chain1, chain2),
    level = 0.95,
    alpha = 0.6,
    M = 5
  )

  sm <- summary(ci)
  expect_s3_class(sm, "summary.ngme_batch_ci")
  expect_true(is.data.frame(sm$table))
  expect_true(all(c("parameter", "mean", "std_error", "lower", "upper") %in% names(sm$table)))
  expect_equal(nrow(sm$table), length(ci$estimates))
  expect_true(is.matrix(sm$covariance))
  expect_equal(dim(sm$covariance), dim(ci$covariance))
  expect_equal(sm$alpha, ci$alpha)

  expect_no_error(print(sm))
})
