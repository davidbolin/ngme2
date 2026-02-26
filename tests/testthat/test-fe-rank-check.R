test_that("grouped fe() terms fail early when design matrix is rank deficient", {
  df <- data.frame(
    lon = c(11.2, 11.2, 13.1, 13.1, 15, 15),
    lat = c(46.7, 46.7, 46.7, 46.7, 46.7, 46.7),
    direction = c("u_wind", "v_wind", "u_wind", "v_wind", "u_wind", "v_wind"),
    wind = c(0.1, -0.2, 0.3, 0.1, 0.8, 0.3)
  )

  expect_error(
    ngme(
      formula = wind ~
        0 +
        fe(~ 1 + lon, which = "u_wind") +
        fe(~ 1 + lat, which = "v_wind"),
      data = df,
      group = df$direction,
      family = noise_normal(),
      control_opt = control_opt(estimation = FALSE)
    ),
    "not full rank.*lat_v_wind"
  )
})

test_that("grouped fe() regression matches lm.fit and is invariant to SVD standardization", {
  withr::local_seed(20260225)

  n_per_group <- 200
  n <- 2 * n_per_group

  direction <- rep(c("u_wind", "v_wind"), each = n_per_group)
  lon <- runif(n, min = 8, max = 16)
  lat <- runif(n, min = 42, max = 50)

  true_coef <- c(
    "(Intercept)_u_wind" = 1.5,
    "lon_u_wind" = 0.8,
    "(Intercept)_v_wind" = -2.0,
    "lat_v_wind" = 0.35
  )

  wind <- ifelse(
    direction == "u_wind",
    true_coef["(Intercept)_u_wind"] + true_coef["lon_u_wind"] * lon,
    true_coef["(Intercept)_v_wind"] + true_coef["lat_v_wind"] * lat
  )

  df <- data.frame(
    lon = lon,
    lat = lat,
    direction = direction,
    wind = wind
  )

  fit_std <- ngme(
    formula = wind ~
      0 +
      fe(~ 1 + lon, which = "u_wind") +
      fe(~ 1 + lat, which = "v_wind"),
    data = df,
    group = df$direction,
    family = noise_normal(),
    control_opt = control_opt(
      estimation = FALSE,
      standardize_fixed = TRUE
    )
  )

  fit_nosvd <- ngme(
    formula = wind ~
      0 +
      fe(~ 1 + lon, which = "u_wind") +
      fe(~ 1 + lat, which = "v_wind"),
    data = df,
    group = df$direction,
    family = noise_normal(),
    control_opt = control_opt(
      estimation = FALSE,
      standardize_fixed = FALSE
    )
  )

  ngme_rep_std <- fit_std$replicates[[1]]
  ngme_rep_nosvd <- fit_nosvd$replicates[[1]]

  coef_ngme_std <- ngme_rep_std$feff
  coef_ngme_nosvd <- ngme_rep_nosvd$feff
  coef_lm <- stats::lm.fit(ngme_rep_std$X, ngme_rep_std$Y)$coefficients

  pred_ngme_std <- as.numeric(ngme_rep_std$X %*% coef_ngme_std)
  pred_ngme_nosvd <- as.numeric(ngme_rep_nosvd$X %*% coef_ngme_nosvd)
  pred_lm <- as.numeric(ngme_rep_std$X %*% coef_lm)

  # Structural zeros from fe() should remain zeros outside each group.
  expect_equal(
    unname(ngme_rep_std$X[direction != "u_wind", "lon_u_wind"]),
    rep(0, sum(direction != "u_wind")),
    tolerance = 1e-10
  )
  expect_equal(
    unname(ngme_rep_std$X[direction != "v_wind", "lat_v_wind"]),
    rep(0, sum(direction != "v_wind")),
    tolerance = 1e-10
  )

  expect_equal(max(abs(coef_ngme_std - coef_lm)), 0, tolerance = 1e-10)
  expect_equal(max(abs(pred_ngme_std - pred_lm)), 0, tolerance = 1e-10)
  expect_equal(max(abs(pred_ngme_std - wind)), 0, tolerance = 1e-10)

  # Compare with no-SVD fixed-effect standardization.
  expect_equal(coef_ngme_std, coef_ngme_nosvd, tolerance = 1e-10)
  expect_equal(pred_ngme_std, pred_ngme_nosvd, tolerance = 1e-10)

  # Slopes should be recovered exactly in this noiseless simulation.
  expect_equal(coef_ngme_std["lon_u_wind"], true_coef["lon_u_wind"], tolerance = 1e-10)
  expect_equal(coef_ngme_std["lat_v_wind"], true_coef["lat_v_wind"], tolerance = 1e-10)
})

test_that("demeaned fixed effects are mapped back to raw parameterization", {
  withr::local_seed(123)

  df <- data.frame(
    x1 = rnorm(80),
    x2 = rnorm(80)
  )
  df$y <- 2 - 0.5 * df$x1 + 1.3 * df$x2 + rnorm(80, sd = 0.2)

  fit <- ngme(
    y ~ x1 + x2,
    data = df,
    control_opt = control_opt(
      estimation = FALSE,
      standardize_fixed = TRUE
    )
  )

  coef_ngme <- fit$replicates[[1]]$feff
  coef_lm <- stats::lm.fit(model.matrix(~ x1 + x2, data = df), df$y)$coefficients
  expect_equal(coef_ngme, coef_lm, tolerance = 1e-10)

  beta_init <- c("(Intercept)" = 2, "x1" = -0.5, "x2" = 1.3)
  fit_init <- ngme(
    y ~ x1 + x2,
    data = df,
    control_ngme = control_ngme(beta = beta_init),
    control_opt = control_opt(
      estimation = FALSE,
      standardize_fixed = TRUE
    )
  )
  expect_equal(fit_init$replicates[[1]]$feff, beta_init, tolerance = 1e-10)
})

test_that("no-intercept fixed effects are not centered and keep original semantics", {
  withr::local_seed(456)

  n <- 100
  df <- data.frame(
    x1 = runif(n),
    x2 = rnorm(n)
  )
  beta_true <- c("x1" = 1.8, "x2" = -0.7)
  X_raw <- model.matrix(~ 0 + x1 + x2, data = df)
  df$y <- as.numeric(X_raw %*% beta_true) + rnorm(n, sd = 0.1)

  fit <- ngme(
    y ~ 0 + x1 + x2,
    data = df,
    control_opt = control_opt(
      estimation = FALSE,
      standardize_fixed = TRUE
    )
  )

  rep1 <- fit$replicates[[1]]
  coef_lm <- stats::lm.fit(X_raw, df$y)$coefficients

  expect_equal(as.numeric(rep1$X), as.numeric(X_raw), tolerance = 1e-10)
  expect_equal(dim(rep1$X), dim(X_raw))
  expect_equal(colnames(rep1$X), colnames(X_raw))
  expect_equal(rep1$feff, coef_lm, tolerance = 1e-10)
})

test_that("fixed-effect standardization excludes intercept columns", {
  withr::local_seed(789)

  n <- 120
  df <- data.frame(
    x1 = rnorm(n),
    x2 = runif(n)
  )
  df$y <- 1.5 + 0.7 * df$x1 - 0.3 * df$x2 + rnorm(n, sd = 0.2)

  fit <- ngme(
    y ~ x1 + x2,
    data = df,
    control_opt = control_opt(
      estimation = FALSE,
      standardize_fixed = TRUE
    )
  )

  rep1 <- fit$replicates[[1]]
  X_raw <- model.matrix(~ x1 + x2, data = df)
  expect_true(rep1$standardize)
  expect_equal(rep1$svd$idx, which(!startsWith(colnames(X_raw), "(Intercept)")))
  expect_equal(colnames(rep1$svd$u), colnames(X_raw)[rep1$svd$idx])
  expect_equal(as.numeric(rep1$X), as.numeric(X_raw), tolerance = 1e-10)
})

test_that("restore uses replicate indices when svd$u is global", {
  x <- c(0, 1, 2, 10, 11, 12)
  X_raw <- cbind("(Intercept)" = 1, "x" = x)
  centered <- ngme2:::ngme_center_non_intercept_cols(X_raw)
  X_centered <- centered$X

  svd_full <- svd(X_centered)
  colnames(svd_full$u) <- colnames(X_raw)

  beta_raw <- c("(Intercept)" = 3, "x" = 0.5)
  beta_centered <- ngme2:::ngme_beta_raw_to_center(beta_raw, centered$center)
  beta_svd <- as.numeric(svd_full$d * drop(t(svd_full$v) %*% beta_centered))
  names(beta_svd) <- names(beta_raw)

  idx_rep2 <- 4:6
  ngme_rep <- list(
    X = svd_full$u[idx_rep2, , drop = FALSE],
    feff = beta_svd,
    X_center = centered$center,
    X_center_target = centered$center_target,
    standardize = TRUE,
    svd = svd_full,
    data_idx = idx_rep2
  )

  restored <- ngme2:::ngme_restore_fixed_effect_scale(ngme_rep)

  expect_equal(restored$X, X_raw[idx_rep2, , drop = FALSE], tolerance = 1e-10)
  expect_equal(restored$feff, beta_raw, tolerance = 1e-10)
})

test_that("isotropic normal beta prior is invariant to fixed-effect standardization", {
  withr::local_seed(20260226)

  n <- 80
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  df$y <- 1.2 + 0.8 * df$x1 - 1.5 * df$x2 + rnorm(n, sd = 0.4)

  beta_prior <- priors(beta = prior_normal(mean = 0, sd = 0.3))

  fit_std <- ngme(
    y ~ x1 + x2,
    data = df,
    prior_beta = beta_prior,
    control_opt = control_opt(
      optimizer = adam(),
      burnin = 100,
      iterations = 400,
      n_parallel_chain = 1,
      seed = 101,
      standardize_fixed = TRUE,
      verbose = FALSE,
      print_check_info = FALSE
    )
  )

  fit_raw <- ngme(
    y ~ x1 + x2,
    data = df,
    prior_beta = beta_prior,
    control_opt = control_opt(
      optimizer = adam(),
      burnin = 100,
      iterations = 400,
      n_parallel_chain = 1,
      seed = 101,
      standardize_fixed = FALSE,
      verbose = FALSE,
      print_check_info = FALSE
    )
  )

  coef_std <- fit_std$replicates[[1]]$feff
  coef_raw <- fit_raw$replicates[[1]]$feff

  expect_equal(coef_std, coef_raw, tolerance = 1e-5)
})

test_that("incompatible custom beta prior disables fixed-effect standardization", {
  withr::local_seed(20260226)

  n <- 40
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  df$y <- 0.5 + 1.1 * df$x1 - 0.7 * df$x2 + rnorm(n, sd = 0.2)

  # Different variances across fixed effects => not isotropic normal.
  prior_bad <- priors(
    "(Intercept)" = prior_normal(mean = 0, sd = 1.0),
    x1 = prior_normal(mean = 0, sd = 0.3),
    x2 = prior_normal(mean = 0, sd = 0.8)
  )

  expect_warning(
    fit <- ngme(
      y ~ x1 + x2,
      data = df,
      prior_beta = prior_bad,
      control_opt = control_opt(
        estimation = FALSE,
        standardize_fixed = TRUE
      )
    ),
    "standardize_fixed = TRUE\\) was disabled"
  )

  rep1 <- fit$replicates[[1]]
  X_raw <- model.matrix(~ x1 + x2, data = df)

  expect_false(rep1$standardize)
  expect_equal(as.numeric(rep1$X), as.numeric(X_raw), tolerance = 1e-10)
  expect_equal(dim(rep1$X), dim(X_raw))
  expect_equal(colnames(rep1$X), colnames(X_raw))
})

test_that("start works when previous fit disabled fixed-effect standardization", {
  withr::local_seed(20260226)

  n <- 40
  df <- data.frame(
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
  df$y <- 0.8 + 1.2 * df$x1 - 0.4 * df$x2 + rnorm(n, sd = 0.2)

  prior_bad <- priors(
    "(Intercept)" = prior_normal(mean = 0, sd = 1.0),
    x1 = prior_normal(mean = 0, sd = 0.3),
    x2 = prior_normal(mean = 0, sd = 0.8)
  )

  expect_warning(
    fit_prev <- ngme(
      y ~ x1 + x2,
      data = df,
      prior_beta = prior_bad,
      control_opt = control_opt(
        estimation = FALSE,
        standardize_fixed = TRUE
      )
    ),
    "standardize_fixed = TRUE\\) was disabled"
  )
  expect_false(fit_prev$replicates[[1]]$standardize)

  expect_error(
    fit_next <- ngme(
      y ~ x1 + x2,
      data = df,
      start = fit_prev,
      control_opt = control_opt(
        estimation = FALSE,
        standardize_fixed = TRUE
      )
    ),
    NA
  )

  expect_true(fit_next$replicates[[1]]$standardize)
})
