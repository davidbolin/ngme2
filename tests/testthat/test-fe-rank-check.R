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
