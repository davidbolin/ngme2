# The optimizer records one trajectory row per *free* parameter. Slicing those
# rows with `a:b` counts backwards when a component owns no rows at all, which
# handed the measurement noise a non-existent row plus a row belonging to the
# last latent parameter, and blew up traceplot() with "subscript out of bounds".

fit_with_fixed_sigma_e <- function(fixed) {
  set.seed(2)
  n <- 300
  y <- as.numeric(arima.sim(list(ar = 0.5), n)) + rnorm(n, sd = 1)
  ngme(
    y ~ 0 + f(x, model = ar1(), noise = noise_normal()),
    family = if (fixed) noise_normal(theta_sigma = log(1), fix_theta_sigma = TRUE)
             else noise_normal(),
    data = data.frame(x = 1:n, y = y),
    control_opt = control_opt(seed = 1, burnin = 20, iterations = 60,
                              n_parallel_chain = 2, warn_no_convergence = FALSE)
  )
}

test_that("a fully fixed measurement noise records no trajectory rows", {
  rep1 <- fit_with_fixed_sigma_e(TRUE)$replicates[[1]]
  traj <- attr(rep1, "block_traj")

  expect_equal(nrow(traj[[1]]), 0)
  expect_false(anyNA(traj[[1]]))

  # the latent rows must not have been handed to the block
  lat <- attr(rep1$models[[1]], "lat_traj")
  expect_equal(nrow(lat[[1]]), rep1$models[[1]]$n_params)
  expect_equal(nrow(lat[[1]]) + nrow(traj[[1]]), rep1$n_params)
})

test_that("a free measurement noise is unaffected", {
  rep1 <- fit_with_fixed_sigma_e(FALSE)$replicates[[1]]
  traj <- attr(rep1, "block_traj")
  lat <- attr(rep1$models[[1]], "lat_traj")

  expect_equal(nrow(traj[[1]]), 1)
  expect_false(anyNA(traj[[1]]))
  expect_equal(nrow(lat[[1]]) + nrow(traj[[1]]), rep1$n_params)
})

test_that("traceplot works with a fixed measurement noise", {
  fit <- fit_with_fixed_sigma_e(TRUE)

  # the combined plot simply omits the fixed parameter
  expect_s3_class(traceplot(fit), "ggplot")
  expect_s3_class(traceplot(fit, "field1"), "ggplot")

  # asking for the fully fixed component alone is a clear error, not a crash
  expect_error(traceplot(fit, "noise"), "all of its parameters are fixed")
})
