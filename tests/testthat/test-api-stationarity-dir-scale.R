# The flat-direction test compares dir = |sum d| / sum|d| against
# stationarity_dir_lim. Under a stationary null E[dir] = 1/sqrt(n), so an
# UNSCALED threshold means something different at every history length: the
# same 0.5 called a stationary parameter flat ~64% of the time at n = 6 and
# ~83% at the old cap of n = 12. Scaling by sqrt(n) holds the test's size
# constant and makes the history length a free choice.
#
# The default limit is therefore 0.5 * sqrt(6), chosen so the scaled test
# accepts exactly what the old unscaled 0.5 accepted at the point the test
# first becomes eligible. That equivalence is the whole reason for the
# constant, so it is pinned here: changing either half alone silently moves
# the criterion's false-positive rate with nothing else to catch it.

test_that("scaled directionality default reproduces the historical threshold", {
  ctl <- control_opt()

  expect_true(isTRUE(ctl$stationarity_dir_scaled))
  expect_equal(ctl$stationarity_dir_lim, 0.5 * sqrt(6))

  # With scaling on, the history is held to stationarity_min_checks and the
  # test needs that many points, so n is exactly min_checks when it fires.
  n <- ctl$stationarity_min_checks
  expect_equal(n, 6L)

  # Scaled test at the default limit <=> unscaled test at the historical 0.5.
  historical_lim <- 0.5
  expect_equal(ctl$stationarity_dir_lim / sqrt(n), historical_lim)

  # and the equivalence is an "accepts the same dir" statement, not just
  # arithmetic: values either side of the boundary must classify identically.
  dir_vals <- c(0.30, 0.49, 0.50, 0.51, 0.70)
  scaled_flat   <- dir_vals * sqrt(n) < ctl$stationarity_dir_lim
  unscaled_flat <- dir_vals < historical_lim
  expect_identical(scaled_flat, unscaled_flat)
})

test_that("the scaled limit is passed through to the fit controls", {
  ctl <- control_opt(stationarity_dir_scaled = FALSE, stationarity_dir_lim = 0.5)
  expect_false(isTRUE(ctl$stationarity_dir_scaled))
  expect_equal(ctl$stationarity_dir_lim, 0.5)
})
