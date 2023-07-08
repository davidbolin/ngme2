# tests regarding create noise

test_that("noise_normal() works", {
  # basic things
  tmp <- noise_normal(sigma = 3, n = 10)
  expect_equal(tmp$theta_sigma, log(3))
  # expect_equal(tmp$n_noise, 100)

  # test non-stationary mu
  expect_error(
    noise_normal(theta_sigma = c(1, 2))
  )

  expect_no_error(
    noise_normal(theta_sigma = c(1, 2),
      B_sigma = matrix(c(rep(1, 10), rep(2, 10)), ncol=2))
  )

  # test auto-complete using n works
  noise_normal(
    theta_sigma = c(1.2, 1.2),
    B_sigma = matrix(c(1,2), ncol=2),
    n = 10
  )
})
