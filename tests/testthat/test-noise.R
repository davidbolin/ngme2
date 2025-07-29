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

test_that("fix_theta_sigma works in noise_normal()", {
  B_sigma <- matrix(c(1,0,0,1), nrow=2)
  
  # Test vector fixing
  n1 <- noise_normal(
    theta_sigma = c(0, 1),
    B_sigma = B_sigma,
    fix_theta_sigma = c(TRUE, FALSE)
  )
  expect_equal(n1$n_theta_sigma, 1)
  expect_equal(n1$fix_theta_sigma, c(TRUE, FALSE))
  
  # Test single value expansion
  n2 <- noise_normal(
    theta_sigma = c(0, 1),
    B_sigma = B_sigma,
    fix_theta_sigma = FALSE
  )
  expect_equal(n2$n_theta_sigma, 2)
  expect_equal(n2$fix_theta_sigma, c(FALSE, FALSE))
})
