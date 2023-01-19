# tests for nig distribution

test_that("test basic", {
  # load_all()
  # plot nig
  plot(
    noise_nig(mu=1, sigma=2, nu=0.2),
    noise_nig(mu=1, sigma=2, nu=2)
  )

})
