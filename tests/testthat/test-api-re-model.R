test_that("re model works with f", {
  out <- f(~ 1 + x,
    model = re(),
    noise = noise_nig(),
    data = data.frame(x = 1:7)
  )
  expect_true(inherits(out, "ngme_model"))
  expect_equal(out$operator$model, "re")
})
