test_that("control_opt rejects model starts passed via partial matching", {
  expect_error(
    control_opt(start = list(dummy = TRUE)),
    "`start` is not a `control_opt\\(\\)` argument"
  )
})

test_that("control_opt validates start_sd eagerly", {
  expect_error(
    control_opt(start_sd = list(dummy = TRUE)),
    "start_sd must be a numeric scalar"
  )
})
