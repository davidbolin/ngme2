test_that("posterior_plot returns ggplot for SGLD sample data", {
  withr::local_seed(11)

  samples <- data.frame(
    .chain = rep(1:2, each = 40),
    .draw = rep(1:40, times = 2),
    .iter = rep(seq_len(40), times = 2),
    alpha = rnorm(80, mean = 0.5, sd = 0.2),
    beta = rnorm(80, mean = -1, sd = 0.4)
  )
  attr(samples, "name") <- "general"

  p <- posterior_plot(
    samples,
    parameters = c("alpha", "beta"),
    by_chain = TRUE
  )

  expect_s3_class(p, "ggplot")
  expect_true(inherits(plot.ngme_sgld_ci(ngme_sgld_ci(samples)), "ggplot"))
})

test_that("posterior_plot validates parameter selection and layout inputs", {
  samples <- data.frame(
    .chain = 1,
    .draw = 1:5,
    .iter = 1:5,
    alpha = rnorm(5)
  )

  expect_error(
    posterior_plot(samples, parameters = "missing"),
    "Unknown parameter"
  )

  expect_error(
    posterior_plot(samples, ncol = 0),
    "`ncol` must be a positive integer"
  )
})
