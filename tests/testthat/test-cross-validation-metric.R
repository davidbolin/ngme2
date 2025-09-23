test_that("prepare_metric_data applies custom metric across groups", {
  y_data <- c(1, 2, 10, 20)
  group_data <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))

  Y_N_1_thin <- matrix(
    c(1.1, 0.9,
      2.1, 1.9,
      10.1, 9.9,
      20.2, 19.8),
    ncol = 2,
    byrow = TRUE
  )

  Y_N_2_thin <- matrix(
    c(1.3, 1.1,
      2.3, 2.1,
      10.4, 10.0,
      20.4, 20.0),
    ncol = 2,
    byrow = TRUE
  )

  metric_fn <- function(data) {
    val <- 2 * data$y["A"] + data$y["B"]
    names(val) <- "combo"
    val
  }

  res <- ngme2:::prepare_metric_data(metric_fn, y_data, Y_N_1_thin, Y_N_2_thin, group_data)

  expect_equal(res$y, c(12, 24))
  expect_equal(levels(res$group), "combo")
  expect_equal(nrow(res$samples1), 2)
  expect_equal(ncol(res$samples1), 2)

  expected_samples1 <- matrix(
    c(2 * 1.1 + 10.1, 2 * 0.9 + 9.9,
      2 * 2.1 + 20.2, 2 * 1.9 + 19.8),
    nrow = 2,
    byrow = TRUE
  )

  expected_samples2 <- matrix(
    c(2 * 1.3 + 10.4, 2 * 1.1 + 10.0,
      2 * 2.3 + 20.4, 2 * 2.1 + 20.0),
    nrow = 2,
    byrow = TRUE
  )

  expect_equal(res$samples1, expected_samples1)
  expect_equal(res$samples2, expected_samples2)

  expected_pred <- 0.5 * (rowMeans(expected_samples1) + rowMeans(expected_samples2))
  expect_equal(res$pred, expected_pred)
})


test_that("prepare_metric_data identity fallback matches input", {
  y_data <- c(5, 6)
  group_data <- factor(c("single", "single"))
  Y_N_1_thin <- matrix(c(5.1, 4.9), ncol = 1)
  Y_N_2_thin <- matrix(c(5.2, 4.8), ncol = 1)

  res <- ngme2:::prepare_metric_data(NULL, y_data, Y_N_1_thin, Y_N_2_thin, group_data)

  expect_equal(res$y, y_data)
  expect_equal(res$group, droplevels(group_data))
  expect_equal(res$samples1, as.matrix(Y_N_1_thin))
  expect_equal(res$samples2, as.matrix(Y_N_2_thin))
})


test_that("prepare_metric_data errors when metric references unknown groups", {
  y_data <- c(1, 2, 10, 20)
  group_data <- factor(c("A", "A", "B", "B"), levels = c("A", "B"))
  Y_N_1_thin <- matrix(1, nrow = 4, ncol = 1)
  Y_N_2_thin <- matrix(1, nrow = 4, ncol = 1)

  bad_metric <- function(data) {
    val <- 2 * data$y["C"] + data$y["B"]
    names(val) <- "combo"
    val
  }

  expect_error(
    ngme2:::prepare_metric_data(bad_metric, y_data, Y_N_1_thin, Y_N_2_thin, group_data),
    "Available groups: A, B",
    class = "invalid_metric_error"
  )
})


test_that("prepare_metric_data errors when group counts differ", {
  y_data <- c(1, 2, 10, 20)
  group_data <- factor(c("A", "A", "A", "B"), levels = c("A", "B"))
  Y_N_1_thin <- matrix(1, nrow = 4, ncol = 1)
  Y_N_2_thin <- matrix(1, nrow = 4, ncol = 1)

  metric_fn <- function(data) sum(data$y)

  expect_error(
    ngme2:::prepare_metric_data(metric_fn, y_data, Y_N_1_thin, Y_N_2_thin, group_data),
    "same number of observations"
  )
})
