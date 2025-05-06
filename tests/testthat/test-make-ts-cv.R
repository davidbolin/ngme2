# Test that the basic expanding window works

test_that("make_time_series_cv_index: basic expanding window works", {
  # Simple sequence 1:5
  cv <- make_time_series_cv_index(1:5)
  
  # Should have 4 folds for expanding window with 5 points
  expect_equal(length(cv$train), 4)
  expect_equal(length(cv$test), 4)
  
  # Check contents of each fold
  expect_equal(cv$train[[1]], 1)
  expect_equal(cv$test[[1]], 2)
  
  expect_equal(cv$train[[2]], c(1, 2))
  expect_equal(cv$test[[2]], 3)
  
  expect_equal(cv$train[[3]], c(1, 2, 3))
  expect_equal(cv$test[[3]], 4)
  
  expect_equal(cv$train[[4]], c(1, 2, 3, 4))
  expect_equal(cv$test[[4]], 5)
})

test_that("make_time_series_cv_index: fixed window size works", {
  # Using train_length = 2
  cv <- make_time_series_cv_index(1:5, train_length = 2)
  cv$train
  
  # Should have 3 folds for window of size 2
  expect_equal(length(cv$train), 3)
  expect_equal(length(cv$test), 3)
  
  # Check contents of each fold
  expect_equal(cv$train[[1]], c(1, 2))
  expect_equal(cv$test[[1]], 3)
  
  expect_equal(cv$train[[2]], c(2, 3))
  expect_equal(cv$test[[2]], 4)
})

test_that("make_time_series_cv_index: test_length parameter works", {
  # Using test_length = 2
  cv <- make_time_series_cv_index(1:5, test_length = 2)
  
  # Should have 3 folds for expanding window with 5 points and test_length=2
  expect_equal(length(cv$train), 3)
  expect_equal(length(cv$test), 3)
  
  # Check contents of each fold
  expect_equal(cv$train[[1]], 1)
  expect_equal(cv$test[[1]], c(2, 3))
  
  expect_equal(cv$train[[2]], c(1, 2))
  expect_equal(cv$test[[2]], c(3, 4))
  
  expect_equal(cv$train[[3]], c(1, 2, 3))
  expect_equal(cv$test[[3]], c(4, 5))
})

test_that("make_time_series_cv_index: replicates parameter works", {
  # Create time series with replicates
  time_idx <- c(1, 1, 2, 2, 3, 3)
  replicates <- c("A", "A", "B", "B", "C", "C")
  
  cv <- make_time_series_cv_index(time_idx, replicate = replicates)
  
  # Should have 2 folds with 3 replicate groups
  expect_equal(length(cv$train), 2)
  expect_equal(length(cv$test), 2)
  
  # First fold: train on A, test on B
  # Since we're selecting indices by group, the specific indices will depend on how they're ordered
  expect_setequal(cv$train[[1]], c(1, 2))  # indices of group A
  expect_setequal(cv$test[[1]], c(3, 4))   # indices of group B
  
  # Second fold: train on A and B, test on C
  expect_setequal(cv$train[[2]], c(1, 2, 3, 4))  # indices of groups A and B
  expect_setequal(cv$test[[2]], c(5, 6))         # indices of group C
})

test_that("make_time_series_cv_index: gap parameter works", {
  # Using gap = 1 (skip one observation between train and test)
  cv <- make_time_series_cv_index(1:6, gap = 1)
  
  # Should have 4 folds for expanding window with 6 points and gap=1
  expect_equal(length(cv$train), 4)
  expect_equal(length(cv$test), 4)
  
  # Check contents of each fold
  expect_equal(cv$train[[1]], 1)
  expect_equal(cv$test[[1]], 3)  # Skip 2
  
  expect_equal(cv$train[[2]], c(1, 2))
  expect_equal(cv$test[[2]], 4)  # Skip 3
  
  expect_equal(cv$train[[3]], c(1, 2, 3))
  expect_equal(cv$test[[3]], 5)  # Skip 4
  
  expect_equal(cv$train[[4]], c(1, 2, 3, 4))
  expect_equal(cv$test[[4]], 6)  # Skip 5
})

test_that("make_time_series_cv_index: gap and train_length work together", {
  # Using train_length = 2, gap = 1
  cv <- make_time_series_cv_index(1:6, train_length = 2, gap = 1)
  
  # Check contents of each fold
  expect_equal(cv$train[[1]], c(1, 2))
  expect_equal(cv$test[[1]], 4)  # Skip 3
  
  expect_equal(cv$train[[2]], c(2, 3))
  expect_equal(cv$test[[2]], 5)  # Skip 4
})

test_that("make_time_series_cv_index: replicates and gap work together", {
  # Create time series with replicates
  time_idx <- c(1, 1, 2, 2, 3, 3, 4, 4)
  replicates <- c("A", "A", "B", "B", "C", "C", "D", "D")
  
  cv <- make_time_series_cv_index(time_idx, replicate = replicates, gap = 1)
  
  # Should have 2 folds with gap=1
  expect_equal(length(cv$train), 2)
  expect_equal(length(cv$test), 2)
  
  # First fold: train on A, skip B, test on C
  expect_setequal(cv$train[[1]], c(1, 1))  # indices of group A
  expect_setequal(cv$test[[1]], c(3, 3))   # indices of group B (skipping A)
  
  # Second fold: train on A and B, skip C, test on D
  expect_setequal(cv$train[[2]], c(1, 1, 2, 2))  # indices of groups A and B
  expect_setequal(cv$test[[2]], c(4, 4))         # indices of group D (skipping C)
})

test_that("make_time_series_cv_index: error handling works", {
  # Invalid time_idx
  expect_error(make_time_series_cv_index(NULL))
  expect_error(make_time_series_cv_index(1))  # Need at least 2 points
  
  # Invalid test_length
  expect_error(make_time_series_cv_index(1:5, test_length = 0))
  expect_error(make_time_series_cv_index(1:5, test_length = -1))
  
  # Invalid train_length
  expect_error(make_time_series_cv_index(1:5, train_length = 0))
  expect_error(make_time_series_cv_index(1:5, train_length = 5))  # Must be less than length
  
  # Invalid gap
  expect_error(make_time_series_cv_index(1:5, gap = -1))
  
  # Not enough data for train_length + gap + test_length
  expect_error(make_time_series_cv_index(1:5, train_length = 2, gap = 2, test_length = 2))
  
  # Replicate length mismatch
  expect_error(make_time_series_cv_index(1:5, replicate = 1:4))
})

test_that("make_time_series_cv_index: unsorted time indices are handled correctly", {
  # Unsorted indices
  cv <- make_time_series_cv_index(c(3, 1, 5, 2, 4))
  
  # Should have 4 folds
  expect_equal(length(cv$train), 4)
  
  # Check that indices were sorted
  expect_equal(cv$train[[1]], 1)
  expect_equal(cv$test[[1]], 2)
  
  expect_equal(cv$train[[4]], c(1, 2, 3, 4))
  expect_equal(cv$test[[4]], 5)
})

test_that("make_time_series_cv_index: full combinations of parameters work", {
  # Complex case with all parameters
  time_idx <- c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5)
  
  cv <- make_time_series_cv_index(
    time_idx = time_idx,
    train_length = 2,  # Use 2 replicate groups for training
    test_length = 2,   # Test on 2 replicate groups
    gap = 1            # Skip 1 replicate group between train and test
  )
  
  # Should have 1 fold
  expect_equal(length(cv$train), 1)
  expect_equal(length(cv$test), 1)
  
  # Train on A and B, skip C, test on D and E
  expect_setequal(cv$train[[1]], c(1, 1, 2, 2))  # A and B
  expect_setequal(cv$test[[1]], c(4, 4, 5, 5))  # D and E (skipping C)
})
