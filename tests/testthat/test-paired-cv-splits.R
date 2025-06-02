test_that("create_paired_cv_splits basic functionality", {
  # Create test data
  data_long <- data.frame(
    lon = rep(c(1, 2, 3, 4), 2),
    lat = rep(c(1, 1, 2, 2), 2),
    direction = rep(c("u_wind", "v_wind"), each = 4),
    wind = rnorm(8)
  )
  
  # Test basic functionality
  splits <- create_paired_cv_splits(
    data = data_long,
    loc_col = c("lon", "lat"),
    group = "direction",
    k = 2,
    seed = 123
  )
  
  # Check return structure
  expect_type(splits, "list")
  expect_named(splits, c("test_idx", "train_idx"))
  expect_length(splits$test_idx, 2)
  expect_length(splits$train_idx, 2)
  
  # Check that all indices are covered
  all_test_idx <- sort(unlist(splits$test_idx))
  all_train_idx <- sort(unlist(splits$train_idx))
  expect_equal(length(all_test_idx), nrow(data_long))
  expect_equal(length(all_train_idx), nrow(data_long))
  
  # Check that test and train are complementary
  for (i in 1:2) {
    expect_equal(
      sort(c(splits$test_idx[[i]], splits$train_idx[[i]])),
      1:nrow(data_long)
    )
    expect_length(intersect(splits$test_idx[[i]], splits$train_idx[[i]]), 0)
  }
  
  # Check attributes
  attrs <- attributes(splits$test_idx)
  expect_equal(attrs$n_locations, 4)
  expect_equal(attrs$n_groups, 2)
  expect_equal(attrs$unique_groups, c("u_wind", "v_wind"))
  expect_equal(attrs$loc_col, c("lon", "lat"))
  expect_equal(attrs$group_col, "direction")
})

test_that("create_paired_cv_splits preserves pairs", {
  # Create test data
  data_long <- data.frame(
    x = rep(c(1, 2, 3), 2),
    y = rep(c(1, 2, 3), 2),
    group = rep(c("A", "B"), each = 3),
    value = 1:6
  )
  
  splits <- create_paired_cv_splits(
    data = data_long,
    loc_col = c("x", "y"),
    group = "group",
    k = 2,
    seed = 456
  )
  
  # Check that each fold contains complete pairs
  for (i in 1:2) {
    test_data <- data_long[splits$test_idx[[i]], ]
    
    # Get unique locations in test set
    test_locations <- unique(paste(test_data$x, test_data$y, sep = "_"))
    
    # Check that each location has both groups
    for (loc in test_locations) {
      loc_data <- test_data[paste(test_data$x, test_data$y, sep = "_") == loc, ]
      expect_equal(length(unique(loc_data$group)), 2)
      expect_setequal(unique(loc_data$group), c("A", "B"))
    }
  }
})

test_that("create_paired_cv_splits with different k values", {
  # Test with k=1 (all data in one fold)
  data_small <- data.frame(
    lon = rep(c(1, 2), 2),
    lat = rep(c(1, 2), 2),
    direction = rep(c("u", "v"), each = 2),
    value = 1:4
  )
  
  splits_k1 <- create_paired_cv_splits(
    data = data_small,
    loc_col = c("lon", "lat"),
    group = "direction",
    k = 1,
    seed = 789
  )
  
  expect_length(splits_k1$test_idx, 1)
  expect_length(splits_k1$train_idx, 1)
  expect_equal(length(splits_k1$test_idx[[1]]), 4)
  expect_equal(length(splits_k1$train_idx[[1]]), 0)
  
  # Test with k=3
  data_larger <- data.frame(
    x = rep(1:6, 2),
    y = rep(1:6, 2),
    group = rep(c("A", "B"), each = 6),
    value = 1:12
  )
  
  splits_k3 <- create_paired_cv_splits(
    data = data_larger,
    loc_col = c("x", "y"),
    group = "group",
    k = 3,
    seed = 101
  )
  
  expect_length(splits_k3$test_idx, 3)
  expect_length(splits_k3$train_idx, 3)
  
  # Check that all observations are covered exactly once in test sets
  all_test_obs <- sort(unlist(splits_k3$test_idx))
  expect_equal(all_test_obs, 1:12)
})

test_that("create_paired_cv_splits input validation", {
  # Test missing data frame
  expect_error(
    create_paired_cv_splits(
      data = NULL,
      loc_col = c("x", "y"),
      group = "group",
      k = 2
    ),
    "data must be a data frame"
  )
  
  # Test non-character loc_col
  data_test <- data.frame(x = 1, y = 1, group = "A", value = 1)
  expect_error(
    create_paired_cv_splits(
      data = data_test,
      loc_col = c(1, 2),
      group = "group",
      k = 2
    ),
    "loc_col must be character vector"
  )
  
  # Test non-character group
  expect_error(
    create_paired_cv_splits(
      data = data_test,
      loc_col = c("x", "y"),
      group = 1,
      k = 2
    ),
    "group must be character"
  )
  
  # Test invalid k
  expect_error(
    create_paired_cv_splits(
      data = data_test,
      loc_col = c("x", "y"),
      group = "group",
      k = 0
    ),
    "k must be positive integer"
  )
  
  # Test missing location columns
  expect_error(
    create_paired_cv_splits(
      data = data_test,
      loc_col = c("missing_col", "y"),
      group = "group",
      k = 2
    ),
    "location columns must exist in data"
  )
  
  # Test missing group column
  expect_error(
    create_paired_cv_splits(
      data = data_test,
      loc_col = c("x", "y"),
      group = "missing_group",
      k = 2
    ),
    "group column must exist in data"
  )
})

test_that("create_paired_cv_splits handles edge cases", {
  # Test with k > number of locations
  data_small <- data.frame(
    x = rep(c(1, 2), 2),
    y = rep(c(1, 2), 2),
    group = rep(c("A", "B"), each = 2),
    value = 1:4
  )
  
  expect_warning(
    splits <- create_paired_cv_splits(
      data = data_small,
      loc_col = c("x", "y"),
      group = "group",
      k = 5,  # More folds than locations
      seed = 222
    ),
    "Number of unique locations.*is less than k"
  )
  
  # Should automatically reduce k to number of locations
  expect_length(splits$test_idx, 2)
  expect_length(splits$train_idx, 2)
})

test_that("create_paired_cv_splits with incomplete data", {
  # Create data where some locations don't have all groups
  data_incomplete <- data.frame(
    x = c(1, 1, 2, 3, 3),
    y = c(1, 1, 2, 3, 3),
    group = c("A", "B", "A", "A", "B"),  # Location 2_2 missing group B
    value = 1:5
  )
  
  expect_warning(
    splits <- create_paired_cv_splits(
      data = data_incomplete,
      loc_col = c("x", "y"),
      group = "group",
      k = 2,
      seed = 333
    ),
    "Some locations don't have observations for all groups"
  )
  
  # Function should still work but with warning
  expect_type(splits, "list")
  expect_named(splits, c("test_idx", "train_idx"))
})

test_that("create_paired_cv_splits with single location column", {
  # Test with only one location column
  data_1d <- data.frame(
    time = rep(1:4, 2),
    group = rep(c("morning", "evening"), each = 4),
    temp = rnorm(8)
  )
  
  splits <- create_paired_cv_splits(
    data = data_1d,
    loc_col = "time",
    group = "group",
    k = 2,
    seed = 444
  )
  
  expect_type(splits, "list")
  expect_equal(attributes(splits$test_idx)$n_locations, 4)
  expect_equal(attributes(splits$test_idx)$loc_col, "time")
})

test_that("create_paired_cv_splits reproducibility", {
  # Test that same seed produces same results
  data_test <- data.frame(
    lon = rep(c(1, 2, 3, 4), 2),
    lat = rep(c(1, 1, 2, 2), 2),
    direction = rep(c("u_wind", "v_wind"), each = 4),
    wind = rnorm(8)
  )
  
  splits1 <- create_paired_cv_splits(
    data = data_test,
    loc_col = c("lon", "lat"),
    group = "direction",
    k = 2,
    seed = 555
  )
  
  splits2 <- create_paired_cv_splits(
    data = data_test,
    loc_col = c("lon", "lat"),
    group = "direction",
    k = 2,
    seed = 555
  )
  
  expect_identical(splits1$test_idx, splits2$test_idx)
  expect_identical(splits1$train_idx, splits2$train_idx)
  
  # Test that different seeds produce different results
  splits3 <- create_paired_cv_splits(
    data = data_test,
    loc_col = c("lon", "lat"),
    group = "direction",
    k = 2,
    seed = 999  # Use a much different seed
  )
  
  expect_false(identical(splits1$test_idx, splits3$test_idx))
})

test_that("create_paired_cv_splits with more than 2 groups", {
  # Test with 3 groups per location
  data_3groups <- data.frame(
    x = rep(1:3, 3),
    y = rep(1:3, 3),
    group = rep(c("A", "B", "C"), each = 3),
    value = 1:9
  )
  
  splits <- create_paired_cv_splits(
    data = data_3groups,
    loc_col = c("x", "y"),
    group = "group",
    k = 2,
    seed = 666
  )
  
  expect_equal(attributes(splits$test_idx)$n_groups, 3)
  expect_setequal(attributes(splits$test_idx)$unique_groups, c("A", "B", "C"))
  
  # Check that each test fold has complete triplets
  for (i in 1:2) {
    test_data <- data_3groups[splits$test_idx[[i]], ]
    test_locations <- unique(paste(test_data$x, test_data$y, sep = "_"))
    
    for (loc in test_locations) {
      loc_data <- test_data[paste(test_data$x, test_data$y, sep = "_") == loc, ]
      expect_equal(length(unique(loc_data$group)), 3)
      expect_setequal(unique(loc_data$group), c("A", "B", "C"))
    }
  }
})

test_that("create_paired_cv_splits with large dataset", {
  # Test performance and correctness with larger dataset
  set.seed(777)
  n_locations <- 50
  locs <- data.frame(
    lon = runif(n_locations, 0, 10),
    lat = runif(n_locations, 0, 10)
  )
  
  data_large <- data.frame(
    lon = rep(locs$lon, 2),
    lat = rep(locs$lat, 2),
    group = rep(c("group1", "group2"), each = n_locations),
    value = rnorm(n_locations * 2)
  )
  
  splits <- create_paired_cv_splits(
    data = data_large,
    loc_col = c("lon", "lat"),
    group = "group",
    k = 5,
    seed = 888
  )
  
  expect_length(splits$test_idx, 5)
  expect_equal(attributes(splits$test_idx)$n_locations, n_locations)
  
  # Verify all observations are used exactly once in test sets
  all_test_obs <- sort(unlist(splits$test_idx))
  expect_equal(all_test_obs, 1:(n_locations * 2))
  expect_equal(length(all_test_obs), n_locations * 2)
}) 