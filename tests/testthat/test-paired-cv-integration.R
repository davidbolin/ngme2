test_that("create_paired_cv_splits integration with cross_validation", {
  # Create sample data for testing integration
  set.seed(42)
  data_long <- data.frame(
    lon = rep(c(1, 2, 3, 4), 2),
    lat = rep(c(1, 1, 2, 2), 2),
    direction = rep(c("u_wind", "v_wind"), each = 4),
    wind = rnorm(8, mean = 5, sd = 2)
  )
  
  # Create paired CV splits
  cv_splits <- create_paired_cv_splits(
    data = data_long,
    loc_col = c("lon", "lat"),
    group = "direction",
    k = 2,
    seed = 42
  )
  
  # Test that the splits have the correct structure for cross_validation
  expect_type(cv_splits, "list")
  expect_named(cv_splits, c("test_idx", "train_idx"))
  expect_length(cv_splits$test_idx, 2)
  expect_length(cv_splits$train_idx, 2)
  
  # Test that we can use these splits with cross_validation parameters
  # This is a structural test - we verify the format is correct for the CV function
  for (i in 1:2) {
    # Each test_idx should be a vector of integers
    expect_type(cv_splits$test_idx[[i]], "integer")
    expect_type(cv_splits$train_idx[[i]], "integer")
    
    # Test and train should be disjoint
    expect_length(intersect(cv_splits$test_idx[[i]], cv_splits$train_idx[[i]]), 0)
    
    # Test and train should cover all data
    expect_setequal(
      c(cv_splits$test_idx[[i]], cv_splits$train_idx[[i]]),
      1:nrow(data_long)
    )
    
    # Test that paired observations are in the same fold
    test_data <- data_long[cv_splits$test_idx[[i]], ]
    test_locations <- unique(paste(test_data$lon, test_data$lat, sep = "_"))
    
    for (loc in test_locations) {
      loc_data <- test_data[paste(test_data$lon, test_data$lat, sep = "_") == loc, ]
      # Each location should have both wind components
      expect_setequal(unique(loc_data$direction), c("u_wind", "v_wind"))
    }
  }
  
  # Verify the format matches what cross_validation expects
  # These parameters would be passed to cross_validation like this:
  # cross_validation(
  #   models,
  #   type = "custom",
  #   test_idx = cv_splits$test_idx,
  #   train_idx = cv_splits$train_idx,
  #   merge_groups = TRUE,
  #   merged_group_name = "wind_vector"
  # )
  
  # Test with merge_groups context - verify pairing is preserved
  all_pairs_preserved <- TRUE
  for (i in 1:2) {
    test_data <- data_long[cv_splits$test_idx[[i]], ]
    # Group by location and check each has both components
    location_groups <- split(test_data, paste(test_data$lon, test_data$lat, sep = "_"))
    for (loc_group in location_groups) {
      if (length(unique(loc_group$direction)) != 2) {
        all_pairs_preserved <- FALSE
        break
      }
    }
    if (!all_pairs_preserved) break
  }
  expect_true(all_pairs_preserved, "All location pairs should be preserved in CV splits")
})

test_that("create_paired_cv_splits validation for different data structures", {
  # Test with time series data (1D locations)
  data_time <- data.frame(
    time = rep(1:6, 2),
    component = rep(c("real", "imaginary"), each = 6),
    signal = rnorm(12)
  )
  
  splits_time <- create_paired_cv_splits(
    data = data_time,
    loc_col = "time",
    group = "component",
    k = 3,
    seed = 123
  )
  
  expect_length(splits_time$test_idx, 3)
  expect_equal(attributes(splits_time$test_idx)$n_locations, 6)
  expect_equal(attributes(splits_time$test_idx)$n_groups, 2)
  
  # Test with geographic data (multiple location columns)
  data_geo <- data.frame(
    longitude = rep(seq(-10, 10, length.out = 5), each = 3),  # each location repeated 3 times
    latitude = rep(seq(40, 50, length.out = 5), each = 3),   # each location repeated 3 times
    measurement_type = rep(c("temperature", "humidity", "pressure"), times = 5),  # 3 types for each location
    value = rnorm(15)
  )
  
  splits_geo <- create_paired_cv_splits(
    data = data_geo,
    loc_col = c("longitude", "latitude"),
    group = "measurement_type",
    k = 2,
    seed = 456
  )
  
  expect_length(splits_geo$test_idx, 2)
  expect_equal(attributes(splits_geo$test_idx)$n_locations, 5)
  expect_equal(attributes(splits_geo$test_idx)$n_groups, 3)
  expect_equal(attributes(splits_geo$test_idx)$loc_col, c("longitude", "latitude"))
})

test_that("create_paired_cv_splits handles real-world wind data format", {
  # Simulate real wind data format
  set.seed(789)
  n_stations <- 10
  
  # Create station locations
  stations <- data.frame(
    station_id = 1:n_stations,
    lon = runif(n_stations, -180, 180),
    lat = runif(n_stations, -90, 90)
  )
  
  # Create wind observations for each station
  wind_data <- data.frame(
    station_id = rep(stations$station_id, 2),
    lon = rep(stations$lon, 2),
    lat = rep(stations$lat, 2),
    component = rep(c("u_wind", "v_wind"), each = n_stations),
    wind_speed = rnorm(n_stations * 2, mean = 5, sd = 3),
    time = as.POSIXct("2023-01-01") + rep(1:n_stations, 2) * 3600
  )
  
  # Test creating splits
  splits_wind <- create_paired_cv_splits(
    data = wind_data,
    loc_col = c("lon", "lat"),
    group = "component",
    k = 3,
    seed = 789
  )
  
  # Verify structure
  expect_length(splits_wind$test_idx, 3)
  expect_equal(attributes(splits_wind$test_idx)$n_locations, n_stations)
  expect_equal(attributes(splits_wind$test_idx)$n_groups, 2)
  
  # Verify that wind components stay paired by station
  for (fold in 1:3) {
    test_data <- wind_data[splits_wind$test_idx[[fold]], ]
    
    # Check that each station has both wind components
    station_components <- aggregate(
      component ~ station_id, 
      data = test_data, 
      FUN = function(x) length(unique(x))
    )
    
    # All stations in test set should have both components
    expect_true(all(station_components$component == 2),
                paste("Fold", fold, "should have paired wind components for all stations"))
  }
  
  # Verify no station appears in multiple test folds
  all_test_stations <- list()
  for (fold in 1:3) {
    test_data <- wind_data[splits_wind$test_idx[[fold]], ]
    all_test_stations[[fold]] <- unique(test_data$station_id)
  }
  
  # Check that test sets are disjoint by station
  for (i in 1:2) {
    for (j in (i+1):3) {
      overlap <- intersect(all_test_stations[[i]], all_test_stations[[j]])
      expect_length(
        overlap,
        0
      )
    }
  }
}) 