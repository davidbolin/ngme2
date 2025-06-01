# Test mesh list functionality for different replicates
library(testthat)

test_that("Basic mesh list functionality with 2 replicates", {
  # Setup data
  n_obs_per_rep <- 10
  n_reps <- 2
  total_obs <- n_obs_per_rep * n_reps
  
  # Generate locations
  rep1_locs <- seq(0, 1, length.out = n_obs_per_rep)
  rep2_locs <- seq(0, 2, length.out = n_obs_per_rep)
  locs <- c(rep1_locs, rep2_locs)
  
  # Create replicate factor
  replicates <- factor(rep(1:n_reps, each = n_obs_per_rep))
  
  # Generate response data
  set.seed(123)
  Y <- rnorm(total_obs)
  
  # Create data frame
  data <- data.frame(
    Y = Y,
    locs = locs,
    replicate = replicates
  )
  
  # Create different meshes for each replicate
  mesh1 <- fmesher::fm_mesh_1d(loc = seq(0, 1, length.out = 5))
  mesh2 <- fmesher::fm_mesh_1d(loc = seq(0, 2, length.out = 8))
  mesh_list <- list(mesh1, mesh2)
  
  # Test the functionality
  result <- ngme(
    formula = Y ~ 1 + f(locs, model = "matern", mesh = mesh_list),
    data = data,
    replicate = replicates,
    control_opt = control_opt(estimation = FALSE)
  )
  
  # Verify results
  expect_equal(result$n_repls, 2)
  expect_s3_class(result, "ngme")
  
  # Check that each replicate has the correct mesh
  for (i in 1:result$n_repls) {
    mesh_used <- result$replicates[[i]]$models[[1]]$mesh
    expected_n_nodes <- mesh_list[[i]]$n
    actual_n_nodes <- mesh_used$n
    
    expect_equal(actual_n_nodes, expected_n_nodes, 
                 info = paste("Replicate", i, "should use correct mesh"))
  }
})

test_that("Three replicates with different meshes", {
  # Setup data for 3 replicates
  n_obs_per_rep <- 10
  n_reps_3 <- 3
  total_obs_3 <- n_obs_per_rep * n_reps_3
  
  # Generate locations for 3 replicates
  rep1_locs <- seq(0, 1, length.out = n_obs_per_rep)
  rep2_locs <- seq(0, 2, length.out = n_obs_per_rep)  
  rep3_locs <- seq(0, 3, length.out = n_obs_per_rep)
  locs_3 <- c(rep1_locs, rep2_locs, rep3_locs)
  
  replicates_3 <- factor(rep(1:n_reps_3, each = n_obs_per_rep))
  set.seed(456)
  Y_3 <- rnorm(total_obs_3)
  
  data_3 <- data.frame(
    Y = Y_3,
    locs = locs_3,
    replicate = replicates_3
  )
  
  # Create three different meshes
  mesh1_3 <- fmesher::fm_mesh_1d(loc = seq(0, 1, length.out = 4))
  mesh2_3 <- fmesher::fm_mesh_1d(loc = seq(0, 2, length.out = 6))
  mesh3_3 <- fmesher::fm_mesh_1d(loc = seq(0, 3, length.out = 10))
  mesh_list_3 <- list(mesh1_3, mesh2_3, mesh3_3)
  
  # Test functionality
  result_3 <- ngme(
    formula = Y ~ 1 + f(locs, model = "matern", mesh = mesh_list_3),
    data = data_3,
    replicate = replicates_3,
    control_opt = control_opt(estimation = FALSE)
  )
  
  # Verify results
  expect_equal(result_3$n_repls, 3)
  expect_s3_class(result_3, "ngme")
  
  # Check that each replicate has the correct mesh
  for (i in 1:result_3$n_repls) {
    mesh_used <- result_3$replicates[[i]]$models[[1]]$mesh
    expected_n_nodes <- mesh_list_3[[i]]$n
    actual_n_nodes <- mesh_used$n
    
    expect_equal(actual_n_nodes, expected_n_nodes,
                 info = paste("Replicate", i, "should use correct mesh"))
  }
})

test_that("Error handling for insufficient meshes", {
  # Setup data for 3 replicates
  n_obs_per_rep <- 10
  n_reps_3 <- 3
  total_obs_3 <- n_obs_per_rep * n_reps_3
  
  rep1_locs <- seq(0, 1, length.out = n_obs_per_rep)
  rep2_locs <- seq(0, 2, length.out = n_obs_per_rep)  
  rep3_locs <- seq(0, 3, length.out = n_obs_per_rep)
  locs_3 <- c(rep1_locs, rep2_locs, rep3_locs)
  
  replicates_3 <- factor(rep(1:n_reps_3, each = n_obs_per_rep))
  set.seed(789)
  Y_3 <- rnorm(total_obs_3)
  
  data_3 <- data.frame(
    Y = Y_3,
    locs = locs_3,
    replicate = replicates_3
  )
  
  # Create only 2 meshes for 3 replicates (should fail)
  mesh1 <- fmesher::fm_mesh_1d(loc = seq(0, 1, length.out = 5))
  mesh2 <- fmesher::fm_mesh_1d(loc = seq(0, 2, length.out = 8))
  insufficient_mesh_list <- list(mesh1, mesh2)
  
  # Should throw an error about insufficient meshes
  expect_error(
    ngme(
      formula = Y ~ 1 + f(locs, model = "matern", mesh = insufficient_mesh_list),
      data = data_3,
      replicate = replicates_3,
      control_opt = control_opt(estimation = FALSE)
    ),
    "Not enough meshes provided"
  )
})

test_that("Single mesh compatibility maintained", {
  # Setup data
  n_obs_per_rep <- 10
  n_reps <- 2
  total_obs <- n_obs_per_rep * n_reps
  
  rep1_locs <- seq(0, 1, length.out = n_obs_per_rep)
  rep2_locs <- seq(0, 2, length.out = n_obs_per_rep)
  locs <- c(rep1_locs, rep2_locs)
  
  replicates <- factor(rep(1:n_reps, each = n_obs_per_rep))
  set.seed(101)
  Y <- rnorm(total_obs)
  
  data <- data.frame(
    Y = Y,
    locs = locs,
    replicate = replicates
  )
  
  # Test with single mesh (should work as before)
  single_mesh <- fmesher::fm_mesh_1d(loc = seq(0, 2, length.out = 12))
  
  result_single <- ngme(
    formula = Y ~ 1 + f(locs, model = "matern", mesh = single_mesh),
    data = data,
    replicate = replicates,
    control_opt = control_opt(estimation = FALSE)
  )
  
  # Verify results
  expect_equal(result_single$n_repls, 2)
  expect_s3_class(result_single, "ngme")
  
  # Both replicates should use the same mesh
  for (i in 1:result_single$n_repls) {
    mesh_used <- result_single$replicates[[i]]$models[[1]]$mesh
    expect_equal(mesh_used$n, single_mesh$n,
                 info = paste("Replicate", i, "should use single mesh correctly"))
  }
})

test_that("Different model types work with mesh lists (RW1)", {
  # Setup data
  n_obs_per_rep <- 10
  n_reps <- 2
  total_obs <- n_obs_per_rep * n_reps
  
  rep1_locs <- seq(0, 1, length.out = n_obs_per_rep)
  rep2_locs <- seq(0, 2, length.out = n_obs_per_rep)
  locs <- c(rep1_locs, rep2_locs)
  
  replicates <- factor(rep(1:n_reps, each = n_obs_per_rep))
  set.seed(202)
  Y <- rnorm(total_obs)
  
  data <- data.frame(
    Y = Y,
    locs = locs,
    replicate = replicates
  )
  
  # Create different meshes for RW1 model
  mesh1_rw1 <- fmesher::fm_mesh_1d(loc = seq(0, 1, length.out = 7))
  mesh2_rw1 <- fmesher::fm_mesh_1d(loc = seq(0, 2, length.out = 12))
  mesh_list_rw1 <- list(mesh1_rw1, mesh2_rw1)
  
  # Test RW1 model with mesh list
  result_rw1 <- ngme(
    formula = Y ~ 1 + f(locs, model = "rw1", mesh = mesh_list_rw1),
    data = data,
    replicate = replicates,
    control_opt = control_opt(estimation = FALSE)
  )
  
  # Verify results
  expect_equal(result_rw1$n_repls, 2)
  expect_s3_class(result_rw1, "ngme")
  
  # Check that each replicate has the correct mesh
  for (i in 1:result_rw1$n_repls) {
    mesh_used <- result_rw1$replicates[[i]]$models[[1]]$mesh
    expected_n_nodes <- mesh_list_rw1[[i]]$n
    actual_n_nodes <- mesh_used$n
    
    expect_equal(actual_n_nodes, expected_n_nodes,
                 info = paste("Replicate", i, "should use correct mesh for RW1"))
  }
})

test_that("Different model types work with mesh lists (OU)", {
  # Setup data
  n_obs_per_rep <- 10
  n_reps <- 2
  total_obs <- n_obs_per_rep * n_reps
  
  rep1_locs <- seq(0, 1, length.out = n_obs_per_rep)
  rep2_locs <- seq(0, 2, length.out = n_obs_per_rep)
  locs <- c(rep1_locs, rep2_locs)
  
  replicates <- factor(rep(1:n_reps, each = n_obs_per_rep))
  set.seed(303)
  Y <- rnorm(total_obs)
  
  data <- data.frame(
    Y = Y,
    locs = locs,
    replicate = replicates
  )
  
  # Create different meshes for OU model
  mesh1_ou <- fmesher::fm_mesh_1d(loc = seq(0, 1, length.out = 6))
  mesh2_ou <- fmesher::fm_mesh_1d(loc = seq(0, 2, length.out = 9))
  mesh_list_ou <- list(mesh1_ou, mesh2_ou)
  
  # Test OU model with mesh list
  result_ou <- ngme(
    formula = Y ~ 1 + f(locs, model = "ou", mesh = mesh_list_ou),
    data = data,
    replicate = replicates,
    control_opt = control_opt(estimation = FALSE)
  )
  
  # Verify results
  expect_equal(result_ou$n_repls, 2)
  expect_s3_class(result_ou, "ngme")
  
  # Check that each replicate has the correct mesh
  for (i in 1:result_ou$n_repls) {
    mesh_used <- result_ou$replicates[[i]]$models[[1]]$mesh
    expected_n_nodes <- mesh_list_ou[[i]]$n
    actual_n_nodes <- mesh_used$n
    
    expect_equal(actual_n_nodes, expected_n_nodes,
                 info = paste("Replicate", i, "should use correct mesh for OU"))
  }
})

test_that("Mesh list with non-Gaussian noise", {
  # Setup data
  n_obs_per_rep <- 10
  n_reps <- 2
  total_obs <- n_obs_per_rep * n_reps
  
  rep1_locs <- seq(0, 1, length.out = n_obs_per_rep)
  rep2_locs <- seq(0, 2, length.out = n_obs_per_rep)
  locs <- c(rep1_locs, rep2_locs)
  
  replicates <- factor(rep(1:n_reps, each = n_obs_per_rep))
  set.seed(404)
  Y <- rnorm(total_obs)
  
  data <- data.frame(
    Y = Y,
    locs = locs,
    replicate = replicates
  )
  
  # Create different meshes
  mesh1 <- fmesher::fm_mesh_1d(loc = seq(0, 1, length.out = 5))
  mesh2 <- fmesher::fm_mesh_1d(loc = seq(0, 2, length.out = 8))
  mesh_list <- list(mesh1, mesh2)
  
  # Test with NIG noise
  result_nig <- ngme(
    formula = Y ~ 1 + f(locs, model = "matern", mesh = mesh_list, noise = noise_nig()),
    data = data,
    replicate = replicates,
    control_opt = control_opt(estimation = FALSE)
  )
  
  # Verify results
  expect_equal(result_nig$n_repls, 2)
  expect_s3_class(result_nig, "ngme")
  
  # Check that each replicate has the correct mesh
  for (i in 1:result_nig$n_repls) {
    mesh_used <- result_nig$replicates[[i]]$models[[1]]$mesh
    expected_n_nodes <- mesh_list[[i]]$n
    actual_n_nodes <- mesh_used$n
    
    expect_equal(actual_n_nodes, expected_n_nodes,
                 info = paste("Replicate", i, "should use correct mesh with NIG noise"))
  }
}) 