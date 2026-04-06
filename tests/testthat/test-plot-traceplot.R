test_that("trace transformations for stationary nu include nu_lower_bound", {
  noise <- noise_nig(
    theta_nu = 0,
    B_nu = matrix(1, 1, 1),
    nu_lower_bound = 0.7,
    fix_theta_mu = TRUE,
    fix_theta_sigma = TRUE
  )

  rep1 <- list(noise = noise, feff = numeric(0))
  attr(rep1, "block_traj") <- list(matrix(0, nrow = 1, ncol = 3))
  fake_ngme <- structure(list(replicates = list(rep1)), class = "ngme")

  tr <- get_trace_trajectories(fake_ngme, name = "data", apply_transform = TRUE)

  expect_equal(tr$parameter_names, "nu")
  expect_equal(as.numeric(tr$trajectories$nu[, 1]), rep(1.7, 3))
})
