library(ngme2)
library(testthat)

test_that("bv model with c1, c2 parameters works", {
  n_obs <- 1000
  n_each <- n_obs / 2
  group <- rep(c("W1", "W2"), n_each)
  reorder_loc <- 1:n_each

  # Simulate with c1=2, c2=0.5
  sim_fields <- simulate(
    f(
      c(reorder_loc, reorder_loc),
      model = bv(
        mesh = reorder_loc,
        theta = 0.5, rho = 0.5,
        c1 = 2, c2 = 0.5,
        use_c_param = TRUE,
        sub_models = list(W1 = ar1(rho = 0.5), W2 = ar1(rho = 0.5))
      ),
      group = group,
      noise = list(W1 = noise_normal(), W2 = noise_normal())
    )
  )[[1]]

  # Fit
  res <- ngme(
    Y ~ 0 + f(
      c(reorder_loc, reorder_loc),
      model = bv(
        mesh = reorder_loc,
        use_c_param = TRUE,
        sub_models = list(W1 = ar1(), W2 = ar1())
      ),
      group = group,
      noise = list(W1 = noise_normal(), W2 = noise_normal())
    ),
    data = data.frame(Y = sim_fields),
    control_opt = control_opt(
      estimation = TRUE,
      iterations = 10,
      print_check_info = FALSE
    )
  )

  expect_true(inherits(res, "ngme"))
  # Check parameter names
  print(res$replicates[[1]]$par_names)
  expect_true("c1" %in% res$replicates[[1]]$par_names)
  expect_true("c2" %in% res$replicates[[1]]$par_names)
})
