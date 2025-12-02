test_that("ngme_optimizers works", {
  expect_true(all(c("sgd", "precond_sgd", "momentum", "adagrad", "rmsprop", "adam", "adamW", "bfgs", "adaptive_gd") %in% ngme_optimizers()))
})
