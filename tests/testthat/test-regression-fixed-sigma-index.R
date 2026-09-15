# theta_sigma keeps every column of B_sigma, but the optimiser carries only the
# unfixed components, so everything indexed by the parameter vector -- the
# gradient, the Rao-Blackwell traces and the preconditioner blocks -- has to be
# mapped through the free columns rather than the leading ones. The two agree
# whenever the free components come first, which is the only arrangement the
# rest of the suite exercises, so this builds the case that separates them.
#
# Swapping the columns of B_sigma together with the fix flags and the values
# describes the SAME model, so the free component must come back identical.

library(ngme2)
library(testthat)

test_that("a fixed theta_sigma before a free one does not change the fit", {
  skip_on_cran()
  n <- 600L
  grp <- rep(c(0, 1), length.out = n)
  idx <- seq_len(n)
  d <- withr::with_seed(4, {
    W <- as.numeric(simulate(f(map = idx, model = ar1(mesh = idx, rho = 0.6),
                               noise = noise_normal(sigma = 1)), seed = 4)[[1]])
    data.frame(Y = W + rnorm(n, 0, exp(0.3 * (1 - grp) - 0.2 * grp)))
  })
  ctrl <- control_opt(seed = 7, iterations = 150, n_parallel_chain = 2,
                      warn_no_convergence = FALSE, optimizer = precond_sgd())

  free_first <- unlist(ngme_result(ngme(
    Y ~ 0 + f(map = idx, model = ar1(mesh = idx), name = "ar"), data = d,
    family = noise_normal(B_sigma = cbind(1 - grp, grp),
                          theta_sigma = c(0.3, -0.2),
                          fix_theta_sigma = c(FALSE, TRUE)),
    control_opt = ctrl)))

  fixed_first <- unlist(ngme_result(ngme(
    Y ~ 0 + f(map = idx, model = ar1(mesh = idx), name = "ar"), data = d,
    family = noise_normal(B_sigma = cbind(grp, 1 - grp),
                          theta_sigma = c(-0.2, 0.3),
                          fix_theta_sigma = c(TRUE, FALSE)),
    control_opt = ctrl)))

  # The reported names are positional, so the free component is the first slot
  # in one arrangement and the second in the other.
  expect_equal(unname(free_first[["data.theta_sigma1"]]),
               unname(fixed_first[["data.theta_sigma2"]]))
  expect_equal(unname(free_first[["ar.rho"]]), unname(fixed_first[["ar.rho"]]))
  expect_equal(unname(free_first[["ar.sigma"]]),
               unname(fixed_first[["ar.sigma"]]))

  # The fixed component must be left exactly where it was put.
  expect_equal(unname(free_first[["data.theta_sigma2"]]), -0.2)
  expect_equal(unname(fixed_first[["data.theta_sigma1"]]), -0.2)
})
