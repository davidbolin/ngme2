# The Rao-Blackwell traces are tr(QQ^-1 A' diag(d) B). Forming that triple
# product is pointless -- it is as dense as QQ and is built once per latent
# parameter per iteration, only to be handed to the probe estimator, which needs
# it applied to the probe block. It is applied as A' (d * (B QU)) instead.
#
# Same quantity, different association, so the two agree to round-off rather
# than bitwise. Setting NGME2_NO_FACTORED_TRACE restores the assembled form,
# which is what lets the two be compared here.

library(ngme2)
library(testthat)

test_that("factored and assembled Rao-Blackwell traces agree", {
  skip_on_cran()
  skip_if_not_installed("fmesher")
  # spacetime is the case that matters: non-symmetric, no exact-trace shortcut,
  # and a 2-d mesh, so the selected inverse is ruled out and the traces really
  # do go through the probe estimator.
  st <- withr::with_seed(1, {
    n <- 200
    loc <- matrix(runif(n * 2), n, 2)
    list(mesh_s = fmesher::fm_mesh_2d(loc = loc, cutoff = 0.25,
                                      max.edge = c(0.5, 1)),
         mesh_t = fmesher::fm_mesh_1d(1:5),
         data = data.frame(y = rnorm(n), c1 = loc[, 1], c2 = loc[, 2],
                           time = sample(1:5, n, TRUE)))
  })
  mesh_st <- list(st$mesh_t, st$mesh_s)
  st_data <- st$data
  ctrl <- control_opt(seed = 99, burnin = 3, iterations = 12, n_batch = 1,
                      n_parallel_chain = 1, max_num_threads = 1,
                      warn_no_convergence = FALSE, polish_iterations = 0,
                      R_hat_conv_check = FALSE, trend_std_conv_check = FALSE)
  # ngme() resolves its formula in the calling frame, so these are inlined
  # rather than wrapped in a closure.
  factored <- ngme_numeric_digest(ngme(
    y ~ f(list(time, cbind(c1, c2)), model = spacetime(mesh = mesh_st, alpha = 2),
          noise = noise_normal()),
    data = st_data, family = "normal", control_opt = ctrl))
  assembled <- with_assembled_traces(ngme_numeric_digest(ngme(
    y ~ f(list(time, cbind(c1, c2)), model = spacetime(mesh = mesh_st, alpha = 2),
          noise = noise_normal()),
    data = st_data, family = "normal", control_opt = ctrl)))
  expect_identical(names(factored), names(assembled))
  expect_equal(ngme_digest_max_diff(factored, assembled), 0, tolerance = 1e-6)
})
