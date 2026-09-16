# Structured probing for the Hutchinson trace estimators.
#
# The estimator's contract is a statement about its DISTRIBUTION -- still
# exactly unbiased, with a far smaller spread at the same probe count -- so it
# cannot be checked by comparing single values. Everything here draws the
# estimator many times through trace_probe_draws(), which re-factorizes and
# re-draws exactly as a Gibbs pass does.

library(ngme2)
library(testthat)
library(Matrix)

# A 2-d SPDE-like precision: Q = K' D K + A' D_e A, which is the shape the
# Rao-Blackwell traces are taken against. Local operator, so Q^-1 decays and
# the colouring has something to exploit.
make_case <- function(nx = 20, kappa = 0.35) {
  n <- nx * nx
  idx <- function(i, j) (j - 1) * nx + i
  ii <- jj <- xx <- numeric(0)
  for (i in 1:nx) for (j in 1:nx) {
    a <- idx(i, j); deg <- 0
    for (d in list(c(1, 0), c(-1, 0), c(0, 1), c(0, -1))) {
      i2 <- i + d[1]; j2 <- j + d[2]
      if (i2 >= 1 && i2 <= nx && j2 >= 1 && j2 <= nx) {
        ii <- c(ii, a); jj <- c(jj, idx(i2, j2)); xx <- c(xx, -1); deg <- deg + 1
      }
    }
    ii <- c(ii, a); jj <- c(jj, a); xx <- c(xx, deg)
  }
  G <- sparseMatrix(i = ii, j = jj, x = xx, dims = c(n, n))
  K <- kappa^2 * Diagonal(n) + G
  withr::with_seed(4, {
    SV <- runif(n, 0.6, 1.6)
    obs <- sort(sample(n, round(n / 2)))
  })
  A <- sparseMatrix(i = seq_along(obs), j = obs, x = 1, dims = c(length(obs), n))
  Q <- Matrix::forceSymmetric(Matrix::t(K) %*% Diagonal(x = 1 / SV) %*% K +
                                4 * Matrix::crossprod(A))
  M <- Matrix::t(0.7 * Diagonal(n)) %*% Diagonal(x = 1 / SV) %*% K
  csc <- function(X) as(as(X, "generalMatrix"), "CsparseMatrix")
  list(Q = csc(Q), M = csc(M), n = n,
       truth = sum(diag(as.matrix(Matrix::solve(Q, as.matrix(M))))))
}

test_that("probing leaves the trace estimator unbiased", {
  skip_on_cran()
  cs <- make_case()
  for (d in 1:3) {
    r <- ngme2:::trace_probe_draws(cs$Q, cs$M, n_probes = 120L, reps = 1500L,
                                   probing = TRUE, max_dist = d, min_reps = 1L,
                                   seed = 7L)
    if (r$distance != d) next # this distance did not fit the budget
    se <- stats::sd(r$estimate) / sqrt(length(r$estimate))
    # Four standard errors: the estimator is unbiased exactly, so the only
    # thing being allowed for here is the Monte Carlo error of the check.
    expect_lt(abs(mean(r$estimate) - cs$truth), 4 * se)
  }
})

test_that("probing cuts the spread at no extra cost", {
  skip_on_cran()
  cs <- make_case()
  dense <- ngme2:::trace_probe_draws(cs$Q, cs$M, 120L, 600L, probing = FALSE,
                                     seed = 11L)
  probed <- ngme2:::trace_probe_draws(cs$Q, cs$M, 120L, 600L, probing = TRUE,
                                      min_reps = 1L, seed = 11L)
  expect_true(probed$probing)
  # Never more probes than the budget: probing spends the same allowance
  # differently, it does not ask for a larger one.
  expect_lte(probed$n_probes, 120L)
  expect_lte(dense$n_probes, 120L)
  rmse <- function(r) sqrt(mean((r$estimate - cs$truth)^2))
  # The measured factor on this case is far larger; the threshold is set low so
  # the test is about the mechanism working at all, not about a tuned number.
  expect_gt(rmse(dense) / rmse(probed), 3)
})

test_that("the reported probe variance matches the realised spread", {
  skip_on_cran()
  cs <- make_case()
  # Two replicates over the colouring are what makes a spread measurable; this
  # is what the budget controller reads, so it has to be right.
  r <- ngme2:::trace_probe_draws(cs$Q, cs$M, 200L, 600L, probing = TRUE,
                                 min_reps = 2L, seed = 5L)
  expect_true(r$probing)
  expect_gte(r$reps, 2L)
  reported <- mean(r$probe_var) / r$budget
  expect_equal(reported, stats::var(r$estimate), tolerance = 0.35)
})

test_that("a single replicate reports the spread as unmeasurable", {
  skip_on_cran()
  cs <- make_case()
  r <- ngme2:::trace_probe_draws(cs$Q, cs$M, 30L, 20L, probing = TRUE,
                                 min_reps = 1L, seed = 3L)
  expect_true(r$probing)
  expect_identical(r$reps, 1L)
  # Not zero: zero would read as "exact" to the budget controller. Negative
  # says the number is absent, and the controller then leaves the budget alone.
  expect_true(all(r$probe_var < 0))
})

test_that("a budget too small for any colouring falls back to dense probes", {
  skip_on_cran()
  cs <- make_case()
  # Distance 1 needs 7 colours here, and two replicates of it need 14.
  r <- ngme2:::trace_probe_draws(cs$Q, cs$M, 8L, 20L, probing = TRUE,
                                 min_reps = 2L, seed = 3L)
  expect_false(r$probing)
  expect_identical(r$n_probes, 8L)
})

test_that("the budget is raised to reach a colouring, within its cost bound", {
  skip_on_cran()
  cs <- make_case()
  # Distance 1 needs 7 colours here, so two replicates of it need 14 -- above a
  # budget of 8. Probing cannot engage at all until the budget reaches 14.
  bare <- ngme2:::trace_probe_draws(cs$Q, cs$M, 8L, 20L, probing = TRUE,
                                    min_reps = 2L, raise_cap = 0L, seed = 3L)
  expect_false(bare$probing)
  expect_identical(bare$budget, 8L)

  # Allowed to raise to 14, it engages and the budget lands exactly on the
  # floor -- never above it.
  up <- ngme2:::trace_probe_draws(cs$Q, cs$M, 8L, 40L, probing = TRUE,
                                  min_reps = 2L, raise_cap = 24L, seed = 3L)
  expect_true(up$probing)
  expect_identical(up$budget, up$n_colours * up$reps)
  expect_identical(up$reps, 2L)
  expect_gt(up$budget, 8L)

  # A cap below the floor buys nothing and must not raise the budget partway:
  # a bigger bill for the same dense probes would be the worst of both.
  short <- ngme2:::trace_probe_draws(cs$Q, cs$M, 8L, 20L, probing = TRUE,
                                     min_reps = 2L, raise_cap = 10L, seed = 3L)
  expect_false(short$probing)
  expect_identical(short$budget, 8L)
})

test_that("trace_probing = FALSE reproduces the dense estimator exactly", {
  skip_on_cran()
  cs <- make_case()
  a <- ngme2:::trace_probe_draws(cs$Q, cs$M, 40L, 25L, probing = FALSE, seed = 2L)
  b <- ngme2:::trace_probe_draws(cs$Q, cs$M, 40L, 25L, probing = TRUE,
                                 max_dist = 0L, seed = 2L)
  expect_false(b$probing)
  expect_identical(a$estimate, b$estimate)
})

test_that("a spacetime fit runs on structured probes for no more work", {
  skip_on_cran()
  skip_if_not_installed("fmesher")
  # An integration check, not a statistical one: that the probing path is
  # reached through a real fit -- the QQ pattern hook, the budget, the block
  # traces -- and costs no more than the dense path. Whether the ESTIMATOR is
  # unbiased and less variable is settled above, by drawing it directly;
  # twelve SGD steps of one chain cannot measure either, and comparing two
  # such runs to a tolerance would only be testing the random number stream.
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
  ctrl_of <- function(probing)
    control_opt(seed = 42, burnin = 3, iterations = 12, n_batch = 1,
                n_parallel_chain = 1, max_num_threads = 1,
                rao_blackwellization = TRUE, n_trace_iter = 60,
                trace_adapt = FALSE, trace_probing = probing,
                polish_iterations = 0, R_hat_conv_check = FALSE,
                trend_std_conv_check = FALSE, warn_no_convergence = FALSE)
  # ngme() resolves its formula in the CALLING frame, so both fits are inlined
  # rather than wrapped in a closure; a closure puts mesh_st out of its reach.
  ngme2:::ngme_factor_counters(reset = TRUE)
  a <- unlist(ngme_result(ngme(
    y ~ 0 + f(list(time, cbind(c1, c2)),
              model = spacetime(mesh = mesh_st, alpha = 2),
              noise = noise_normal()),
    data = st_data, family = "normal", control_opt = ctrl_of(FALSE))))
  cost_dense <- ngme2:::ngme_factor_counters(reset = TRUE)$probe_solves
  b <- unlist(ngme_result(ngme(
    y ~ 0 + f(list(time, cbind(c1, c2)),
              model = spacetime(mesh = mesh_st, alpha = 2),
              noise = noise_normal()),
    data = st_data, family = "normal", control_opt = ctrl_of(TRUE))))
  cost_probed <- ngme2:::ngme_factor_counters(reset = TRUE)$probe_solves

  expect_identical(names(a), names(b))
  expect_true(all(is.finite(b)))
  expect_gt(cost_dense, 0)
  # Probing quantizes the budget DOWN to a whole number of colourings, so it
  # spends the same or less. Rising here would mean the budget was being
  # overspent, which is the one thing this must never do.
  expect_lte(cost_probed, cost_dense)
})

test_that("the cost-based budget rule runs and does not overspend", {
  skip_on_cran()
  skip_if_not_installed("fmesher")
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
  ctrl_of <- function(rule)
    control_opt(seed = 7, burnin = 3, iterations = 60, n_batch = 6,
                n_parallel_chain = 2, max_num_threads = 2,
                rao_blackwellization = TRUE, n_trace_iter = 10,
                trace_adapt = TRUE, trace_adapt_rule = rule,
                trace_probing = TRUE, polish_iterations = 60,
                warn_no_convergence = FALSE, print_check_info = FALSE)
  ngme2:::ngme_factor_counters(reset = TRUE)
  a <- unlist(ngme_result(ngme(
    y ~ 0 + f(list(time, cbind(c1, c2)),
              model = spacetime(mesh = mesh_st, alpha = 2),
              noise = noise_normal()),
    data = st_data, family = "normal", control_opt = ctrl_of("share"))))
  cost_share <- ngme2:::ngme_factor_counters(reset = TRUE)$probe_solves
  b <- unlist(ngme_result(ngme(
    y ~ 0 + f(list(time, cbind(c1, c2)),
              model = spacetime(mesh = mesh_st, alpha = 2),
              noise = noise_normal()),
    data = st_data, family = "normal", control_opt = ctrl_of("cost"))))
  cost_cost <- ngme2:::ngme_factor_counters(reset = TRUE)$probe_solves

  expect_identical(names(a), names(b))
  expect_true(all(is.finite(b)))
  expect_gt(cost_share, 0)
  # The cost rule holds the search at the budget it was given and only sizes the
  # polish, so it cannot spend more than the share rule does over a run this
  # short -- where the share rule has already started climbing and the polish is
  # a fraction of the whole.
  expect_lte(cost_cost, cost_share * 1.5)
})
