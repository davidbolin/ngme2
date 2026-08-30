# Regression tests for the QQ / K factorization caches.
#
# Background: QQ = K' diag(1/SV) K + (AZ)' D (AZ) is the conditional precision
# of the latent field. It depends on the current parameters (through K, Z, the
# latent sigmas and the measurement precision) and on the variance vectors V --
# and on nothing else. In particular it does not depend on W, so for a model
# whose V never moves (Gaussian latent noise and Gaussian measurement noise) it
# is constant across all Gibbs draws of an optimizer iteration.
#
# The code therefore caches QQ, its Cholesky factor and the AZ block matrix,
# invalidating them from every mutator that can change them, and reruns the
# symbolic (analyze) phase only when the sparsity pattern actually changes --
# which happens for rational approximations when the smoothness crosses an
# integer.
#
# These tests check two things:
#   1. the caches change no numbers, for every model family, by fitting each
#      model twice -- once normally and once with NGME2_DISABLE_FACTOR_CACHE
#      set, which restores the uncached "rebuild on every use" behaviour;
#   2. the caches actually skip work, via the instrumentation counters.

library(ngme2)
library(testthat)

# ---------------------------------------------------------------- fixtures --

make_ts_data <- function(n = 200, seed = 1) {
  withr::with_seed(seed, {
    x1 <- runif(n)
    x2 <- rnorm(n)
    w <- as.numeric(arima.sim(n = n, list(ar = 0.5)))
    data.frame(y = 1.5 + 2 * x1 - x2 + w + rnorm(n, sd = 0.4),
               x1 = x1, x2 = x2, idx = seq_len(n))
  })
}

make_spatial_fixture <- function(n = 150, seed = 2) {
  withr::with_seed(seed, {
    loc <- matrix(runif(n * 2), n, 2)
    z <- rnorm(n)
    list(
      mesh = fmesher::fm_mesh_2d(loc = loc, cutoff = 0.1, max.edge = c(0.25, 0.7)),
      data = data.frame(y = 1 + 0.5 * z + rnorm(n, sd = 0.3), z = z,
                        c1 = loc[, 1], c2 = loc[, 2])
    )
  })
}

# Evaluate `expr` (a fit, a prediction, a set of posterior samples) twice --
# once with the caches active and once with them disabled -- and require the
# two results to agree in every single number.
#
# `expr` is evaluated in the caller's frame rather than in a closure, because
# ngme() resolves the symbols in its formula in the frame that called it.
expect_cache_invariant <- function(expr, label = NULL) {
  q <- substitute(expr)
  pf <- parent.frame()
  cached <- ngme_numeric_digest(eval(q, pf))
  uncached <- with_factor_cache_disabled(ngme_numeric_digest(eval(q, pf)))
  expect_gt(length(cached), 0)
  expect_identical(names(cached), names(uncached), info = label)
  expect_identical(ngme_digest_max_diff(cached, uncached), 0, info = label)
}

# ------------------------------------------------- 1. numerical invariance --

test_that("caching does not change results for time series models", {
  skip_on_cran()
  d <- make_ts_data()
  n <- nrow(d)
  ctl <- ngme_reproducible_control()

  expect_cache_invariant(ngme(
    y ~ x1 + x2 + f(idx, model = ar1(), noise = noise_normal()),
    data = d, family = "normal", control_opt = ctl), "ar1, gaussian throughout")
  expect_cache_invariant(ngme(
    y ~ x1 + x2 + f(idx, model = ar1(), noise = noise_nig()),
    data = d, family = "normal", control_opt = ctl), "ar1, NIG latent noise")
  expect_cache_invariant(ngme(
    y ~ x1 + x2 + f(idx, model = ar1(), noise = noise_normal()),
    data = d, family = noise_nig(), control_opt = ctl), "ar1, NIG measurement noise")
  expect_cache_invariant(ngme(
    y ~ x1 + f(idx, model = ar1(), noise = noise_gal()),
    data = d, family = "normal", control_opt = ctl), "ar1, GAL latent noise")
  expect_cache_invariant(ngme(
    y ~ x1 + f(idx, model = ar1(), noise = noise_t()),
    data = d, family = "normal", control_opt = ctl), "ar1, t latent noise")
  expect_cache_invariant(ngme(
    y ~ x1 + f(idx, model = rw1(seq_len(n)), noise = noise_normal()),
    data = d, family = "normal", control_opt = ctl), "rw1")
  expect_cache_invariant(ngme(
    y ~ x1 + f(idx, model = rw2(seq_len(n)), noise = noise_nig()),
    data = d, family = "normal", control_opt = ctl), "rw2, NIG latent noise")
  expect_cache_invariant(ngme(
    y ~ x1 + f(idx, model = ou(seq_len(n)), noise = noise_normal()),
    data = d, family = "normal", control_opt = ctl), "ou")
  # Gaussian ARMA(1,1) against free Gaussian measurement noise is not
  # identifiable and ngme() says so; irrelevant here, this only checks that
  # caching does not change the numbers.
  expect_cache_invariant(suppressWarnings(ngme(
    y ~ x1 + f(idx, model = arma11(seq_len(n)), noise = noise_normal()),
    data = d, family = "normal", control_opt = ctl)), "arma(1,1)")
  expect_cache_invariant(ngme(
    y ~ x1 + f(idx, model = iid(seq_len(n)), noise = noise_normal()),
    data = d, family = "normal", control_opt = ctl), "iid")
})

test_that("caching does not change results without a latent field or with several", {
  skip_on_cran()
  d <- make_ts_data()
  n <- nrow(d)
  ctl <- ngme_reproducible_control()

  expect_cache_invariant(ngme(y ~ x1 + x2, data = d, family = "normal",
    control_opt = ctl), "fixed effects only, no latent field")
  expect_cache_invariant(ngme(
    y ~ x1 + f(idx, model = ar1(), name = "a", noise = noise_normal()) +
        f(idx, model = rw1(seq_len(n)), name = "b", noise = noise_normal()),
    data = d, family = "normal", control_opt = ctl), "sum of two latent fields")
})

test_that("caching does not change results under robust, Rao-Blackwell or replicates", {
  skip_on_cran()
  d <- make_ts_data()
  d_rep <- rbind(cbind(d, rep = 1L), cbind(d, rep = 2L))

  # robust = TRUE symmetrizes QQ and adds a jitter to its diagonal before
  # factorizing; the cache must reproduce that matrix, not the plainer one the
  # constructor assembles.
  expect_cache_invariant(ngme(
    y ~ x1 + x2 + f(idx, model = ar1(), noise = noise_normal()),
    data = d, family = "normal",
    control_opt = ngme_reproducible_control(robust = TRUE)), "robust = TRUE")
  expect_cache_invariant(ngme(
    y ~ x1 + x2 + f(idx, model = ar1(), noise = noise_normal()),
    data = d, family = "normal",
    control_opt = ngme_reproducible_control(rao_blackwellization = TRUE)),
    "Rao-Blackwellization, all gaussian")
  expect_cache_invariant(ngme(
    y ~ x1 + f(idx, model = ar1(), noise = noise_nig()),
    data = d, family = noise_nig(),
    control_opt = ngme_reproducible_control(rao_blackwellization = TRUE)),
    "Rao-Blackwellization, NIG throughout")
  expect_cache_invariant(ngme(
    y ~ x1 + f(idx, model = ar1(), noise = noise_normal()),
    data = d_rep, replicate = d_rep$rep, family = "normal",
    control_opt = ngme_reproducible_control()), "two replicates")
})

test_that("caching does not change results with replicates", {
  skip_on_cran()
  # Every replicate is its own BlockModel and therefore owns its caches; the
  # risk to check for is a replicate whose parameters were updated lazily
  # keeping a factor that belongs to another replicate's state.
  n <- 150
  d3 <- withr::with_seed(21, do.call(rbind, lapply(1:3, function(i) {
    x <- runif(n)
    w <- as.numeric(arima.sim(n = n, list(ar = 0.5)))
    data.frame(y = 1 + 2 * x + w + rnorm(n, sd = 0.4), x = x,
               idx = seq_len(n), rep = i)
  })))
  # unequal replicate sizes exercise a different W_size per BlockModel
  du <- withr::with_seed(22, do.call(rbind, lapply(1:3, function(i) {
    ni <- c(80, 160, 240)[i]
    x <- runif(ni)
    w <- as.numeric(arima.sim(n = ni, list(ar = 0.4)))
    data.frame(y = 1 + 2 * x + w + rnorm(ni, sd = 0.4), x = x,
               idx = seq_len(ni), rep = i)
  })))

  expect_cache_invariant(ngme(
    y ~ x + f(idx, model = ar1(), noise = noise_normal()), data = d3,
    replicate = d3$rep, family = "normal",
    control_opt = ngme_reproducible_control()), "3 replicates, gaussian")
  expect_cache_invariant(ngme(
    y ~ x + f(idx, model = ar1(), noise = noise_nig()), data = d3,
    replicate = d3$rep, family = "normal",
    control_opt = ngme_reproducible_control()), "3 replicates, NIG latent")
  expect_cache_invariant(ngme(
    y ~ x + f(idx, model = ar1(), noise = noise_normal()), data = d3,
    replicate = d3$rep, family = noise_nig(),
    control_opt = ngme_reproducible_control()), "3 replicates, NIG measurement")
  expect_cache_invariant(ngme(
    y ~ x + f(idx, model = ar1(), noise = noise_normal()), data = du,
    replicate = du$rep, family = "normal",
    control_opt = ngme_reproducible_control()), "unequal replicate sizes")
  expect_cache_invariant(ngme(
    y ~ x + f(idx, model = ar1(), noise = noise_normal()), data = d3,
    replicate = d3$rep, family = "normal",
    control_opt = ngme_reproducible_control(robust = TRUE)),
    "3 replicates, robust")
  expect_cache_invariant(ngme(
    y ~ x + f(idx, model = ar1(), noise = noise_normal()), data = d3,
    replicate = d3$rep, family = "normal",
    control_opt = ngme_reproducible_control(rao_blackwellization = TRUE)),
    "3 replicates, Rao-Blackwell")

  # sampling_strategy = "ws" only syncs the replicate it draws, so the others
  # keep a cache built from older parameters. That is correct precisely because
  # their parameters are older too -- this checks it.
  expect_cache_invariant(ngme(
    y ~ x + f(idx, model = ar1(), noise = noise_normal()), data = d3,
    replicate = d3$rep, family = "normal",
    control_opt = ngme_reproducible_control(sampling_strategy = "ws")),
    "3 replicates, weighted-sampling strategy")
  expect_cache_invariant(ngme(
    y ~ x + f(idx, model = ar1(), noise = noise_nig()), data = d3,
    replicate = d3$rep, family = "normal",
    control_opt = ngme_reproducible_control(sampling_strategy = "ws")),
    "3 replicates, NIG, weighted-sampling strategy")
  expect_cache_invariant(ngme(
    y ~ x + f(idx, model = ar1(), noise = noise_normal()), data = du,
    replicate = du$rep, family = "normal",
    control_opt = ngme_reproducible_control(sampling_strategy = "ws")),
    "unequal replicates, weighted-sampling strategy")

  # per-replicate posterior samples and prediction
  fit <- ngme(y ~ x + f(idx, model = ar1(), noise = noise_nig()), data = d3,
              replicate = d3$rep, family = "normal",
              control_opt = ngme_reproducible_control())
  for (r in 1:3) {
    expect_cache_invariant(ngme_post_samples(fit, 1, "W", replicate = r))
    expect_cache_invariant(ngme_post_samples(fit, 1, "V", replicate = r))
  }
  expect_cache_invariant(predict(fit, map = list(field1 = seq_len(n)),
    data = d3[d3$rep == 1, ], sampling_size = 25, burnin_size = 5, seed = 7))
})

test_that("caching does not change space-time or tensor product results with replicates", {
  skip_on_cran()
  skip_if_not_installed("fmesher")
  nst <- 150
  st <- withr::with_seed(25, {
    loc <- matrix(runif(nst * 3 * 2), nst * 3, 2)
    list(mesh_s = fmesher::fm_mesh_2d(loc = loc, cutoff = 0.16,
                                      max.edge = c(0.4, 1.0)),
         data = data.frame(y = rnorm(nst * 3), c1 = loc[, 1], c2 = loc[, 2],
                           time = sample(1:5, nst * 3, replace = TRUE),
                           rep = rep(1:3, each = nst)))
  })
  st_data <- st$data
  st_mesh_s <- st$mesh_s
  st_mesh_t <- fmesher::fm_mesh_1d(1:5)
  st_map <- list(st_data$time, cbind(st_data$c1, st_data$c2))
  ctl <- ngme_reproducible_control(12)

  expect_cache_invariant(ngme(
    y ~ f(st_map, model = spacetime(mesh = list(st_mesh_t, st_mesh_s), alpha = 2),
      noise = noise_normal()),
    data = st_data, replicate = st_data$rep, family = "normal",
    control_opt = ctl), "space-time, 3 replicates, gaussian")
  expect_cache_invariant(ngme(
    y ~ f(st_map, model = spacetime(mesh = list(st_mesh_t, st_mesh_s), alpha = 2),
      noise = noise_nig()),
    data = st_data, replicate = st_data$rep, family = "normal",
    control_opt = ctl), "space-time, 3 replicates, NIG")
  expect_cache_invariant(ngme(
    y ~ f(st_map, model = spacetime(mesh = list(st_mesh_t, st_mesh_s), alpha = 2),
      noise = noise_nig()),
    data = st_data, replicate = st_data$rep, family = "normal",
    control_opt = ngme_reproducible_control(12, sampling_strategy = "ws")),
    "space-time, 3 replicates, NIG, weighted sampling")
  expect_cache_invariant(ngme(
    y ~ f(st_map, model = tp(first = ar1(mesh = 1:5),
      second = matern(mesh = st_mesh_s)), noise = noise_normal()),
    data = st_data, replicate = st_data$rep, family = "normal",
    control_opt = ctl), "tensor product, 3 replicates, gaussian")
  expect_cache_invariant(ngme(
    y ~ f(st_map, model = tp(first = ar1(mesh = 1:5),
      second = matern(mesh = st_mesh_s)), noise = noise_nig()),
    data = st_data, replicate = st_data$rep, family = "normal",
    control_opt = ctl), "tensor product, 3 replicates, NIG")
})

test_that("caching does not change results for correlated noise or random effects", {
  skip_on_cran()
  d <- make_ts_data(n = 120, seed = 5)
  pairs <- rep(seq_len(nrow(d) / 2), each = 2)
  # The correlated-measurement path builds QQ from Q_eps rather than from a
  # diagonal precision, and Q_eps is refreshed whenever rho moves.
  expect_cache_invariant(ngme(
    y ~ x1 + f(idx, model = ar1(), noise = noise_normal()), data = d,
    family = noise_normal(corr_measurement = TRUE, index_corr = pairs),
    control_opt = ngme_reproducible_control()), "correlated measurement noise")

  re_d <- withr::with_seed(6, {
    ng <- 30
    ni <- 5
    g <- rep(seq_len(ng), each = ni)
    x <- runif(ng * ni)
    b0 <- rnorm(ng, sd = 0.7)
    data.frame(y = 1 + 2 * x + b0[g] + rnorm(ng * ni, sd = 0.4), x = x, g = g)
  })
  expect_cache_invariant(ngme(
    y ~ x + f(g, model = re(~1, data = re_d), noise = noise_normal()),
    data = re_d, family = "normal",
    control_opt = ngme_reproducible_control()), "random effects")
})

test_that("caching does not change results for bivariate models", {
  skip_on_cran()
  nb <- 120
  bv_d <- withr::with_seed(7, data.frame(
    y = rnorm(2 * nb), idx = rep(seq_len(nb), 2),
    type = rep(c("a", "b"), each = nb), x = runif(2 * nb)))
  expect_cache_invariant(ngme(
    y ~ x + f(idx,
      model = bv(mesh = seq_len(nb), rho = 0.3,
                 sub_models = list(a = ar1(rho = 0.4), b = ar1(rho = -0.2))),
      group = bv_d$type,
      noise = list(a = noise_normal(), b = noise_normal())),
    data = bv_d, family = "normal",
    control_opt = ngme_reproducible_control()), "bivariate ar1 / ar1")
})

test_that("caching does not change results for spatial models", {
  skip_on_cran()
  skip_if_not_installed("fmesher")
  sp <- make_spatial_fixture()
  sp_mesh <- sp$mesh
  sp_data <- sp$data
  ctl <- ngme_reproducible_control(15)

  expect_cache_invariant(ngme(
    y ~ z + f(~ c1 + c2, model = matern(mesh = sp_mesh, alpha = 2),
      noise = noise_normal()),
    data = sp_data, family = "normal", control_opt = ctl), "matern alpha = 2")
  expect_cache_invariant(ngme(
    y ~ z + f(~ c1 + c2, model = matern(mesh = sp_mesh, alpha = 2),
      noise = noise_nig()),
    data = sp_data, family = "normal", control_opt = ctl),
    "matern alpha = 2, NIG latent noise")
  expect_cache_invariant(ngme(
    y ~ z + f(~ c1 + c2, model = matern(mesh = sp_mesh, alpha = 4),
      noise = noise_normal()),
    data = sp_data, family = "normal", control_opt = ctl), "matern alpha = 4")
})

test_that("caching does not change results for space-time models", {
  skip_on_cran()
  skip_if_not_installed("fmesher")
  n <- 200
  st <- withr::with_seed(3, {
    loc <- matrix(runif(n * 2), n, 2)
    list(mesh_s = fmesher::fm_mesh_2d(loc = loc, cutoff = 0.15,
                                      max.edge = c(0.35, 0.9)),
         data = data.frame(y = rnorm(n), c1 = loc[, 1], c2 = loc[, 2],
                           time = sample(1:5, n, replace = TRUE)))
  })
  st_data <- st$data
  st_mesh_s <- st$mesh_s
  st_mesh_t <- fmesher::fm_mesh_1d(1:5)
  st_map <- list(st_data$time, cbind(st_data$c1, st_data$c2))
  ctl <- ngme_reproducible_control(12)

  # This is the family the performance regression was reported against.
  expect_cache_invariant(ngme(
    y ~ f(st_map, model = spacetime(mesh = list(st_mesh_t, st_mesh_s), alpha = 2),
      noise = noise_normal()),
    data = st_data, family = "normal", control_opt = ctl),
    "non-separable space-time, gaussian noise")
  expect_cache_invariant(ngme(
    y ~ f(st_map, model = spacetime(mesh = list(st_mesh_t, st_mesh_s), alpha = 2),
      noise = noise_nig()),
    data = st_data, family = "normal", control_opt = ctl),
    "non-separable space-time, NIG noise")
  expect_cache_invariant(ngme(
    y ~ f(st_map, model = tp(first = ar1(mesh = 1:5),
      second = matern(mesh = st_mesh_s)), noise = noise_normal()),
    data = st_data, family = "normal", control_opt = ctl),
    "tensor product ar1 x matern")
})

test_that("caching does not change results for rational (fractional) models", {
  skip_on_cran()
  skip_if_not_installed("fmesher")
  skip_if_not_installed("rSPDE")
  sp <- make_spatial_fixture()
  sp_mesh <- sp$mesh
  sp_data <- sp$data
  ctl <- ngme_reproducible_control(15)

  expect_cache_invariant(ngme(
    y ~ z + f(~ c1 + c2, model = matern(mesh = sp_mesh, alpha = 2.3,
      rational_order = 1, fix_alpha = TRUE), noise = noise_normal()),
    data = sp_data, family = "normal", control_opt = ctl),
    "rational matern, alpha fixed")

  # alpha estimated: this is the one family whose operator sparsity can change
  # during the optimization, so it is the only one that ever needs a second
  # symbolic factorization.
  expect_cache_invariant(ngme(
    y ~ z + f(~ c1 + c2, model = matern(mesh = sp_mesh, alpha = 2.3,
      rational_order = 1, fix_alpha = FALSE), noise = noise_normal()),
    data = sp_data, family = "normal", control_opt = ctl),
    "rational matern, alpha estimated")
  expect_cache_invariant(ngme(
    y ~ z + f(~ c1 + c2, model = matern(mesh = sp_mesh, alpha = 2.3,
      rational_order = 1, fix_alpha = FALSE), noise = noise_normal()),
    data = sp_data, family = "normal",
    control_opt = ngme_reproducible_control(15, robust = TRUE)),
    "rational matern, alpha estimated, robust = TRUE")
  expect_cache_invariant(ngme(
    y ~ z + f(~ c1 + c2, model = matern(mesh = sp_mesh, alpha = 2.3,
      rational_order = 1, fix_alpha = FALSE), noise = noise_nig()),
    data = sp_data, family = "normal", control_opt = ctl),
    "rational matern, alpha estimated, NIG latent noise")
})

test_that("caching does not change posterior sampling or prediction", {
  skip_on_cran()
  d <- make_ts_data()
  n <- nrow(d)
  ctl <- ngme_reproducible_control()

  for (latent_noise in list(noise_normal(), noise_nig())) {
    fit <- ngme(y ~ x1 + f(idx, model = ar1(), noise = latent_noise),
                data = d, family = "normal", control_opt = ctl)
    for (type in c("W", "V")) {
      expect_cache_invariant(ngme_post_samples(fit, 1, type))
    }
    expect_cache_invariant(predict(fit, map = list(field1 = seq_len(n)),
      data = d, sampling_size = 25, burnin_size = 5, seed = 7))
  }
})

test_that("caching does not change results on the platform default backend", {
  skip_on_cran()
  # Covers the default solver backend as well. On macOS that is Accelerate,
  # which is not bitwise reproducible even between two identical fits, so this
  # is a tolerance check rather than an exact one.
  d <- make_ts_data()
  ctl <- control_opt(seed = 4321, burnin = 8, iterations = 20, n_batch = 1,
    n_parallel_chain = 2, max_num_threads = 2, verbose = FALSE,
    R_hat_conv_check = FALSE,
    trend_std_conv_check = FALSE)
  fit_expr <- quote(ngme(
    y ~ x1 + x2 + f(idx, model = ar1(), noise = noise_normal()),
    data = d, family = "normal", control_opt = ctl))
  a <- ngme_numeric_digest(eval(fit_expr))
  b <- with_factor_cache_disabled(ngme_numeric_digest(eval(fit_expr)))
  expect_identical(names(a), names(b))
  expect_lt(ngme_digest_max_diff(a, b), 1e-6)
})

# ---------------------------------------------------- 2. the caches do work --

test_that("a gaussian model assembles and factorizes QQ once per iteration", {
  skip_on_cran()
  d <- make_ts_data()
  iterations <- 20
  burnin <- 8

  ngme_read_factor_counters()
  ngme(y ~ x1 + x2 + f(idx, model = ar1(), noise = noise_normal()), data = d,
       family = "normal",
       control_opt = ngme_reproducible_control(iterations, burnin = burnin))
  cached <- ngme_read_factor_counters()

  # Nothing entering QQ moves between the Gibbs draws of a gaussian model, so
  # QQ is assembled once per optimizer iteration -- when the parameters change
  # -- plus once for the burn-in sweeps, which all share the same QQ.
  expect_equal(unname(cached[["QQ_builds"]]), iterations + 1)

  uncached <- with_factor_cache_disabled({
    ngme_read_factor_counters()
    ngme(y ~ x1 + x2 + f(idx, model = ar1(), noise = noise_normal()), data = d,
         family = "normal",
         control_opt = ngme_reproducible_control(iterations, burnin = burnin))
    ngme_read_factor_counters()
  })
  # The uncached path rebuilds on every burn-in sweep, on every Gibbs draw of
  # every iteration, and on every draw of the posterior sample pass that ngme()
  # runs at the end. Assert the order of magnitude rather than that schedule,
  # which is an internal detail.
  expect_gt(uncached[["QQ_builds"]], 5 * cached[["QQ_builds"]])
})

test_that("non-gaussian models still refresh QQ on every draw", {
  skip_on_cran()
  d <- make_ts_data()
  iterations <- 20
  burnin <- 8

  count_nig_latent <- function() {
    ngme_read_factor_counters()
    ngme(y ~ x1 + f(idx, model = ar1(), noise = noise_nig()), data = d,
         family = "normal",
         control_opt = ngme_reproducible_control(iterations, burnin = burnin))
    ngme_read_factor_counters()
  }
  count_nig_measure <- function() {
    ngme_read_factor_counters()
    ngme(y ~ x1 + f(idx, model = ar1(), noise = noise_normal()), data = d,
         family = noise_nig(),
         control_opt = ngme_reproducible_control(iterations, burnin = burnin))
    ngme_read_factor_counters()
  }

  # V is redrawn every Gibbs sweep for NIG latent noise, so QQ genuinely
  # changes and the cache must not suppress a single rebuild: the cached and
  # uncached build counts have to agree exactly.
  latent_cached <- count_nig_latent()
  latent_uncached <- with_factor_cache_disabled(count_nig_latent())
  expect_identical(unname(latent_cached[["QQ_builds"]]),
                   unname(latent_uncached[["QQ_builds"]]))
  expect_gt(latent_cached[["QQ_builds"]], 5 * iterations)

  # Likewise when only the measurement noise is non-gaussian: noise_V moves
  # every sweep and feeds the measurement block of QQ.
  measure_cached <- count_nig_measure()
  measure_uncached <- with_factor_cache_disabled(count_nig_measure())
  expect_identical(unname(measure_cached[["QQ_builds"]]),
                   unname(measure_uncached[["QQ_builds"]]))
  expect_gt(measure_cached[["QQ_builds"]], 5 * iterations)
})

test_that("the symbolic factorization is not repeated for a fixed sparsity pattern", {
  skip_on_cran()
  d <- make_ts_data()
  ngme_read_factor_counters()
  ngme(y ~ x1 + x2 + f(idx, model = ar1(), noise = noise_normal()), data = d,
       family = "normal",
       control_opt = ngme_reproducible_control(20, robust = TRUE))
  counters <- ngme_read_factor_counters()
  # ar1 has a fixed bidiagonal K and hence a fixed QQ pattern, so the
  # fill-reducing ordering is computed once per model object and never again --
  # even under robust = TRUE, which used to force a re-analysis on every
  # update of QQ and of K.
  expect_equal(unname(counters[["QQ_analyzes"]]), 0)
  expect_lte(counters[["K_analyzes"]], 4)
})
