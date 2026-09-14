# group_cv() computes the same leave-group-out predictive that
# cross_validation() reaches by refitting each fold, but from one full-data
# chain and a rank-|I| downdate. These tests pin it against that path.

fit_small_model <- function(noise_type = c("normal", "nig"), n = 300, seed = 7) {
  noise_type <- match.arg(noise_type)
  set.seed(seed)
  h <- rep(1, n)
  K <- ngme2::ar1(seq_len(n), rho = 0.7)$K
  if (noise_type == "nig") {
    V <- ngme2::rig(n, 0.5, 0.5 * h^2, seed = seed)
    e <- -1.5 * h + 1.5 * V + sqrt(V) * rnorm(n)
  } else {
    e <- rnorm(n)
  }
  W <- as.numeric(Matrix::solve(K, e))
  dat <- data.frame(y = 1 + W + rnorm(n, 0, 0.5), t = seq_len(n))
  latent_noise <- if (noise_type == "nig") noise_nig() else noise_normal()
  ngme(
    y ~ 1 + f(t, model = ar1(seq_len(n)), noise = latent_noise),
    data = dat, family = "normal",
    control_opt = control_opt(
      iterations = 60, n_parallel_chain = 1, verbose = FALSE, seed = seed
    )
  )
}


test_that("group_cv reproduces cross_validation's leave-group-out predictive", {
  skip_on_cran()
  fit <- fit_small_model("normal")
  test_groups <- list(20L, 75L, 140L, 200L:202L)

  gcv <- group_cv(fit, test_idx = test_groups, n_gibbs_samples = 400,
                  n_burnin = 100, seed = 11)

  expect_s3_class(gcv, "ngme_group_cv")
  expect_equal(nrow(gcv), length(test_groups))
  expect_equal(gcv$size, lengths(test_groups))
  expect_true(all(is.finite(gcv$mean)))
  expect_true(all(gcv$sd > 0))

  # Reference: refit each group the way cross_validation() does.
  n_data <- attr(fit, "fit")$n_data
  ref <- vapply(test_groups, function(g) {
    cv <- cross_validation(
      fit, type = "custom",
      test_idx = list(as.integer(g)),
      train_idx = list(setdiff(seq_len(n_data), g)),
      N_sim = 2, n_gibbs_samples = 400, n_burnin = 100,
      print = FALSE, seed = 11, keep_pred = TRUE
    )
    # keep_pred hands these back as attributes, not list elements
    mean(unlist(attr(cv, "Y_1")))
  }, numeric(1))

  # A Gaussian model makes group_cv exact, so the whole discrepancy is the
  # reference's own Monte Carlo error.
  expect_equal(gcv$mean, ref, tolerance = 0.08)
})


test_that("group_cv handles a NIG field and reports diagnostics", {
  skip_on_cran()
  fit <- fit_small_model("nig")
  gcv <- group_cv(fit, test_idx = list(30L, 90L, 150L), n_gibbs_samples = 500,
                  n_burnin = 150, seed = 3)

  expect_equal(nrow(gcv), 3L)
  expect_true(all(is.finite(gcv$log_score)))
  expect_true(all(gcv$`neg.CRPS` > 0))
  expect_true(all(gcv$ess > 0))
  expect_true(all(gcv$r_eff > 0 & gcv$r_eff <= 1))
  expect_true(all(c("low_ess", "high_khat", "unreliable") %in% names(gcv)))
  # mse is the squared version of mae for singleton groups
  expect_equal(gcv$MSE, gcv$MAE^2, tolerance = 1e-10)
})


test_that("group_cv handles non-Gaussian MEASUREMENT noise on singleton groups", {
  skip_on_cran()
  # The k == 1 path marginalises the held-out noise_V over its prior; with a
  # non-normal family that prior draw is a k x n_q matrix, and it is easy to let
  # it collapse to a vector.
  set.seed(11); n <- 200
  K <- ngme2::ar1(seq_len(n), rho = 0.75)$K
  V <- ngme2::rig(n, 0.4, 0.4, seed = 11)
  W <- as.numeric(Matrix::solve(K, -2 + 2 * V + sqrt(V) * rnorm(n)))
  nV <- ngme2::rig(n, 1, 1, seed = 12)
  dat <- data.frame(y = 1 + W + 0.6 * (nV - 1) + 0.5 * sqrt(nV) * rnorm(n),
                    t = seq_len(n))
  fit <- ngme(y ~ 1 + f(t, model = ar1(seq_len(n)), noise = noise_nig()),
              data = dat, family = "nig",
              control_opt = control_opt(iterations = 60, n_parallel_chain = 1,
                                        verbose = FALSE, seed = 11))
  g <- group_cv(fit, test_idx = list(25L, 80L), n_gibbs_samples = 400,
                n_burnin = 100, seed = 2)
  expect_equal(nrow(g), 2L)
  expect_true(all(is.finite(g$mean)))
  expect_true(all(is.finite(g$log_score)))
  expect_true(all(g$sd > 0))

  # and the same through the leave-one-out entry point
  gl <- group_cv(fit, type = "loo", n_gibbs_samples = 200, n_burnin = 50,
                 seed = 2)
  expect_equal(nrow(gl), n)
  expect_false(anyNA(gl$mean))
})


test_that("group_cv runs leave-one-out over the whole data set", {
  skip_on_cran()
  fit <- fit_small_model("normal", n = 120)
  gcv <- group_cv(fit, type = "loo", n_gibbs_samples = 200, n_burnin = 50,
                  seed = 5)
  expect_equal(nrow(gcv), 120L)
  expect_equal(gcv$group, seq_len(120L))
  expect_true(all(gcv$size == 1L))
  expect_false(anyNA(gcv$log_score))
})


test_that("group_cv validates its groups", {
  skip_on_cran()
  fit <- fit_small_model("normal", n = 120)
  expect_error(group_cv(fit, test_idx = list(integer(0))), "empty group")
  expect_error(group_cv(fit, test_idx = list(500L)), "outside 1:120")
  expect_error(group_cv(fit, test_idx = list(c(3L, 3L))), "repeated index")
})


test_that("psis_weights and mcmc_ess behave on known input", {
  set.seed(1)
  # uniform weights: ESS is the sample size
  w <- ngme2:::psis_weights(rep(0, 500))
  expect_equal(w$ess, 500, tolerance = 1e-8)
  expect_equal(sum(w$w), 1, tolerance = 1e-12)

  # one dominating weight: ESS collapses towards 1
  lw <- c(50, rep(0, 499))
  expect_lt(ngme2:::psis_weights(lw)$ess, 50)

  # r_eff scales the reported ESS
  expect_equal(ngme2:::psis_weights(rep(0, 500), r_eff = 0.1)$ess, 50,
               tolerance = 1e-8)

  # an iid sequence has ESS close to its length; a slow AR(1) has far less
  expect_gt(ngme2:::mcmc_ess(rnorm(2000)), 1000)
  ar <- as.numeric(stats::filter(rnorm(2000), 0.95, method = "recursive"))
  expect_lt(ngme2:::mcmc_ess(ar), 400)
})


test_that("mixture_abs_moments matches the closed forms for one component", {
  # For a single Gaussian both moments are known exactly, and CRPS built from
  # them must equal the standard closed form.
  y <- 0.7; m <- 0.2; s <- 1.3
  z <- (y - m) / s
  mm <- ngme2:::mixture_abs_moments(y, m, s, 1)
  expect_equal(mm$e_xy, (y - m) * (2 * pnorm(z) - 1) + 2 * s * dnorm(z),
               tolerance = 1e-12)
  expect_equal(mm$e_xx, 2 * s / sqrt(pi), tolerance = 1e-12)   # E|X - X'|
  expect_equal(mm$e_xy - 0.5 * mm$e_xx,
               s * (z * (2 * pnorm(z) - 1) + 2 * dnorm(z) - 1 / sqrt(pi)),
               tolerance = 1e-10)

  # and the one-component fast path must agree with the general mixture code
  mult <- ngme2:::mixture_abs_moments(y, c(m, m), c(s, s), c(0.5, 0.5))
  expect_equal(mult$e_xy, mm$e_xy, tolerance = 1e-10)
  expect_equal(mult$e_xx, mm$e_xx, tolerance = 1e-10)
})


test_that("group_cv dispatches on model class", {
  skip_on_cran()
  g <- group_cv(fit_small_model("normal", n = 150), test_idx = list(20L, 60L),
                n_gibbs_samples = 500, n_burnin = 100, seed = 4)
  expect_equal(attr(g, "method"), "lowrank")   # Gaussian latent
  expect_true(all(g$exact))                    # ... and fully Gaussian: exact
  expect_equal(attr(g, "n_draws"), 1L)         # forced to a single draw
  expect_false(any(g$unreliable))

  gn <- group_cv(fit_small_model("nig", n = 150), test_idx = list(20L, 60L),
                 n_gibbs_samples = 400, n_burnin = 100, seed = 4)
  expect_equal(attr(gn, "method"), "exact")    # NIG latent -> exact chains
  expect_false(any(gn$exact))
  # ess is the honest MCMC effective sample size of the scored trace, not the
  # draw count: there are no importance weights, but the draws are still an
  # autocorrelated chain. Reporting ess == S claimed perfect mixing.
  expect_true(all(gn$ess > 0 & gn$ess <= 400))
  expect_true(all(is.na(gn$khat)))
})


test_that("the exact C++ path agrees with cross_validation", {
  skip_on_cran()
  fit <- fit_small_model("nig", n = 200)
  idx <- c(30L, 90L, 150L)
  g <- group_cv(fit, test_idx = as.list(idx), n_gibbs_samples = 3000,
                n_burnin = 500, seed = 6)
  n_data <- attr(fit, "fit")$n_data
  ref <- vapply(idx, function(i) {
    cv <- cross_validation(fit, type = "custom", test_idx = list(i),
                           train_idx = list(setdiff(seq_len(n_data), i)),
                           N_sim = 2, n_gibbs_samples = 3000, n_burnin = 500,
                           print = FALSE, seed = 6, keep_pred = TRUE)
    mean(unlist(attr(cv, "Y_1")))
  }, numeric(1))
  # Both are the same estimator, so the gap is only Monte Carlo error -- but
  # a relative tolerance on the MEAN is the wrong yardstick for a heavy-tailed
  # latent field. Where a held-out y is an outlier, BOTH estimators move by
  # several units from seed to seed, so a tight relative tolerance there
  # passes or fails by luck. Judge the gap against the predictive spread,
  # which is what sets that error.
  expect_lt(max(abs(g$mean - ref) / g$sd), 0.75)
})


test_that("the observation mask removes a fold exactly", {
  skip_on_cran()
  # Zeroing an observation's measurement precision must reproduce a model that
  # never saw it. Compare the exact path on n-1 observations against the same
  # model fitted without the held-out row in the first place.
  fit <- fit_small_model("normal", n = 120)
  rep1 <- fit$replicates[[1]]
  drop <- 60L
  keep <- setdiff(seq_len(120), drop)

  masked <- ngme2:::group_cv_exact_cpp(list(rep1), list(list(drop)),
                                       n = 200L, n_burnin = 100L,
                                       seed = 3L, num_threads = 1L,
                                       n_chains = 1L, chain_starts = list())
  eta_masked <- as.numeric(masked[[1]]$eta[[1]])

  sub <- rep1
  sub$X <- sub$X[keep, , drop = FALSE]
  sub$Y <- sub$Y[keep]
  sub$noise <- ngme2:::subset_noise(sub$noise, sub_idx = keep, compute_corr = TRUE)
  A_pred <- rep1$models[[1]]$A[drop, , drop = FALSE]
  for (j in seq_along(sub$models))
    sub$models[[j]]$A <- sub$models[[j]]$A[keep, , drop = FALSE]
  Ws <- ngme2:::sampling_cpp(sub, n = 200L, n_burnin = 100L,
                             posterior = TRUE, seed = 3L)[["W"]]
  eta_sub <- vapply(Ws, function(w) as.numeric(A_pred %*% w), numeric(1))

  # Different RNG streams, so compare distributions rather than draws -- and
  # compare the means in units of their own MC error. eta sits near zero here,
  # so a RELATIVE tolerance on the mean is meaningless: two short
  # autocorrelated chains straddling zero fail it however well they agree.
  se <- sqrt(stats::var(eta_masked) + stats::var(eta_sub)) / sqrt(200 / 10)
  expect_lt(abs(mean(eta_masked) - mean(eta_sub)), 3 * se)
  expect_equal(sd(eta_masked), sd(eta_sub), tolerance = 0.25)
})



test_that("group_cv reports the same scores as cross_validation on the same folds", {
  skip_on_cran()
  fit <- fit_small_model("normal", n = 300)
  n <- 300
  folds <- ngme2:::resolve_cv_type("k-fold", NULL, NULL, n, 5, 0.2, 10, 1)$test
  expect_equal(length(folds), 5L)
  expect_equal(sort(unlist(folds)), seq_len(n))       # a partition

  g <- group_cv(fit, type = "custom", test_idx = folds, n_gibbs_samples = 2000,
                n_burnin = 300, seed = 1)
  cv <- cross_validation(fit, type = "custom", test_idx = folds,
                         train_idx = lapply(folds, function(f) setdiff(seq_len(n), f)),
                         N_sim = 1, n_gibbs_samples = 2000, n_burnin = 300,
                         print = FALSE, seed = 1)
  # mean.scores is a data.frame with one row per `group` level, not a list.
  old <- as.numeric(cv$mean.scores[1, ])
  new <- c(mean(g$MAE), mean(g$MSE), mean(g$`neg.CRPS`), mean(g$`neg.sCRPS`))
  expect_equal(new, old, tolerance = 0.02)
})


test_that("group_cv supports the cross_validation fold types", {
  skip_on_cran()
  fit <- fit_small_model("normal", n = 200)
  gk <- group_cv(fit, type = "k-fold", k = 4, seed = 1)
  expect_equal(nrow(gk), 4L)
  expect_equal(sum(gk$size), 200L)

  gl <- group_cv(fit, type = "lpo", percent = 0.1, times = 5, seed = 1)
  expect_equal(nrow(gl), 5L)
  expect_true(all(gl$size == 20L))

  # folds are reproducible from the seed and leave the caller's stream alone
  set.seed(99); before <- runif(1)
  f1 <- ngme2:::resolve_cv_type("k-fold", NULL, NULL, 200, 4, 0.2, 10, 7)$test
  set.seed(99); after <- runif(1)
  expect_equal(before, after)
  f2 <- ngme2:::resolve_cv_type("k-fold", NULL, NULL, 200, 4, 0.2, 10, 7)$test
  expect_equal(f1, f2)

  expect_error(group_cv(fit, type = "custom"), "requires `test_idx`")
})


test_that("transform rescales the scores exactly as it should", {
  skip_on_cran()
  fit <- fit_small_model("normal", n = 200)
  grp <- as.list(1:40)
  gi <- group_cv(fit, type = "custom", test_idx = grp, seed = 2)
  gt <- group_cv(fit, type = "custom", test_idx = grp,
                 transform = function(x) 2 * x + 1, seed = 2)
  # A linear scale multiplies MAE and CRPS by |a| and MSE by a^2.
  expect_equal(mean(gt$MAE) / mean(gi$MAE), 2, tolerance = 0.05)
  expect_equal(mean(gt$MSE) / mean(gi$MSE), 4, tolerance = 0.08)
  expect_equal(mean(gt$`neg.CRPS`) / mean(gi$`neg.CRPS`), 2, tolerance = 0.05)
  # no Jacobian available, so the density score is not defined on a new scale
  expect_true(all(is.na(gt$log_score)))

  # a transform whose domain excludes part of the real line fails loudly
  expect_error(group_cv(fit, type = "custom", test_idx = as.list(1:5),
                        transform = log, seed = 2), "non-finite")

  # several scales at once
  gs <- group_cv(fit, type = "custom", test_idx = grp, seed = 2,
                 transform = list(raw = identity, doubled = function(x) 2 * x + 1))
  expect_named(gs, c("raw", "doubled"))
  expect_equal(mean(gs$doubled$MAE) / mean(gs$raw$MAE), 2, tolerance = 0.05)
})


test_that("several models are scored on shared folds", {
  skip_on_cran()
  fit <- fit_small_model("normal", n = 200)
  gm <- group_cv(list(a = fit, b = fit), type = "k-fold", k = 4, seed = 1)
  expect_s3_class(gm, "ngme_group_cv_list")
  expect_named(gm, c("a", "b"))
  # identical models on identical folds must give identical scores
  expect_equal(gm$a$MAE, gm$b$MAE)
  expect_equal(gm$a$group, gm$b$group)
  expect_error(group_cv(list(fit, "not a fit")), "list of them")
})


test_that("method selection can be forced and warns when it should", {
  skip_on_cran()
  fit <- fit_small_model("nig", n = 200)
  expect_equal(attr(group_cv(fit, test_idx = list(10L), n_gibbs_samples = 200,
                             n_burnin = 50, seed = 1), "method"), "exact")
  expect_warning(
    group_cv(fit, test_idx = list(10L), method = "lowrank",
             n_gibbs_samples = 200, n_burnin = 50, seed = 1),
    "degenerate")
  # a Gaussian fit keeps the low-rank path even for a big group, where it is
  # still exact, while a large group on a NIG latent field goes to the chains
  fg <- fit_small_model("normal", n = 200)
  expect_equal(attr(group_cv(fg, type = "k-fold", k = 2, seed = 1), "method"),
               "lowrank")
})


test_that("an exclusion buffer holds out more than it scores", {
  skip_on_cran()
  # INLA's group CV conditions on y_{-G_i} but scores only y_i. Without a buffer
  # a spatial leave-one-out is defeated by the neighbours, so this is the
  # configuration real analyses use -- and the one where a score averaged over
  # ALL held-out rows silently returns NA, because the buffer rows are never
  # scored. Both methods are checked: the bug was in the exact path only.
  for (nt in c("normal", "nig")) {
    fit <- fit_small_model(nt, n = 200)
    grp <- list(50L, 120L)
    tr <- lapply(grp, function(g) setdiff(seq_len(200), (g - 3):(g + 3)))
    g <- group_cv(fit, type = "custom", test_idx = grp, train_idx = tr,
                  n_gibbs_samples = 300, n_burnin = 50, seed = 1)
    expect_equal(g$size, c(1L, 1L))         # scored
    expect_equal(g$n_held_out, c(7L, 7L))   # held out
    expect_false(anyNA(g$MAE))
    expect_false(anyNA(g$`neg.CRPS`))
    expect_false(anyNA(g$`neg.sCRPS`))
    expect_true(all(g$MAE >= 0))
  }

  # a buffer must widen the predictive: less data to borrow from
  fit <- fit_small_model("normal", n = 200)
  no_buf <- group_cv(fit, type = "custom", test_idx = list(50L), seed = 1)
  with_buf <- group_cv(fit, type = "custom", test_idx = list(50L),
                       train_idx = list(setdiff(seq_len(200), 40:60)), seed = 1)
  expect_gt(with_buf$sd, no_buf$sd)

  # omitting train_idx means "train on everything not scored" -- no buffer
  expect_equal(no_buf$n_held_out, 1L)
  expect_equal(with_buf$n_held_out, 21L)

  # a scored observation cannot also be trained on
  expect_error(group_cv(fit, type = "custom", test_idx = list(50L),
                        train_idx = list(seq_len(200)), seed = 1),
               "cannot also be conditioned on")
  expect_error(group_cv(fit, type = "custom", test_idx = list(50L),
                        train_idx = list(1:10, 11:20)), "same length")
})


test_that("a position-dependent transform sees the right observation", {
  skip_on_cran()
  fit <- fit_small_model("normal", n = 150)
  grp <- as.list(c(30L, 60L, 90L))
  # scale observation i by i itself: MAE must scale by the same factor, which
  # only holds if the transform is handed the ORIGINAL data index.
  lev <- seq_len(150) / 10
  gi <- group_cv(fit, type = "custom", test_idx = grp, seed = 3)
  gt <- group_cv(fit, type = "custom", test_idx = grp,
                 transform = function(x, i) x * lev[i], seed = 3)
  expect_equal(gt$MAE, gi$MAE * lev[unlist(grp)], tolerance = 0.05)

  # and it still works with a buffer, where the index is not the group position
  ge <- group_cv(fit, type = "custom", test_idx = grp,
                 train_idx = lapply(grp, function(g)
                   setdiff(seq_len(150), (g - 2):(g + 2))),
                 transform = function(x, i) x * lev[i], seed = 3)
  expect_false(anyNA(ge$MAE))
  expect_equal(length(ge$MAE), 3L)
})


test_that("a replicate holding no groups is not an error", {
  skip_on_cran()
  # With a single fold only one replicate contains it; the others get an empty
  # group list. The low-rank C++ path used to throw there, which never showed
  # up because every full run had folds in every replicate.
  fit <- fit_small_model("normal", n = 150)
  expect_silent(g <- group_cv(fit, type = "custom", test_idx = list(40L),
                              n_gibbs_samples = 200, n_burnin = 50, seed = 1))
  expect_equal(nrow(g), 1L)
  expect_true(is.finite(g$mean))

  gn <- group_cv(fit_small_model("nig", n = 150), type = "custom",
                 test_idx = list(40L), n_gibbs_samples = 200, n_burnin = 50,
                 seed = 1)
  expect_equal(nrow(gn), 1L)
  expect_true(is.finite(gn$mean))
})


test_that("n_draw_out returns draws for the scored rows only", {
  skip_on_cran()
  fit <- fit_small_model("normal", n = 150)
  g <- group_cv(fit, type = "custom", test_idx = list(50L, 90L),
                train_idx = list(setdiff(seq_len(150), 47:53),
                                 setdiff(seq_len(150), 87:93)),
                n_gibbs_samples = 200, n_burnin = 50, seed = 1,
                n_draw_out = 100L)
  d <- attr(g, "draws")
  expect_equal(length(d), 2L)
  # one row per SCORED observation, not per held-out one
  expect_equal(unname(vapply(d, nrow, 1L)), c(1L, 1L))
  expect_equal(unname(vapply(d, ncol, 1L)), c(100L, 100L))
  expect_false(anyNA(unlist(d)))
  # and the draws must be consistent with the reported summary
  expect_equal(mean(d[[1]]), g$mean[1], tolerance = 0.15)
})


test_that("the buffered downdate matches an analytic leave-H-out at k > 1", {
  skip_on_cran()
  # The rank-1 case was pinned early, but the buffered case -- hold out H,
  # score a subset of it -- was only ever checked on the predictive MEAN. A
  # wrong k > 1 covariance would sail through that while every reported sd,
  # CRPS and sCRPS was off. Here the whole leave-H-out posterior is rebuilt
  # from scratch in R and compared to what the downdate produces.
  fit <- fit_small_model("normal", n = 300)
  m <- fit$replicates[[1]]
  op <- m$models[[1]]$operator
  # NOTE: operator$K is stale in a fitted object -- it holds K at the INITIAL
  # theta_K. Rebuilding is what C++ does internally, and skipping it silently
  # scores a rho = 0 model.
  Kf <- op$update_K(op$theta_K)
  h <- op$h
  AZ <- m$models[[1]]$A %*% op$Z
  sig <- as.numeric(exp(m$models[[1]]$noise$B_sigma %*%
                          m$models[[1]]$noise$theta_sigma))
  sig_e <- as.numeric(exp(m$noise$B_sigma %*% m$noise$theta_sigma))
  fe <- as.numeric(m$X %*% m$feff)
  prec <- 1 / sig_e^2
  Q <- Matrix::t(Kf) %*% Matrix::Diagonal(x = 1 / (sig^2 * h)) %*% Kf
  n <- nrow(AZ)

  analytic <- function(H, sc) {
    keep <- setdiff(seq_len(n), H)
    A_k <- AZ[keep, , drop = FALSE]
    QQ <- Q + Matrix::t(A_k) %*% Matrix::Diagonal(x = prec[keep]) %*% A_k
    b <- as.numeric(Matrix::t(A_k) %*% (prec[keep] * (m$Y[keep] - fe[keep])))
    a <- as.numeric(AZ[sc, ])
    c(mean = fe[sc] + sum(a * as.numeric(Matrix::solve(QQ, b))),
      sd = sqrt(sum(a * as.numeric(Matrix::solve(QQ, a))) + sig_e[sc]^2))
  }

  for (H in list(150L, 148:152, 143:157)) {
    an <- analytic(H, 150L)
    g <- group_cv(fit, type = "custom", test_idx = list(150L),
                  train_idx = list(setdiff(seq_len(300), H)), seed = 1)
    expect_equal(g$n_held_out, length(H))
    expect_equal(g$mean, an[["mean"]], tolerance = 1e-8)
    expect_equal(g$sd, an[["sd"]], tolerance = 1e-8)
  }

  # Score every member of one held-out set in turn. Reading a scored row's
  # variance means indexing a packed lower triangle, and an index that is
  # wrong from the third position on -- exactly what a dropped pair of
  # parentheses produced -- still passes every check above, because they all
  # score the member that happens to sort first. Holding out the same five
  # observations while scoring a different one each time walks the scored
  # position across the set. A fully Gaussian fit is exact from a single
  # draw, so each row's CRPS has a closed form.
  Hs <- 140:144
  gg <- group_cv(fit, type = "custom", test_idx = as.list(Hs),
                 train_idx = lapply(Hs, function(h) setdiff(seq_len(300), Hs)),
                 seed = 1)
  expect_equal(unique(gg$n_held_out), length(Hs))
  for (jj in seq_along(Hs)) {
    an <- analytic(Hs, Hs[jj])
    mu_j <- an[["mean"]]; sd_j <- an[["sd"]]
    yy <- m$Y[Hs[jj]]
    z <- (yy - mu_j) / sd_j
    e_xy <- (yy - mu_j) * (2 * stats::pnorm(z) - 1) + 2 * sd_j * stats::dnorm(z)
    e_xx <- 2 * sd_j / sqrt(pi)
    expect_equal(gg$mean[jj], mu_j, tolerance = 1e-8)
    expect_equal(gg$sd[jj], sd_j, tolerance = 1e-8)
    expect_equal(gg$MAE[jj], abs(mu_j - yy), tolerance = 1e-8)
    expect_equal(gg$`neg.CRPS`[jj], e_xy - 0.5 * e_xx, tolerance = 1e-8)
  }
})


test_that("the E|X - X'| estimator is unbiased and budget-insensitive", {
  skip_on_cran()
  # CRPS and sCRPS need E|X - X'| over the predictive mixture, an O(J^2) double
  # sum. It used to be made affordable by THINNING the mixture to a fixed
  # number of equal-weight components, picked by a systematic walk over the
  # cumulative weight in whatever order the components happened to be in. That
  # error is fixed once the cap is chosen -- common to every seed -- so
  # independent runs agreed with each other while all of them sat the same
  # distance from the truth. Sampling pairs is unbiased instead.
  set.seed(3)
  J <- 4000L
  m <- rnorm(J, 0, 2); s <- runif(J, 0.2, 1.5)
  w <- runif(J); w <- w / sum(w)
  y <- 0.7

  A <- function(d, sd) d * (2 * pnorm(d / sd) - 1) + 2 * sd * dnorm(d / sd)
  truth_xx <- sum(outer(w, w) * A(outer(m, m, "-"), sqrt(outer(s^2, s^2, "+"))))
  truth_xy <- sum(w * A(y - m, s))

  # E|X - y| is a single sum and must be exact whatever the budget.
  for (b in c(1000L, 200000L)) {
    mm <- ngme2:::mixture_abs_moments(y, m, s, w, max_comp = b, seed = 1L)
    expect_equal(mm$e_xy, truth_xy, tolerance = 1e-12)
  }
  # Below the budget the double sum is taken exactly.
  big <- ngme2:::mixture_abs_moments(y, m, s, w, max_comp = J^2 + 1L, seed = 1L)
  expect_equal(big$e_xx, truth_xx, tolerance = 1e-12)
  # Above it, sampled pairs: close, and closer as the budget grows.
  err <- vapply(c(2000L, 500000L), function(b)
    abs(ngme2:::mixture_abs_moments(y, m, s, w, max_comp = b, seed = 1L)$e_xx -
        truth_xx), 0)
  expect_lt(err[2], err[1])
  expect_lt(err[2] / truth_xx, 0.01)
  # and the error must be honest MC error -- it moves with the seed, rather
  # than being the same wrong number every time.
  v <- vapply(1:6, function(sd_)
    ngme2:::mixture_abs_moments(y, m, s, w, max_comp = 2000L, seed = sd_)$e_xx, 0)
  expect_gt(stats::sd(v), 0)
  expect_lt(abs(mean(v) - truth_xx) / truth_xx, 0.01)

  # A single-component mixture keeps its closed form.
  one <- ngme2:::mixture_abs_moments(y, 0.3, 0.8, 1, max_comp = 10L, seed = 1L)
  expect_equal(one$e_xx, 2 * 0.8 / sqrt(pi), tolerance = 1e-12)
})


test_that("the scored-row centre and variance match the per-draw factorisation", {
  skip_on_cran()
  # The k > 1 scorer used to build S x nq x k arrays of centres and variances
  # inside its Cholesky loop and then read a single slice out of them. Those
  # two quantities are scalar arithmetic -- Sig[i, i] is Vp[i, i] + sg_i^2 V_i,
  # no factor needed -- so they are now built vectorised for the scored rows
  # only. That turns on two things a typo would break: the index of (i, i) in
  # a column-major packed lower triangle, and the row-s/column-q recycling.
  k <- 7L
  S <- 5L
  nq <- 3L
  set.seed(4)
  # packed lower triangles of k x k SPD matrices, one column per draw
  Cp <- vapply(seq_len(S), function(s) {
    A <- matrix(rnorm(k * k), k, k)
    M <- crossprod(A) + diag(k)
    M[lower.tri(M, diag = TRUE)]
  }, numeric(k * (k + 1) / 2))
  Vq <- matrix(rgamma(k * S * nq, 2, 2), nrow = k, ncol = S * nq)
  fe <- rnorm(k); mu <- rnorm(k); sg <- runif(k, 0.3, 1.2)
  M <- matrix(rnorm(k * S), k, S)

  # what the loop produced, element by element
  cen_ref <- vv_ref <- array(NA_real_, c(S, nq, k))
  for (s in seq_len(S)) {
    Vp <- ngme2:::unpack_sym(Cp[, s], k)
    for (q in seq_len(nq)) {
      Vsq <- Vq[, (s - 1L) * nq + q]
      cen_ref[s, q, ] <- fe + M[, s] + mu * (Vsq - 1)
      vv_ref[s, q, ] <- diag(Vp + diag(sg^2 * Vsq, k))
    }
  }
  # what the vectorised form produces
  pos_ii <- function(i) ((i - 1L) * (2L * k - i + 2L)) %/% 2L + 1L
  for (i in seq_len(k)) {
    # the packed index really is the diagonal
    expect_equal(Cp[pos_ii(i), 3L],
                 ngme2:::unpack_sym(Cp[, 3L], k)[i, i], tolerance = 1e-12)
    Vm <- matrix(Vq[i, ], nrow = S, ncol = nq, byrow = TRUE)
    expect_equal(as.numeric(fe[i] + M[i, ] + mu[i] * (Vm - 1)),
                 as.numeric(cen_ref[, , i]), tolerance = 1e-12)
    expect_equal(as.numeric(Cp[pos_ii(i), ] + sg[i]^2 * Vm),
                 as.numeric(vv_ref[, , i]), tolerance = 1e-12)
  }
})


test_that("a cross_validation() call runs unchanged under group_cv()", {
  skip_on_cran()
  # The point of the argument list: an existing call should need nothing but
  # the function name changed. Every argument cross_validation() takes must be
  # accepted, and the ones with no counterpart here must say so rather than be
  # silently dropped.
  fit <- fit_small_model("normal", n = 150)
  n <- 150

  expect_true(all(setdiff(names(formals(cross_validation)), "...") %in%
                  names(formals(group_cv))))

  # k-fold, exactly as it would be written for cross_validation()
  g <- suppressWarnings(group_cv(fit, type = "k-fold", k = 4, seed = 2,
                                 n_gibbs_samples = 200, n_burnin = 50,
                                 print = FALSE))
  expect_equal(nrow(g), 4L)
  expect_false(anyNA(g$MAE))

  # a custom design written for cross_validation(), handed over verbatim
  te <- list(30L, 90L)
  tr <- lapply(te, function(i) setdiff(seq_len(n), (i - 2):(i + 2)))
  gc <- group_cv(fit, type = "custom", test_idx = te, train_idx = tr,
                 n_gibbs_samples = 200, n_burnin = 50, seed = 1)
  expect_equal(gc$size, c(1L, 1L))
  expect_equal(gc$n_held_out, c(5L, 5L))

  # N_sim is not "ignored": it is cross_validation()'s name for the number of
  # independent runs, which is what n_chains counts here. The routes that give
  # each fold its own chain honour it; the low-rank route, where one full-data
  # chain serves every fold, says so instead of dropping it.
  nig <- fit_small_model("nig", n = 200)
  te2 <- list(50L, 100L, 150L)
  g1 <- group_cv(nig, type = "custom", test_idx = te2, N_sim = 1,
                 n_gibbs_samples = 200, n_burnin = 50, seed = 1)
  g3 <- group_cv(nig, type = "custom", test_idx = te2, N_sim = 3,
                 n_gibbs_samples = 200, n_burnin = 50, seed = 1)
  expect_equal(attr(g1, "method"), "exact")
  expect_true(all(is.na(g1$rhat)))            # one chain: nothing to compare
  expect_false(anyNA(g3$rhat))                # three: a real between-chain spread
  expect_gt(sum(g3$ess), sum(g1$ess))         # and more draws behind the scores
  # an explicit n_chains wins over N_sim
  expect_true(all(is.na(group_cv(nig, type = "custom", test_idx = te2, N_sim = 3,
                                 n_chains = 1L, n_gibbs_samples = 100,
                                 n_burnin = 20, seed = 1,
                                 print = FALSE)$rhat)))
  # on the low-rank route it is reported, not silently dropped
  expect_warning(group_cv(fit, type = "custom", test_idx = te, N_sim = 3,
                          n_gibbs_samples = 100, n_burnin = 20, seed = 1),
                 "ignored by method")

  # Threads are spent over (fold, chain) pairs, and the answer must not depend
  # on how that loop was scheduled. It did once: reseeding the block alone left
  # each latent carrying state from whatever chain the worker ran before, so a
  # fold's draws moved with the thread count.
  s1 <- group_cv(nig, type = "custom", test_idx = te2, n_chains = 3L,
                 n_gibbs_samples = 200, n_burnin = 50, seed = 1,
                 num_threads = 1, print = FALSE)
  s4 <- group_cv(nig, type = "custom", test_idx = te2, n_chains = 3L,
                 n_gibbs_samples = 200, n_burnin = 50, seed = 1,
                 num_threads = 4, print = FALSE)
  expect_equal(s1$mean, s4$mean)
  expect_equal(s1$sd, s4$sd)
  expect_equal(s1$ess, s4$ess)

  # arguments describing machinery this function lacks
  expect_warning(group_cv(fit, type = "custom", test_idx = te, parallel = TRUE,
                          n_gibbs_samples = 100, n_burnin = 20, seed = 1),
                 "parallel")
  # ... and print = FALSE silences exactly that, without changing the answer
  q <- expect_silent(group_cv(fit, type = "custom", test_idx = te, N_sim = 3,
                              print = FALSE, n_gibbs_samples = 100,
                              n_burnin = 20, seed = 1))
  r <- suppressWarnings(group_cv(fit, type = "custom", test_idx = te, N_sim = 3,
                                 n_gibbs_samples = 100, n_burnin = 20, seed = 1))
  expect_equal(q$mean, r$mean)

  # refused rather than ignored
  expect_error(group_cv(fit, metric = function(x) x), "does not take a custom")
  expect_error(group_cv(fit, merge_groups = TRUE), "merge_groups")

  # max_num_threads is cross_validation()'s name for num_threads
  gt <- group_cv(fit, type = "custom", test_idx = te, max_num_threads = 2,
                 n_gibbs_samples = 200, n_burnin = 50, seed = 1)
  expect_equal(gt$mean, group_cv(fit, type = "custom", test_idx = te,
                                 num_threads = 2, n_gibbs_samples = 200,
                                 n_burnin = 50, seed = 1)$mean)

  # the RESULT has to be reachable the same way too: callers pull scores out
  # with `r$mean.scores[["abs.MAE"]]`, so those have to exist, be data frames
  # (a matrix does not answer `[[` by column name), and carry the same
  # transform-prefixed column names -- for one scale and for several.
  for (r in list(group_cv(fit, type = "custom", test_idx = te,
                          n_gibbs_samples = 200, n_burnin = 50, seed = 1),
                 group_cv(fit, type = "custom", test_idx = te,
                          transform = list(abs = function(z) abs(z)),
                          n_gibbs_samples = 200, n_burnin = 50, seed = 1))) {
    pre <- if (inherits(r, "ngme_group_cv")) "" else "abs."
    for (nm in paste0(pre, c("MAE", "MSE", "neg.CRPS", "neg.sCRPS"))) {
      expect_true(is.data.frame(r$mean.scores))
      expect_true(nm %in% names(r$mean.scores))
      expect_true(is.finite(as.numeric(r$mean.scores[[nm]])))
      expect_true(is.finite(as.numeric(r$sd.scores[[nm]])))
    }
  }
  # and the score table is still reachable as before
  gs <- group_cv(fit, type = "custom", test_idx = te, n_gibbs_samples = 200,
                 n_burnin = 50, seed = 1)
  expect_false(is.null(gs$ess))
  expect_false(is.null(gs$mean))
  expect_false("mean.scores" %in% names(gs))   # computed, not a column

  # keep_pred is the cross_validation() spelling of n_draw_out
  gk <- group_cv(fit, type = "custom", test_idx = te, keep_pred = TRUE,
                 n_gibbs_samples = 200, n_burnin = 50, seed = 1)
  expect_equal(length(attr(gk, "draws")), 2L)
})
