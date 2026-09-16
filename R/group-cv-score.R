# Scoring for group_cv(): turn the per-draw leave-group-out Gaussians returned
# by group_cv_cpp() into predictive scores.
#
# The C++ side returns, for each group and draw, the mean and covariance of
# eta_I = (AZ W)_I under p(W | y_-I, V). Both are free of the group's own
# measurement mixing variable, so the predictive is
#
#   y_I | y_-I, V  ~  N( X_I beta + m_s + mu_I (nV_I - 1),  Vp_s + sigma_I^2 nV_I )
#
# marginalised over the prior of nV_I -- once y_I is removed there is nothing
# left to inform it. Averaging over draws with weights w_s proportional to
# 1 / p(y_I | y_-I, V_s) gives the leave-group-out predictive.

`%||%` <- function(a, b) if (is.null(a)) b else a

log_sum_exp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}


# Initial-positive-sequence estimate of the effective sample size of an MCMC
# chain. Weight-based ESS assumes independent draws; these are not, and at an
# influential observation the two differ by more than an order of magnitude.
mcmc_ess <- function(x) {
  n <- length(x)
  if (n < 8 || stats::sd(x) == 0) return(as.numeric(n))
  a <- stats::acf(x - mean(x), lag.max = min(n - 1L, 500L),
                  plot = FALSE, demean = FALSE)$acf
  a <- as.numeric(a)
  s <- 0
  k <- 2L
  while (k + 1L <= length(a) && a[k] + a[k + 1L] > 0) {
    s <- s + a[k] + a[k + 1L]
    k <- k + 2L
  }
  max(1, n / (1 + 2 * s))
}


# Generalised Pareto fit to the tail of the weights (Zhang & Stephens
# empirical-Bayes estimator, as used by PSIS).
gpd_fit <- function(x) {
  n <- length(x)
  m <- 30L + floor(sqrt(n))
  bs <- 1 - sqrt(m / (seq_len(m) - 0.5))
  bs <- bs / (3 * x[max(1L, floor(n / 4 + 0.5))]) + 1 / x[n]
  ks <- vapply(bs, function(b) mean(log1p(-b * x)), numeric(1))
  ok <- is.finite(ks) & ks < 0 & is.finite(bs)
  if (!any(ok)) return(list(k = NA_real_, sigma = NA_real_))
  bs <- bs[ok]; ks <- ks[ok]
  L <- n * (log(-bs / ks) - ks - 1)
  w <- exp(L - max(L)); w <- w / sum(w)
  b <- sum(bs * w)
  k <- mean(log1p(-b * x))
  if (!is.finite(k) || k >= 1) return(list(k = k, sigma = NA_real_))
  list(k = k, sigma = -k / b)
}

qgpd <- function(p, k, sigma) {
  if (abs(k) < 1e-10) -sigma * log1p(-p) else sigma * expm1(-k * log1p(-p)) / k
}


#' Pareto-smoothed importance weights
#'
#' @param lw unnormalised log weights.
#' @param r_eff relative efficiency of the draws (chain ESS / number of draws).
#'   The tail length and the reported ESS both scale with it.
#' @return list with normalised weights `w`, the Pareto shape `k`, and the
#'   effective sample size `ess`.
#' @noRd
psis_weights <- function(lw, r_eff = 1) {
  S <- length(lw)
  ess_of <- function(w) r_eff / sum(w^2)
  plain <- function(k) {
    z <- lw - max(lw); w <- exp(z); w <- w / sum(w)
    list(w = w, k = k, ess = ess_of(w))
  }
  M <- min(floor(S / 5), ceiling(3 * sqrt(S / max(r_eff, 1e-8))))
  if (!is.finite(M) || M < 5 || M >= S) return(plain(NA_real_))

  ord <- order(lw)
  cutoff <- lw[ord][S - M]
  tail_idx <- ord[(S - M + 1L):S]
  x <- exp(lw[tail_idx] - cutoff) - 1
  o <- order(x); xs <- x[o]
  if (!all(is.finite(xs)) || xs[M] < 1e-10) return(plain(NA_real_))

  fit <- gpd_fit(xs)
  if (!is.finite(fit$sigma)) return(plain(fit$k))
  lw[tail_idx[o]] <- pmin(
    log1p(qgpd((seq_len(M) - 0.5) / M, fit$k, fit$sigma)) + cutoff, max(lw)
  )
  z <- lw - max(lw); w <- exp(z); w <- w / sum(w)
  list(w = w, k = fit$k, ess = ess_of(w))
}



# Exact moments of the measurement mixing variable, where the family has them.
#
# The predictive variance of a held-out observation is dominated by
# sigma^2 E[V], so estimating E[V] by averaging n_q prior draws puts the
# noisiest part of the answer on the smallest sample. At small nu the mixing
# prior is heavy tailed enough that a handful of draws misses its mean by
# enough to move the reported sd, and the error is not a Monte Carlo error
# that more Gibbs draws remove. Both NIG and GAL draw V with mean h and
# variance h/nu, so the mixture moments can be written down instead:
#
#   E[y_i | s]   = fe_i + m_s,i          (the mu (V - h) term is mean zero)
#   Var[y_i | s] = Vp_ii,s + mu_i^2 Var(V) + sigma_i^2 E[V]
#
# The t families mix on an inverse gamma whose variance is infinite for
# nu <= 4, so they keep the sampled mixture. n_q still controls the SCORES,
# which need the whole mixture and not just its first two moments.
noise_V_moments <- function(noise_type, nu, h = 1) {
  if (identical(noise_type, "normal")) return(list(EV = rep(1, length(nu)), VarV = rep(0, length(nu))))
  if (noise_type %in% c("nig", "normal_nig", "gal"))
    return(list(EV = rep(h, length(nu)), VarV = h / nu))
  NULL
}

# Prior draws of the held-out group's measurement mixing variable, FRESH for
# every Gibbs draw.
#
# Reusing one set of n_q draws across all S Gibbs draws (as this did) makes the
# marginalisation error a fixed offset that never averages away: raising S does
# nothing and only raising n_q helps, and at small nu the mixing prior is heavy
# tailed enough for that offset to be large. Drawing per Gibbs draw gives
# S * n_q distinct values and the error shrinks with S, matching what
# cross_validation() does. Note that the predictive moments no longer depend
# on these draws at all where the family has closed-form moments; n_q governs
# the scores and the weights.
prior_noise_V <- function(noise_type, k, mu, sigma, nu, n_q, single_V, seed) {
  if (identical(noise_type, "normal")) {
    return(matrix(1, nrow = k, ncol = 1L))
  }
  # One call for the whole k x n_q block rather than one per column.
  if (single_V) {
    # One shared V per column, so draw n_q scalars and broadcast down the rows.
    e <- simulate_noise(noise_type, rep(1, n_q), rep(mu[1], n_q),
                        rep(sigma[1], n_q), rep(nu[1], n_q),
                        seed = seed, single_V = FALSE)
    return(matrix(rep(as.numeric(attr(e, "V")), each = k), nrow = k, ncol = n_q))
  }
  e <- simulate_noise(noise_type, rep(1, k * n_q), rep(mu, n_q),
                      rep(sigma, n_q), rep(nu, n_q),
                      seed = seed, single_V = FALSE)
  matrix(as.numeric(attr(e, "V")), nrow = k, ncol = n_q)
}


unpack_sym <- function(v, k) {
  M <- matrix(0, k, k)
  M[lower.tri(M, diag = TRUE)] <- v
  M[upper.tri(M)] <- t(M)[upper.tri(M)]
  M
}


# Draw from the leave-group-out predictive, which is the weighted mixture of
# Gaussians the scoring already has in hand: pick a component with probability
# wq, then a normal from it. Deterministic given `seed`.
sample_lgo_mixture <- function(cen, vv, wq, n_draw, seed) {
  set.seed(seed)
  j <- sample.int(length(wq), n_draw, replace = TRUE, prob = wq)
  stats::rnorm(n_draw, cen[j], sqrt(vv[j]))
}


# The two absolute moments every proper score here is built from, in closed form
# for a weighted Gaussian mixture:
#     e_xy = E|X - y|        e_xx = E|X - X'|
# The cross term is O(J^2), so the mixture is first thinned to at most
# `max_comp` components by systematic resampling -- deterministic, and unbiased
# in weight.
mixture_abs_moments <- function(y, m, s, w, max_comp = 200000L, seed = 1L) {
  keep <- w > 0
  m <- m[keep]; s <- s[keep]; w <- w[keep]
  if (length(m) == 1L) {
    # One component (any fully Gaussian fit): both moments are closed forms and
    # the O(J^2) cross term collapses. E|X - X'| for iid N(m, s^2) is 2s/sqrt(pi).
    z <- (y - m) / max(s, .Machine$double.eps)
    return(list(e_xy = (y - m) * (2 * stats::pnorm(z) - 1) +
                       2 * s * stats::dnorm(z),
                e_xx = 2 * s / sqrt(pi)))
  }
  w <- w / sum(w)
  J <- length(m)
  A <- function(d, sd) {
    sd <- pmax(sd, .Machine$double.eps)
    z <- d / sd
    d * (2 * stats::pnorm(z) - 1) + 2 * sd * stats::dnorm(z)
  }
  e_xy <- sum(w * A(y - m, s))
  # E|X - X'| is a double sum. Exact while that is affordable; beyond the cap,
  # sample pairs (i, j) ~ w x w instead of thinning the mixture.
  e_xx <- if (J <= sqrt(max_comp)) {          # J^2 overflows integer at J > 46341
    sum(outer(w, w) * A(outer(m, m, "-"), sqrt(outer(s^2, s^2, "+"))))
  } else {
    withr_seed(seed, {
      i <- sample.int(J, max_comp, replace = TRUE, prob = w)
      j <- sample.int(J, max_comp, replace = TRUE, prob = w)
      mean(A(m[i] - m[j], sqrt(s[i]^2 + s[j]^2)))
    })
  }
  list(e_xy = e_xy, e_xx = e_xx)
}

# MAE / MSE / neg.CRPS / neg.sCRPS from the two moments, using exactly the
# definitions in compute_score_given_pred() so the two functions are comparable.
# `neg.CRPS` is the ordinary CRPS (lower is better) despite the name; the sign
# convention is inherited from cross_validation().
scores_from_moments <- function(y, pred, e_xy, e_xx) {
  denom <- max(e_xx, .Machine$double.eps)
  c(MAE = abs(pred - y), MSE = (pred - y)^2,
    neg.CRPS = e_xy - 0.5 * denom,
    neg.sCRPS = e_xy / denom + 0.5 * log(denom))
}

# Split-Rhat (Gelman-Rubin) from a draws x chains matrix. Splitting each chain
# in half catches a chain that is still drifting within itself, which the
# plain statistic misses.
split_rhat <- function(mat) {
  n <- nrow(mat); m <- ncol(mat)
  if (n < 8L) return(NA_real_)
  h <- n %/% 2L
  X <- cbind(mat[seq_len(h), , drop = FALSE],
             mat[(n - h + 1L):n, , drop = FALSE])
  N <- nrow(X); M <- ncol(X)
  if (M < 2L) return(NA_real_)
  cm <- colMeans(X)
  W <- mean(apply(X, 2L, stats::var))
  if (!is.finite(W) || W <= 0) return(NA_real_)
  B <- N * stats::var(cm)
  sqrt(((N - 1) / N * W + B / N) / W)
}


# One data.frame for the whole replicate, built from a numeric matrix.
assemble_score_rows <- function(rows, global_ids, n_draw_out) {
  m <- matrix(unlist(rows, use.names = FALSE), nrow = length(rows), byrow = TRUE,
              dimnames = list(NULL, names(rows[[1]])))
  out <- as.data.frame(m)
  out$exact <- out$exact > 0
  if (n_draw_out > 0L)
    attr(out, "draws") <- stats::setNames(lapply(rows, attr, "draws"),
                                          as.character(global_ids))
  out
}

stop_transform_domain <- function() {
  # The predictive is Gaussian on the original scale, so it has mass wherever
  # the real line goes; a transform with a restricted domain (log, sqrt) will be
  # fed values outside it. cross_validation() behaves the same way -- say so
  # once, clearly, instead of emitting a NaN warning per draw.
  stop("`transform` produced non-finite values: the leave-group-out ",
       "predictive is Gaussian on the original scale, so it puts mass outside ",
       "this transform's domain. Use a transform defined on the whole real ",
       "line, or model the transformed response directly.", call. = FALSE)
}

# 32-node Gauss-Hermite, rescaled for N(0,1): nodes x, probability weights w.
gauss_hermite_32 <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    n <- 32L
    b <- sqrt(seq_len(n - 1L) / 2)
    J <- diag(0, n); J[cbind(2:n, 1:(n-1))] <- b; J[cbind(1:(n-1), 2:n)] <- b
    e <- eigen(J, symmetric = TRUE)
    o <- order(e$values)
    x <- e$values[o] * sqrt(2)
    w <- (e$vectors[1, o]^2)
    cache <<- list(x = x, w = w / sum(w))
    cache
  }
})

# Apply a transform that may depend on WHICH observation is being scored.
# A one-argument function is the usual elementwise scale; a two-argument one
# also receives the observation's index in the original data, which is what a
# per-pixel rescaling like exp(x) * level[i] needs. The old cross_validation()
# spelling of this was a closure holding a counter that had to be called in a
# particular order -- stateless is safer, and order-independent.
apply_tf <- function(transform, x, idx) {
  if (length(formals(transform)) >= 2L) transform(x, idx) else transform(x)
}

# Score one observation. With no transform every quantity is closed form; with
# one, the scale is applied to the observation and to every draw (matching
# compute_score_given_pred()) and the moments are estimated from samples.
group_cv_one_score <- function(y, ci, vi, wq, transform, crps_max,
                               n_score_draws, seed, pred_mean, obs_idx = NA_integer_) {
  if (identical(transform, identity)) {
    mm <- mixture_abs_moments(y, ci, sqrt(vi), wq, crps_max, seed)
    return(scores_from_moments(y, pred_mean, mm$e_xy, mm$e_xx))
  }
  if (!is.function(transform)) stop("`transform` must be a function.")
  keep <- wq > 0
  if (sum(keep) <= 4L) {
    # A small mixture makes the two moments ordinary integrals against normals,
    # so Gauss-Hermite is both cheaper than sampling and deterministic.
    gh <- gauss_hermite_32()
    m <- ci[keep]; sd <- sqrt(vi[keep]); w <- wq[keep] / sum(wq[keep])
    ty <- suppressWarnings(apply_tf(transform, y, obs_idx))
    gx <- lapply(seq_along(m), function(j)
      suppressWarnings(apply_tf(transform, m[j] + sd[j] * gh$x, obs_idx)))
    if (!is.finite(ty) || any(!vapply(gx, function(z) all(is.finite(z)), TRUE)))
      stop_transform_domain()
    e_xy <- sum(vapply(seq_along(m), function(j)
      w[j] * sum(gh$w * abs(gx[[j]] - ty)), 0))
    e_xx <- 0
    for (j in seq_along(m)) for (l in seq_along(m))
      e_xx <- e_xx + w[j] * w[l] *
        sum(outer(gh$w, gh$w) * abs(outer(gx[[j]], gx[[l]], "-")))
    pm <- sum(vapply(seq_along(m), function(j) w[j] * sum(gh$w * gx[[j]]), 0))
    return(scores_from_moments(ty, pm, e_xy, e_xx))
  }
  ty <- suppressWarnings(apply_tf(transform, y, obs_idx))
  x1 <- suppressWarnings(apply_tf(transform,
          sample_lgo_mixture(ci, vi, wq, n_score_draws, seed), obs_idx))
  x2 <- suppressWarnings(apply_tf(transform,
          sample_lgo_mixture(ci, vi, wq, n_score_draws, seed + 1L), obs_idx))
  if (!is.finite(ty) || anyNA(x1) || anyNA(x2)) stop_transform_domain()
  scores_from_moments(ty, mean(x1), mean(abs(x1 - ty)), mean(abs(x1 - x2)))
}


score_group_cv_replicate <- function(rep_model, draws, local_groups, global_ids,
                                     score_pos, score_global,
                                     use_r_eff, n_q = 8L, crps_max = 200000L,
                                     seed = 1L, n_draw_out = 0L,
                                     transform = identity,
                                     n_score_draws = 4000L) {
  noise <- rep_model$noise
  family <- noise$noise_type
  if (length(family) != 1L) {
    stop("group_cv() does not support bivariate measurement noise yet; use ",
         "cross_validation().")
  }
  mu_all <- as.numeric(noise$B_mu %*% noise$theta_mu)
  sigma_all <- as.numeric(exp(noise$B_sigma %*% noise$theta_sigma))
  nu_all <- as.numeric(noise$nu_lower_bound + exp(noise$B_nu %*% noise$theta_nu))
  fe_all <- as.numeric(rep_model$X %*% rep_model$feff)
  y_all <- rep_model$Y
  gaussian_noise <- identical(family, "normal")
  if (gaussian_noise) n_q <- 1L
  # A fully Gaussian model has no mixing variables at all: QQ is constant, the
  # importance weights are uniform, and the leave-group-out predictive is exact
  # in closed form. Nothing here is an approximation, so nothing is gated.
  gaussian_exact <- gaussian_noise && isTRUE(rep_model$all_gaussian)
  identity_tf <- identical(transform, identity)

  # For normal measurement noise the "prior draw" of the mixing variable is the
  # constant 1, so build it once rather than per group.
  Vq_const <- if (gaussian_noise) matrix(1, nrow = 1L, ncol = 1L) else NULL

  rows <- lapply(seq_along(local_groups), function(gi) {
    idx <- local_groups[[gi]]
    k <- length(idx)
    sp <- score_pos[[gi]]          # which members of the held-out set are scored
    sidx <- score_global[[gi]]     # and their indices in the original data
    M <- draws$mean[[gi]]          # k x S
    Cp <- draws$cov[[gi]]          # k(k+1)/2 x S
    S <- ncol(M)
    y <- y_all[idx]; fe <- fe_all[idx]
    mu <- mu_all[idx]; sg <- sigma_all[idx]; nu <- nu_all[idx]

    nq <- if (gaussian_noise) 1L else n_q
    # k x (S * nq): fresh draws, blocked by Gibbs draw.
    Vq <- if (gaussian_noise) matrix(1, nrow = k, ncol = S * nq)
          else prior_noise_V(family, k, mu, sg, nu, S * nq,
                             isTRUE(noise$single_V), seed + 977L * gi)

    if (k == 1L) {
      # The leave-one-out case, and the one that has to be fast. Row s of the
      # S x nq grid gets draw s's own nV values, so a length-S vector recycles
      # down the columns exactly as needed.
      Vm <- matrix(Vq[1, ], nrow = S, ncol = nq, byrow = TRUE)
      cen <- as.numeric(fe + M[1, ]) + mu * (Vm - 1)
      vv <- as.numeric(Cp[1, ]) + sg^2 * Vm
      # dnorm copies dim from its longest argument, and with S == nq == 1 that
      # is the scalar `y`, so the result comes back without dim and apply()
      # below fails. Restore the shape explicitly.
      ld <- matrix(stats::dnorm(y, cen, sqrt(vv), log = TRUE), nrow = S, ncol = nq)
    } else {
      # The Cholesky below is needed ONLY for ld, the k-variate density behind
      # the weights. The per-observation centre and variance the scores use are
      # scalar arithmetic -- Sig[i, i] is Vp[i, i] + sg_i^2 V_i, no factor
      # required -- so they are built vectorised, for the SCORED rows only,
      # once the loop is done. Filling S x nq x k arrays inside the loop and
      # then reading one slice of them was the bulk of its cost.
      ld <- matrix(NA_real_, S, nq)
      cen <- vv <- NULL
      for (s in seq_len(S)) {
        Vp <- unpack_sym(Cp[, s], k)
        for (q in seq_len(nq)) {
          Vsq <- Vq[, (s - 1L) * nq + q]        # this draw's own nV
          ctr <- fe + M[, s] + mu * (Vsq - 1)
          Sig <- Vp + diag(sg^2 * Vsq, k)
          ch <- chol(Sig)
          z <- backsolve(ch, y - ctr, transpose = TRUE)
          ld[s, q] <- -0.5 * (k * log(2 * pi) + 2 * sum(log(diag(ch))) + sum(z^2))
        }
      }
      # Index of (i, i) in a column-major packed lower triangle. The outer
      # parentheses matter: %/% binds tighter than * in R.
      pos_ii <- function(i) ((i - 1L) * (2L * k - i + 2L)) %/% 2L + 1L
      # Row s of the S x nq grid holds draw s's own nV values, so a length-S
      # vector recycles down the columns exactly as needed -- the same layout
      # the k == 1 branch uses, and the one `wq` below assumes.
      cen <- lapply(seq_len(k), function(i) NULL)
      vv <- cen
      for (i in sp) {
        Vm <- matrix(Vq[i, ], nrow = S, ncol = nq, byrow = TRUE)
        cen[[i]] <- as.numeric(fe[i] + M[i, ] + mu[i] * (Vm - 1))
        vv[[i]] <- as.numeric(Cp[pos_ii(i), ] + sg[i]^2 * Vm)
      }
    }

    # log p(y_I | y_-I, V_s), the held-out mixing variable integrated out.
    # apply() over a 1-row matrix is pure overhead, and the exact Gaussian path
    # always has S == 1.
    lp <- if (S == 1L) log_sum_exp(as.numeric(ld)) - log(nq)
          else apply(ld, 1L, log_sum_exp) - log(nq)
    if (gaussian_exact) {
      # QQ does not move across draws, so every lp is identical, the weights are
      # exactly uniform and the answer is exact from a single draw. Short-circuit
      # rather than running PSIS on a constant vector.
      reff <- 1
      w <- rep(1 / S, S)
      ps <- list(w = w, k = NA_real_, ess = S)
    } else {
      reff <- if (use_r_eff) min(1, mcmc_ess(lp) / S) else 1
      ps <- psis_weights(-lp, r_eff = reff)
      w <- ps$w
    }

    log_score <- log_sum_exp(log(pmax(w, 0)) + lp)

    wq <- as.numeric(outer(w, rep(1 / nq, nq)))   # weight of each (s, q) cell
    pred_mean <- numeric(k); pred_sd <- numeric(k)
    sc <- matrix(NA_real_, k, 4L)
    # One row per SCORED observation, not per held-out one: with a buffer the
    # unscored rows would stay NA and poison anything the caller does with them.
    draw_out <- if (n_draw_out > 0L) matrix(NA_real_, length(sp), n_draw_out) else NULL
    # Exact where the family allows it; NULL falls back to the sampled mixture.
    vmom <- noise_V_moments(family, nu)
    Vp_ii <- if (!is.null(vmom)) {
      if (k == 1L) matrix(Cp[1, ], nrow = 1L)
      else vapply(seq_len(S), function(s) diag(unpack_sym(Cp[, s], k)), numeric(k))
    } else NULL
    for (i in sp) {
      ci <- if (k == 1L) as.numeric(cen) else cen[[i]]
      vi <- if (k == 1L) as.numeric(vv) else vv[[i]]
      if (!is.null(vmom)) {
        # One value per Gibbs draw, weighted by w -- no n_q average anywhere.
        ct <- fe[i] + M[i, ] + mu[i] * (vmom$EV[i] - 1)
        vt <- Vp_ii[i, ] + mu[i]^2 * vmom$VarV[i] + sg[i]^2 * vmom$EV[i]
        pred_mean[i] <- sum(w * ct)
        pred_sd[i] <- sqrt(max(0, sum(w * (vt + ct^2)) - pred_mean[i]^2))
      } else {
        pred_mean[i] <- sum(wq * ci)
        pred_sd[i] <- sqrt(max(0, sum(wq * (vi + ci^2)) - pred_mean[i]^2))
      }
      sc[i, ] <- group_cv_one_score(y[i], ci, vi, wq, transform, crps_max,
                                    n_score_draws, seed + 613L * gi + i,
                                    pred_mean[i], sidx[match(i, sp)])
      if (n_draw_out > 0L) {
        draw_out[match(i, sp), ] <- sample_lgo_mixture(ci, vi, wq, n_draw_out,
                                                       seed + 613L * gi + i)
      }
    }

    # A plain numeric row: building a one-row data.frame per group and
    # rbind-ing thousands of them was the single largest cost in the exact
    # Gaussian path (69% of runtime at 3000 groups), dwarfing the arithmetic.
    row <- c(group = global_ids[gi], size = length(sp), n_held_out = k,
             mean = mean(pred_mean[sp]), sd = mean(pred_sd[sp]),
             log_score = if (identity_tf) log_score else NA_real_,
             MAE = mean(sc[sp, 1]), MSE = mean(sc[sp, 2]),
             neg.CRPS = mean(sc[sp, 3]), neg.sCRPS = mean(sc[sp, 4]),
             khat = ps$k,
             # psis_weights() already returns the effective sample size in
             # draws (r_eff / sum(w^2)); bounded by S, never rescaled.
             ess = ps$ess, r_eff = reff, rhat = NA_real_,
             exact = as.numeric(gaussian_exact))
    if (n_draw_out > 0L) attr(row, "draws") <- draw_out
    row
  })

  assemble_score_rows(rows, global_ids, n_draw_out)
}


#' @export
# Subsetting drops the attributes the print method reads (method, summary,
# n_draws), and the result then printed as a fully Gaussian exact fit whatever
# it actually was. A subset is a plain data frame.
#' @export
`[.ngme_group_cv` <- function(x, ...) {
  y <- NextMethod()
  if (is.data.frame(y)) class(y) <- "data.frame"
  y
}


# cross_validation() hands back a list with `mean.scores` and `sd.scores`, and
# callers reach for them by name. Computing them on demand keeps that spelling
# working without turning them into columns of the score table, which would
# show up in names() and in every print.
#
# Both are a one-row data frame, columns named as cross_validation() names
# them: bare for a single scale, prefixed by the transform's name when several
# were scored. A data frame rather than a matrix because callers reach in with
# `ms[["abs.MAE"]]`, which a matrix does not answer.
score_summary <- function(x, fun, label = "model_1") {
  cols <- c("MAE", "MSE", "neg.CRPS", "neg.sCRPS")
  one <- function(g, prefix) {
    v <- vapply(cols, function(cc) fun(as.data.frame(g)[[cc]]), 0)
    stats::setNames(v, if (nzchar(prefix)) paste0(prefix, ".", cols) else cols)
  }
  v <- if (inherits(x, "ngme_group_cv")) one(x, "")
       else unlist(lapply(names(x), function(nm) one(unclass(x)[[nm]], nm)))
  out <- as.data.frame(as.list(v), check.names = FALSE)
  rownames(out) <- label
  out
}

#' @export
`$.ngme_group_cv` <- function(x, name) {
  if (identical(name, "mean.scores")) return(score_summary(x, mean))
  if (identical(name, "sd.scores")) return(score_summary(x, stats::sd))
  NextMethod()
}

#' @export
`$.ngme_group_cv_list` <- function(x, name) {
  if (identical(name, "mean.scores")) return(score_summary(x, mean))
  if (identical(name, "sd.scores")) return(score_summary(x, stats::sd))
  NextMethod()
}


#' @export
print.ngme_group_cv <- function(x, ...) {
  m <- attr(x, "method")
  cat("Leave-group-out cross-validation\n")
  cat("  method: ", switch(m %||% "lowrank",
        lowrank = if (isTRUE(all(x$exact)))
          "exact closed form (fully Gaussian model)"
        else "low-rank updates on one full-data chain",
        exact = "exact per-fold chains (C++, parallel)"), "\n", sep = "")
  cat("  groups:", nrow(x), " draws:", attr(x, "n_draws"), "\n\n")
  print(round(attr(x, "summary"), 5))
  if (isTRUE(attr(x, "transform")))
    cat("\n  scores are on the transformed scale; log_score is not defined there\n")

  n_low <- attr(x, "n_low_ess")
  n_khat <- attr(x, "n_high_khat")
  bad <- attr(x, "n_unreliable")
  if (isTRUE(bad > 0)) {
    cat("\n", bad, " of ", nrow(x), " group(s) unreliable",
        " (ESS < ", attr(x, "min_ess"), ": ", n_low,
        "; Pareto k > ", attr(x, "khat_threshold"), ": ", n_khat, ")\n",
        "  Re-score those with cross_validation(type = \"custom\", test_idx = ...).\n",
        "  Their scores are reported but are not trustworthy.\n",
        sep = "")
    cat("  Raising n_gibbs_samples is NOT a reliable remedy: where the weights are\n",
        "  degenerate the effective sample size does not grow with the draw count,\n",
        "  and the estimate can converge tightly onto a wrong value -- the spread\n",
        "  shrinks while the error grows. Use the exact refit.\n", sep = "")
    if (n_low > 0 && n_khat == 0) {
      cat("  Note: every flag came from the effective sample size, not Pareto k.\n",
          "  That is the expected pattern -- k diagnoses the weight tail, while what\n",
          "  fails here is chain mixing. Do not raise min_ess to silence these.\n",
          sep = "")
    }
  }
  invisible(x)
}


# Scoring for the exact path: group_cv_exact_cpp() returns draws of
# eta_I = (AZ W)_I from a genuine leave-group-out chain, so there are no
# importance weights at all -- every draw carries weight 1/S. The measurement
# noise is layered on exactly as cross_validation() does, by drawing its mixing
# variable from the prior.
score_group_cv_exact <- function(rep_model, draws, local_groups, global_ids,
                                 score_pos, score_global, n_q = 8L, crps_max = 200000L, seed = 1L,
                                 n_draw_out = 0L, transform = identity,
                                 n_score_draws = 4000L, n_chains = 1L) {
  noise <- rep_model$noise
  family <- noise$noise_type
  if (length(family) != 1L) {
    stop("group_cv() does not support bivariate measurement noise yet; use ",
         "cross_validation().")
  }
  mu_all <- as.numeric(noise$B_mu %*% noise$theta_mu)
  sigma_all <- as.numeric(exp(noise$B_sigma %*% noise$theta_sigma))
  nu_all <- as.numeric(noise$nu_lower_bound + exp(noise$B_nu %*% noise$theta_nu))
  fe_all <- as.numeric(rep_model$X %*% rep_model$feff)
  y_all <- rep_model$Y
  if (identical(family, "normal")) n_q <- 1L

  rows <- lapply(seq_along(local_groups), function(gi) {
    idx <- local_groups[[gi]]
    k <- length(idx)
    sp <- score_pos[[gi]]; sidx <- score_global[[gi]]
    E <- draws$eta[[gi]]                      # k x S draws of eta_I
    S <- ncol(E)
    y <- y_all[idx]; fe <- fe_all[idx]
    mu <- mu_all[idx]; sg <- sigma_all[idx]; nu <- nu_all[idx]
    nq <- if (identical(family, "normal")) 1L else n_q
    Vq <- if (identical(family, "normal")) matrix(1, nrow = k, ncol = S * nq)
          else prior_noise_V(family, k, mu, sg, nu, S * nq,
                             isTRUE(noise$single_V), seed + 977L * gi)

    vmom <- noise_V_moments(family, nu)
    pred_mean <- numeric(k); pred_sd <- numeric(k)
    sc <- matrix(NA_real_, k, 4L)
    ld <- matrix(0, S, nq)
    draw_out <- if (n_draw_out > 0L) matrix(NA_real_, length(sp), n_draw_out) else NULL
    for (i in sp) {
      Vm <- matrix(Vq[i, ], nrow = S, ncol = nq, byrow = TRUE)
      cen <- (fe[i] + E[i, ]) + mu[i] * (Vm - 1)
      vv <- sg[i]^2 * Vm
      ld <- ld + matrix(stats::dnorm(y[i], cen, sqrt(vv), log = TRUE),
                        nrow = S, ncol = nq)
      wq <- rep(1 / (S * nq), S * nq)
      ci <- as.numeric(cen); vi <- as.numeric(vv)
      if (!is.null(vmom)) {
        # Draws of eta carry no weights here, so the mixture is a plain mean.
        ct <- fe[i] + E[i, ] + mu[i] * (vmom$EV[i] - 1)
        vt <- mu[i]^2 * vmom$VarV[i] + sg[i]^2 * vmom$EV[i]
        pred_mean[i] <- mean(ct)
        pred_sd[i] <- sqrt(max(0, mean(vt + ct^2) - pred_mean[i]^2))
      } else {
        pred_mean[i] <- sum(wq * ci)
        pred_sd[i] <- sqrt(max(0, sum(wq * (vi + ci^2)) - pred_mean[i]^2))
      }
      sc[i, ] <- group_cv_one_score(y[i], ci, vi, wq, transform, crps_max,
                                    n_score_draws, seed + 613L * gi + i,
                                    pred_mean[i], sidx[match(i, sp)])
      if (n_draw_out > 0L)
        draw_out[match(i, sp), ] <- sample_lgo_mixture(ci, vi, wq, n_draw_out,
                                                       seed + 613L * gi + i)
    }
    # Honest effective sample size: the autocorrelation of the scored
    # observations' own trace, which is the quantity the scores are built from.
    ess_exact <- if (length(sp) == 1L) mcmc_ess(E[sp, ])
                 else mean(vapply(sp, function(i) mcmc_ess(E[i, ]), 0))
    # Convergence across the independent chains, on the quantity being scored.
    # C++ concatenates the chains, so S = n_chains * (draws per chain) and the
    # matrix below puts one chain per column.
    n_ch <- as.integer(n_chains); S_ch <- S %/% max(1L, n_ch)
    rh <- if (n_ch > 1L && S_ch * n_ch == S)
      mean(vapply(sp, function(i)
        split_rhat(matrix(E[i, ], nrow = S_ch, ncol = n_ch)), 0), na.rm = TRUE)
      else NA_real_

    # Average over the scored rows only: `sc` has a row per held-out
    # observation and the buffer rows are never filled, so mean(sc[, j]) over
    # all of them is NA whenever an exclusion buffer is in use.
    row <- c(group = global_ids[gi], size = length(sp), n_held_out = k,
             mean = mean(pred_mean[sp]), sd = mean(pred_sd[sp]),
             log_score = if (identical(transform, identity))
               log_sum_exp(as.numeric(ld)) - log(S * nq) else NA_real_,
             MAE = mean(sc[sp, 1]), MSE = mean(sc[sp, 2]),
             neg.CRPS = mean(sc[sp, 3]), neg.sCRPS = mean(sc[sp, 4]),
             # There are no importance weights here, so there is no Pareto
             # tail to diagnose -- but the draws are still an autocorrelated
             # Markov chain, so the effective sample size is not S. Reporting
             # ess = S claimed perfect mixing and was simply wrong.
             khat = NA_real_, ess = ess_exact, r_eff = ess_exact / S,
             rhat = rh, exact = 0)
    if (n_draw_out > 0L) attr(row, "draws") <- draw_out
    row
  })
  assemble_score_rows(rows, global_ids, n_draw_out)
}
