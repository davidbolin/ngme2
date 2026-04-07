#!/usr/bin/env Rscript

# Demo: verify analytical Hessian of log f(W|V) w.r.t. theta_K
# Model: K(a,b) = exp(a) M_a + 3 b^2 M_b
#
# Notation (matching C++ implementation and derivation):
# - A = diag(1 / (sigma^2 * V))
# - m = mu * (V - h)
# - e = W - K^{-1} m
# - KWm = K W - m
#
# Hessian formula for entries H_ij with i,j in {a,b}:
#   H_ij = - tr(K^{-1} K_j K^{-1} K_i) + tr(K^{-1} K_ij)
#          - [ (K_j W)^T A K_i e + (K W - m)^T A K_ij e
#              + (K W - m)^T A K_i K^{-1} K_j K^{-1} m ]
#
suppressPackageStartupMessages({
  # No extra packages strictly needed; keep it base R compatible
})

make_spd <- function(n, jitter = 1e-1) {
  M <- matrix(rnorm(n * n), n, n)
  S <- crossprod(M)
  S <- 0.5 * (S + t(S)) + diag(jitter, n)
  S
}

ell_logf <- function(par, env) {
  a <- par[1]; b <- par[2]
  with(env, {
    K <- exp(a) * Ma + 3 * b^2 * Mb
    # logdet(K)
    # Use base::determinant for robust log |K|
    ld <- as.numeric(determinant(K, logarithm = TRUE)$modulus)
    KWm <- K %*% W - m
    qf <- sum((KWm * A_diag) * KWm) # (KWm)^T A (KWm), A is diagonal
    ld - 0.5 * qf
  })
}

analytic_hessian <- function(par, env) {
  a <- par[1]; b <- par[2]
  with(env, {
    K   <- exp(a) * Ma + 3 * b^2 * Mb
    Ka  <- exp(a) * Ma
    Kb  <- 6 * b * Mb
    Kaa <- exp(a) * Ma
    Kab <- 0 * Ma
    Kbb <- 6 * Mb

    # solvers
    Kinvm <- solve(K, m)
    e     <- W - Kinvm
    KWm   <- K %*% W - m

    # Helper trace(K^{-1} X)
    tr_Kinv_X <- function(X) {
      # tr(K^{-1} X) = sum_i e_i^T K^{-1} X e_i ~ exact via solve of each column (small n)
      # For small n this is fine; for large n use randomized trace.
      sum(diag(solve(K, X)))
    }

    # - tr(K^{-1} K_j K^{-1} K_i)
    tr_Kinv_Kj_Kinv_Ki <- function(Ki, Kj) {
      # compute tr(K^{-1} Kj K^{-1} Ki) exactly (small n)
      X <- solve(K, Kj)
      Y <- solve(K, Ki)
      sum(diag(X %*% Y))
    }

    # Common pieces
    A_dot <- A_diag

    # Terms builder
    Hij_fun <- function(Ki, Kj, Kij) {
      term_tr  <- - tr_Kinv_Kj_Kinv_Ki(Ki, Kj) + tr_Kinv_X(Kij)
      term1    <- sum((Kj %*% W) * A_dot * (Ki %*% e))
      term2    <- sum(KWm * A_dot * (Kij %*% e))
      yj1      <- Kj %*% Kinvm
      yj2      <- solve(K, yj1)
      yi3      <- Ki %*% yj2
      term3    <- sum(KWm * A_dot * yi3)
      term_tr - (term1 + term2 + term3)
    }

    H <- matrix(0, 2, 2)
    # (a,a)
    H[1,1] <- Hij_fun(Ka, Ka, Kaa)
    # (a,b) and (b,a)
    H[1,2] <- Hij_fun(Ka, Kb, Kab)
    H[2,1] <- H[1,2]
    # (b,b)
    H[2,2] <- Hij_fun(Kb, Kb, Kbb)
    H
  })
}

numeric_hessian <- function(par, fn, env, eps = 1e-5) {
  d <- length(par)
  H <- matrix(0, d, d)
  f00 <- fn(par, env)
  # central difference second derivatives
  for (i in 1:d) {
    for (j in i:d) {
      e_i <- rep(0, d); e_i[i] <- eps
      e_j <- rep(0, d); e_j[j] <- eps
      fpp <- fn(par + e_i + e_j, env)
      fpm <- fn(par + e_i - e_j, env)
      fmp <- fn(par - e_i + e_j, env)
      fmm <- fn(par - e_i - e_j, env)
      H[i, j] <- (fpp - fpm - fmp + fmm) / (4 * eps * eps)
      if (j != i) H[j, i] <- H[i, j]
    }
  }
  # symmetrize for numerical noise
  0.5 * (H + t(H))
}

run_once <- function(n = 10, seed = 1L, a = 0.1, b = -0.2) {
  set.seed(seed)
  Ma <- make_spd(n)
  Mb <- make_spd(n)
  sigma <- exp(rnorm(n, sd = 0.3))
  V <- rexp(n, rate = 1) + 0.5
  A_diag <- 1 / (sigma^2 * V)
  mu <- rnorm(n)
  h <- rep(1, n)
  m <- mu * (V - h)
  W <- rnorm(n)

  env <- list(Ma = Ma, Mb = Mb, A_diag = A_diag, mu = mu, V = V, h = h, m = m, W = W)
  par <- c(a, b)

  H_ana <- analytic_hessian(par, env)
  H_num <- numeric_hessian(par, ell_logf, env, eps = 1e-5)
  list(analytic = H_ana, numeric = H_num, diff = H_ana - H_num,
       max_abs_err = max(abs(H_ana - H_num)))
}

main <- function() {
  cat("Hessian check for K(a,b) = exp(a) M_a + 3 b^2 M_b\n")
  res <- run_once(n = 12, seed = 42, a = 0.3, b = -0.4)
  cat("Analytical Hessian:\n"); print(signif(res$analytic, 5))
  cat("Numeric Hessian:\n"); print(signif(res$numeric, 5))
  cat("Difference (ana - num):\n"); print(signif(res$diff, 5))
  cat(sprintf("Max abs error: %.3e\n", res$max_abs_err))

  # Run a few random trials
  errs <- numeric(5)
  for (k in 1:5) {
    a <- rnorm(1)
    b <- rnorm(1)
    tmp <- run_once(n = 10, seed = 100 + k, a = a, b = b)
    errs[k] <- tmp$max_abs_err
  }
  cat("Random trials max abs errors:\n")
  print(signif(errs, 5))
  cat(sprintf("Median max abs error across trials: %.3e\n", median(errs)))
}

if (identical(environment(), globalenv())) {
  main()
}

