#' @title Calibrate Inverse-Exponential Prior for NIG Driven Noise
#' @description
#' Calibrate \eqn{\lambda} in the prior
#' \eqn{\kappa = 1/\nu \sim \mathrm{Exp}(\lambda)}
#' using a tail-inflation target for driven NIG noise.
#'
#' The calibration target is
#' \deqn{\Pr(R_c(\nu) > r_\text{target}) = \alpha,}
#' where
#' \deqn{R_c(\nu) = \frac{\Pr(|U| > c \mid \nu)}{\Pr(|Z| > c)},\quad Z\sim N(0,1)}
#' and
#' \deqn{U = \frac{\mu(V-h) + \sigma\sqrt{V}Z}{\sigma\sqrt{h}},\quad V\sim\mathrm{GIG}(-1/2,\nu,\nu h^2).}
#'
#' The function solves \eqn{R_c(\nu_r) = r_\text{target}} using log-scale
#' bisection, then returns
#' \deqn{\lambda = -\nu_r \log(\alpha).}
#'
#' @param r_target target tail inflation level \eqn{r_\text{target} > 1}
#' @param alpha target prior probability in \eqn{(0,1)}
#' @param c tail threshold for \eqn{|U|>c}, typically `2.5` or `3`
#' @param mu NIG drift parameter in the driven noise term
#' @param sigma NIG scale parameter (> 0)
#' @param h positive increment scaling (> 0)
#' @param n_samples Monte Carlo sample size used per evaluation of \eqn{R_c(\nu)}
#' @param nu_lower initial lower bracket for \eqn{\nu}
#' @param nu_upper initial upper bracket for \eqn{\nu}
#' @param tol relative tolerance for solving \eqn{R_c(\nu_r)=r_\text{target}}
#' @param max_iter maximum number of bisection iterations
#' @param max_expand maximum bracket expansion steps on each side
#' @param seed optional integer seed for reproducible calibration
#' @param ... arguments forwarded to \code{calibrate_inv_exp_lambda_driven_nig()}
#'
#' @return A list with calibrated `lambda`, solved `nu_r`, achieved `rc_nu_r`,
#'   and diagnostics.
#' @export
calibrate_inv_exp_lambda_driven_nig <- function(
    r_target = 2,
    alpha = 0.1,
    c = 3,
    mu = 0,
    sigma = 1,
    h = 1,
    n_samples = 1e5,
    nu_lower = 0.1,
    nu_upper = 100,
    tol = 0.05,
    max_iter = 30,
    max_expand = 30,
    seed = NULL) {
  stopifnot(
    is.numeric(r_target), length(r_target) == 1, is.finite(r_target), r_target > 1,
    is.numeric(alpha), length(alpha) == 1, is.finite(alpha), alpha > 0, alpha < 1,
    is.numeric(c), length(c) == 1, is.finite(c), c > 0,
    is.numeric(mu), length(mu) == 1, is.finite(mu),
    is.numeric(sigma), length(sigma) == 1, is.finite(sigma), sigma > 0,
    is.numeric(h), length(h) == 1, is.finite(h), h > 0,
    is.numeric(n_samples), length(n_samples) == 1, is.finite(n_samples), n_samples > 0,
    is.numeric(nu_lower), length(nu_lower) == 1, is.finite(nu_lower), nu_lower > 0,
    is.numeric(nu_upper), length(nu_upper) == 1, is.finite(nu_upper), nu_upper > nu_lower,
    is.numeric(tol), length(tol) == 1, is.finite(tol), tol > 0,
    is.numeric(max_iter), length(max_iter) == 1, is.finite(max_iter), max_iter >= 1,
    is.numeric(max_expand), length(max_expand) == 1, is.finite(max_expand), max_expand >= 1
  )

  n_samples <- as.integer(n_samples)
  max_iter <- as.integer(max_iter)
  max_expand <- as.integer(max_expand)

  gaussian_tail <- 2 * stats::pnorm(-c)
  if (gaussian_tail <= 0) {
    stop("Gaussian tail probability underflows for this c; use a smaller c.")
  }

  normalize_seed <- function(x, offset = 0L) {
    as.integer((abs(as.double(x)) + as.double(offset)) %%
      (.Machine$integer.max - 1L) + 1L)
  }

  rc_cache <- new.env(parent = emptyenv())

  eval_rc <- function(nu_value) {
    key <- formatC(nu_value, digits = 17, format = "fg", flag = "#")
    if (exists(key, envir = rc_cache, inherits = FALSE)) {
      return(get(key, envir = rc_cache, inherits = FALSE))
    }

    if (!is.null(seed)) {
      seed_v <- normalize_seed(seed, 101L)
      seed_z <- normalize_seed(seed, 202L)
    } else {
      seed_v <- 0L
      seed_z <- NULL
    }

    V <- rgig(n = n_samples, p = -0.5, a = nu_value, b = nu_value * h * h, seed = seed_v)
    Z <- if (!is.null(seed_z)) {
      withr::with_seed(seed_z, stats::rnorm(n_samples))
    } else {
      stats::rnorm(n_samples)
    }

    U <- (mu * (V - h) + sigma * sqrt(V) * Z) / (sigma * sqrt(h))
    rc <- mean(abs(U) > c) / gaussian_tail
    if (!is.finite(rc)) rc <- NA_real_

    assign(key, rc, envir = rc_cache)
    rc
  }

  nu_lo <- nu_lower
  nu_hi <- nu_upper
  bracket_found <- FALSE
  rc_lo <- NA_real_
  rc_hi <- NA_real_
  rc_min <- NA_real_
  rc_max <- NA_real_

  for (k in seq_len(max_expand + 1L)) {
    nu_grid <- exp(seq(log(nu_lo), log(nu_hi), length.out = 41))
    rc_grid <- vapply(nu_grid, eval_rc, numeric(1))
    ok <- is.finite(rc_grid)
    if (!any(ok)) {
      stop("Numerical failure while evaluating R_c(nu); no finite values were obtained.")
    }

    nu_grid <- nu_grid[ok]
    rc_grid <- rc_grid[ok]
    f_grid <- rc_grid - r_target
    rc_min <- min(rc_grid)
    rc_max <- max(rc_grid)

    crossings <- which(f_grid[-length(f_grid)] * f_grid[-1] <= 0)
    if (length(crossings) > 0) {
      # Use the rightmost crossing (closest to Gaussian regime).
      i <- crossings[length(crossings)]
      nu_lo <- nu_grid[i]
      nu_hi <- nu_grid[i + 1]
      rc_lo <- rc_grid[i]
      rc_hi <- rc_grid[i + 1]
      bracket_found <- TRUE
      break
    }

    if (k > max_expand) break

    if (rc_max < r_target) {
      nu_lo <- nu_lo / 2
    } else if (rc_min > r_target) {
      nu_hi <- nu_hi * 2
    } else {
      # No clean crossing due non-monotonic/noisy curve; expand both sides.
      nu_lo <- nu_lo / 2
      nu_hi <- nu_hi * 2
    }
  }

  if (!bracket_found) {
    stop(
      "Could not bracket nu_r for R_c(nu)=r_target. ",
      "Observed scanned range was R_c in [", signif(rc_min, 4), ", ", signif(rc_max, 4), "]. ",
      "Target may be unattainable for the chosen (c, mu, sigma, h). ",
      "Try larger c (e.g. 2.5 or 3), smaller r_target, or larger n_samples."
    )
  }

  converged <- FALSE
  f_lo <- rc_lo - r_target
  f_hi <- rc_hi - r_target
  nu_mid <- exp((log(nu_lo) + log(nu_hi)) / 2)
  rc_mid <- eval_rc(nu_mid)

  for (iter in seq_len(max_iter)) {
    nu_mid <- exp((log(nu_lo) + log(nu_hi)) / 2)
    rc_mid <- eval_rc(nu_mid)
    f_mid <- rc_mid - r_target

    rel_err <- abs(rc_mid - r_target) / max(1, r_target)
    if (rel_err <= tol) {
      converged <- TRUE
      break
    }

    if (f_lo * f_mid <= 0) {
      nu_hi <- nu_mid
      rc_hi <- rc_mid
      f_hi <- f_mid
    } else {
      nu_lo <- nu_mid
      rc_lo <- rc_mid
      f_lo <- f_mid
    }
  }

  nu_r <- nu_mid
  lambda <- -nu_r * log(alpha)

  structure(
    list(
      lambda = lambda,
      nu_r = nu_r,
      rc_nu_r = rc_mid,
      r_target = r_target,
      alpha = alpha,
      alpha_implied = exp(-lambda / nu_r),
      c = c,
      mu = mu,
      sigma = sigma,
      h = h,
      n_samples = n_samples,
      bracket_nu = c(lower = nu_lo, upper = nu_hi),
      bracket_rc = c(lower = rc_lo, upper = rc_hi),
      rc_scan_range = c(min = rc_min, max = rc_max),
      converged = converged,
      max_iter = max_iter,
      tol = tol,
      seed = seed,
      root_signs = c(lower = f_lo, upper = f_hi)
    ),
    class = "ngme_inv_exp_calibration"
  )
}

#' @rdname calibrate_inv_exp_lambda_driven_nig
#' @export
calibrate_inv_exp_lambda <- function(...) {
  calibrate_inv_exp_lambda_driven_nig(...)
}
