# Helpers for "did this refactor change any number?" style regression tests.
#
# ngme_numeric_digest() walks an arbitrary R object (a fitted ngme object, a
# prediction, a list of samples) and returns every numeric leaf it contains as
# a named list. Comparing two digests therefore compares parameter estimates,
# optimizer trajectories, sampled latent fields, posterior samples and the
# fitted noise parameters all at once, without having to enumerate them.

ngme_numeric_digest_env <- function(x, prefix = "", depth = 0,
                                    acc = new.env(parent = emptyenv())) {
  if (depth > 12L) return(invisible(acc))
  if (is.function(x) || is.environment(x)) return(invisible(acc))
  if (isS4(x)) {
    # sparse matrices and friends: keep their numeric slots (e.g. @x, @i, @p)
    for (s in methods::slotNames(x)) {
      v <- methods::slot(x, s)
      if (is.numeric(v) || is.logical(v)) {
        ngme_numeric_digest_env(v, paste0(prefix, "@", s), depth + 1L, acc)
      }
    }
    return(invisible(acc))
  }
  if (is.numeric(x) || is.logical(x)) {
    v <- as.numeric(x)
    if (length(v)) assign(prefix, v, envir = acc)
    return(invisible(acc))
  }
  if (is.list(x)) {
    nm <- names(x)
    if (is.null(nm)) nm <- rep("", length(x))
    for (i in seq_along(x)) {
      key <- if (nzchar(nm[i])) nm[i] else as.character(i)
      ngme_numeric_digest_env(x[[i]], paste0(prefix, "$", key), depth + 1L, acc)
    }
    # results also hang off attributes (e.g. `chain_outputs` on an ngme fit)
    a <- attributes(x)
    skip <- c("names", "class", "row.names", "dim", "dimnames", "levels")
    for (an in setdiff(names(a), skip)) {
      ngme_numeric_digest_env(a[[an]], paste0(prefix, "@@", an), depth + 1L, acc)
    }
  }
  invisible(acc)
}

ngme_numeric_digest <- function(x) {
  e <- ngme_numeric_digest_env(x, "")
  keys <- sort(ls(e))
  out <- lapply(keys, function(k) get(k, envir = e))
  names(out) <- keys
  out
}

# Largest absolute difference between two digests; Inf if they do not even
# have the same shape.
ngme_digest_max_diff <- function(a, b) {
  if (!identical(names(a), names(b))) return(Inf)
  worst <- 0
  for (k in names(a)) {
    x <- a[[k]]; y <- b[[k]]
    if (length(x) != length(y)) return(Inf)
    d <- suppressWarnings(max(abs(x - y)))
    if (!is.finite(d)) d <- if (identical(x, y)) 0 else Inf
    worst <- max(worst, d)
  }
  worst
}

# Run `expr` with the factorization caches disabled, i.e. reassembling and
# refactorizing QQ (and redoing the symbolic phase) on every single use, the
# way the code behaved before the caches were introduced.
with_factor_cache_disabled <- function(expr) {
  withr::with_envvar(c(NGME2_DISABLE_FACTOR_CACHE = "1"), expr)
}

# Read the factorization counters, resetting them for the next measurement.
ngme_read_factor_counters <- function(reset = TRUE) {
  unlist(ngme2:::ngme_factor_counters(reset = reset))
}

# Counting the factorizations only works if every draw goes through one
# process; parallel chains make the totals thread-interleaved but still exact,
# so no skipping is needed. Deterministic comparisons, on the other hand, need
# a reproducible backend: Apple's Accelerate sparse solver is not bitwise
# reproducible even between two identical fits, so the equivalence tests pin
# the backend to the bundled CHOLMOD.
ngme_reproducible_control <- function(iterations = 20, burnin = 8, ...) {
  ngme2::control_opt(
    seed = 4321,
    burnin = burnin,
    iterations = iterations,
    n_batch = 1,
    n_parallel_chain = 1,
    max_num_threads = 1,
    solver_backend = "cholmod",
    verbose = FALSE,
    print_check_info = FALSE,
    R_hat_conv_check = FALSE,
    trend_std_conv_check = FALSE,
    # These fits pin exact iteration and factorization counts, so the
    # post-convergence polish phase (which adds iterations of its own) is off.
    polish_iterations = 0L,
    ...
  )
}
