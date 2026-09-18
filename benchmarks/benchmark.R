#!/usr/bin/env Rscript
# Timing and correctness harness for the model suite in models.R.
#
#   Rscript benchmarks/benchmark.R                      run everything, print a table
#   Rscript benchmarks/benchmark.R --save ref.rds       record a reference
#   Rscript benchmarks/benchmark.R --compare ref.rds    time, and check estimates
#                                                       against that reference
#   Rscript benchmarks/benchmark.R --models ar1,tp_nig  restrict the set
#   Rscript benchmarks/benchmark.R --reps 3             median of n runs each
#   Rscript benchmarks/benchmark.R --solver-order auto  pick the ordering by fill
#
# The point of --save/--compare is that a performance change must not move the
# numbers. Record a reference before the change, compare after; anything beyond
# the tolerance is a behaviour change, not an optimization.
#
# Timings are single-threaded and single-chain on purpose: parallel chains make
# wall clock depend on machine load, and the convergence machinery needs more
# than one chain, which would make run length vary between configurations and
# stop the timings comparing like with like. Iterations are fixed for the same
# reason.

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1]
}
save_to    <- getarg("--save")
compare_to <- getarg("--compare")
only       <- getarg("--models")
no_large   <- "--no-large" %in% args
reps       <- as.integer(getarg("--reps", "1"))
iters      <- as.integer(getarg("--iterations", "100"))
tol        <- as.numeric(getarg("--tol", "1e-8"))
# Apple's Accelerate sparse solver is NOT bitwise reproducible between two
# identical fits, so an estimate comparison against a reference is meaningless
# under it -- two runs of the same code differ by ~1e-6 on a spacetime model.
# Correctness comparisons therefore pin CHOLMOD. Pass --backend accelerate for
# timings representative of the macOS default, but do not trust --compare then.
backend    <- getarg("--backend", "cholmod")
# Fill-reducing ordering for the sparse factorizations; see control_opt().
solver_ord <- getarg("--solver-order", "default")
# Dropped when the installed build predates the option, so a reference can be
# recorded against an older library with the same harness.
solver_arg <- function() {
  if (!"solver_order" %in% names(formals(ngme2::control_opt))) list()
  else list(solver_order = solver_ord)
}

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "benchmarks"
source(file.path(here, "models.R"))

models <- BM_MODELS
if (!is.null(only)) models <- models[strsplit(only, ",")[[1]]]
if (no_large) models <- models[!vapply(models, function(m) isTRUE(m$large), TRUE)]

# Every parameter the fit estimated, flattened and named, so two runs can be
# compared entry by entry.
estimates <- function(fit) {
  r <- fit$replicates[[1]]
  out <- c()
  for (m in r$models) out <- c(out, m$theta_K, m$noise$theta_sigma,
                               m$noise$theta_mu, m$noise$theta_nu)
  c(out, r$feff, r$noise$theta_sigma, r$noise$theta_mu, r$noise$theta_nu)
}

run_one <- function(nm) {
  truth <- models[[nm]]$setup()
  fitfun <- models[[nm]]$fit
  it <- if (is.null(models[[nm]]$iters)) iters else models[[nm]]$iters
  el <- numeric(reps); est <- NULL; dim <- NA_integer_
  for (k in seq_len(reps)) {
    t0 <- Sys.time()
    fit <- suppressWarnings(do.call(fitfun, c(list(
      seed = 4321, iterations = it, burnin = min(20L, it %/% 2L),
      # control_opt requires iterations to be a whole number of checkpoints.
      n_parallel_chain = 1,
      solver_backend = backend,
      max_num_threads = 1, print_check_info = FALSE, verbose = FALSE,
      polish_iterations = 0L), solver_arg())))
    el[k] <- as.numeric(Sys.time() - t0, units = "secs")
    if (is.null(est)) {
      est <- estimates(fit)
      dim <- sum(vapply(fit$replicates[[1]]$models, function(m) nrow(m$operator$K), 0L))
    }
  }
  list(time = median(el), est = est, truth = truth, iters = it, dim = dim)
}

# One untimed fit before anything is measured. The first fits in a fresh R
# process pay lazy loading, allocator growth and BLAS thread setup, which on the
# fastest models is larger than the fit itself -- enough that a median over a
# couple of reps can land on a warm-up run and report several times the true
# time. Warming once per process is enough; it is not per model.
warmup <- function() {
  nm <- names(models)[which.min(vapply(models, function(m)
    if (isTRUE(m$large)) 1e9 else 1, 0))]
  m <- models[[nm]]
  m$setup()
  invisible(try(suppressWarnings(do.call(m$fit, c(list(
    seed = 4321, iterations = 10L, burnin = 5L,
    n_parallel_chain = 1, solver_backend = backend,
    max_num_threads = 1, print_check_info = FALSE, verbose = FALSE,
    polish_iterations = 0L), solver_arg()))), silent = TRUE))
}
invisible(capture.output(suppressMessages(warmup())))

res <- list()
cat(sprintf("backend %s, %d rep(s)\n", backend, reps))
cat(sprintf("%-15s %8s %6s %7s %10s\n", "model", "seconds", "iters", "dim", "ms/iter"))
for (nm in names(models)) {
  r <- try(run_one(nm), silent = TRUE)
  if (inherits(r, "try-error")) {
    cat(sprintf("%-15s %8s  %s\n", nm, "ERROR", conditionMessage(attr(r, "condition"))))
    next
  }
  res[[nm]] <- r
  cat(sprintf("%-15s %8.2f %6d %7d %10.1f\n", nm, r$time, r$iters, r$dim,
              1000 * r$time / r$iters))
}
cat(sprintf("%-15s %8.2f\n", "TOTAL", sum(vapply(res, function(x) x$time, 0))))

if (!is.null(save_to)) {
  saveRDS(res, save_to); cat("\nreference written to ", save_to, "\n", sep = "")
}

if (!is.null(compare_to)) {
  ref <- readRDS(compare_to)
  cat("\nestimates against ", compare_to, " (tolerance ", tol, ")\n", sep = "")
  bad <- 0
  for (nm in names(res)) {
    if (is.null(ref[[nm]])) { cat(sprintf("  %-15s not in reference\n", nm)); next }
    a <- res[[nm]]$est; b <- ref[[nm]]$est
    if (length(a) != length(b)) {
      cat(sprintf("  %-15s LENGTH %d vs %d\n", nm, length(a), length(b))); bad <- bad + 1; next
    }
    d <- max(abs(a - b) / pmax(abs(b), 1e-8))
    ok <- d <= tol
    if (!ok) bad <- bad + 1
    cat(sprintf("  %-15s %-4s max rel diff %.3e  %7.2f s -> %7.2f s  (%.2fx)\n",
                nm, if (ok) "OK" else "DIFF", d, ref[[nm]]$time, res[[nm]]$time,
                ref[[nm]]$time / res[[nm]]$time))
  }
  cat(if (bad == 0) "\nall estimates unchanged\n" else
      sprintf("\n%d model(s) CHANGED beyond tolerance\n", bad))
  quit(status = if (bad == 0) 0 else 1)
}
