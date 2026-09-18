#!/usr/bin/env Rscript
# Phase breakdown for the benchmark models.
#
#   Rscript benchmarks/profile.R [--models a,b] [--backend accelerate] [--top 8]
#
# The phase counters are compiled out unless the package was built with
# -DNGME_PHASE_TIMING, so this reports nothing useful otherwise and says so.
# To build with them:
#
#   add -DNGME_PHASE_TIMING to PKG_CPPFLAGS in src/Makevars, remove src/*.o,
#   R CMD INSTALL .            (and revert afterwards -- it is not free)
#
# Shares are of wall clock for that model. They are NOT a partition: several
# phases nest inside others, and `opt_step` encloses almost everything. The
# leaves are what to read.

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1]
}
only    <- getarg("--models")
backend <- getarg("--backend", "accelerate")
top     <- as.integer(getarg("--top", "8"))
iters   <- as.integer(getarg("--iterations", "100"))

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
if (is.na(here) || !nzchar(here)) here <- "benchmarks"
source(file.path(here, "models.R"))

if (!isTRUE(ngme2:::factorization_timing()[["enabled"]] == 1)) {
  cat("phase timing is compiled out; rebuild with -DNGME_PHASE_TIMING (see header)\n")
  quit(status = 2)
}

models <- BM_MODELS
if (!is.null(only)) models <- models[strsplit(only, ",")[[1]]]

# Leaves worth reading, in the order they occur in a pass.
LEAVES <- c(
  op_build = "build K, Z", op_dK = "dK / dZ", op_trace = "tr(K^-1 dK)",
  k_numeric = "factorize K", qq_prod = "K' D K", qq_add = "QQ = Q + meas",
  qq_measure = "rebuild meas block", qq_numeric = "factorize QQ",
  rmvn = "draw W + cond mean", sw_M = "rhs M", sw_G = "G", sw_H = "H",
  rb_qu_solve = "probe block solve", rb_product = "trace products",
  grad_V = "sample V", grad_score = "score s", grad_assemble = "gradient assembly",
  grad_prec_lat = "precond latent", grad_prec_ZGN = "precond Z/GN",
  grad_prec_merr = "precond meas")

for (nm in names(models)) {
  models[[nm]]$setup()
  it <- if (is.null(models[[nm]]$iters)) iters else models[[nm]]$iters
  ngme2:::factorization_timing(reset = TRUE)
  t0 <- Sys.time()
  invisible(suppressWarnings(models[[nm]]$fit(
    seed = 4321, iterations = it, burnin = min(20L, it %/% 2L),
    n_parallel_chain = 1,
    max_num_threads = 1, solver_backend = backend, print_check_info = FALSE,
    verbose = FALSE, polish_iterations = 0L)))
  W <- as.numeric(Sys.time() - t0, units = "secs")
  tm <- ngme2:::factorization_timing(reset = TRUE)
  v <- vapply(names(LEAVES), function(k) as.numeric(tm[[k]]), 0)
  named <- sum(v)
  cat(sprintf("\n%s  (%.2f s, %d iterations)\n", nm, W, it))
  ord <- order(v, decreasing = TRUE)
  for (i in head(ord, top)) if (v[i] > 0.002 * W)
    cat(sprintf("   %-22s %6.2f%%  %7.3f s\n", LEAVES[i], 100 * v[i] / W, v[i]))
  cat(sprintf("   %-22s %6.2f%%  %7.3f s\n", "(everything else)",
              100 * (W - named) / W, W - named))
}
