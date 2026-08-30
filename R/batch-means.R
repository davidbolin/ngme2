#' Batch-Means Covariance Estimator for SGD Trajectories
#'
#' Estimate the asymptotic covariance from a single SGD trajectory using the
#' increasing-batch construction from Xi et al. (2020).
#'
#' @param trajectory numeric matrix with rows = iterations and columns = parameters.
#' @param alpha stepsize decay exponent in
#'   \eqn{\eta_i = \eta i^{-\alpha}}. Must satisfy \eqn{1/2 < \alpha < 1}.
#' @param M number of retained batches (excluding burn-in batch 0). If `NULL`,
#'   use \eqn{\lfloor n^{(1-\alpha)/2} \rfloor}.
#' @param N decorrelation constant in
#'   \eqn{e_k = \lfloor ((k+1)N)^{1/(1-\alpha)} \rfloor}. If `NULL`, use
#'   \eqn{N = n^{1-\alpha}/(M+1)}.
#' @param drop_burnin logical; if `TRUE`, discard batch 0.
#' @param burnin_iter non-negative integer. Explicitly discard the first
#'   `burnin_iter` iterations before building Xi-style batches. After trimming,
#'   batch boundaries are rebuilt from iteration 1 of the retained trajectory.
#'
#' @return A list containing the covariance estimate, pooled mean, batch sizes,
#' batch boundaries, and chain-level metadata.
#' @export
batch_means_estimator <- function(
    trajectory,
    alpha = 0.501,
    M = NULL,
    N = NULL,
    drop_burnin = TRUE,
    burnin_iter = 0) {
  stopifnot(
    is.matrix(trajectory) || is.data.frame(trajectory),
    is.numeric(alpha),
    length(alpha) == 1,
    alpha > 0.5,
    alpha < 1,
    is.logical(drop_burnin),
    length(drop_burnin) == 1,
    is.numeric(burnin_iter),
    length(burnin_iter) == 1,
    is.finite(burnin_iter),
    burnin_iter >= 0,
    abs(burnin_iter - round(burnin_iter)) < sqrt(.Machine$double.eps)
  )

  trajectory <- as.matrix(trajectory)
  storage.mode(trajectory) <- "double"

  burnin_iter <- as.integer(round(burnin_iter))
  n_iter_total <- nrow(trajectory)
  n_param <- ncol(trajectory)
  stopifnot(
    n_iter_total >= 2,
    n_param >= 1
  )

  if (burnin_iter > 0) {
    stopifnot(
      "burnin_iter is too large; at least 2 iterations must remain after trimming." =
        n_iter_total - burnin_iter >= 2
    )
    trajectory <- trajectory[seq.int(burnin_iter + 1L, n_iter_total), , drop = FALSE]
  }
  n_iter <- nrow(trajectory)

  boundary <- .batch_means_boundaries(
    n_iter = n_iter,
    alpha = alpha,
    M = M,
    N = N
  )

  keep_idx <- if (drop_burnin) {
    seq.int(2L, length(boundary$end))
  } else {
    seq_along(boundary$end)
  }

  stopifnot(
    "No valid batch remains after burn-in removal; increase iterations or reduce M." =
      length(keep_idx) > 0
  )

  if (length(keep_idx) < 2) {
    warning(
      "Only one batch is retained; covariance estimate may be unstable.",
      call. = FALSE
    )
  }

  batch_means <- t(vapply(keep_idx, function(k) {
    idx <- boundary$start[k]:boundary$end[k]
    colMeans(trajectory[idx, , drop = FALSE])
  }, numeric(n_param)))

  batch_sizes <- boundary$size[keep_idx]
  global_mean <- colSums(batch_means * batch_sizes) / sum(batch_sizes)

  centered <- sweep(batch_means, 2, global_mean, "-")
  covariance <- matrix(0, n_param, n_param)
  for (i in seq_len(nrow(centered))) {
    covariance <- covariance + batch_sizes[i] * tcrossprod(centered[i, ])
  }
  covariance <- covariance / length(keep_idx)

  param_names <- colnames(trajectory)
  if (is.null(param_names)) {
    param_names <- paste0("param_", seq_len(n_param))
  }
  names(global_mean) <- param_names
  dimnames(covariance) <- list(param_names, param_names)

  ret <- list(
    mean = global_mean,
    covariance = covariance,
    n_eff = sum(batch_sizes),
    batch_sizes = batch_sizes,
    batch_start = boundary$start[keep_idx],
    batch_end = boundary$end[keep_idx],
    batch_means = batch_means,
    alpha = alpha,
    M = length(keep_idx),
    N = boundary$N,
    dropped_burnin = drop_burnin,
    burnin_iter = burnin_iter,
    n_iterations = n_iter,
    n_iterations_total = n_iter_total
  )

  class(ret) <- "ngme_batch_means"
  ret
}

#' Pooled Batch-Means Confidence Intervals from Multiple Chains
#'
#' Compute pooled point estimates, covariance, and Wald confidence intervals from
#' multiple chain trajectories.
#'
#' @param chain_trajectories list of numeric matrices, each with rows = iterations
#'   and columns = parameters.
#' @param level confidence level.
#' @param alpha stepsize decay exponent in \eqn{(1/2, 1)}.
#' @param M number of retained batches per chain (excluding burn-in batch 0).
#' @param N decorrelation constant used in the batch boundary formula.
#' @param drop_burnin logical; if `TRUE`, discard batch 0.
#' @param burnin_iter non-negative integer. Explicitly discard the first
#'   `burnin_iter` iterations of each chain before Xi-style batching.
#'
#' @return A list with pooled estimates, covariance, standard errors, and
#' confidence intervals.
#' @export
batch_means_ci <- function(
    chain_trajectories,
    level = 0.95,
    alpha = 0.501,
    M = NULL,
    N = NULL,
    drop_burnin = TRUE,
    burnin_iter = 0) {
  if (is.matrix(chain_trajectories) || is.data.frame(chain_trajectories)) {
    chain_trajectories <- list(as.matrix(chain_trajectories))
  }

  stopifnot(
    is.list(chain_trajectories),
    length(chain_trajectories) >= 1,
    is.numeric(level),
    length(level) == 1,
    level > 0,
    level < 1
  )

  chain_trajectories <- lapply(chain_trajectories, function(x) {
    stopifnot(is.matrix(x) || is.data.frame(x))
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    x
  })

  n_param <- unique(vapply(chain_trajectories, ncol, integer(1)))
  stopifnot(
    "All chains must have the same number of parameters." = length(n_param) == 1
  )

  n_iter_vec <- vapply(chain_trajectories, nrow, integer(1))
  if (length(unique(n_iter_vec)) > 1) {
    min_iter <- min(n_iter_vec)
    warning(
      "Chains have different lengths; truncating all chains to the minimum length.",
      call. = FALSE
    )
    chain_trajectories <- lapply(chain_trajectories, function(x) {
      x[seq_len(min_iter), , drop = FALSE]
    })
  }

  param_names <- colnames(chain_trajectories[[1]])
  if (is.null(param_names)) {
    param_names <- paste0("param_", seq_len(n_param))
  }

  chain_results <- lapply(chain_trajectories, function(x) {
    if (is.null(colnames(x))) {
      colnames(x) <- param_names
    }
    batch_means_estimator(
      trajectory = x,
      alpha = alpha,
      M = M,
      N = N,
      drop_burnin = drop_burnin,
      burnin_iter = burnin_iter
    )
  })

  n_chain <- length(chain_results)
  chain_means <- do.call(rbind, lapply(chain_results, `[[`, "mean"))
  pooled_mean <- colMeans(chain_means)

  pooled_cov <- Reduce(`+`, lapply(chain_results, `[[`, "covariance")) / n_chain
  dimnames(pooled_cov) <- list(param_names, param_names)

  n_eff_vec <- vapply(chain_results, `[[`, numeric(1), "n_eff")
  n_eff <- min(n_eff_vec)
  if (length(unique(n_eff_vec)) > 1) {
    warning(
      "Effective sample sizes differ across chains; using the minimum for CI width.",
      call. = FALSE
    )
  }

  z_value <- stats::qnorm(1 - (1 - level) / 2)
  std_error <- sqrt(pmax(diag(pooled_cov), 0) / (n_eff * n_chain))
  conf_int <- cbind(
    lower = pooled_mean - z_value * std_error,
    upper = pooled_mean + z_value * std_error
  )

  pooled_mean <- stats::setNames(pooled_mean, param_names)
  std_error <- stats::setNames(std_error, param_names)
  rownames(conf_int) <- param_names

  ret <- list(
    estimates = pooled_mean,
    std_error = std_error,
    conf_int = conf_int,
    covariance = pooled_cov,
    level = level,
    z_value = z_value,
    alpha = alpha,
    M = chain_results[[1]]$M,
    N = chain_results[[1]]$N,
    burnin_iter = burnin_iter,
    n_eff = n_eff,
    n_chains = n_chain,
    n_iterations = chain_results[[1]]$n_iterations,
    n_iterations_total = chain_results[[1]]$n_iterations_total,
    chain_results = chain_results
  )

  class(ret) <- "ngme_batch_ci"
  ret
}

#' Batch-Means Confidence Intervals from an ngme Fit
#'
#' Compute Xi-style batch-means confidence intervals directly from trajectories
#' stored in an `ngme` object.
#'
#' @param ngme fitted `ngme` object with `store_traj = TRUE` during optimization.
#' @param name either `"all"` (default) to use all parameters jointly, a latent
#'   model name, or `"general"` for measurement-noise/fixed-effect block.
#' @param level confidence level.
#' @param alpha stepsize decay exponent in \eqn{(1/2, 1)}.
#' @param M number of retained batches (excluding burn-in batch 0).
#' @param N decorrelation constant in the batch boundary formula.
#' @param apply_transform logical; if `TRUE`, inference is first computed in the
#'   raw optimization space and then mapped to user scale via post-hoc Delta
#'   method (Post-delta).
#' @param drop_burnin logical; if `TRUE`, discard batch 0.
#' @param burnin_iter non-negative integer. Explicitly discard the first
#'   `burnin_iter` optimization iterations before Xi-style batching.
#'
#' @return Same structure as [batch_means_ci()] with extra metadata fields
#' `name` and `apply_transform`.
#' @export
ngme_batch_ci <- function(
    ngme,
    name = "all",
    level = 0.95,
    alpha = 0.501,
    M = NULL,
    N = NULL,
    apply_transform = TRUE,
    drop_burnin = TRUE,
    burnin_iter = 0) {
  stopifnot(
    inherits(ngme, "ngme"),
    is.character(name),
    length(name) == 1,
    is.logical(apply_transform),
    length(apply_transform) == 1
  )

  chain_data <- .ngme_chain_trajectories(
    ngme = ngme,
    name = name,
    apply_transform = FALSE
  )

  raw_ret <- batch_means_ci(
    chain_trajectories = chain_data$chains,
    level = level,
    alpha = alpha,
    M = M,
    N = N,
    drop_burnin = drop_burnin,
    burnin_iter = burnin_iter
  )

  ret <- if (apply_transform) {
    .delta_transform_batch_ci(raw_ret, chain_data$transforms)
  } else {
    raw_ret
  }

  ret$name <- name
  ret$apply_transform <- apply_transform
  if (apply_transform) {
    ret$raw_estimates <- raw_ret$estimates
    ret$raw_covariance <- raw_ret$covariance
    ret$raw_std_error <- raw_ret$std_error
    ret$raw_conf_int <- raw_ret$conf_int
  }
  ret
}

#' Refit an Existing ngme Object with SGD and Compute Batch-Means CI
#'
#' Run one additional SGD stage (warm-started from an existing `ngme` fit),
#' then compute Xi-style batch-means confidence intervals from the newly stored
#' trajectories.
#'
#' @param fit existing fitted `ngme` object (can be obtained with any optimizer).
#' @param iterations optimization iterations for the SGD CI stage.
#' @param alpha Xi-style exponent used for both polynomial schedule and
#'   batch-means inference.
#' @param optimizer optimizer for the CI stage; must be `sgd(...)` for Xi theory.
#' @param burnin burn-in iterations before SGD optimization.
#' @param n_batch number of optimization checkpoints.
#' @param n_parallel_chain number of parallel chains for SGD.
#' @param t0 non-negative schedule offset in
#'   \eqn{\eta_t = \eta_0 (t + t0)^{-\alpha}}.
#' @param start_sd standard deviation for randomized chain initialization.
#' @param seed random seed for CI-stage optimization. Defaults to a seed drawn
#'   from the current R random number stream, so \code{set.seed()} makes the
#'   result reproducible.
#' @param verbose logical; print optimization progress.
#' @param level confidence level for CI.
#' @param name parameter block for CI (`"all"`, latent name, or `"general"`).
#' @param M number of retained Xi batches (excluding batch 0).
#' @param N decorrelation constant in Xi batch boundaries.
#' @param apply_transform logical; if `TRUE`, apply Post-delta transform after
#'   raw-space inference.
#' @param drop_burnin logical; if `TRUE`, discard Xi batch 0.
#' @param burnin_iter non-negative integer used both as optimizer schedule
#'   warmup (`stepsize_schedule_burnin_iter`) and as explicit CI trajectory
#'   trimming before Xi-style batching.
#' @param control_opt optional pre-built `control_opt` object for the CI stage.
#'   If supplied, it is used directly (with `store_traj` forced to `TRUE`).
#' @param ... additional arguments forwarded to [control_opt_batch_ci()] when
#'   `control_opt` is `NULL`.
#'
#' @return An object of class `ngme_batch_ci` with extra fields:
#'   `refit` (the SGD-refit model), `refit_control_opt`, and `refit_source`.
#' @export
compute_ngme_ci <- function(
    fit,
    iterations = 4000,
    alpha = 0.501,
    optimizer = sgd(stepsize = 0.03),
    burnin = 100,
    n_batch = 20,
    n_parallel_chain = 4,
    t0 = 1,
    start_sd = 0.2,
    seed = ngme_random_seed(),
    verbose = FALSE,
    level = 0.95,
    name = "all",
    M = NULL,
    N = NULL,
    apply_transform = TRUE,
    drop_burnin = TRUE,
    burnin_iter = 0,
    control_opt = NULL,
    ...) {
  stopifnot(inherits(fit, "ngme"))

  fit_meta <- attr(fit, "fit")
  stopifnot(
    "Missing fit metadata in `fit`; please pass a fitted ngme object from ngme()." =
      !is.null(fit_meta),
    "Missing formula in fit metadata." = !is.null(fit_meta$formula),
    "Missing data in fit metadata." = !is.null(fit_meta$data)
  )

  dots <- list(...)
  if ("schedule_burnin_iter" %in% names(dots)) {
    stop(
      "In compute_ngme_ci(), `schedule_burnin_iter` is unified with `burnin_iter`. ",
      "Please set only `burnin_iter`.",
      call. = FALSE
    )
  }

  if (is.null(control_opt)) {
    control_opt <- do.call(
      control_opt_batch_ci,
      utils::modifyList(
        list(
          optimizer = optimizer,
          burnin = burnin,
          iterations = iterations,
          n_batch = n_batch,
          n_parallel_chain = n_parallel_chain,
          alpha = alpha,
          t0 = t0,
          schedule_burnin_iter = burnin_iter,
          start_sd = start_sd,
          seed = seed,
          verbose = verbose
        ),
        dots
      )
    )
  } else {
    stopifnot(
      inherits(control_opt, "control_opt"),
      "When `control_opt` is supplied, extra `...` arguments are not used." =
        length(dots) == 0
    )
  }

  stopifnot(
    "compute_ngme_ci currently requires optimizer = sgd(...) for Xi-style guarantees." =
      identical(control_opt$sgd_method, "sgd")
  )

  if (isTRUE(control_opt$trend_std_conv_check) ||
    isTRUE(control_opt$R_hat_conv_check)) {
    warning(
      "Convergence diagnostics are enabled in `control_opt`; early stopping may reduce fixed-iteration Xi-style validity.",
      call. = FALSE
    )
  }

  control_opt$stepsize_schedule_burnin_iter <- as.integer(round(burnin_iter))
  control_opt$store_traj <- TRUE

  fit_env <- build_cv_rebuild_env(
    formula = fit_meta$formula,
    ngme_obj = fit,
    data = fit_meta$data
  )

  control_ngme_refit <- fit$replicates[[1]]$control_ngme
  if (is.null(control_ngme_refit)) {
    control_ngme_refit <- control_ngme()
  }

  refit <- do.call(
    "ngme",
    list(
      formula = fit_meta$formula,
      data = fit_meta$data,
      family = fit_meta$family,
      control_opt = control_opt,
      control_ngme = control_ngme_refit,
      group = fit_meta$group,
      replicate = fit_meta$replicate,
      start = fit
    ),
    envir = fit_env
  )

  ci <- ngme_batch_ci(
    ngme = refit,
    name = name,
    level = level,
    alpha = alpha,
    M = M,
    N = N,
    apply_transform = apply_transform,
    drop_burnin = drop_burnin,
    burnin_iter = burnin_iter
  )

  ci$refit <- refit
  ci$refit_control_opt <- control_opt
  ci$refit_source <- fit
  ci
}

#' @rdname compute_ngme_ci
#' @export
compute_ngme_CI <- compute_ngme_ci

#' Refit an Existing ngme Object with SGLD and Extract Samples
#'
#' Run one additional SGLD stage (warm-started from an existing `ngme` fit),
#' then extract posterior-like samples from stored trajectories.
#'
#' @param fit existing fitted `ngme` object (can be obtained with any optimizer).
#' @param iterations optimization iterations for the SGLD stage.
#' @param optimizer optimizer for the sampling stage; must be `sgld(...)`.
#' @param burnin burn-in iterations before optimization.
#' @param n_batch number of optimization checkpoints.
#' @param n_parallel_chain number of parallel chains.
#' @param alpha polynomial schedule exponent used by `poly_decay(alpha, t0)`.
#' @param t0 non-negative schedule offset.
#' @param start_sd standard deviation for randomized chain initialization.
#' @param seed random seed for SGLD stage. Defaults to a seed drawn from the
#'   current R random number stream, so \code{set.seed()} makes the result
#'   reproducible.
#' @param verbose logical; print optimization progress.
#' @param name parameter block to extract: `"all"` (default), latent model name,
#'   or `"general"`.
#' @param burnin_iter non-negative integer used both as optimizer schedule
#'   warmup (`stepsize_schedule_burnin_iter`) and as explicit trajectory
#'   trimming before sampling.
#' @param thinning positive integer thinning interval.
#' @param n_gibbs_samples optional positive integer. If supplied, override
#'   `control_ngme$n_gibbs_samples` for the SGLD refit stage only; otherwise
#'   inherit the value from `fit`.
#' @param apply_transform logical; apply parameter transforms to user scale.
#' @param combine_chains logical; if `TRUE`, return one combined data.frame,
#'   otherwise return one data.frame per chain.
#' @param control_opt optional pre-built `control_opt` object for the SGLD
#'   stage. If supplied, it is used directly (with `store_traj` forced to
#'   `TRUE`).
#' @param ... additional arguments forwarded to [control_opt()] when
#'   `control_opt` is `NULL`.
#'
#' @return A data.frame (or list of data.frames) from [ngme_sgld_samples()]
#' with extra attributes `refit`, `refit_control_opt`, and `refit_source`.
#' @export
compute_ngme_sgld_samples <- function(
    fit,
    iterations = 4000,
    optimizer = sgld(stepsize = 0.01, temperature = 1),
    burnin = 100,
    n_batch = 20,
    n_parallel_chain = 4,
    alpha = 0.6,
    t0 = 10,
    start_sd = 0.2,
    seed = ngme_random_seed(),
    verbose = FALSE,
    name = "all",
    burnin_iter = 0,
    thinning = 1,
    n_gibbs_samples = NULL,
    apply_transform = TRUE,
    combine_chains = TRUE,
    control_opt = NULL,
    ...) {
  stopifnot(inherits(fit, "ngme"))
  stopifnot(
    is.null(n_gibbs_samples) || (
      is.numeric(n_gibbs_samples) &&
        length(n_gibbs_samples) == 1 &&
        is.finite(n_gibbs_samples) &&
        n_gibbs_samples > 0 &&
        abs(n_gibbs_samples - round(n_gibbs_samples)) < sqrt(.Machine$double.eps)
    )
  )

  fit_meta <- attr(fit, "fit")
  stopifnot(
    "Missing fit metadata in `fit`; please pass a fitted ngme object from ngme()." =
      !is.null(fit_meta),
    "Missing formula in fit metadata." = !is.null(fit_meta$formula),
    "Missing data in fit metadata." = !is.null(fit_meta$data)
  )

  dots <- list(...)
  if ("schedule_burnin_iter" %in% names(dots)) {
    stop(
      "In compute_ngme_sgld_samples(), `schedule_burnin_iter` is unified with `burnin_iter`. ",
      "Please set only `burnin_iter`.",
      call. = FALSE
    )
  }

  if (is.null(control_opt)) {
    control_opt <- do.call(
      "control_opt",
      utils::modifyList(
        list(
          optimizer = optimizer,
          burnin = burnin,
          iterations = iterations,
          n_batch = n_batch,
          n_parallel_chain = n_parallel_chain,
          store_traj = TRUE,
          trend_std_conv_check = FALSE,
          R_hat_conv_check = FALSE,
          stepsize_control = poly_decay(
            alpha = alpha,
            t0 = t0,
            burnin_iter = burnin_iter
          ),
          start_sd = start_sd,
          seed = seed,
          verbose = verbose
        ),
        dots
      )
    )
  } else {
    stopifnot(
      inherits(control_opt, "control_opt"),
      "When `control_opt` is supplied, extra `...` arguments are not used." =
        length(dots) == 0
    )
  }

  stopifnot(
    "compute_ngme_sgld_samples currently requires optimizer = sgld(...)." =
      identical(control_opt$sgd_method, "sgld")
  )

  if (isTRUE(control_opt$trend_std_conv_check) ||
    isTRUE(control_opt$R_hat_conv_check)) {
    warning(
      "Convergence diagnostics are enabled in `control_opt`; early stopping may reduce fixed-iteration sampling consistency.",
      call. = FALSE
    )
  }

  control_opt$stepsize_schedule_burnin_iter <- as.integer(round(burnin_iter))
  control_opt$store_traj <- TRUE

  fit_env <- build_cv_rebuild_env(
    formula = fit_meta$formula,
    ngme_obj = fit,
    data = fit_meta$data
  )

  control_ngme_refit <- fit$replicates[[1]]$control_ngme
  if (is.null(control_ngme_refit)) {
    control_ngme_refit <- control_ngme()
  }
  if (!is.null(n_gibbs_samples)) {
    control_ngme_refit$n_gibbs_samples <- as.integer(round(n_gibbs_samples))
  }

  refit <- do.call(
    "ngme",
    list(
      formula = fit_meta$formula,
      data = fit_meta$data,
      family = fit_meta$family,
      control_opt = control_opt,
      control_ngme = control_ngme_refit,
      group = fit_meta$group,
      replicate = fit_meta$replicate,
      start = fit
    ),
    envir = fit_env
  )

  ret <- ngme_sgld_samples(
    ngme = refit,
    name = name,
    burnin_iter = burnin_iter,
    thinning = thinning,
    apply_transform = apply_transform,
    combine_chains = combine_chains
  )

  attr(ret, "refit") <- refit
  attr(ret, "refit_control_opt") <- control_opt
  attr(ret, "refit_source") <- fit
  ret
}

#' Extract Posterior-like Samples from Stored SGLD Trajectories
#'
#' Build posterior-like samples from optimizer trajectories by dropping an
#' initial burn-in segment and applying thinning.
#'
#' @param ngme fitted `ngme` object with `store_traj = TRUE`.
#' @param name parameter block to extract: `"all"` (default), latent model name,
#'   or `"general"`.
#' @param burnin_iter non-negative integer. Number of initial iterations to
#'   discard before sampling.
#' @param thinning positive integer thinning interval.
#' @param apply_transform logical; apply parameter transforms to user scale.
#' @param combine_chains logical; if `TRUE`, return one combined data.frame,
#'   otherwise return one data.frame per chain.
#'
#' @return A data.frame (or list of data.frames when `combine_chains = FALSE`)
#' with columns `.chain`, `.draw`, `.iter`, and one column per parameter.
#' @export
ngme_sgld_samples <- function(
    ngme,
    name = "all",
    burnin_iter = 0,
    thinning = 1,
    apply_transform = TRUE,
    combine_chains = TRUE) {
  stopifnot(
    inherits(ngme, "ngme"),
    is.character(name),
    length(name) == 1,
    is.numeric(burnin_iter),
    length(burnin_iter) == 1,
    is.finite(burnin_iter),
    burnin_iter >= 0,
    abs(burnin_iter - round(burnin_iter)) < sqrt(.Machine$double.eps),
    is.numeric(thinning),
    length(thinning) == 1,
    is.finite(thinning),
    thinning >= 1,
    abs(thinning - round(thinning)) < sqrt(.Machine$double.eps),
    is.logical(apply_transform),
    length(apply_transform) == 1,
    is.logical(combine_chains),
    length(combine_chains) == 1
  )

  chain_data <- .ngme_chain_trajectories(
    ngme = ngme,
    name = name,
    apply_transform = apply_transform
  )

  burnin_iter <- as.integer(round(burnin_iter))
  thinning <- as.integer(round(thinning))

  n_iter_vec <- vapply(chain_data$chains, nrow, integer(1))
  n_iter <- min(n_iter_vec)
  if (length(unique(n_iter_vec)) > 1) {
    warning(
      "Chains have different lengths; truncating to the minimum length.",
      call. = FALSE
    )
  }

  stopifnot(
    "burnin_iter must be smaller than the available trajectory length." =
      burnin_iter < n_iter
  )

  keep_iter <- seq.int(burnin_iter + 1L, n_iter, by = thinning)
  stopifnot(length(keep_iter) >= 1)

  per_chain <- lapply(seq_along(chain_data$chains), function(chain_idx) {
    mat <- chain_data$chains[[chain_idx]][seq_len(n_iter), , drop = FALSE]
    mat <- mat[keep_iter, , drop = FALSE]
    df <- as.data.frame(mat, check.names = FALSE)
    df$.chain <- chain_idx
    df$.draw <- seq_len(nrow(df))
    df$.iter <- keep_iter
    df <- df[, c(".chain", ".draw", ".iter", colnames(mat)), drop = FALSE]
    rownames(df) <- NULL
    df
  })

  if (isTRUE(combine_chains)) {
    ret <- do.call(rbind, per_chain)
    rownames(ret) <- NULL
  } else {
    ret <- per_chain
  }

  attr(ret, "name") <- name
  attr(ret, "burnin_iter") <- burnin_iter
  attr(ret, "thinning") <- thinning
  attr(ret, "apply_transform") <- apply_transform
  attr(ret, "n_chains") <- length(per_chain)
  attr(ret, "n_samples_per_chain") <- length(keep_iter)
  ret
}

#' Quantile Confidence Intervals from SGLD Samples
#'
#' Compute parameter-wise quantile confidence intervals from posterior-like
#' SGLD samples returned by [ngme_sgld_samples()].
#'
#' @param samples data.frame (or list of data.frames) returned by
#'   [ngme_sgld_samples()].
#' @param lower lower quantile probability (e.g. 0.025).
#' @param upper upper quantile probability (e.g. 0.975).
#'
#' @return A list with posterior-like point estimates (`estimates`), quantile
#' intervals (`ci`), and sample covariance matrix (`covariance`).
#' Class: `ngme_sgld_ci`.
#' @export
ngme_sgld_ci <- function(samples, lower = 0.025, upper = 0.975) {
  stopifnot(
    is.data.frame(samples) ||
      (is.list(samples) &&
        length(samples) >= 1 &&
        all(vapply(samples, is.data.frame, logical(1)))),
    is.numeric(lower),
    length(lower) == 1,
    is.finite(lower),
    lower > 0,
    lower < 1,
    is.numeric(upper),
    length(upper) == 1,
    is.finite(upper),
    upper > 0,
    upper < 1,
    lower < upper
  )

  if (is.list(samples) && !is.data.frame(samples)) {
    samples <- do.call(rbind, samples)
    rownames(samples) <- NULL
  }

  burnin_iter <- attr(samples, "burnin_iter")
  thinning <- attr(samples, "thinning")
  apply_transform <- attr(samples, "apply_transform")
  name <- attr(samples, "name")

  param_cols <- setdiff(colnames(samples), c(".chain", ".draw", ".iter"))
  stopifnot(length(param_cols) >= 1)
  stopifnot(
    "All parameter columns in `samples` must be numeric." =
      all(vapply(samples[, param_cols, drop = FALSE], is.numeric, logical(1)))
  )

  sample_mat <- as.matrix(samples[, param_cols, drop = FALSE])
  storage.mode(sample_mat) <- "double"
  colnames(sample_mat) <- param_cols

  estimates <- colMeans(sample_mat)
  std_error <- if (nrow(sample_mat) > 1) {
    apply(sample_mat, 2, stats::sd) / sqrt(nrow(sample_mat))
  } else {
    setNames(rep(0, ncol(sample_mat)), colnames(sample_mat))
  }

  qmat <- t(vapply(seq_len(ncol(sample_mat)), function(j) {
    stats::quantile(
      sample_mat[, j],
      probs = c(lower, upper),
      names = FALSE,
      na.rm = TRUE
    )
  }, numeric(2)))
  colnames(qmat) <- c("lower", "upper")
  rownames(qmat) <- colnames(sample_mat)

  cov_mat <- if (nrow(sample_mat) > 1) {
    stats::cov(sample_mat)
  } else {
    matrix(0, ncol(sample_mat), ncol(sample_mat))
  }
  dimnames(cov_mat) <- list(colnames(sample_mat), colnames(sample_mat))

  n_chains <- attr(samples, "n_chains")
  n_samples_per_chain <- attr(samples, "n_samples_per_chain")
  if (is.null(n_chains)) {
    n_chains <- if (".chain" %in% colnames(samples)) {
      length(unique(samples$.chain))
    } else {
      1L
    }
  }
  if (is.null(n_samples_per_chain)) {
    n_samples_per_chain <- if (".chain" %in% colnames(samples)) {
      as.integer(min(table(samples$.chain)))
    } else {
      nrow(samples)
    }
  }

  ret <- list(
    estimates = stats::setNames(estimates, colnames(sample_mat)),
    std_error = stats::setNames(as.numeric(std_error), colnames(sample_mat)),
    ci = qmat,
    covariance = cov_mat,
    lower = lower,
    upper = upper,
    burnin_iter = burnin_iter,
    thinning = thinning,
    n_chains = n_chains,
    n_samples_per_chain = n_samples_per_chain,
    n_samples = nrow(sample_mat),
    name = name,
    apply_transform = apply_transform,
    samples = samples
  )

  class(ret) <- "ngme_sgld_ci"
  ret
}

#' Summary for Batch-Means CI Results
#'
#' Summarize an object of class `ngme_batch_ci`, including parameter means,
#' standard errors, confidence intervals, and covariance matrix.
#'
#' @param object an object of class `ngme_batch_ci`.
#' @param digits number of digits used for printing.
#' @param ... currently unused.
#'
#' @return A list with summary table and covariance matrix.
#' @export
summary.ngme_batch_ci <- function(object, digits = 4, ...) {
  stopifnot(inherits(object, "ngme_batch_ci"))

  conf_int <- as.matrix(object$conf_int)
  est <- as.numeric(object$estimates)
  se <- as.numeric(object$std_error)
  pnames <- names(object$estimates)
  if (is.null(pnames)) {
    pnames <- rownames(conf_int)
  }
  if (is.null(pnames)) {
    pnames <- paste0("param_", seq_along(est))
  }

  table <- data.frame(
    parameter = pnames,
    mean = est,
    std_error = se,
    lower = conf_int[, "lower"],
    upper = conf_int[, "upper"],
    check.names = FALSE,
    row.names = NULL
  )

  ret <- list(
    table = table,
    covariance = object$covariance,
    level = object$level,
    alpha = object$alpha,
    M = object$M,
    N = object$N,
    burnin_iter = object$burnin_iter,
    n_eff = object$n_eff,
    n_chains = object$n_chains,
    n_iterations = object$n_iterations,
    n_iterations_total = object$n_iterations_total,
    name = object$name,
    apply_transform = object$apply_transform,
    digits = digits
  )

  if (!is.null(object$raw_estimates)) {
    raw_conf <- as.matrix(object$raw_conf_int)
    raw_names <- names(object$raw_estimates)
    if (is.null(raw_names)) {
      raw_names <- rownames(raw_conf)
    }
    if (is.null(raw_names)) {
      raw_names <- paste0("param_", seq_along(object$raw_estimates))
    }

    ret$raw_table <- data.frame(
      parameter = raw_names,
      mean = as.numeric(object$raw_estimates),
      std_error = as.numeric(object$raw_std_error),
      lower = raw_conf[, "lower"],
      upper = raw_conf[, "upper"],
      check.names = FALSE,
      row.names = NULL
    )
    ret$raw_covariance <- object$raw_covariance
  }

  class(ret) <- "summary.ngme_batch_ci"
  ret
}

#' @export
print.summary.ngme_batch_ci <- function(x, ...) {
  .format_table <- function(tbl, digits) {
    out <- tbl
    is_num <- vapply(out, is.numeric, logical(1))
    out[is_num] <- lapply(out[is_num], signif, digits = digits)
    out
  }

  cat("Batch-means CI summary\n")
  cat(
    "  level =", x$level,
    "| alpha =", x$alpha,
    "| M =", x$M,
    "| n_eff =", x$n_eff,
    "| n_chains =", x$n_chains, "\n"
  )
  if (!is.null(x$name)) {
    cat("  block =", x$name, "\n")
  }
  if (!is.null(x$apply_transform)) {
    cat("  apply_transform =", x$apply_transform, "\n")
  }

  cat("\nParameter summary:\n")
  print(.format_table(x$table, x$digits), row.names = FALSE)

  cat("\nCovariance matrix:\n")
  print(signif(x$covariance, x$digits))

  if (!is.null(x$raw_table)) {
    cat("\nRaw-space parameter summary (pre-transform):\n")
    print(.format_table(x$raw_table, x$digits), row.names = FALSE)
    cat("\nRaw-space covariance matrix:\n")
    print(signif(x$raw_covariance, x$digits))
  }

  invisible(x)
}

.batch_means_boundaries <- function(n_iter, alpha, M = NULL, N = NULL) {
  stopifnot(
    is.numeric(n_iter),
    length(n_iter) == 1,
    n_iter >= 2,
    is.numeric(alpha),
    length(alpha) == 1,
    alpha > 0.5,
    alpha < 1
  )

  if (is.null(M)) {
    M <- floor(n_iter^((1 - alpha) / 2))
  }
  M <- as.integer(M)

  stopifnot(
    is.numeric(M),
    length(M) == 1,
    M >= 1
  )

  if (is.null(N)) {
    N <- (n_iter^(1 - alpha)) / (M + 1)
  }

  stopifnot(
    is.numeric(N),
    length(N) == 1,
    N > 0
  )

  end <- floor(((seq_len(M + 1L)) * N)^(1 / (1 - alpha)))
  end <- pmin(pmax(as.integer(end), 1L), as.integer(n_iter))
  end <- cummax(end)
  end[length(end)] <- as.integer(n_iter)

  # Remove zero-size batches induced by integer rounding
  keep <- c(TRUE, diff(end) > 0)
  end <- end[keep]

  start <- c(1L, head(end, -1) + 1L)
  size <- end - start + 1L

  list(
    start = start,
    end = end,
    size = size,
    N = N
  )
}

.ngme_chain_trajectories <- function(ngme, name = "all", apply_transform = TRUE) {
  ngme_rep <- ngme$replicates[[1]]

  if (!identical(name, "all")) {
    tr <- get_trace_trajectories(
      ngme = ngme,
      name = name,
      apply_transform = apply_transform
    )

    tr_list <- unname(tr$trajectories)
    parameter_names <- make.unique(tr$parameter_names)

    chains <- lapply(seq_len(tr$n_chains), function(chain_idx) {
      mat <- do.call(cbind, lapply(tr_list, function(x) {
        x[, chain_idx, drop = FALSE]
      }))
      mat <- as.matrix(mat)
      colnames(mat) <- parameter_names
      mat
    })

    transforms <- tr$transformations
    stopifnot(length(transforms) == length(parameter_names))
    return(list(chains = chains, parameter_names = parameter_names, transforms = transforms))
  }

  latent_names <- names(ngme_rep$models)
  parts <- lapply(latent_names, function(lat_name) {
    latent <- ngme_rep$models[[lat_name]]
    traj <- attr(latent, "lat_traj")
    stopifnot(
      "Please run ngme() with store_traj=TRUE before using ngme_batch_ci()." =
        !is.null(traj)
    )

    ts <- get_latent_info(latent)
    ts$name <- paste0(ts$name, " (", lat_name, ")")

    list(traj = traj, ts = ts)
  })

  block_traj <- attr(ngme_rep, "block_traj")
  stopifnot(
    "Please run ngme() with store_traj=TRUE before using ngme_batch_ci()." =
      !is.null(block_traj)
  )

  ts_block <- get_noise_info(ngme_rep$noise)
  fe_names <- if (length(ngme_rep$feff) == 0) {
    NULL
  } else if (!is.null(names(ngme_rep$feff)) &&
    all(!is.na(names(ngme_rep$feff))) &&
    all(nzchar(names(ngme_rep$feff)))) {
    names(ngme_rep$feff)
  } else {
    paste("fixed effect", seq_len(length(ngme_rep$feff)))
  }
  fe_trans <- rep(list(identity), length(ngme_rep$feff))
  ts_block$name <- c(ts_block$name, fe_names)
  ts_block$trans <- c(ts_block$trans, fe_trans)

  parts <- c(parts, list(list(traj = block_traj, ts = ts_block)))

  n_chain <- length(parts[[1]]$traj)
  stopifnot(
    "All trajectory components must have the same number of chains." =
      all(vapply(parts, function(x) length(x$traj) == n_chain, logical(1)))
  )

  parameter_names <- make.unique(unlist(lapply(parts, function(x) x$ts$name), use.names = FALSE))
  parameter_transforms <- unlist(lapply(parts, function(x) x$ts$trans), recursive = FALSE, use.names = FALSE)
  stopifnot(length(parameter_transforms) == length(parameter_names))

  chains <- lapply(seq_len(n_chain), function(chain_idx) {
    mats <- lapply(parts, function(part) {
      mat <- as.matrix(part$traj[[chain_idx]])
      if (apply_transform) {
        mat <- .apply_parameter_transforms(mat, part$ts$trans)
      }
      mat
    })

    min_iter <- min(vapply(mats, ncol, integer(1)))
    mats <- lapply(mats, function(x) x[, seq_len(min_iter), drop = FALSE])

    chain_mat <- t(do.call(rbind, mats))
    colnames(chain_mat) <- parameter_names
    chain_mat
  })

  list(chains = chains, parameter_names = parameter_names, transforms = parameter_transforms)
}

.apply_parameter_transforms <- function(mat, transforms) {
  stopifnot(length(transforms) == nrow(mat))

  for (i in seq_len(nrow(mat))) {
    ff <- transforms[[i]]
    if (is.null(ff)) {
      next
    }

    vals <- ff(as.numeric(mat[i, ]))
    stopifnot(length(vals) == ncol(mat))
    mat[i, ] <- vals
  }

  mat
}

.delta_transform_batch_ci <- function(ci, transforms) {
  stopifnot(
    is.list(ci),
    is.numeric(ci$estimates),
    is.matrix(ci$covariance),
    is.list(transforms),
    length(transforms) == length(ci$estimates)
  )

  param_names <- names(ci$estimates)
  mean_raw <- as.numeric(ci$estimates)

  mean_trans <- numeric(length(mean_raw))
  deriv_trans <- numeric(length(mean_raw))

  for (i in seq_along(mean_raw)) {
    ff <- transforms[[i]]
    if (is.null(ff)) {
      ff <- identity
    }
    mean_trans[i] <- ff(mean_raw[i])
    deriv_trans[i] <- .delta_derivative(ff, mean_raw[i])
  }

  jacobian <- diag(deriv_trans, nrow = length(deriv_trans), ncol = length(deriv_trans))
  cov_trans <- jacobian %*% ci$covariance %*% jacobian

  std_error <- sqrt(pmax(diag(cov_trans), 0) / (ci$n_eff * ci$n_chains))
  conf_int <- cbind(
    lower = mean_trans - ci$z_value * std_error,
    upper = mean_trans + ci$z_value * std_error
  )

  names(mean_trans) <- param_names
  names(std_error) <- param_names
  dimnames(cov_trans) <- list(param_names, param_names)
  rownames(conf_int) <- param_names

  ci$estimates <- mean_trans
  ci$covariance <- cov_trans
  ci$std_error <- std_error
  ci$conf_int <- conf_int

  if (!is.null(ci$chain_results)) {
    ci$chain_results <- lapply(ci$chain_results, function(cr) {
      cr_mean <- as.numeric(cr$mean)
      cr_mean_trans <- numeric(length(cr_mean))
      cr_deriv <- numeric(length(cr_mean))
      for (j in seq_along(cr_mean)) {
        ff <- transforms[[j]]
        if (is.null(ff)) {
          ff <- identity
        }
        cr_mean_trans[j] <- ff(cr_mean[j])
        cr_deriv[j] <- .delta_derivative(ff, cr_mean[j])
      }
      cr_jacobian <- diag(cr_deriv, nrow = length(cr_deriv), ncol = length(cr_deriv))
      cr$mean <- stats::setNames(cr_mean_trans, param_names)
      cr$covariance <- cr_jacobian %*% cr$covariance %*% cr_jacobian
      dimnames(cr$covariance) <- list(param_names, param_names)
      if (!is.null(cr$batch_means) && is.matrix(cr$batch_means)) {
        for (j in seq_along(transforms)) {
          ff <- transforms[[j]]
          if (is.null(ff)) {
            ff <- identity
          }
          cr$batch_means[, j] <- ff(as.numeric(cr$batch_means[, j]))
        }
      }
      cr
    })
  }

  ci
}

.delta_derivative <- function(ff, x) {
  h <- 1e-6 * max(1, abs(x))
  if (!is.finite(h) || h <= 0) {
    h <- 1e-6
  }

  x_plus <- x + h
  x_minus <- x - h
  y_plus <- ff(x_plus)
  y_minus <- ff(x_minus)

  if (is.finite(y_plus) && is.finite(y_minus)) {
    return((y_plus - y_minus) / (2 * h))
  }

  # fallback to one-sided differences for domain-constrained transforms
  y0 <- ff(x)
  if (is.finite(y_plus) && is.finite(y0)) {
    return((y_plus - y0) / h)
  }
  if (is.finite(y0) && is.finite(y_minus)) {
    return((y0 - y_minus) / h)
  }

  stop("Failed to compute Delta-method derivative for parameter transform.")
}
