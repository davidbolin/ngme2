#' Generate control specifications for \code{ngme()} function.
#'
#' These are configurations for \code{ngme}
#' optimization process.
#'
#' @details
#' Convergence diagnostics (multi-chain):
#' * R-hat: per-parameter Gelman–Rubin statistic; passes if \code{R_hat <= max_R_hat}.
#' * Trend/Std: uses the last \code{n_slope_check} checkpoints after at least \code{n_min_batch} batches.
#'   Passes when both the relative std (\code{sqrt(var)/|mean| <= std_lim}) and linear trend
#'   of the means (\code{|slope| <= trend_lim}) satisfy their thresholds.
#' * Pflug: per-chain criterion \code{pflug_sum < pflug_alpha * max_pflug_sum} in the latest batch;
#'   if all chains satisfy it, overall convergence is declared.
#' Checks are evaluated every \code{iters_per_check = iterations / n_batch}. A parameter is marked
#' converged if any enabled parameter-level diagnostic (R-hat or Trend/Std) passes; the run stops
#' when all parameters converge or when the Pflug diagnostic triggers.
#' @param seed  set the seed for pesudo random number generator
#' @param burnin          interations for burn-in periods (before optimization)
#' @param iterations      optimization iterations
#' @param estimation      run the estimation process (call C++ in backend)
#' @param standardize_fixed  whether or not standardize the fixed effect. When
#'   \code{TRUE} (default), the design matrix is SVD-standardized internally and
#'   any user-supplied \code{control_ngme(beta_init = ...)} is automatically
#'   mapped from the original design scale to that standardized basis. Set to
#'   \code{FALSE} to keep both the design matrix and any provided
#'   \code{beta_init} on
#'   their original scale.
#'
#' @param n_parallel_chain number of parallel chains
#' @param n_batch     number of checkpoints; optimization is split into \code{n_batch} equal batches
#' @param iters_per_check run how many iterations between each check point (or specify \code{n_batch})
#' @param n_min_batch   minimum number of checkpoints before any convergence diagnostic is attempted
#' @param n_slope_check number of checkpoints used as the regression window for the trend test
#' @param std_lim         maximum allowed standard deviation
#' @param trend_lim       maximum allowed slope
#' @param print_check_info print the convergence information
#' @param start deprecated guard argument. Do not pass model starts through
#'   \code{control_opt()}; use \code{ngme(..., start = previous_fit)} instead.
#' @param start_sd        standard deviation of the initial parameter (1st chain fixed, other chains random), set 0 to be fixed for all chains
#' @param optimizer choose different sgd optimization method,
#' currently support "sgd", "sgld", "precond_sgd", "momentum", "adagrad", "rmsprop", "adam", "adamW"
#' see ?sgd, ?precond_sgd, ?momentum, ?adagrad, ?rmsprop, ?adam, ?adamW
#'
#' @param max_num_threads maximum number of threads used for parallel computing, by default will be set same as n_parallel_chain.
#' If it is more than n_parallel_chain, the rest will be used to parallel different replicates of the model.
#' @param max_relative_step   max relative step allowed in 1 iteration
#' @param max_absolute_step   max absolute step allowed in 1 iteration
#' @param trend_std_conv_check enable the trend/std diagnostic (uses \code{std_lim}, \code{trend_lim}, \code{n_slope_check})
#' @param solver_backend backend in ("eigen", "cholmod", "accelerate", "pardiso")
#' @param solver_type factorization type: "llt" or "ldlt"
#' @param rao_blackwellization  use rao_blackwellization
#' @param n_trace_iter  use how many iterations to approximate the trace (Hutchinson’s trick)
#'
#' @param verbose print estimation
#' @param store_traj store the optimizer trajectory for diagnostics (set FALSE to reduce memory)
#' @param sampling_strategy subsampling method of replicates of model, c("all", "ws")
#' "all" means using all replicates in each iteration,
#' "ws" means weighted sampling (each iteration use 1 replicate to compute the gradient, the sample probability is proption to its number of observations)
#' @param stepsize_control unified stepsize configuration created by
#'   \code{stepsize_control()}, \code{poly_decay()}, or \code{batch_decay()}.
#'   For polynomial schedule, \code{poly_decay(..., burnin_iter = B)} keeps
#'   schedule scale at 1 for the first \code{B} iterations, then starts decay
#'   with reset local time index.
#' @param robust use robust mode in the backend optimizer/model updates
#' @param R_hat_conv_check use the R-hat diagnostic for convergence checking
#' @param pflug_conv_check use Pflug diagnostic for convergence check
#' @param pflug_alpha scaling factor (0-1] for Pflug criterion: require \code{pflug_sum < pflug_alpha * max_pflug_sum}
#' @param max_R_hat maximum allowed R_hat
#' @return list of control variables
#' @export
control_opt <- function(
    seed = Sys.time(),
    burnin = 100,
    iterations = 500,
    estimation = TRUE,
    standardize_fixed = TRUE,
    n_batch = 10,
    iters_per_check = iterations / n_batch,
    optimizer = adam(),
    start = NULL,
    start_sd = 0.5,
    # parallel options
    n_parallel_chain = 4,
    max_num_threads = n_parallel_chain,
    print_check_info = FALSE,
    max_relative_step = 0.5,
    max_absolute_step = 0.5,
    rao_blackwellization = FALSE,
    n_trace_iter = 10,
    sampling_strategy = "all",
    solver_backend = if (Sys.info()["sysname"] == "Darwin") "accelerate" else "cholmod",
    solver_type = "llt",
    # opt print
    verbose = FALSE,
    store_traj = TRUE,
    robust = FALSE,
    stepsize_control = NULL,
    n_min_batch = min(n_batch, 3),
    n_slope_check = min(n_batch, 3),
    trend_std_conv_check = TRUE,
    std_lim = 0.01,
    trend_lim = 0.01,
    R_hat_conv_check = TRUE,
    max_R_hat = 1.1,
    pflug_conv_check = TRUE,
    pflug_alpha = 0.9) {
  strategy_list <- c("all", "ws")
  preconditioner_list <- c("none", "fast", "full")
  solver_backend_list <- c("eigen", "cholmod", "accelerate", "pardiso")
  solver_factor_list <- c("llt", "ldlt")
  stepsize_decay_list <- c("none", "grad_norm_plateau")
  stepsize_schedule_list <- c("constant", "poly")

  if (!is.null(start)) {
    stop(
      "`start` is not a `control_opt()` argument. ",
      "Pass previous fits via `ngme(..., start = previous_fit)` instead.",
      call. = FALSE
    )
  }

  if (is.null(stepsize_control)) {
    stepsize_control <- .default_stepsize_control()
  }
  stopifnot(inherits(stepsize_control, "ngme_stepsize_control"))
  stepsize_decay_method <- stepsize_control$decay$method
  stepsize_decay_patience <- stepsize_control$decay$patience
  stepsize_decay_gamma <- stepsize_control$decay$gamma
  stepsize_decay_min_delta <- stepsize_control$decay$min_delta
  stepsize_decay_warmup <- stepsize_control$decay$warmup
  stepsize_decay_min_stepsize <- stepsize_control$decay$min_stepsize
  stepsize_schedule_method <- stepsize_control$schedule$method
  stepsize_schedule_alpha <- stepsize_control$schedule$alpha
  stepsize_schedule_t0 <- stepsize_control$schedule$t0
  stepsize_schedule_burnin_iter <- stepsize_control$schedule$burnin_iter

  # read preconditioner from optimizer
  preconditioner <- "none"
  numerical_eps <- 1e-5
  precond_by_diff_chain <- FALSE
  compute_precond_each_iter <- FALSE
  if (optimizer$method == "precond_sgd") {
    preconditioner <- optimizer$preconditioner
    numerical_eps <- optimizer$numerical_eps
    precond_by_diff_chain <- optimizer$precond_by_diff_chain
    compute_precond_each_iter <- optimizer$compute_precond_each_iter
  }

  # if user inputs iters_per_check
  if (!missing(iters_per_check) && !missing(n_batch)) {
    stop("Specify only one of iters_per_check and n_batch")
  } else if (!missing(iters_per_check)) {
    stopifnot(
      "iterations should be multiple of iters_per_check" = iterations %% iters_per_check == 0
    )
    n_batch <- iterations / iters_per_check
  }

  # resolve solver backend + factorization; send both to C++ and let it map
  solver_backend <- match.arg(solver_backend, solver_backend_list)
  solver_factor <- match.arg(solver_type, solver_factor_list)
  stepsize_decay_method <- match.arg(stepsize_decay_method, stepsize_decay_list)
  stepsize_schedule_method <- match.arg(stepsize_schedule_method, stepsize_schedule_list)
  solver_backend_idx <- match(solver_backend, solver_backend_list) - 1L
  solver_factor_idx <- match(solver_factor, solver_factor_list) - 1L

  if (identical(optimizer$method, "sgld") &&
      identical(stepsize_schedule_method, "constant") &&
      identical(stepsize_decay_method, "none")) {
    warning(
      "optimizer = sgld() is using constant stepsize (no decay). ",
      "For asymptotically exact SGLD sampling, use diminishing stepsize ",
      "(e.g., stepsize_control = poly_decay(alpha = 0.6, t0 = 10)).",
      call. = FALSE
    )
  }

  # platform guard: accelerate only on macOS; pardiso disabled on macOS builds without MKL
  is_mac <- Sys.info()["sysname"] == "Darwin"
  if (solver_backend == "pardiso" && !has_pardiso()) {
    stop(
      "solver_backend 'pardiso' requires MKL (compiled with USEMKL). ",
      "Reinstall ngme2 with MKLROOT set so MKL/Pardiso can be enabled."
    )
  }
  if (is_mac && solver_backend == "pardiso" && !has_pardiso()) {
    stop("solver_backend 'pardiso' is not available on this macOS build; reinstall with MKLROOT to enable Pardiso.")
  }
  if (!is_mac && solver_backend == "accelerate") {
    stop("solver_backend 'accelerate' is only available on macOS")
  }

  stopifnot(
    sampling_strategy %in% strategy_list,
    preconditioner %in% preconditioner_list,
    "start_sd must be a numeric scalar" =
      is.numeric(start_sd) && length(start_sd) == 1 && is.finite(start_sd),
    is.numeric(max_num_threads) && length(max_num_threads) == 1,
    iterations > 0 && n_batch > 0,
    "iterations should be multiple of n_batch" = iterations %% n_batch == 0,
    "n_min_batch must be numeric" = is.numeric(n_min_batch),
    "n_min_batch must be a single value" = length(n_min_batch) == 1,
    "n_min_batch must be greater than 0" = n_min_batch > 0,
    "n_min_batch cannot be greater than n_batch" = n_min_batch <= n_batch,
    is.numeric(n_slope_check) && length(n_slope_check) == 1 &&
      n_slope_check > 0 && n_slope_check <= n_batch,
    inherits(optimizer, "ngme_optimizer"),
    "solver backend must map to 0:3" = solver_backend_idx %in% 0:3,
    "solver factor must be 0 (llt) or 1 (ldlt)" = solver_factor_idx %in% 0:1,
    is.numeric(pflug_alpha) && length(pflug_alpha) == 1 && pflug_alpha > 0 && pflug_alpha <= 1,
    "stepsize_decay must be one of 'none' or 'grad_norm_plateau'" =
      stepsize_decay_method %in% stepsize_decay_list,
    "stepsize_decay_patience must be >= 1" =
      stepsize_decay_method == "none" ||
        (is.numeric(stepsize_decay_patience) && length(stepsize_decay_patience) == 1 &&
          stepsize_decay_patience >= 1),
    "stepsize_decay_gamma must be in (0, 1)" =
      stepsize_decay_method == "none" ||
        (is.numeric(stepsize_decay_gamma) && length(stepsize_decay_gamma) == 1 &&
          stepsize_decay_gamma > 0 && stepsize_decay_gamma < 1),
    "stepsize_decay_min_delta must be >= 0" =
      stepsize_decay_method == "none" ||
        (is.numeric(stepsize_decay_min_delta) &&
          length(stepsize_decay_min_delta) == 1 && stepsize_decay_min_delta >= 0),
    "stepsize_decay_warmup must be >= 0" =
      stepsize_decay_method == "none" ||
        (is.numeric(stepsize_decay_warmup) && length(stepsize_decay_warmup) == 1 &&
          stepsize_decay_warmup >= 0),
    "stepsize_decay_min_stepsize must be >= 0" =
      stepsize_decay_method == "none" ||
        (is.numeric(stepsize_decay_min_stepsize) &&
          length(stepsize_decay_min_stepsize) == 1 && stepsize_decay_min_stepsize >= 0),
    "stepsize_schedule must be one of 'constant' or 'poly'" =
      stepsize_schedule_method %in% stepsize_schedule_list,
    "stepsize_schedule_alpha must be numeric scalar" =
      is.numeric(stepsize_schedule_alpha) && length(stepsize_schedule_alpha) == 1,
    "stepsize_schedule_t0 must be a non-negative numeric scalar" =
      is.numeric(stepsize_schedule_t0) && length(stepsize_schedule_t0) == 1 &&
        stepsize_schedule_t0 >= 0,
    "stepsize_schedule_burnin_iter must be a non-negative integer scalar" =
      is.numeric(stepsize_schedule_burnin_iter) &&
        length(stepsize_schedule_burnin_iter) == 1 &&
        is.finite(stepsize_schedule_burnin_iter) &&
        stepsize_schedule_burnin_iter >= 0 &&
        abs(stepsize_schedule_burnin_iter - round(stepsize_schedule_burnin_iter)) <
          sqrt(.Machine$double.eps),
    "for stepsize_schedule='poly', alpha must satisfy 1/2 < alpha < 1" =
      stepsize_schedule_method != "poly" ||
        (stepsize_schedule_alpha > 0.5 && stepsize_schedule_alpha < 1)
  )

  if (n_parallel_chain == 1) {
    precond_by_diff_chain <- FALSE
  }

  # variance reduction techniques (not used for now)
  {
    reduce_var <- FALSE
    reduce_power <- 0.75
    threshold <- 1e-5
    window_size <- 1
    stopifnot(
      "reduceVar should be in (0.5,1]" =
        (reduce_power > 0.5) && (reduce_power <= 1)
    )
  }

  control <- list(
    seed = seed,
    start_sd = start_sd,
    burnin = burnin,
    iterations = iterations,
    estimation = estimation,
    standardize_fixed = standardize_fixed,
    n_parallel_chain = n_parallel_chain,
    n_batch = n_batch,
    n_min_batch = n_min_batch, # minimum batches before checking
    n_slope_check = n_slope_check, # window for trend regression
    std_lim = std_lim,
    trend_lim = trend_lim,
    num_threads = c(
      max(n_parallel_chain, 1),
      max(floor(max_num_threads / n_parallel_chain), 1)
    ),
    rao_blackwellization = rao_blackwellization,
    n_trace_iter = n_trace_iter,
    print_check_info = print_check_info,
    verbose = verbose,
    store_traj = store_traj,
    sampling_strategy = which(strategy_list == sampling_strategy) - 1, # start from 0,

    max_relative_step = max_relative_step,
    max_absolute_step = max_absolute_step,
    trend_std_conv_check = trend_std_conv_check,

    # preconditioner related
    numerical_eps = numerical_eps,
    precond_by_diff_chain = precond_by_diff_chain,
    compute_precond_each_iter = compute_precond_each_iter,
    precond_strategy = which(preconditioner_list == preconditioner) - 1, # start from 0

    # optimization method related
    stepsize = optimizer$stepsize,
    sgd_method = optimizer$method,
    sgd_parameters = optimizer$sgd_parameters,
    line_search = optimizer$line_search,

    # solver related
    solver_backend = solver_backend_idx,
    solver_factor = solver_factor_idx,

    # stepsize decay
    stepsize_decay = stepsize_decay_method,
    stepsize_decay_patience = stepsize_decay_patience,
    stepsize_decay_gamma = stepsize_decay_gamma,
    stepsize_decay_min_delta = stepsize_decay_min_delta,
    stepsize_decay_warmup = stepsize_decay_warmup,
    stepsize_decay_min_stepsize = stepsize_decay_min_stepsize,
    stepsize_schedule = stepsize_schedule_method,
    stepsize_schedule_alpha = stepsize_schedule_alpha,
    stepsize_schedule_t0 = stepsize_schedule_t0,
    stepsize_schedule_burnin_iter = as.integer(round(stepsize_schedule_burnin_iter)),

    # variance reduction (not used for now)
    reduce_var = reduce_var,
    reduce_power = reduce_power,
    threshold = threshold,
    window_size = window_size,
    robust = robust,
    R_hat_conv_check = R_hat_conv_check,
    max_R_hat = max_R_hat,
    pflug_conv_check = pflug_conv_check,
    pflug_alpha = pflug_alpha
  )

  class(control) <- "control_opt"
  control
}

#' Generate CI-focused control settings for batch-means inference
#'
#' @description
#' Convenience wrapper for \code{control_opt()} with defaults tailored for
#' Xi-style batch-means confidence intervals.
#'
#' @details
#' This helper enforces trajectory-friendly defaults:
#' \itemize{
#'   \item \code{store_traj = TRUE}
#'   \item \code{trend_std_conv_check = FALSE}
#'   \item \code{R_hat_conv_check = FALSE}
#'   \item \code{pflug_conv_check = FALSE}
#'   \item \code{stepsize_control = poly_decay(alpha, t0, schedule_burnin_iter)}
#' }
#' Any of these can still be overridden through \code{...}.
#'
#' @param optimizer optimizer object, default \code{sgd(stepsize = 0.03)}.
#' @param burnin burn-in iterations before optimization.
#' @param iterations optimization iterations.
#' @param n_batch number of checkpoints.
#' @param n_parallel_chain number of parallel chains.
#' @param alpha polynomial stepsize exponent for \code{poly_decay(alpha, t0)}.
#' @param t0 non-negative schedule offset.
#' @param schedule_burnin_iter non-negative integer. Initial optimization
#'   iterations where polynomial schedule scaling is disabled.
#' @param ... additional arguments forwarded to \code{control_opt()} to override
#'   defaults.
#'
#' @return object of class \code{control_opt}.
#' @export
control_opt_batch_ci <- function(
    optimizer = sgd(stepsize = 0.03),
    burnin = 100,
    iterations = 2000,
    n_batch = 20,
    n_parallel_chain = 4,
    alpha = 0.501,
    t0 = 1,
    schedule_burnin_iter = 0,
    ...) {
  defaults <- list(
    optimizer = optimizer,
    burnin = burnin,
    iterations = iterations,
    n_batch = n_batch,
    n_parallel_chain = n_parallel_chain,
    store_traj = TRUE,
    trend_std_conv_check = FALSE,
    R_hat_conv_check = FALSE,
    pflug_conv_check = FALSE,
    stepsize_control = poly_decay(
      alpha = alpha,
      t0 = t0,
      burnin_iter = schedule_burnin_iter
    )
  )

  args <- utils::modifyList(defaults, list(...))
  do.call(control_opt, args)
}


#' Generate control specifications for the ngme model
#'
#' @param n_gibbs_samples    number of gibbs samples at each iteration
#' @param fix_beta       logical, fix fixed effect estimation
#' @param beta_init      fixed effect initial value on the original design
#'   scale. If \code{control_opt(standardize_fixed = TRUE)} (the default), the
#'   vector is internally rotated/scaled to match the SVD-standardized design
#'   matrix before being passed to the optimizer. Provide values on the raw
#'   design scale; no manual rescaling is required.
#' @param n_post_samples number of posterior samples, see ?ngme_post_samples()
#' @param debug          debug mode
#' @param ... additional arguments. Legacy aliases \code{feff} and
#'   \code{fix_feff} are still recognized and mapped to \code{beta_init} and
#'   \code{fix_beta}.
#' @return a list of control variables for ngme
#' @export
control_ngme <- function(
    n_gibbs_samples = 5,
    fix_beta = FALSE,
    n_post_samples = 100,
    beta_init = NULL,
    debug = FALSE,
    ...) {
  dots <- list(...)
  # backward compatibility for legacy arguments
  if (!is.null(dots$feff)) beta_init <- dots$feff
  if (!is.null(dots$fix_feff)) fix_beta <- dots$fix_feff

  control <- list(
    n_gibbs_samples = n_gibbs_samples,
    fix_beta = fix_beta,
    beta_init = beta_init,
    # legacy names preserved for downstream code until fully migrated
    fix_feff = fix_beta,
    feff = beta_init,
    n_post_samples = n_post_samples,
    debug = debug
  )

  class(control) <- "control_ngme"
  control
}

update_control_ngme <- function(control_ngme, control_opt) {
  control_ngme$rao_blackwellization <- control_opt$rao_blackwellization
  control_ngme$n_trace_iter <- control_opt$n_trace_iter
  control_ngme$stepsize <- control_opt$stepsize
  control_ngme$solver_backend <- control_opt$solver_backend
  control_ngme$solver_factor <- control_opt$solver_factor
  control_ngme$robust <- control_opt$robust

  control_ngme
}
