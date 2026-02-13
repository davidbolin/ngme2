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
#' @param start_sd        standard deviation of the initial parameter (1st chain fixed, other chains random), set 0 to be fixed for all chains
#' @param optimizer choose different sgd optimization method,
#' currently support "sgd", "precond_sgd", "momentum", "adagrad", "rmsprop", "adam", "adamW"
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
#' @param stepsize_decay stepsize decay strategy. Either a character string
#'   ("none" or "grad_norm_plateau") or an object created by
#'   \code{stepsize_decay()}.
#'   decays the stepsize when the mean grad.norm() across chains has not decreased for
#'   \code{stepsize_decay_patience} consecutive epochs (checked at each checkpoint, i.e.
#'   every \code{iterations / n_batch} iterations). All chains decay together after trigger
#'   (after \code{stepsize_decay_warmup} epochs).
#' @param stepsize_decay_patience number of consecutive epochs without grad.norm() improvement
#'   before decaying stepsize.
#' @param stepsize_decay_gamma decay factor applied when triggered (0 < gamma < 1).
#' @param stepsize_decay_min_delta minimum required decrease in grad.norm() to be counted as improvement.
#' @param stepsize_decay_warmup number of initial epochs to skip decay checks.
#' @param stepsize_decay_min_stepsize lower bound for stepsize after decay (absolute value).
#' @param pflug_conv_check use Pflug diagnostic for convergence check
#' @param pflug_alpha scaling factor (0-1] for Pflug criterion: require \code{pflug_sum < pflug_alpha * max_pflug_sum}
#' @param max_R_hat_conv_check use max_R_hat for convergence check
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
    stepsize_decay = c("none", "grad_norm_plateau"),
    stepsize_decay_patience = 3,
    stepsize_decay_gamma = 0.5,
    stepsize_decay_min_delta = 0,
    stepsize_decay_warmup = 0,
    stepsize_decay_min_stepsize = 0,
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

  if (inherits(stepsize_decay, "ngme_stepsize_decay")) {
    if (!missing(stepsize_decay_patience) || !missing(stepsize_decay_gamma) ||
      !missing(stepsize_decay_min_delta) || !missing(stepsize_decay_warmup) ||
      !missing(stepsize_decay_min_stepsize)) {
      warning(
        "stepsize_decay_* arguments are ignored when stepsize_decay() is supplied."
      )
    }
    stepsize_decay_method <- stepsize_decay$method
    stepsize_decay_patience <- stepsize_decay$patience
    stepsize_decay_gamma <- stepsize_decay$gamma
    stepsize_decay_min_delta <- stepsize_decay$min_delta
    stepsize_decay_warmup <- stepsize_decay$warmup
    stepsize_decay_min_stepsize <- stepsize_decay$min_stepsize
  } else {
    stepsize_decay_method <- stepsize_decay
  }

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
  solver_backend_idx <- match(solver_backend, solver_backend_list) - 1L
  solver_factor_idx <- match(solver_factor, solver_factor_list) - 1L

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
          length(stepsize_decay_min_stepsize) == 1 && stepsize_decay_min_stepsize >= 0)
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
