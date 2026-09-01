#' Generate control specifications for \code{ngme()} function.
#'
#' These are configurations for \code{ngme}
#' optimization process.
#'
#' @details
#' Convergence diagnostics (multi-chain):
#' * R-hat: per-parameter Gelman–Rubin statistic; passes if \code{R_hat <= max_R_hat}.
#' * Trend: fits a weighted linear trend over a window of \code{n_slope_check} points and
#'   passes unless the drift is both distinguishable from zero (\code{|slope| / se(slope) >
#'   trend_lim}) and materially fast (relative movement per 100 iterations \code{>
#'   trend_rel_lim}). This is scale-free and independent of \code{n_batch}. Setting
#'   \code{trend_use_tstat = FALSE} recovers the absolute test \code{|slope| <= trend_lim},
#'   and \code{use_std_check = TRUE} additionally requires \code{sqrt(var)/|mean| <= std_lim}.
#'   The window points are sub-batch means spread evenly over the SECOND HALF of the run so
#'   far, so the test is available from the first checkpoint and never regresses through the
#'   optimiser's initial transient.
#' Checks are evaluated every \code{iters_per_check = iterations / n_batch}. A parameter is marked
#' converged only if every enabled parameter-level diagnostic (R-hat and Trend/Std) passes, so a
#' single diagnostic cannot declare convergence on its own; the run stops when all parameters
#' converge. Disable a diagnostic
#' (\code{R_hat_conv_check = FALSE} or \code{trend_std_conv_check = FALSE}) to drop it from the
#' requirement. The criteria must hold on \code{n_conv_batch} consecutive checkpoints, since a
#' single passing checkpoint is weak evidence.
#' @param seed  random seed. Defaults to a seed drawn from the current R
#'   random number stream, so \code{set.seed()} makes the result reproducible.
#' @param burnin          interations for burn-in periods (before optimization)
#' @param iterations      optimization iterations
#' @param estimation      run the estimation process (call C++ in backend)
#' @param standardize_fixed  whether or not standardize the fixed effect. When
#'   \code{TRUE} (default), the design matrix is SVD-standardized internally and
#'   any user-supplied \code{control_ngme(beta_init = ...)} is automatically
#'   mapped from the original design scale to that standardized basis. Set to
#'   \code{FALSE} to keep both the design matrix and any provided
#'   \code{beta_init} on their original scale.
#'
#' @param n_parallel_chain number of parallel chains
#' @param n_batch     number of checkpoints; optimization is split into \code{n_batch} equal batches
#' @param iters_per_check run how many iterations between each check point (or specify \code{n_batch})
#' @param n_min_batch   minimum number of checkpoints before any convergence
#'   diagnostic is attempted (default 1).
#' @param n_slope_check number of checkpoints used as the regression window for the trend test
#' @param trend_use_tstat compare the fitted slope to its own standard error
#'   (\code{|slope| / se(slope) <= trend_lim}) instead of to an absolute bound. Scale-free
#'   and independent of \code{n_batch}; with this on, \code{trend_lim} is a t-value (~2-3).
#' @param use_std_check include the relative-standard-deviation part of the trend/std
#'   diagnostic. The raw coefficient of variation conflates chain disagreement with
#'   parameter imprecision, so weakly identified parameters can never pass it.
#' @param n_conv_batch number of consecutive checkpoints that must satisfy the criteria
#'   before convergence is declared (default 2). Guards against a single lucky
#'   checkpoint.
#' @param warn_no_convergence emit a warning when the iteration budget is exhausted without
#'   the convergence criteria being met. Set \code{FALSE} for short runs where convergence
#'   is not expected (e.g. fast unit tests).
#' @param std_lim         maximum allowed standard deviation
#' @param trend_lim       maximum allowed drift. With \code{trend_use_tstat = TRUE}
#'   (default) this is a t-value on the fitted slope (2 ~ "not distinguishable from no
#'   drift"); with \code{trend_use_tstat = FALSE} it is an absolute bound on the slope.
#' @param trend_rel_lim how fast a parameter must actually be MOVING before a
#'   detectable drift counts as non-convergence, as a fraction of the parameter's
#'   own scale per 100 iterations The default 0.01 reads as "1\% of itself per
#'   100 iterations". A parameter fails the trend test only when its drift is
#'   both statistically detectable and faster than this.
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
#' @param solver_backend backend in ("eigen", "cholmod", "accelerate", "pardiso").
#'   Defaults to "accelerate" on macOS and "cholmod" elsewhere.
#'   Note that the speed of "cholmod" is set by the BLAS R itself is linked against,
#'   and with R's bundled reference BLAS those kernels run one to two orders of
#'   magnitude below a tuned one. Thus, link R against OpenBLAS or MKL before fitting
#'   large models.
#' @param solver_type factorization type: "llt" or "ldlt"
#' @param nonsym_solver how the operator matrix \code{K} of a non-symmetric
#'   model is factorized when estimating \code{tr(K^-1 dK)}.
#'   "normal_equations" (default) takes a Cholesky of \code{t(K) K} and applies
#'   \code{K^-1} as \code{(t(K) K)^-1 t(K)}; "lu" takes a sparse LU of
#'   \code{K}. The Cholesky is the cheaper of the two for every non-symmetric
#'   operator here. Its one cost is that it squares the
#'   condition number of \code{K}, so switch to "lu" for an operator so
#'   ill-conditioned that the Cholesky of \code{t(K) K} fails.
#' @param rao_blackwellization  replace the sampled latent field by its conditional
#'   expectation \code{E[W | Y, V]} in the gradient (default \code{TRUE}).
#' @param n_trace_iter  use how many iterations to approximate the trace (Hutchinson’s trick).
#'   The starting probe budget; with \code{trace_adapt = TRUE} it is retuned during the run.
#' @param polish_iterations iterations of a post-convergence polish phase, run once the
#'   stopping rule has fired (or the budget is spent) with the step size scaled by
#'   \code{polish_stepsize_factor}. The reported estimate is the average of the iterates over
#'   this window (Polyak-Ruppert), which is what reduces its variance; the smaller step only
#'   shrinks the ball being averaged over. Capped at the number of iterations the fit actually
#'   ran, so a very short fit is not handed a disproportionate tail.
#' @param polish_stepsize_factor multiplier applied to the step size during the polish phase.
#'   Measured over six optimizer seeds at fixed data, 0.3 cut the estimator's standard
#'   deviation by about a quarter; too small a factor freezes the chains so there is little
#'   left to average over, and no reduction at all lets them keep wandering.
#' @param selinv_max_fill use the exact selected (Takahashi) inverse for the
#'   Rao-Blackwell traces when the Cholesky factor of QQ has \code{nnz(L)/n} at or below
#'   this, and Hutchinson probes otherwise. Triangular operators (ar1, ou, arma) and 1-d
#'   meshes give a banded QQ with a ratio near 2-3, where the exact route is both cheaper
#'   and free of probe noise; a 2-d mesh is nearer 35, where it is far slower. Set to 0 to
#'   always probe.
#' @param trace_adapt size the Hutchinson probe count automatically (default \code{TRUE}).
#'   The estimator's own variance falls as 1/N and is measured from the spread of the probes,
#'   while the Gibbs sampling noise is measured across iterations; the budget is set so the
#'   former is \code{trace_adapt_frac} of the latter. Only active when Rao-Blackwellisation
#'   is on, since that is where the trace estimators feed the gradient.
#' @param trace_adapt_frac target share of the total gradient variance carried by the trace
#'   estimator (default 0.1, i.e. it adds about 5\% to the gradient standard deviation). A
#'   deadband leaves the budget alone while the measured share is within a factor of two of
#'   this, so the achieved share settles between roughly \code{0.5 * trace_adapt_frac} and
#'   \code{2 * trace_adapt_frac}. Updates are damped and capped at 25\% per step, so the
#'   budget approaches its target over several updates rather than jumping.
#' @param trace_adapt_every how many iterations between probe-budget updates.
#' @param trace_adapt_min,trace_adapt_max bounds on the adapted probe count.
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
#' @param nig_param_std coordinates the optimiser uses for stationary NIG
#'   noise, following Cabral, Bolin and Rue (2023). In the native parameters the
#'   variance \code{h(sigma^2 + mu^2/nu)} is shared by all three, giving a long
#'   flat ridge; these coordinates turn that ridge into an axis. Only the
#'   optimiser's coordinates change and the objective, the priors and the
#'   reported estimates stay native, so fits remain comparable across settings.
#'
#'   \describe{
#'     \item{0}{native \code{(theta_mu, log sigma, log nu)}}
#'     \item{1}{standardised \code{(log sigma_marg, zeta,
#'       log eta)} with \code{eta = 1/nu}, \code{zeta = mu/sigma} and
#'       \code{sigma_marg = sqrt(sigma^2 + mu^2/nu)} the marginal SD}
#'     \item{2}{additionally orthogonalised,
#'       \code{zeta* = zeta sqrt(eta)}, \code{eta* = eta / xi^2} with
#'       \code{xi = 1 + zeta*^2 - |zeta*| sqrt(1 + zeta*^2)}, which makes the
#'       kurtosis invariant to skewness}
#'     \item{3}{as 2 with \code{zeta*} carried as \code{asinh(zeta*)}, which
#'       keeps the skewness coordinate unbounded and better scaled (default)}
#'   }
#' @param robust use robust mode in the backend optimizer/model updates
#' @param continue_chains make \code{ngme(start = previous_fit)} a true
#'   continuation (default \code{TRUE}). Three things follow from it:
#'   \code{start_sd} is not applied,  each chain resumes
#'   from its OWN final parameters and its own latent state (\code{W} and the
#'   mixing variables \code{V}); and
#'   the stored optimisation trajectory is concatenated across the restarts
#'   instead of replaced. A run split into chunks is then the same optimisation
#'   as one long run.
#'
#'   Set \code{FALSE} to recover the old behaviour, which is what you want if
#'   you are deliberately re-dispersing chains from a fitted point to test
#'   whether they return to it. The per-chain part is skipped, with a warning,
#'   when the previous fit has a different parameterization; the no-jitter part
#'   still applies, since it needs nothing from the previous fit but the fact
#'   that it is one.
#' @param R_hat_conv_check use the R-hat diagnostic for convergence checking
#' @param max_R_hat maximum allowed R_hat
#' @return list of control variables
#' @export
control_opt <- function(
    seed = ngme_random_seed(),
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
    rao_blackwellization = TRUE,
    n_trace_iter = 10,
    polish_iterations = 50L,
    polish_stepsize_factor = 0.3,
    selinv_max_fill = 4,
    trace_adapt = TRUE,
    trace_adapt_frac = 0.1,
    trace_adapt_every = 100L,
    trace_adapt_min = 5L,
    trace_adapt_max = 200L,
    sampling_strategy = "all",
    solver_backend = if (Sys.info()["sysname"] == "Darwin") "accelerate" else "cholmod",
    solver_type = "llt",
    nonsym_solver = "normal_equations",
    # opt print
    verbose = FALSE,
    store_traj = TRUE,
    robust = FALSE,
    nig_param_std = 3,
    stepsize_control = NULL,
    n_min_batch = 1,
    n_slope_check = min(n_batch, 3),
    trend_std_conv_check = TRUE,
    std_lim = 0.01,
    trend_lim = 2,
    trend_rel_lim = 0.01,
    trend_use_tstat = TRUE,
    use_std_check = FALSE,
    n_conv_batch = 2,
    warn_no_convergence = TRUE,
    continue_chains = TRUE,
    R_hat_conv_check = TRUE,
    max_R_hat = 1.1) {
  strategy_list <- c("all", "ws")
  solver_backend_list <- c("eigen", "cholmod", "accelerate", "pardiso")
  solver_factor_list <- c("llt", "ldlt")
  nonsym_solver_list <- c("lu", "normal_equations")
  stepsize_decay_list <- c("none", "grad_norm_plateau", "trend")
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

  numerical_eps <- 1e-5
  if (optimizer$method == "precond_sgd") {
    numerical_eps <- optimizer$numerical_eps
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
  nonsym_solver <- match.arg(nonsym_solver, nonsym_solver_list)
  stepsize_decay_method <- match.arg(stepsize_decay_method, stepsize_decay_list)
  stepsize_schedule_method <- match.arg(stepsize_schedule_method, stepsize_schedule_list)
  solver_backend_idx <- match(solver_backend, solver_backend_list) - 1L
  solver_factor_idx <- match(solver_factor, solver_factor_list) - 1L
  nonsym_solver_idx <- match(nonsym_solver, nonsym_solver_list) - 1L

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
    trend_rel_lim = trend_rel_lim,
    trend_use_tstat = trend_use_tstat,
    use_std_check = use_std_check,
    n_conv_batch = n_conv_batch,
    warn_no_convergence = warn_no_convergence,
    continue_chains = continue_chains,
    num_threads = c(
      max(n_parallel_chain, 1),
      max(floor(max_num_threads / n_parallel_chain), 1)
    ),
    rao_blackwellization = rao_blackwellization,
    n_trace_iter = n_trace_iter,
    polish_iterations = polish_iterations,
    polish_stepsize_factor = polish_stepsize_factor,
    selinv_max_fill = selinv_max_fill,
    trace_adapt = trace_adapt,
    trace_adapt_frac = trace_adapt_frac,
    trace_adapt_every = trace_adapt_every,
    trace_adapt_min = trace_adapt_min,
    trace_adapt_max = trace_adapt_max,
    print_check_info = print_check_info,
    verbose = verbose,
    store_traj = store_traj,
    sampling_strategy = which(strategy_list == sampling_strategy) - 1, # start from 0,

    max_relative_step = max_relative_step,
    max_absolute_step = max_absolute_step,
    trend_std_conv_check = trend_std_conv_check,

    # preconditioner related
    numerical_eps = numerical_eps,

    # optimization method related
    stepsize = optimizer$stepsize,
    sgd_method = optimizer$method,
    sgd_parameters = optimizer$sgd_parameters,
    line_search = optimizer$line_search,

    # solver related
    solver_backend = solver_backend_idx,
    solver_factor = solver_factor_idx,
    nonsym_solver = nonsym_solver_idx,

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
    nig_param_std = nig_param_std,
    R_hat_conv_check = R_hat_conv_check,
    max_R_hat = max_R_hat
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
#' @param n_gibbs_samples    number of gibbs samples at each iteration.
#'   Raising it lowers the variance of each gradient but costs proportionally more per
#'   iteration, so the number of iterations needed to satisfy the
#'   convergence checks is largely unchanged and the wall-clock time to
#'   converge goes up. Raise it when the estimates are noisy at
#'   convergence.
#' @param post_burnin number of Gibbs sweeps to discard before recording the
#'   posterior samples of W and V returned by \code{ngme()}. A burn-in is
#'   needed because the chain is restarted from the stored W and V, which are
#'   averaged over the parallel chains and so are not themselves a draw from
#'   the posterior. Set to 0 only if you know the stored state is a valid draw.
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
#'   \code{fix_beta}. Anything else is ignored with a warning naming it.
#' @return a list of control variables for ngme
#' @export
control_ngme <- function(
    n_gibbs_samples = 5,
    post_burnin = 100,
    fix_beta = FALSE,
    n_post_samples = 100,
    beta_init = NULL,
    debug = FALSE,
    ...) {
  dots <- list(...)
  # backward compatibility for legacy arguments
  if (!is.null(dots$feff)) beta_init <- dots$feff
  if (!is.null(dots$fix_feff)) fix_beta <- dots$fix_feff

  # Anything else in ... used to be dropped without a word, which quietly hides
  # real mistakes: control_ngme(rao_blackwellization = FALSE) looks like it
  # works but does nothing, since that argument belongs to control_opt(). Warn
  # rather than stop, so existing scripts keep running.
  legacy <- c("feff", "fix_feff")
  nms <- names(dots)
  if (is.null(nms)) nms <- rep("", length(dots))
  unknown <- setdiff(nms[nzchar(nms)], legacy)
  n_unnamed <- sum(!nzchar(nms))
  if (length(unknown) || n_unnamed) {
    bits <- character(0)
    if (length(unknown))
      bits <- c(bits, paste0("unrecognised argument(s) ",
                             paste0("`", unknown, "`", collapse = ", ")))
    if (n_unnamed)
      bits <- c(bits, sprintf("%d unnamed argument(s)", n_unnamed))
    msg <- paste0("control_ngme(): ", paste(bits, collapse = " and "),
                  " ignored.")
    # The most likely mistake is reaching for a control_opt() argument.
    opt_args <- setdiff(names(formals(control_opt)), "...")
    in_opt <- intersect(unknown, opt_args)
    if (length(in_opt))
      msg <- paste0(msg, "\n  ", paste0("`", in_opt, "`", collapse = ", "),
                    if (length(in_opt) == 1L) " is an argument of"
                    else " are arguments of",
                    " control_opt(), not control_ngme().")
    # Otherwise offer the nearest formal, which catches plain typos.
    own <- setdiff(names(formals(control_ngme)), "...")
    rest <- setdiff(unknown, in_opt)
    near <- unlist(lapply(rest, function(u) {
      cand <- agrep(u, own, max.distance = 0.3, value = TRUE)
      if (length(cand)) sprintf("`%s` -> did you mean `%s`?", u, cand[1]) else NULL
    }))
    if (length(near)) msg <- paste0(msg, "\n  ", paste(near, collapse = "\n  "))
    warning(msg, call. = FALSE)
  }

  control <- list(
    n_gibbs_samples = n_gibbs_samples,
    post_burnin = post_burnin,
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
  control_ngme$selinv_max_fill <- control_opt$selinv_max_fill
  control_ngme$trace_adapt <- control_opt$trace_adapt
  control_ngme$trace_adapt_frac <- control_opt$trace_adapt_frac
  control_ngme$trace_adapt_every <- control_opt$trace_adapt_every
  control_ngme$trace_adapt_min <- control_opt$trace_adapt_min
  control_ngme$trace_adapt_max <- control_opt$trace_adapt_max
  control_ngme$stepsize <- control_opt$stepsize
  control_ngme$solver_backend <- control_opt$solver_backend
  control_ngme$solver_factor <- control_opt$solver_factor
  control_ngme$nonsym_solver <- control_opt$nonsym_solver
  control_ngme$robust <- control_opt$robust
  control_ngme$nig_param_std <- control_opt$nig_param_std

  control_ngme
}
