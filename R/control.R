
#' Generate control specifications for \code{ngme()} function.
#'
#' These are configurations for \code{ngme}
#' optimization process.
#'
#' @details
#'  To enable convergence check, we need multiple chains running.
#'  We compare the trend of the estimated parameter of length
#'  \code{n_slope_check} (linear regression) with \code{trend_lim}.
#'  We compare the standard devation of estimated parameters (in different chains)
#'  with std_lim.
#' @param seed  set the seed for pesudo random number generator
#' @param burnin          interations for burn-in periods (before optimization)
#' @param iterations      optimization iterations
#' @param estimation      run the estimation process (call C++ in backend)
#' @param standardize_fixed  whether or not standardize the fixed effect
#'
#' @param n_parallel_chain number of parallel chains
#' @param stop_points     number of stop points for convergence check (or specify iters_per_check)
#' @param iters_per_check run how many iterations between each check point (or specify stop_points)
#' @param exchange_VW     exchange last V and W in each chian
#' @param n_slope_check   number of stop points for regression
#' @param std_lim         maximum allowed standard deviation
#' @param trend_lim       maximum allowed slope
#' @param print_check_info print the convergence information
#' @param optimizer choose different sgd optimization method, 
#' currently support "vanilla", "momentum", "adagrad", "rmsprop", "adam", "adamW"
#' see ?vanilla, ?momentum, ?adagrad, ?rmsprop, ?adam, ?adamW
#'
#'
#' @param preconditioner  preconditioner, can be c("none", "fast", "full")
#' "none" means no preconditioner, "fast" means precondition everything except for the parameter of K matrix (for speed reason), "full" means precondition everything
#' @param precond_eps   numerical, the gap used for estimate preconditioner, default is 1e-5
#' @param precond_by_diff_chain logical, if TRUE, use different chains to estimate preconditioner (only computed at check points), if FALSE, use the same chain to estimate preconditioner (computed at each iteration)
#' @param compute_precond_each_iter logical, if TRUE, compute preconditioner at each iteration, if FALSE, only compute preconditioner at check points (if has only 1 chain running, it will be set TRUE)
#'
#' @param max_num_threads maximum number of threads used for parallel computing, by default will be set same as n_parallel_chain.
#' If it is more than n_parallel_chain, the rest will be used to parallel different replicates of the model.
#' @param max_relative_step   max relative step allowed in 1 iteration
#' @param max_absolute_step   max absolute step allowed in 1 iteration
#' @param rao_blackwellization  use rao_blackwellization
#' @param n_trace_iter  use how many iterations to approximate the trace (Hutchinson’s trick)
#'
#' @param reduce_var      logical, reduce variace
#' @param reduce_power    numerical the power of reduce level
#' @param threshold       till when start to reduce the variance
#' @param window_size     numerical, length of window for final estimates
#'
#' @param verbose print estimation
#' @param sampling_strategy subsampling method of replicates of model, c("all", "is")
#' "all" means using all replicates in each iteration,
#' "ws" means weighted sampling (each iteration use 1 replicate to compute the gradient, the sample probability is proption to its number of observations)
#' @return list of control variables
#' @export
control_opt <- function(
  seed              = Sys.time(),
  burnin            = 100,
  iterations        = 500,
  estimation        = TRUE,
  standardize_fixed  = TRUE,
  stop_points       = 10,
  iters_per_check   = iterations / stop_points,

  # parallel options
  n_parallel_chain  = 4,
  max_num_threads   = n_parallel_chain,

  exchange_VW       = TRUE,
  n_slope_check     = 3,
  std_lim           = 0.1,
  trend_lim         = 0.01,
  print_check_info  = FALSE,

  # preconditioner related
  preconditioner    = "none",
  precond_eps       = 1e-5,
  precond_by_diff_chain = FALSE,
  compute_precond_each_iter = FALSE,

  optimizer = adamW(),

  max_relative_step = 0.5,
  max_absolute_step = 0.5,

  # reduce variance after conv. check
  reduce_var        = FALSE,
  reduce_power      = 0.75,
  threshold         = 1e-5,
  window_size       = 1,

  rao_blackwellization = FALSE,
  n_trace_iter      = 10,
  sampling_strategy = "all",

  # opt print
  verbose           = FALSE
) {
  strategy_list <- c("all", "ws")
  preconditioner_list <- c("none", "fast", "full")

  # if user inputs iters_per_check
  if (!missing(iters_per_check) && !missing(stop_points)) {
    stop("Specify only one of iters_per_check and stop_points")
  } else if (!missing(iters_per_check)) {
    stopifnot("iterations should be multiple of iters_per_check"
      = iterations %% iters_per_check == 0)
    stop_points <- iterations / iters_per_check
  }

  stopifnot(
    sampling_strategy %in% strategy_list,
    preconditioner %in% preconditioner_list,
    is.numeric(max_num_threads) && length(max_num_threads) == 1,
    iterations > 0 && stop_points > 0,
    "iterations should be multiple of stop_points"
      = iterations %% stop_points == 0,
    inherits(optimizer, "ngme_optimizer")
  )

  if ((reduce_power <= 0.5) || (reduce_power > 1)) {
    stop("reduceVar should be in (0.5,1]")
  }

  if (n_parallel_chain == 1) {
    # compute_precond_each_iter <- TRUE
    precond_by_diff_chain <- FALSE
  }

  control <- list(
    seed              = seed,
    burnin            = burnin,
    iterations        = iterations,
    estimation        = estimation,
    standardize_fixed  = standardize_fixed,

    n_parallel_chain  = n_parallel_chain,
    stop_points       = stop_points,
    exchange_VW       = exchange_VW,
    n_slope_check     = n_slope_check, # how many on regression check
    std_lim           = std_lim,
    trend_lim         = trend_lim,

    num_threads       = c(
      max(n_parallel_chain, 1),
      max(floor(max_num_threads / n_parallel_chain), 1)
    ),
    rao_blackwellization = rao_blackwellization,
    n_trace_iter      = n_trace_iter,

    print_check_info  = print_check_info,

    # variance reduction
    max_relative_step = max_relative_step,
    max_absolute_step = max_absolute_step,
    reduce_var        = reduce_var,
    reduce_power      = reduce_power,
    threshold         = threshold,
    window_size       = window_size,

    verbose           = verbose,
    precond_eps       = precond_eps,
    precond_by_diff_chain = precond_by_diff_chain,
    compute_precond_each_iter = compute_precond_each_iter,
    precond_strategy  = which(preconditioner_list == preconditioner) - 1, # start from 0
    sampling_strategy = which(strategy_list == sampling_strategy) - 1, # start from 0,

    # optimization method related 
    stepsize          = optimizer$stepsize,
    sgd_method        = optimizer$method,
    sgd_parameters    = optimizer$sgd_parameters
  )

  class(control) <- "control_opt"
  control
}

#' Generate control specifications for \code{f} function
#'
#' @param numer_grad    whether to use numerical gradient
#' @param improve_hessian  improve numerical hessian by using central difference estimation (O(eps^2) error)
#' default is forward difference estimation (O(eps) error)
#'
#' @param eps           eps for computing numerical gradient
#'
#' @return list of control variables
#' @export
control_f <- function(
  numer_grad       = TRUE,
  improve_hessian  = TRUE,
  eps              = 0.0001
  # use_iter_solver = FALSE
  ) {

  control <- list(
    numer_grad       = numer_grad,
    improve_hessian  = improve_hessian,
    eps              = eps,
    use_iter_solver  = FALSE
  )

  class(control) <- "control_f"
  control
}


#' Generate control specifications for the ngme general model
#'
#' @param init_sample_W  sample W|V at the beginning of each chain
#' @param n_gibbs_samples    number of gibbs sampels
#' @param fix_feff       logical, fix fixed effect
#' @param n_post_samples number of posterior samples, see ?ngme_post_samples()
#' @param feff           fixed effect value
#' @param debug          debug mode
#' @return a list of control variables for block model
#' @export
control_ngme <- function(
  init_sample_W = TRUE,
  n_gibbs_samples = 5,
  fix_feff = FALSE,
  n_post_samples = 100,
  feff = NULL,
  debug = FALSE
) {
  control <- list(
    init_sample_W = init_sample_W,
    n_gibbs_samples = n_gibbs_samples,
    fix_feff = fix_feff,
    feff = feff,
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

  control_ngme
}