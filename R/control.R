
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
#' currently support "precond_sgd", "momentum", "adagrad", "rmsprop", "adam", "adamW"
#' see precond_sgd, ?momentum, ?adagrad, ?rmsprop, ?adam, ?adamW
#'
#' @param max_num_threads maximum number of threads used for parallel computing, by default will be set same as n_parallel_chain.
#' If it is more than n_parallel_chain, the rest will be used to parallel different replicates of the model.
#' @param max_relative_step   max relative step allowed in 1 iteration
#' @param max_absolute_step   max absolute step allowed in 1 iteration
#' @param converge_eps        convergence threshold, test if grad.norm() < converge_eps
#'
#' @param rao_blackwellization  use rao_blackwellization
#' @param n_trace_iter  use how many iterations to approximate the trace (Hutchinson’s trick)
#'
#' @param verbose print estimation
#' @param sampling_strategy subsampling method of replicates of model, c("all", "is")
#' "all" means using all replicates in each iteration,
#' "ws" means weighted sampling (each iteration use 1 replicate to compute the gradient, the sample probability is proption to its number of observations)
#' @param solver_type 
#' "eigen" means using eigen solver
#' "cholmod" means using cholmod solver
#' "supernodal" means using supernodal solver
#' "accelerate" means using accelerate solver
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

  optimizer         = adam(),
  
  # parallel options
  n_parallel_chain  = 4,
  max_num_threads   = n_parallel_chain,

  exchange_VW       = TRUE,
  n_slope_check     = 3,
  std_lim           = 0.01,
  trend_lim         = 0.01,
  print_check_info  = FALSE,

  max_relative_step = 0.5,
  max_absolute_step = 0.5,
  converge_eps      = 1e-5,

  rao_blackwellization = FALSE,
  n_trace_iter      = 10,
  sampling_strategy = "all",
  solver_type       = "eigen",

  # opt print
  verbose           = FALSE
) {
  strategy_list <- c("all", "ws")
  preconditioner_list <- c("none", "fast", "full")
  solver_type_list <- c("eigen", "cholmod", "supernodal", "accelerate")

  # read preconditioner from optimizer
  preconditioner            <- "none"
  numerical_eps             <- 1e-5
  precond_by_diff_chain     <- FALSE
  compute_precond_each_iter <- FALSE
  if (optimizer$method == "precond_sgd") {
    preconditioner    = optimizer$preconditioner
    numerical_eps       = optimizer$numerical_eps
    precond_by_diff_chain = optimizer$precond_by_diff_chain
    compute_precond_each_iter = optimizer$compute_precond_each_iter
  }

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
    inherits(optimizer, "ngme_optimizer"),
    solver_type %in% solver_type_list
  )


  if (n_parallel_chain == 1) {
    precond_by_diff_chain <- FALSE
  }

  # variance reduction techniques (not used for now)
  {
    reduce_var        = FALSE
    reduce_power      = 0.75
    threshold         = 1e-5
    window_size       = 1
    stopifnot("reduceVar should be in (0.5,1]" = 
      (reduce_power > 0.5) && (reduce_power <= 1))
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
    verbose           = verbose,
    sampling_strategy = which(strategy_list == sampling_strategy) - 1, # start from 0,

    max_relative_step = max_relative_step,
    max_absolute_step = max_absolute_step,
    converge_eps      = converge_eps,

    # preconditioner related
    numerical_eps       = numerical_eps,
    precond_by_diff_chain = precond_by_diff_chain,
    compute_precond_each_iter = compute_precond_each_iter,
    precond_strategy  = which(preconditioner_list == preconditioner) - 1, # start from 0

    # optimization method related 
    stepsize          = optimizer$stepsize,
    sgd_method        = optimizer$method,
    sgd_parameters    = optimizer$sgd_parameters,
    line_search       = optimizer$line_search,

    # solver related
    solver_type       = which(solver_type_list == solver_type) - 1, # start from 0

    # variance reduction (not used for now)
    reduce_var        = reduce_var,
    reduce_power      = reduce_power,
    threshold         = threshold,
    window_size       = window_size
  )

  class(control) <- "control_opt"
  control
}

#' Generate control specifications for \code{f} function
#'
#' @param numer_grad    whether to use numerical gradient
#' @param improve_hessian  improve numerical hessian by using central difference estimation (O(eps^2) error)
#' default is forward difference estimation (O(eps) error)
#' @param iterative_solver  whether to use iterative solver for the linear system
#' @param eps           eps for computing numerical gradient
#'
#' @return list of control variables
#' @export
control_f <- function(
  numer_grad       = TRUE,
  improve_hessian  = TRUE,
  eps              = 0.0001,
  iterative_solver = FALSE
  ) {

  control <- list(
    numer_grad       = numer_grad,
    improve_hessian  = improve_hessian,
    eps              = eps,
    iterative_solver = iterative_solver
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
#' @param iterative_solver  whether to use iterative solver for the linear system
#' @return a list of control variables for block model
#' @export
control_ngme <- function(
  init_sample_W = TRUE,
  n_gibbs_samples = 5,
  fix_feff = FALSE,
  n_post_samples = 100,
  feff = NULL,
  debug = FALSE,
  iterative_solver = FALSE
) {
  control <- list(
    init_sample_W = init_sample_W,
    n_gibbs_samples = n_gibbs_samples,
    fix_feff = fix_feff,
    feff = feff,
    n_post_samples = n_post_samples,
    debug = debug,
    iterative_solver = iterative_solver
  )

  class(control) <- "control_ngme"
  control
}

update_control_ngme <- function(control_ngme, control_opt) {
  control_ngme$rao_blackwellization <- control_opt$rao_blackwellization
  control_ngme$n_trace_iter <- control_opt$n_trace_iter
  control_ngme$stepsize <- control_opt$stepsize
  control_ngme$solver_type <- control_opt$solver_type

  control_ngme
}