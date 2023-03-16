#' Generate control specifications for \code{ngme()} function.
#'
#' These are configurations for \code{ngme} estimation and
#' optimization process.
#'
#' @details
#'  To enable convergence check, we need multiple chains running.
#'  We compare the trend of the estimated parameter of length
#'  \code{n_slope_check} (linear regression) with \code{trend_lim}.
#'  We compare the standard devation of estimated parameters (in different chains)
#'  with std_lim.
#' @param seed  set the seed for pesudo random number generator
#' @param burnin          burn-in periods
#' @param iterations      optimizing terations
#' @param stepsize        stepsize
#' @param estimation      estimating the parameters
#'
#' @param n_parallel_chain number of parallel chains
#' @param stop_points     number of stop points for convergence check
#' @param exchange_VW     exchange last V and W in each chian
#' @param n_slope_check   number of stop points for regression
#' @param std_lim         maximum allowed standard deviation
#' @param trend_lim       maximum allowed slope
#' @param print_check_info print the convergence information
#'
#' @param opt_beta        logical, optimize fixed effect
#' @param fix_beta        logical, fix fixed effect
#'
#' @param max_relative_step   max relative step allowed in 1 iteration
#' @param max_absolute_step   max absolute step allowed in 1 iteration
#'
#' @param reduce_var      logical, reduce variace
#' @param reduce_power    numerical the power of reduce level
#' @param threshold       till when start to reduce the variance
#' @param window_size     numerical, length of window for final estimates
#'
#' @return list of control variables
#' @export
control_opt <- function(
  seed              = Sys.time(),
  burnin            = 100,
  iterations        = 500,
  stepsize          = 1,
  estimation        = TRUE,

  # parallel options
  n_parallel_chain  = 2,
  stop_points       = 10,
  exchange_VW       = TRUE,
  n_slope_check     = 3,
  std_lim           = 0.1,
  trend_lim         = 0.05,
  print_check_info  = TRUE,

  max_relative_step = 0.2,
  max_absolute_step = 1,

  # reduce variance after conv. check
  reduce_var        = FALSE,
  reduce_power      = 0.75,
  threshold         = 1e-5,
  window_size       = 1,

  # opt print
  verbose          = FALSE
) {
  if ((reduce_power <= 0.5) || (reduce_power > 1)) {
    stop("reduceVar should be in (0.5,1]")
  }

  if (stop_points > iterations) stop_points <- iterations

  control <- list(
    seed              = seed,
    burnin            = burnin,
    iterations        = iterations,
    stepsize          = stepsize,
    estimation        = estimation,

    n_parallel_chain  = n_parallel_chain,
    stop_points       = stop_points,
    exchange_VW       = exchange_VW,
    n_slope_check     = n_slope_check, # how many on regression check
    std_lim           = std_lim,
    trend_lim         = trend_lim,

    print_check_info  = print_check_info,

    # variance reduction
    max_relative_step = max_relative_step,
    max_absolute_step = max_absolute_step,
    reduce_var        = reduce_var,
    reduce_power      = reduce_power,
    threshold         = threshold,
    window_size       = window_size,

    verbose          = verbose
  )

  class(control) <- "control_opt"
  control
}

#' Generate control specifications for \code{f} function
#'
#' @param numer_grad    whether to use numerical gradient
#' @param use_precond   whether to use preconditioner
#' @param use_num_hess  whether to use numerical hessian
#' @param eps           eps for computing numerical gradient
#'
#' @return list of control variables
#' @export
control_f <- function(
  numer_grad    = FALSE,
  use_precond   = FALSE,
  use_num_hess  = TRUE,
  eps           = 0.005
  # use_iter_solver = FALSE
  ) {

  control <- list(
    numer_grad    = numer_grad,
    use_precond   = use_precond,
    use_num_hess  = use_num_hess,
    eps           = eps,
    use_iter_solver = FALSE
  )

  class(control) <- "control_f"
  control
}


#' Generate control specifications for the ngme general model
#'
#' @param init_sample_W  sample W|V at the beginning of each chain
#' @param n_gibbs_samples    number of gibbs sampels
#' @param fix_beta       logical, fix fixed effect
#' @param post_samples_size number of posterior samples
#' @param beta           fixed effect value
#' @return a list of control variables for block model
#' @export
control_ngme <- function(
  init_sample_W = TRUE,
  n_gibbs_samples = 5,
  fix_beta = FALSE,
  post_samples_size = 100,
  beta = NULL,
  debug = FALSE
) {
  control <- list(
    init_sample_W = init_sample_W,
    n_gibbs_samples = n_gibbs_samples,
    fix_beta = fix_beta,
    beta = beta,
    post_samples_size = post_samples_size,
    stepsize = 1,
    debug = debug
  )

  class(control) <- "control_ngme"
  control
}