#' Generate control specifications for ngme
#'
#' @param burnin          burn-in iterations
#' @param iterations      optimizing terations
#' @param gibbs_sample    number of gibbs sampels
#' @param stepsize        stepsize
#' @param kill_var        whether to kill the variance
#' @param kill_power      the power of variance killing (if kill the variance)
#' @param threshold       till when start to kill the variance
#' @param termination     till when stop optimizing
#' @param estimation      do estimation or not
#' @param opt_beta
#' @param fix_beta
#' @param window.size
#'
#' @return list of control variables
#' @export
#'
#' @examples
ngme.control <- function(
  burnin            = 100,
  iterations        = 100,
  gibbs_sample      = 5,
  stepsize          = 1,
  estimation        = TRUE,
  n_parallel_chain  = 2,
  stop_points       = 10,
  exchange_VW       = FALSE,
  n_slope_check     = 3,

  opt_beta          = TRUE,
  fix_beta          = FALSE,

  kill_var          = FALSE,
  kill_power        = 0.75,
  threshold         = 1e-5,
  termination       = 1e-7,

  window.size       = 1
) {
  if ((kill_power <= 0.5) || (kill_power > 1)) {
    stop("reduceVar should be in (0.5,1]")
  }

  control <- list(
    burnin            = burnin,
    iterations        = iterations,
    gibbs_sample      = gibbs_sample,
    stepsize          = stepsize,
    estimation        = estimation,
    n_parallel_chain  = n_parallel_chain,
    stop_points       = stop_points,
    exchange_VW       = exchange_VW,
    n_slope_check     = n_slope_check, # how many on regression check

    opt_beta          = opt_beta,
    fix_beta          = fix_beta,

    # variance reduction
    kill_var          = kill_var,
    kill_power        = kill_power,
    threshold         = threshold,
    termination       = termination
  )

  class(control) <- "ngme_control"
  control
}

#' Generate control specifications for f function
#'
#' @param fix_operator  whether to fix operator parameters
#' @param numer_grad    whether to use numerical gradient
#' @param use_precond   whether to use preconditioner
#' @param eps           eps for numerical gradient
#' @param theta.K       theta.K for parameter K
#' @param fix_W
#' @param use_iter_solver
#' @param use_num_hess  whether to use numerical hessian
#'
#' @return list of control variables
#' @export
#'
#' @examples
ngme.control.f <- function(
  numer_grad    = FALSE,
  use_precond   = FALSE,
  use_num_hess  = TRUE,
  eps           = 0.01,
  use_iter_solver = FALSE
  ) {

  control <- list(
    numer_grad    = numer_grad,
    use_precond   = use_precond,
    use_num_hess  = use_num_hess,
    eps           = eps,
    use_iter_solver = use_iter_solver
  )

  class(control) <- "ngme_control_f"
  control
}