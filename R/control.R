#' Generate control specifications for ngme
#'
#' @param burnin          burn-in iterations
#' @param iterations      optimizing terations
#' @param gibbs_sample    number of gibbs sampels
#' @param stepsize        stepsize
#' @param opt_fix_effect  whether to optimize fix effects
#' @param kill_var        whether to kill the variance
#' @param kill_power      the power of variance killing (if kill the variance)
#' @param threshold       till when start to kill the variance
#' @param termination     till when stop optimizing
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

  opt_beta          = TRUE,
  # fix estimation
  fix_beta          = FALSE,
  fix_mu            = FALSE,
  fix_sigma         = FALSE,
  fix_var           = FALSE,
  fix_V             = FALSE,
  init_V            = NULL,

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
    opt_beta          = opt_beta,
    fix_beta          = fix_beta,
    fix_mu            = fix_mu,
    fix_sigma         = fix_sigma,
    fix_var           = fix_var,
    fix_V             = fix_V,
    init_V            = init_V,

    kill_var          = FALSE,

    # variance reduction
    kill_var          = kill_var,
    kill_power        = kill_power,
    threshold         = threshold,
    termination       = termination
  )

  class(control) <- "control.ngme"
  control
}

#' Generate control specifications for f function
#'
#' @param fix_operator  whether to fix operator parameters
#' @param fix_mu        whether to fix mu
#' @param fix_sigma     whether to fix sigma
#' @param fix_noise     whether to fix noise parameters
#' @param numer_grad    whether to use numerical gradient
#' @param use_precond   whether to use preconditioner
#' @param eps           eps for numerical gradient
#' @param fix_V         to-do
#' @param theta.K       theta.K for parameter K
#' @param use_num_hess  whether to use numerical hessian
#' @param CG_soler      whether to use conjugate gradient solver
#'
#' @return list of control variables
#' @export
#'
#' @examples
ngme.control.f <- function(
  # fix things
  fix_operator  = FALSE,
  fix_mu        = FALSE,
  fix_sigma     = FALSE,
  fix_var     = FALSE,
  fix_V         = FALSE,
  fix_W         = FALSE,

  # initial things
  init_V        = NULL,
  init_W        = NULL,
  theta.K       = 1, # for spde model, use spde(theta.int=...)
  numer_grad    = FALSE,
  use_precond   = TRUE,
  use_num_hess  = TRUE,
  eps           = 0.01,
  use_iter_solver = FALSE
  ) {

  control = list(
    fix_operator  = fix_operator,
    fix_mu        = fix_mu,
    fix_sigma     = fix_sigma,
    fix_var       = fix_var,
    fix_V         = fix_V,
    fix_W         = fix_W,

    # initial things
    init_V        = init_V,
    init_W        = init_W,
    theta.K       = theta.K,

    numer_grad    = numer_grad,
    use_precond   = use_precond,
    use_num_hess  = use_num_hess,
    eps           = eps,
    use_iter_solver = use_iter_solver
    )

  class(control) <- "control.f"

  return (control)
}

#' Generate debug option
#'
#' @param debug    debug mode
#' @param fix_W    fix W
#' @param fix_feff to-do
#' @param fix_merr fix measurement noise
#' @param not_run  print input without running estimation

#'
#' @return
#' @export
#'
#' @examples
ngme.debug <- function(
  debug     = TRUE,
  not_run   = FALSE
) {

  debug  = list(
    debug     = debug,
    not_run    = not_run
  )

  return (debug)
}
