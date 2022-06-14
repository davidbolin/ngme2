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
control.ngme <- function(burnin            = 100,
                         iterations        = 100,
                         gibbs_sample      = 5,
                         stepsize          = 1,
                         opt_fix_effect    = TRUE,

                         kill_var           = FALSE,
                         kill_power         = 0.75,
                         threshold         = 1e-5,
                         termination       = 1e-7
                         ) {

  if ((kill_power <= 0.5) || (kill_power > 1)) {
    stop("reduceVar should be in (0.5,1]")
  }

  control = list( burnin            = burnin,
                  iterations        = iterations,
                  gibbs_sample      = gibbs_sample,
                  stepsize          = stepsize,
                  opt_fix_effect    = opt_fix_effect,

                  # variance reduction
                  kill_var          = kill_var,
                  kill_power        = kill_power,
                  threshold         = threshold,
                  termination       = termination
                )

  class(control) <- "control.ngme"

  return (control)
}

#' Generate control specifications for f function
#'
#' @param opt_operator  whether to optimize operator parameters
#' @param opt_mu        whether to optimize mu
#' @param opt_sigma     whether to optimize sigma
#' @param opt_var       whether to optimize var parameters
#' @param init_operator initial value for operator parameters
#' @param init_mu       initial value for mu
#' @param init_sigma    initial value for sigma
#' @param init_var      initial value for var parameters
#' @param numer_grad    whether to use numerical gradient
#' @param use_precond   whether to use preconditioner
#' @param eps           eps for numerical gradient
#'
#' @return list of control variables
#' @export
#'
#' @examples
control.f <- function(
  opt_operator  = TRUE,
  opt_mu        = TRUE,
  opt_sigma     = TRUE,
  opt_var       = TRUE,

  init_operator = 1, # for spde model, use spde(theta.int=...)
  init_mu       = 0,
  init_sigma    = 1,
  init_var      = 1,

  numer_grad    = FALSE,
  use_precond   = TRUE,
  eps           = 0.01
  ) {

  control = list(
    opt_operator  = opt_operator,
    opt_mu        = opt_mu,
    opt_sigma     = opt_sigma,
    opt_var       = opt_var,

    init_operator = init_operator,
    init_mu       = init_mu,
    init_sigma    = init_sigma,
    init_var      = init_var,

    numer_grad    = numer_grad,
    use_precond   = use_precond,
    eps           = eps
    )

  class(control) <- "control.f"

  return (control)
}

#' Generate debug option
#'
#' @param debug
#' @param trueW
#' @param trueSV
#' @param fixW
#' @param fixSV
#' @param fixSigEps
#' @param sigEps
#'
#' @return
#' @export
#'
#' @examples
debug.ngme <- function(
  debug     = TRUE,
  fixW      = FALSE,
  fixSV     = FALSE,
  fixSigEps = FALSE,
  sigEps    = NULL,
  trueW     = NULL,
  trueSV    = NULL
) {
  if (fixW  && is.null(trueW))  stop("Please provide W")
  if (fixSV && is.null(trueSV)) stop("Please provide SV")

  control  = list(
    debug     = debug,
    fixW      = fixW,
    fixSV     = fixSV,
    trueW     = trueW,
    trueSV    = trueSV,
    fixSigEps = fixSigEps,
    sigEps    = sigEps
  )

  return (control)
}
