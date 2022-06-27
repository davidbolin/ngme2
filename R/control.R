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
ngme.control <- function(burnin            = 100,
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
#' @param fix_operator  whether to fix operator parameters
#' @param fix_mu        whether to fix mu
#' @param fix_sigma     whether to fix sigma
#' @param fix_noise     whether to fix noise parameters
#' @param numer_grad    whether to use numerical gradient
#' @param use_precond   whether to use preconditioner
#' @param eps           eps for numerical gradient
#'
#' @return list of control variables
#' @export
#'
#' @examples
ngme.control.f <- function(
  fix_operator  = FALSE,
  fix_mu        = FALSE,
  fix_sigma     = FALSE,
  fix_noise     = FALSE,

  theta.K       = 1, # for spde model, use spde(theta.int=...)

  numer_grad    = FALSE,
  use_precond   = TRUE,
  use_num_hess  = TRUE,
  eps           = 0.01
  ) {

  control = list(
    fix_operator  = fix_operator,
    fix_mu        = fix_mu,
    fix_sigma     = fix_sigma,
    fix_noise     = fix_noise,

    theta.K       = theta.K,

    numer_grad    = numer_grad,
    use_precond   = use_precond,
    use_num_hess  = use_num_hess,
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
ngme.debug <- function(
  debug     = TRUE,
  fixW      = FALSE,
  fixSV     = FALSE,
  fixBeta   = FALSE,
  fixSigEps = FALSE,

  trueW     = NULL,
  trueSV    = NULL,
  beta      = NULL,
  sigEps    = NULL
) {
  if (fixW  && is.null(trueW))  stop("Please provide W")
  if (fixSV && is.null(trueSV)) stop("Please provide SV")

  debug  = list(
    debug     = debug,
    fixW      = fixW,
    fixSV     = fixSV,
    trueW     = trueW,
    trueSV    = trueSV,
    fixSigEps = fixSigEps,
    sigEps    = sigEps
  )

  return (debug)
}
