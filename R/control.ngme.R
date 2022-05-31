#' Generate control specifications for ngme
#'
#' @param burnin  burn-in iterations
#' @param iterations
#' @param gibbs_sample
#' @param stepsize
#' @param opt_fix_effect
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
#' @param opt_kappa
#' @param opt_mu
#' @param opt_sigma
#' @param opt_var
#' @param numer_grad
#' @param use_precond
#' @param eps
#'
#' @return list of control variables
#' @export
#'
#' @examples
control.f <- function(opt_kappa     = TRUE,
                      opt_mu        = TRUE,
                      opt_sigma     = TRUE,
                      opt_var       = TRUE,
                      numer_grad    = FALSE,
                      use_precond   = TRUE,
                      init_kappa    = 0.5,
                      init_mu       = 0,
                      init_sigma    = 1,
                      init_nu       = 1,
                      eps           = 0.01) {

  control = list(opt_kappa     = opt_kappa,
                opt_mu        = opt_mu,
                opt_sigma     = opt_sigma,
                opt_var       = opt_var,
                numer_grad    = numer_grad,
                use_precond   = use_precond,
                init_kappa    = init_kappa,
                init_mu       = init_mu,
                init_sigma    = init_sigma,
                init_nu       = init_nu,
                eps           = eps)

  class(control) <- "control.f"

  return (control)
}

#' Generate debug option
#'
#' @param debug
#' @param trueW
#' @param trueSV
#'
#' @return
#' @export
#'
#' @examples
debug.ngme <- function(debug     = TRUE,
                       fixW      = FALSE,
                       fixSV     = FALSE,
                       fixSigEps = FALSE,
                       sigEps    = NULL,
                       trueW     = NULL,
                       trueSV    = NULL
) {
  if (fixW  && is.null(trueW))  stop("Please provide W")
  if (fixSV && is.null(trueSV)) stop("Please provide SV")

  control  = list(debug     = debug,
                  fixW      = fixW,
                  fixSV     = fixSV,
                  trueW     = trueW,
                  trueSV    = trueSV,
                  fixSigEps = fixSigEps,
                  sigEps    = sigEps
  )

  return (control)
}
