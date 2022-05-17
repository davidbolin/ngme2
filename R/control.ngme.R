#' Generate control specifications for ngme
#'
#' @param burnin  burn-in iterations
#' @param iterations
#' @param gibbs_sample
#' @param stepsize
#' @param opt_fix_effect
#' @param fix_trueVW
#' @param trueSV
#' @param trueW
#'
#' @return list of control variables
#' @export
#'
#' @examples
control.ngme <- function(burnin            = 100,
                         iterations        = 100,
                         gibbs_sample      = 5,
                         stepsize          = 0.5,

                         opt_fix_effect    = TRUE,
                         fix_trueVW        = FALSE,
                         trueSV            = NULL,
                         trueW             = NULL) {

  control = list(burnin            = burnin,
                 iterations        = iterations,
                 gibbs_sample      = gibbs_sample,
                 stepsize          = stepsize,

                 fix_trueVW        = fix_trueVW,
                 opt_fix_effect    = opt_fix_effect,
                 trueSV            = trueSV,
                 trueW             = trueW)

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
                      eps           = 0.01) {

  control = list(opt_kappa     = opt_kappa,
                opt_mu        = opt_mu,
                opt_sigma     = opt_sigma,
                opt_var       = opt_var,
                numer_grad    = numer_grad,
                use_precond   = use_precond,
                eps           = eps)

  class(control) <- "control.f"

  return (control)
}
