#' Generate control specifications for ngme
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

  class(control) <- "ngme.control"

  return (control)
}


