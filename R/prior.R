#' @title ngme prior specification
#' @description
#' Prior specification for internal representation of the parameters.
#' We will list all the available priors here.
#' Their PDFs can be found in vignette("ngme2-prior").
#'
#' none prior with no parameter: \eqn{p(\theta) \propto 1}
#'
#' Normal prior with parameter mean \eqn{\mu} and precision \eqn{\tau}: \eqn{\theta \sim N(\mu, 1/ \tau)}
#'
#' Half-Cauchy prior with parameter scale \eqn{\tau}: \eqn{\theta \sim HC(0, \tau)}
#'
#' Penalized complexity prior with parameter \eqn{\lambda}
#'
#'
#' @param type type of prior, say ngme_prior_types()
#' @param param parameters of the prior
#'
#' @return a list of prior specification
#' @export
ngme_prior <- function(type, param=double(0)) {
  stopifnot("Please check if the prior name is in ngme_prior_types()" =
    type %in% ngme2::ngme_prior_types())

  # check if num. of parameter is correct
  switch(type,
    "none" = stopifnot(length(param) == 0),
    "normal" = stopifnot(length(param) == 2),
    # internally log precision
    "pc.sd" = {
      # lambda
      # param = - log(alpha) / u
      stopifnot(length(param) == 1)
    },
    "half.cauchy" = stopifnot(length(param) == 1),
    "pc.cor0" = stopifnot(length(param) == 1)
  )

  structure(
    list(type = type, param = param),
    class = "ngme_prior"
  )
}
