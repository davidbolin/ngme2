#' @title Prior Normal
#' @param mean prior mean
#' @param sd prior standard deviation
#' @param target apply prior on coefficient scale (`"coef"`) or field scale (`"field"`)
#' @return prior specification
#' @export
prior_normal <- function(mean = 0, sd = 1, target = "coef") {
  target <- .validate_prior_target(target)
  stopifnot(
    is.numeric(mean), length(mean) == 1, is.finite(mean),
    is.numeric(sd), length(sd) == 1, is.finite(sd), sd > 0
  )
  structure(
    list(dist = "normal", hyper = c(mean = mean, sd = sd), target = target),
    class = "ngme_prior_spec"
  )
}

#' @title Prior PC-SD
#' @param u tail event threshold
#' @param alpha tail probability
#' @param target apply prior on coefficient scale (`"coef"`) or field scale (`"field"`)
#' @return prior specification
#' @export
prior_pc_sd <- function(u, alpha, target = "coef") {
  target <- .validate_prior_target(target)
  stopifnot(
    is.numeric(u), length(u) == 1, is.finite(u), u > 0,
    is.numeric(alpha), length(alpha) == 1, is.finite(alpha),
    alpha > 0, alpha < 1
  )
  structure(
    list(dist = "pc.sd", hyper = c(u = u, alpha = alpha), target = target),
    class = "ngme_prior_spec"
  )
}

#' @title Prior Half-Cauchy
#' @param scale half-Cauchy scale
#' @param target apply prior on coefficient scale (`"coef"`) or field scale (`"field"`)
#' @return prior specification
#' @export
prior_half_cauchy <- function(scale = 1, target = "coef") {
  target <- .validate_prior_target(target)
  stopifnot(
    is.numeric(scale), length(scale) == 1, is.finite(scale), scale > 0
  )
  structure(
    list(dist = "half.cauchy", hyper = c(scale = scale), target = target),
    class = "ngme_prior_spec"
  )
}

#' @title Prior Inverse-Exponential
#' @description Prior induced by \eqn{\kappa = 1 / \nu \sim \mathrm{Exp}(\lambda)},
#' giving \eqn{p(\nu) = \lambda \exp(-\lambda / \nu)\nu^{-2}} for \eqn{\nu > 0}.
#' Internally this prior is applied to \eqn{\nu = \mathrm{lower} + \exp(\theta)}.
#' @param lambda exponential rate on \eqn{\kappa = 1/\nu}
#' @param lower lower shift used in \eqn{\nu = \mathrm{lower} + \exp(\theta)}
#' @param target apply prior on coefficient scale (`"coef"`) or field scale (`"field"`)
#' @return prior specification
#' @export
prior_inv_exponential <- function(lambda = 1, lower = 0, target = "coef") {
  target <- .validate_prior_target(target)
  stopifnot(
    is.numeric(lambda), length(lambda) == 1, is.finite(lambda), lambda > 0,
    is.numeric(lower), length(lower) == 1, is.finite(lower), lower >= 0
  )
  structure(
    list(
      dist = "inv.exponential",
      hyper = c(lambda = lambda, lower = lower),
      target = target
    ),
    class = "ngme_prior_spec"
  )
}

#' @rdname prior_inv_exponential
#' @export
prior_inv_exp <- function(lambda = 1, lower = 0, target = "coef") {
  prior_inv_exponential(lambda = lambda, lower = lower, target = target)
}

#' @title PC Prior for the NIG Flexibility Parameter
#' @description
#' Penalised-complexity prior for the flexibility parameter
#' \eqn{\eta = 1/\nu} of NIG driving noise, following Cabral, Bolin and Rue
#' (2023). The distance to the Gaussian base model is proportional to
#' \eqn{\eta} (their Corollary 3.1.1), so the PC prior is
#' \eqn{\eta \sim \mathrm{Exp}(\lambda)}, which is exactly
#' \code{\link{prior_inv_exponential}} on \eqn{\nu}.
#'
#' The rate is calibrated as in their Section 3.3, from the statement
#' \eqn{\Pr(\eta > U) = \alpha}, giving \eqn{\lambda = -\log(\alpha)/U}.
#'
#' For symmetric noise the excess kurtosis of a unit-\eqn{h} increment is
#' \eqn{3\eta}, which is the easiest way to read the defaults: \eqn{U = 2.5}
#' is an excess kurtosis of 7.5 (equivalently \eqn{\nu = 0.4}), and
#' \eqn{\alpha = 0.01} puts 1\% of the prior mass above it. That gives
#' \eqn{\lambda \approx 1.84}, a mild contraction: the prior median is
#' \eqn{\nu \approx 2.7} (excess kurtosis 1.1), and the prior is close to
#' neutral over the range of \eqn{\nu} these models typically occupy.
#'
#' The defaults are read from \code{getOption("ngme2.pc_nu_U")} and
#' \code{getOption("ngme2.pc_nu_alpha")}, so the package-wide default used by
#' \code{\link{f}} can be changed without passing a prior to every model.
#'
#' @param U upper value \eqn{U > 0} for \eqn{\eta = 1/\nu} in the calibration
#'   statement \eqn{\Pr(\eta > U) = \alpha}
#' @param alpha prior probability \eqn{\alpha \in (0, 1)} of exceeding \code{U}
#' @param lower lower shift used in \eqn{\nu = \mathrm{lower} + \exp(\theta)}
#' @param target apply prior on coefficient scale (`"coef"`) or field scale (`"field"`)
#' @return prior specification
#' @references
#' Cabral, R., Bolin, D., and Rue, H. (2023). Controlling the Flexibility of
#' Non-Gaussian Processes Through Shrinkage Priors. \emph{Bayesian Analysis},
#' 18(4), 1223-1246. \doi{10.1214/22-BA1342}
#' @seealso \code{\link{prior_inv_exponential}}
#' @export
prior_pc_nu <- function(
    U = getOption("ngme2.pc_nu_U", 2.5),
    alpha = getOption("ngme2.pc_nu_alpha", 0.01),
    lower = 0,
    target = "coef") {
  stopifnot(
    is.numeric(U), length(U) == 1, is.finite(U), U > 0,
    is.numeric(alpha), length(alpha) == 1, is.finite(alpha),
    alpha > 0, alpha < 1
  )
  prior_inv_exponential(
    lambda = -log(alpha) / U, lower = lower, target = target
  )
}

#' @title Prior None
#' @param target apply prior on coefficient scale (`"coef"`) or field scale (`"field"`)
#' @return prior specification
#' @export
prior_none <- function(target = "coef") {
  target <- .validate_prior_target(target)
  structure(
    list(dist = "none", hyper = double(0), target = target),
    class = "ngme_prior_spec"
  )
}

#' @title Prior Container
#' @param ... named prior specifications
#' @return named prior container
#' @export
priors <- function(...) {
  ps <- list(...)
  if (length(ps) == 0) {
    return(structure(ps, class = "ngme_priors"))
  }
  nms <- names(ps)
  stopifnot(
    "priors(...) requires named entries." = !is.null(nms) &&
      all(nzchar(nms)) &&
      length(unique(nms)) == length(nms),
    "All entries in priors(...) must be prior_*() objects." =
      all(vapply(ps, inherits, logical(1), what = "ngme_prior_spec"))
  )
  structure(ps, class = "ngme_priors")
}

.validate_prior_target <- function(target) {
  stopifnot(
    is.character(target),
    length(target) == 1,
    target %in% c("coef", "field")
  )
  target
}

is_prior_spec <- function(x) inherits(x, "ngme_prior_spec")
is_prior_collection <- function(x) inherits(x, "ngme_priors")

as_internal_prior <- function(prior) {
  stopifnot("prior must be a prior_*() object." = is_prior_spec(prior))

  target <- .validate_prior_target(prior$target)
  dist <- prior$dist
  hyper <- prior$hyper
  internal <- switch(dist,
    "none" = list(type = "none", param = double(0), target = target),
    "normal" = {
      mean <- as.numeric(hyper["mean"])
      sd <- as.numeric(hyper["sd"])
      list(type = "normal", param = c(mean, 1 / (sd * sd)), target = target)
    },
    "pc.sd" = {
      u <- as.numeric(hyper["u"])
      alpha <- as.numeric(hyper["alpha"])
      lambda <- -log(alpha) / u
      list(type = "pc.sd", param = c(lambda), target = target)
    },
    "half.cauchy" = {
      list(type = "half.cauchy", param = c(as.numeric(hyper["scale"])), target = target)
    },
    "inv.exponential" = {
      lambda <- as.numeric(hyper["lambda"])
      lower <- as.numeric(hyper["lower"])
      list(type = "inv.exponential", param = c(lambda, lower), target = target)
    },
    stop("Unknown prior dist: ", dist)
  )
  internal
}

default_noise_priors <- function() {
  priors(
    mu = prior_normal(0, 10, target = "coef"),
    sigma = prior_normal(0, 10, target = "coef"),
    nu = prior_normal(0, 10, target = "coef")
  )
}

default_operator_prior <- function() {
  prior_normal(0, 10, target = "coef")
}

default_beta_prior <- function() {
  prior_normal(0, 10, target = "coef")
}

default_beta_prior_for_param <- function(param_name) {
  if (!is.null(param_name) && startsWith(param_name, "(Intercept)")) {
    return(prior_none(target = "coef"))
  }
  default_beta_prior()
}

default_beta_priors <- function(param_names) {
  lapply(param_names, default_beta_prior_for_param)
}

compile_noise_priors <- function(prior) {
  if (is.null(prior)) prior <- default_noise_priors()
  if (is_prior_spec(prior)) {
    prior <- priors(mu = prior, sigma = prior, nu = prior)
  }
  stopifnot("noise prior must be priors(...)." = is_prior_collection(prior))

  allowed <- c("mu", "sigma", "nu")
  bad <- setdiff(names(prior), allowed)
  if (length(bad) > 0) {
    stop("Unknown noise prior names: ", paste(bad, collapse = ", "))
  }

  defaults <- default_noise_priors()
  for (nm in names(prior)) defaults[[nm]] <- prior[[nm]]

  list(
    mu = as_internal_prior(defaults$mu),
    sigma = as_internal_prior(defaults$sigma),
    nu = as_internal_prior(defaults$nu)
  )
}

compile_operator_priors <- function(prior, param_names) {
  if (is.null(param_names) || length(param_names) == 0) {
    return(list())
  }
  if (is.null(prior)) {
    return(rep(list(as_internal_prior(default_operator_prior())), length(param_names)))
  }

  if (is_prior_spec(prior)) {
    prior <- priors(theta = prior)
  }
  stopifnot("operator prior must be prior_*() or priors(...)." = is_prior_collection(prior))

  unknown <- setdiff(names(prior), c("theta", param_names))
  if (length(unknown) > 0) {
    stop("Unknown operator prior names: ", paste(unknown, collapse = ", "))
  }

  per_param <- rep(list(default_operator_prior()), length(param_names))
  if ("theta" %in% names(prior)) {
    per_param <- rep(list(prior$theta), length(param_names))
  }
  for (i in seq_along(param_names)) {
    nm <- param_names[[i]]
    if (nm %in% names(prior)) per_param[[i]] <- prior[[nm]]
  }

  lapply(per_param, as_internal_prior)
}

compile_beta_priors <- function(prior, param_names) {
  if (is.null(param_names) || length(param_names) == 0) {
    return(list())
  }
  if (is.null(prior)) {
    return(lapply(default_beta_priors(param_names), as_internal_prior))
  }

  if (is_prior_spec(prior)) {
    prior <- priors(beta = prior)
  }
  stopifnot("beta prior must be prior_*() or priors(...)." = is_prior_collection(prior))

  index_names <- paste0("feff_", seq_along(param_names))
  unknown <- setdiff(names(prior), c("beta", param_names, index_names))
  if (length(unknown) > 0) {
    stop("Unknown beta prior names: ", paste(unknown, collapse = ", "))
  }

  per_param <- default_beta_priors(param_names)
  if ("beta" %in% names(prior)) {
    per_param <- rep(list(prior$beta), length(param_names))
  }
  for (i in seq_along(param_names)) {
    nm <- param_names[[i]]
    idx_nm <- index_names[[i]]
    if (idx_nm %in% names(prior)) per_param[[i]] <- prior[[idx_nm]]
    if (nm %in% names(prior)) per_param[[i]] <- prior[[nm]]
  }

  internal <- lapply(per_param, as_internal_prior)
  bad_target <- which(vapply(internal, function(x) x$target != "coef", logical(1)))
  if (length(bad_target) > 0) {
    stop("beta prior currently supports target='coef' only.")
  }
  internal
}
