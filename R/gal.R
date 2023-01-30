#############################################################################################
#' @name gal
#' @title  The Generalized Asymmetric Laplace (GAL) Distribution
#' @aliases dgal rgal pgal qgal
#' @description Density, distribution function, quantile function and
#' random generation for the generalized asymmetric Laplace distribution
#'  with parameters \code{mu}, \code{sigma} and \code{nu}, \code{delta}.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n, number of observations.
#' @param delta A numeric value for the location parameter.
#' @param mu    A numeric value for the shift parameter.
#' @param nu    A numeric value for the shape parameter.
#' @param sigma A numeric value for the scaling parameter.
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \eqn{p} are
#' returned as \eqn{log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X\leq x]},
#' otherwise, \eqn{P[X>x]}.
#' @param seed Seed for the random generation.
#' @return
#' dgal gives the density, pgal gives the distribution function,
#' qgal gives the quantile function, and rgal generates random deviates.
#'
#' Invalid arguments will result in return value NaN, with a warning.
#'
#' The length of the result is determined by \code{n} for rgal.
#' @details
#' The generalized asymmetric Laplace distribution has density given
#' by
#'  \deqn{f(x; p, a, b) =
#'  \frac{e^{\nu+\mu(x-\delta)/\sigma^2}\sqrt{\nu\mu^2/\sigma^2+\nu^2}}{\pi\sqrt{\nu\sigma^2+(x-\delta)^2}}
#'  K_1(\sqrt{(\nu\sigma^2+(x-\delta)^2)(\mu^2/\sigma^4+\nu/\sigma^2)}),}
#'  where \eqn{K_p} is modified Bessel function of the second kind of order \eqn{p},
#'  \eqn{x>0}, \eqn{\nu>0} and \eqn{\mu,\delta, \sigma\in\mathbb{R}}.
#'  See Barndorff-Nielsen (1977, 1978 and 1997) for further details.
#'
#' If the mixing variable $V$ follows a Gamma distribution (same parameterization in R):
#' \deqn{V \sim \Gamma(h \nu, \nu),} then the poserior follows the GAL distribution (a special case of GIG distribution):
#' \deqn{
#' -\mu +\mu V + \sigma \sqrt{V} Z \sim GIG(h \nu - 0.5, 2 \nu + (\frac{\mu}{\sigma})^{2}, 0)
#' }
#' @references
#'  Barndorff-Nielsen, O. (1977) Exponentially decreasing distributions for the logarithm of particle size. Proceedings of the Royal Society of London.
#'
#'  Series A, Mathematical and Physical Sciences. The Royal Society. 353, 401–409. \doi{10.1098/rspa.1977.0041}
#'
#'  Barndorff-Nielsen, O. (1978) Hyperbolic Distributions and Distributions on Hyperbolae, Scandinavian Journal of Statistics. 5, 151–157.
#'
#' @seealso
#' \code{\link{dgig}}, \code{\link{dig}}, \code{\link{digam}}
#' @examples
#' rgal(100, delta = 0, mu = 5, sigma = 1, nu = 1)
#' pgal(0.4, delta = 0, mu = 5, sigma = 1, nu = 1)
#' qgal(0.8, delta = 0, mu = 5, sigma = 1, nu = 1)
#' plot(function(x){dgal(x, delta = 0, mu = 5, sigma = 1, nu = 1)}, main =
#' "generalized asymmetric Laplace density", ylab = "Probability density",
#' xlim = c(0,10))
#' @rdname gal
#' @export
dgal <- function(x, delta, mu, nu, sigma, log=FALSE){
  if(missing(delta)){
    stop('argument "delta" missing, with no default')
  }
  if(missing(mu)){
    stop('argument "mu" missing, with no default')
  }
  if(missing(nu)){
    stop('argument "nu" missing, with no default')
  }
  if(missing(sigma)){
    stop('argument "sigma" missing, with no default')
  }
  if (length(delta) != length(mu)){
    if(length(delta)==1){
      delta <- rep(delta, length(mu))
    } else if(length(mu)==1){
      mu <- rep(mu, length(delta))
    } else{
      stop("delta and mu are vectors of different lenghts")
    }
  }
  if (length(mu) != length(nu)){
    if(length(nu)==1){
      nu <- rep(nu, length(mu))
    } else if(length(mu)==1){
      mu <- rep(mu, length(nu))
      if(length(delta)==1){
        delta <- rep(delta, length(nu))
      }
      } else{
      stop("mu and nu are vectors of different lengths")
    }
  }
  if (length(nu) != length(sigma)){
    if(length(sigma)==1){
      sigma <- rep(sigma, length(nu))
    } else if(length(nu)==1){
      nu <- rep(nu, length(sigma))
      if(length(mu)==1){
        mu <- rep(mu,length(sigma))
      }
      if(length(delta)==1){
        delta <- rep(delta, length(sigma))
      }
      } else{
      stop("nu and sigma are vectors of different lengths")
    }
  }
  if (min(nu) < 0)
    stop("vector nu must be  positive")
  n = length(x)
  if(n < length(delta)){
    delta_new = delta[1:n]
    mu_new = mu[1:n]
    nu_new = nu[1:n]
    sigma_new = sigma[1:n]
  } else if(n>length(delta)){
    quot <- n%/%length(delta)
    rem <- n%%length(delta)
    delta_new <- rep(delta,quot)
    if(rem > 0){
      delta_new <- c(delta_new, delta[1:rem])
    }
    mu_new <- rep(mu,quot)
    if(rem > 0){
      mu_new <- c(mu_new, mu[1:rem])
    }
    nu_new <- rep(nu,quot)
    if(rem > 0){
      nu_new  <- c(nu_new, nu[1:rem])
    }
    sigma_new <- rep(sigma,quot)
    if(rem > 0){
      sigma_new <- c(sigma_new, sigma[1:rem])
    }
  } else{
    delta_new = delta
    mu_new = mu
    nu_new = nu
    sigma_new = sigma
  }

  p_vec = nu_new - 0.5
  a_vec = 2 * nu_new + (mu_new / sigma_new)^2

  dgig(x, p_vec, a_vec, b = 1e-14, log)
}

#' @rdname gal
#' @importFrom stats rgamma
#' @export
rgal <- function (n, delta, mu, nu, sigma, seed = 0)
{
  if(missing(delta)){
    stop('argument "delta" missing, with no default')
  }
  if(missing(mu)){
    stop('argument "mu" missing, with no default')
  }
  if(missing(nu)){
    stop('argument "nu" missing, with no default')
  }
  if(missing(sigma)){
    stop('argument "sigma" missing, with no default')
  }
  if (length(delta) != length(mu)){
    if(length(delta)==1){
      delta <- rep(delta, length(mu))
    } else if(length(mu)==1){
      mu <- rep(mu, length(delta))
    } else{
      stop("delta and mu are vectors of different lenghts")
    }
  }
  if (length(mu) != length(nu)){
    if(length(nu)==1){
      nu <- rep(nu, length(mu))
    } else if(length(mu)==1){
      mu <- rep(mu, length(nu))
      if(length(delta)==1){
        delta <- rep(delta, length(nu))
      }
    } else{
      stop("mu and nu are vectors of different lengths")
    }
  }
  if (length(nu) != length(sigma)){
    if(length(sigma)==1){
      sigma <- rep(sigma, length(nu))
    } else if(length(nu)==1){
      nu <- rep(nu, length(sigma))
      if(length(mu)==1){
        mu <- rep(mu,length(sigma))
      }
      if(length(delta)==1){
        delta <- rep(delta, length(sigma))
      }
    } else{
      stop("nu and sigma are vectors of different lengths")
    }
  }
  if (min(nu) < 0)
    stop("vector nu must be  positive")
  if(n < length(delta)){
    delta_new = delta[1:n]
    mu_new = mu[1:n]
    nu_new = nu[1:n]
    sigma_new = sigma[1:n]
  } else if(n>length(delta)){
    quot <- n%/%length(delta)
    rem <- n%%length(delta)
    delta_new <- rep(delta,quot)
    if(rem > 0){
      delta_new <- c(delta_new, delta[1:rem])
    }
    mu_new <- rep(mu,quot)
    if(rem > 0){
      mu_new <- c(mu_new, mu[1:rem])
    }
    nu_new <- rep(nu,quot)
    if(rem > 0){
      nu_new  <- c(nu_new, nu[1:rem])
    }
    sigma_new <- rep(sigma,quot)
    if(rem > 0){
      sigma_new <- c(sigma_new, sigma[1:rem])
    }
  } else{
    delta_new = delta
    mu_new = mu
    nu_new = nu
    sigma_new = sigma
  }

  V = as.vector(stats::rgamma(n, nu_new, nu_new))
  return(delta_new + mu_new*V + sigma_new^2*sqrt(V)*stats::rnorm(n))
}

#' @rdname gal
#' @export
pgal <- function(q, delta, mu, nu, sigma, lower.tail = TRUE, log.p = FALSE){
  if(missing(delta)){
    stop('argument "delta" missing, with no default')
  }
  if(missing(mu)){
    stop('argument "mu" missing, with no default')
  }
  if(missing(nu)){
    stop('argument "nu" missing, with no default')
  }
  if(missing(sigma)){
    stop('argument "sigma" missing, with no default')
  }
  if (length(delta) != length(mu)){
    if(length(delta)==1){
      delta <- rep(delta, length(mu))
    } else if(length(mu)==1){
      mu <- rep(mu, length(delta))
    } else{
      stop("delta and mu are vectors of different lenghts")
    }
  }
  if (length(mu) != length(nu)){
    if(length(nu)==1){
      nu <- rep(nu, length(mu))
    } else if(length(mu)==1){
      mu <- rep(mu, length(nu))
      if(length(delta)==1){
        delta <- rep(delta, length(nu))
      }
    } else{
      stop("mu and nu are vectors of different lengths")
    }
  }
  if (length(nu) != length(sigma)){
    if(length(sigma)==1){
      sigma <- rep(sigma, length(nu))
    } else if(length(nu)==1){
      nu <- rep(nu, length(sigma))
      if(length(mu)==1){
        mu <- rep(mu,length(sigma))
      }
      if(length(delta)==1){
        delta <- rep(delta, length(sigma))
      }
    } else{
      stop("nu and sigma are vectors of different lengths")
    }
  }
  if (min(nu) < 0)
    stop("vector nu must be  positive")
  n = length(q)
  if(n < length(delta)){
    delta_new = delta[1:n]
    mu_new = mu[1:n]
    nu_new = nu[1:n]
    sigma_new = sigma[1:n]
  } else if(n>length(delta)){
    quot <- n%/%length(delta)
    rem <- n%%length(delta)
    delta_new <- rep(delta,quot)
    if(rem > 0){
      delta_new <- c(delta_new, delta[1:rem])
    }
    mu_new <- rep(mu,quot)
    if(rem > 0){
      mu_new <- c(mu_new, mu[1:rem])
    }
    nu_new <- rep(nu,quot)
    if(rem > 0){
      nu_new  <- c(nu_new, nu[1:rem])
    }
    sigma_new <- rep(sigma,quot)
    if(rem > 0){
      sigma_new <- c(sigma_new, sigma[1:rem])
    }
  } else{
    delta_new = delta
    mu_new = mu
    nu_new = nu
    sigma_new = sigma
  }

  # prob_gal <- sapply(1:n, function(i){
  #   if(q[i]<=0){
  #     p_gal <- ifelse(lower.tail, 0, 1)
  #     p_gal <- ifelse(log.p, log(p_gal), p_gal)
  #     return(p_gal)
  #   } else{
  #     p_gal <- stats::integrate(dgal,lower = 0, upper = q[i],
  #                               delta=delta_new[i],mu=mu_new[i],
  #                               nu=nu_new[i],sigma=sigma_new[i])$value
  #     if(p_gal < 10^{-5} & q[i] > 200){
  #       p_gal = 1
  #     }
  #     p_gal <- ifelse(lower.tail, p_gal, 1-p_gal)
  #     p_gal <- ifelse(log.p, log(p_gal), p_gal)
  #     return(p_gal)
  #   }
  # })
  # return(prob_gal)
  p_vec = nu_new - 0.5
  a_vec = 2 * nu_new + (mu_new / sigma_new)^2

  pgig(q, p_vec, a_vec, b = 1e-14, lower.tail, log.p)
}

#' @rdname gal
#' @export
qgal <- function(p, delta, mu, nu, sigma, lower.tail = TRUE, log.p = FALSE){
  if(missing(delta)){
    stop('argument "delta" missing, with no default')
  }
  if(missing(mu)){
    stop('argument "mu" missing, with no default')
  }
  if(missing(nu)){
    stop('argument "nu" missing, with no default')
  }
  if(missing(sigma)){
    stop('argument "sigma" missing, with no default')
  }
  if (length(delta) != length(mu)){
    if(length(delta)==1){
      delta <- rep(delta, length(mu))
    } else if(length(mu)==1){
      mu <- rep(mu, length(delta))
    } else{
      stop("delta and mu are vectors of different lenghts")
    }
  }
  if (length(mu) != length(nu)){
    if(length(nu)==1){
      nu <- rep(nu, length(mu))
    } else if(length(mu)==1){
      mu <- rep(mu, length(nu))
      if(length(delta)==1){
        delta <- rep(delta, length(nu))
      }
    } else{
      stop("mu and nu are vectors of different lengths")
    }
  }
  if (length(nu) != length(sigma)){
    if(length(sigma)==1){
      sigma <- rep(sigma, length(nu))
    } else if(length(nu)==1){
      nu <- rep(nu, length(sigma))
      if(length(mu)==1){
        mu <- rep(mu,length(sigma))
      }
      if(length(delta)==1){
        delta <- rep(delta, length(sigma))
      }
    } else{
      stop("nu and sigma are vectors of different lengths")
    }
  }
  if (min(nu) < 0)
    stop("vector nu must be  positive")
  n = length(p)
  if(n < length(delta)){
    delta_new = delta[1:n]
    mu_new = mu[1:n]
    nu_new = nu[1:n]
    sigma_new = sigma[1:n]
  } else if(n>length(delta)){
    quot <- n%/%length(delta)
    rem <- n%%length(delta)
    delta_new <- rep(delta,quot)
    if(rem > 0){
      delta_new <- c(delta_new, delta[1:rem])
    }
    mu_new <- rep(mu,quot)
    if(rem > 0){
      mu_new <- c(mu_new, mu[1:rem])
    }
    nu_new <- rep(nu,quot)
    if(rem > 0){
      nu_new  <- c(nu_new, nu[1:rem])
    }
    sigma_new <- rep(sigma,quot)
    if(rem > 0){
      sigma_new <- c(sigma_new, sigma[1:rem])
    }
  } else{
    delta_new = delta
    mu_new = mu
    nu_new = nu
    sigma_new = sigma
  }

  # quant_gal <- sapply(1:n, function(i){
  #   if(nu_new[i]<1){
  #     up_bd <- max(1/nu_new[i]^2, 1000)
  #   } else {
  #     up_bd <- 1000
  #   }
  #   if(p[i] < 0 | p[i] > 1){
  #     warn_qgal <- TRUE
  #     return(NaN)} else{
  #       while(pgal(up_bd, delta=delta_new[i],mu=mu_new[i],
  #                  nu=nu_new[i],sigma=sigma_new[i],
  #                  lower.tail=lower.tail) < p[i]){
  #         up_bd = up_bd + 1000
  #       }
  #       q_gal <- stats::uniroot(function(y) {
  #         pgal(y, delta=delta_new[i],mu=mu_new[i],
  #              nu=nu_new[i],sigma=sigma_new[i],
  #              lower.tail=lower.tail) - p[i]
  #       },lower = 0, upper = up_bd)$root
  #       q_gal <- ifelse(log.p, log(q_gal), q_gal)
  #       return(q_gal)
  #     }
  # }
  # )
  # if(any(is.nan(quant_gal))){
  #   warning("NaNs produced", call. = TRUE, domain = "R")
  # }
  # return(quant_gal)

  p_vec = nu_new - 0.5
  a_vec = 2 * nu_new + (mu_new / sigma_new)^2

  qgig(p, p_vec, a_vec, b = 1e-14, lower.tail, log.p)
}
