#############################################################################################
#' @name gig
#' @title  The Generalised Inverse-Gaussian (GIG) Distribution
#' @aliases dgig rgig pgig qgig
#' @description Density, distribution function, quantile function and
#' random generation for the generalised inverse-Gaussian distribution
#'  with parameters \code{p}, \code{a} and \code{b}.
#' @param x,q vector of quantiles.
#' @param prob vector of probabilities.
#' @param n, number of observations.
#' @param a,b parameters \code{a} and \code{b}. Must be positive.
#' @param p parameter \code{p}.
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \eqn{p} are
#' returned as \eqn{log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X\leq x]},
#' otherwise, \eqn{P[X>x]}.
#' @param seed Seed for the random generation.
#' @return
#' dgig gives the density, pgig gives the distribution function,
#' qgig gives the quantile function, and rgig generates random deviates.
#'
#' Invalid arguments will result in return value NaN, with a warning.
#'
#' The length of the result is determined by \code{n} for rgig.
#' @details The generalised inverse-Gaussian distribution has density given
#' by
#'  \deqn{f(x; p, a, b) = ((a/b)^{p/2})/(2K_p(\sqrt{ab})) x^{p-1} \exp\{-(a/2)x - (b/2)/x\},}
#'  where \eqn{K_p} is modified Bessel function of the second kind of order \eqn{p},
#'  \eqn{x>0}, \eqn{a,b>0} and \eqn{p\in\mathbb{R}}.
#'  See Jørgensen (1982) for further details.
#'
#' @references
#' Jørgensen, Bent (1982). Statistical Properties of the Generalized
#' Inverse Gaussian Distribution. Lecture Notes in Statistics. 9.
#' New York–Berlin: Springer-Verlag. \doi{10.1007/978-1-4612-5698-4}
#'
#' @seealso
#' \code{\link{dnig}}, \code{\link{dig}}, \code{\link{digam}}
#' @examples
#' rgig(100, p = 1, a = 1, b = 1)
#' pgig(0.4, p = 1, a = 1, b = 1)
#' qgig(0.8, p = 1, a = 1, b = 1)
#' plot(function(x){dgig(x, p = 1, a = 1, b = 1)}, main =
#' "Generalised inverse-Gaussian density", ylab = "Probability density",
#' xlim = c(0,10))
#' @rdname gig
#' @export
dgig <- function(x, p, a ,b, log=FALSE){
  if(missing(p)){
    stop('argument "p" missing, with no default')
  }
  if(missing(a)){
    stop('argument "a" missing, with no default')
  }
  if(missing(b)){
    stop('argument "b" missing, with no default')
  }
  if (length(a) != length(p)){
    if(length(a) == 1){
      a <- rep(a, length(p))
    } else if(length(p)==1){
      p <- rep(p, length(a))
    } else {
      stop("a and p are vectors of different lengths")
    }
  }
  if (length(a) != length(b)){
    if(length(a) == 1){
      a <- rep(a, length(b))
      if(length(p)==1){
        p <- rep(p, length(b))
      }
    } else if(length(b) == 1){
      b <- rep(b, length(a))
    } else{
      stop("a and b are vectors of different lengths")      
    }
  }
  if (min(a) < 0)
    stop("vector a must be  positive")
  if (min(b) < 0)
    stop("vector b must be  positive")
  n = length(x)
  if(n < length(p)){
    p_new = p[1:n]
    a_new = a[1:n]
    b_new = b[1:n]
  } else if(n>length(p)){
    quot <- n%/%length(p)
    rem <- n%%length(p)
    p_new <- rep(p,quot)
    if(rem > 0){
      p_new  <- c(p_new, p[1:rem])
    }
    a_new <- rep(a,quot)
    if(rem > 0){
      a_new <- c(a_new,a[1:rem])
    }
    b_new <- rep(b,quot)
    if(rem > 0){
      b_new <- c(b_new, b[1:rem])
    }
  } else{
    p_new <- p
    a_new <- a
    b_new <- b
  }

  densgig <- sapply(1:n, function(i){
    if(x[i] <= 0){
      dens <- ifelse(log, -Inf,0)
      return(dens)
    } else{
      l = 0.5 * p_new[i]*(log(a_new[i]) - log(b_new[i]))
      l = l + (p_new[i]-1) * log(x[i])
      l = l - log(2) - log(besselK(sqrt(a_new[i]*b_new[i]),p_new[i]))
      l = l - 0.5*(a_new[i]*x[i] + b_new[i]/x[i])
      dens <- ifelse(log, l, exp(l))
      return(dens)
    }
  })
  return(densgig)
}

#' @rdname gig
#' @export
rgig <- function (n, p, a, b, seed = 0)
{
  if(missing(p)){
    stop('argument "p" missing, with no default')
  }
  if(missing(a)){
    stop('argument "a" missing, with no default')
  }
  if(missing(b)){
    stop('argument "b" missing, with no default')
  }
  if (length(a) != length(p)){
    if(length(a) == 1){
      a <- rep(a, length(p))
    } else if(length(p)==1){
      p <- rep(p, length(a))
    } else {
      stop("a and p are vectors of different lengths")
    }
  }
  if (length(a) != length(b)){
    if(length(a) == 1){
      a <- rep(a, length(b))
      if(length(p)==1){
        p <- rep(p, length(b))
      }
    } else if(length(b) == 1){
      b <- rep(b, length(a))
    } else{
      stop("a and b are vectors of different lengths")      
    }
  }
  if (min(a) < 0)
    stop("vector a must be  positive")
  if (min(b) < 0)
    stop("vector b must be  positive")
  if(n < length(p)){
    p_new = p[1:n]
    a_new = a[1:n]
    b_new = b[1:n]
  } else if(n>length(p)){
    quot <- n%/%length(p)
    rem <- n%%length(p)
    p_new <- rep(p,quot)
    if(rem > 0){
      p_new  <- c(p_new, p[1:rem])
    }
    a_new <- rep(a,quot)
    if(rem > 0){
      a_new <- c(a_new,a[1:rem])
    }
    b_new <- rep(b,quot)
    if(rem > 0){
      b_new <- c(b_new, b[1:rem])
    }
  } else{
    p_new <- p
    a_new <- a
    b_new <- b
  }

  return(rGIG_cpp(p_new, a_new, b_new, seed))
}

#' @rdname gig
#' @export
pgig <- function(q, p, a, b, lower.tail = TRUE, log.p = FALSE){
  if(missing(p)){
    stop('argument "p" missing, with no default')
  }
  if(missing(a)){
    stop('argument "a" missing, with no default')
  }
  if(missing(b)){
    stop('argument "b" missing, with no default')
  }
  if (length(a) != length(p)){
    if(length(a) == 1){
      a <- rep(a, length(p))
    } else if(length(p)==1){
      p <- rep(p, length(a))
    } else {
      stop("a and p are vectors of different lengths")
    }
  }
  if (length(a) != length(b)){
    if(length(a) == 1){
      a <- rep(a, length(b))
      if(length(p)==1){
        p <- rep(p, length(b))
      }
    } else if(length(b) == 1){
      b <- rep(b, length(a))
    } else{
      stop("a and b are vectors of different lengths")      
    }
  }
  if (min(a) < 0)
    stop("vector a must be  positive")
  if (min(b) < 0)
    stop("vector b must be  positive")
  n = length(q)
  if(n < length(p)){
    p_new = p[1:n]
    a_new = a[1:n]
    b_new = b[1:n]
  } else if(n>length(p)){
    quot <- n%/%length(p)
    rem <- n%%length(p)
    p_new <- rep(p,quot)
    if(rem > 0){
      p_new  <- c(p_new, p[1:rem])
    }
    a_new <- rep(a,quot)
    if(rem > 0){
      a_new <- c(a_new,a[1:rem])
    }
    b_new <- rep(b,quot)
    if(rem > 0){
      b_new <- c(b_new, b[1:rem])
    }
  } else{
    p_new <- p
    a_new <- a
    b_new <- b
  }
  prob_gig <- sapply(1:n, function(i){
    if(q[i]<=0){
      p_gig <- ifelse(lower.tail, 0, 1)
      p_gig <- ifelse(log.p, log(p_gig), p_gig)
      return(p_gig)
    } else{
      p_gig <- stats::integrate(dgig,lower = 0, upper = q[i],
                                p=p_new[i],a=a_new[i],b=b_new[i])$value
      if(p_gig < 10^{-5} & q[i] > sqrt(b_new[i]/a_new[i])*
         besselK(sqrt(a_new[i]*b_new[i]),p_new[i]+1)/
         besselK(sqrt(a_new[i]*b_new[i]),p_new[i])){
        p_gig = 1
      }
      p_gig <- ifelse(lower.tail, p_gig, 1-p_gig)
      p_gig <- ifelse(log.p, log(p_gig), p_gig)
      return(p_gig)
    }
  })
  return(prob_gig)
}

#' @rdname gig
#' @export
qgig <- function(prob, p, a, b, lower.tail = TRUE, log.p = FALSE){
  if(missing(p)){
    stop('argument "p" missing, with no default')
  }
  if(missing(a)){
    stop('argument "a" missing, with no default')
  }
  if(missing(b)){
    stop('argument "b" missing, with no default')
  }
  if (length(a) != length(p)){
    if(length(a) == 1){
      a <- rep(a, length(p))
    } else if(length(p)==1){
      p <- rep(p, length(a))
    } else {
      stop("a and p are vectors of different lengths")
    }
  }
  if (length(a) != length(b)){
    if(length(a) == 1){
      a <- rep(a, length(b))
      if(length(p)==1){
        p <- rep(p, length(b))
      }
    } else if(length(b) == 1){
      b <- rep(b, length(a))
    } else{
      stop("a and b are vectors of different lengths")      
    }
  }
  if (min(a) < 0)
    stop("vector a must be  positive")
  if (min(b) < 0)
    stop("vector b must be  positive")
  n = length(prob)
  if(n < length(p)){
    p_new = p[1:n]
    a_new = a[1:n]
    b_new = b[1:n]
  } else if(n>length(p)){
    quot <- n%/%length(p)
    rem <- n%%length(p)
    p_new <- rep(p,quot)
    if(rem > 0){
      p_new  <- c(p_new, p[1:rem])
    }
    a_new <- rep(a,quot)
    if(rem > 0){
      a_new <- c(a_new,a[1:rem])
    }
    b_new <- rep(b,quot)
    if(rem > 0){
      b_new <- c(b_new, b[1:rem])
    }
  } else{
    p_new <- p
    a_new <- a
    b_new <- b
  }
  quant_gig <- sapply(1:n, function(i){
    if(a_new[i]<1 & b_new[i]<1){
      up_bd <- max(1/a_new[i] * 1/b_new[i], 100)
    } else if (a[i]>=1){
      up_bd <- max(a_new[i]/b_new[i], 100)
    } else {
      up_bd <- 100
    }
    if(prob[i] < 0 | prob[i] > 1){
      warn_qgig <- TRUE
      return(NaN)} else{
        while(pgig(up_bd, p=p_new[i],a=a_new[i],
                b=b_new[i],lower.tail=lower.tail) < prob[i]){
          up_bd = up_bd + 10
        }

        q_gig <- stats::uniroot(function(y) {
          pgig(y, p=p_new[i],a=a_new[i],b=b_new[i],lower.tail=lower.tail) - prob[i]
        },lower = 0, upper = up_bd)$root
        q_gig <- ifelse(log.p, log(q_gig), q_gig)
        return(q_gig)
      }
  }
  )
  if(any(is.nan(quant_gig))){
    warning("NaNs produced", call. = TRUE, domain = "R")
  }
  return(quant_gig)
}

#############################################################################################
#' @name ig
#' @title  The Inverse-Gaussian (IG) Distribution
#' @aliases dig rig pig qig
#' @description Density, distribution function, quantile function and
#' random generation for the inverse-Gaussian distribution
#'  with parameters \code{a} and \code{b}.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n, number of observations.
#' @param a,b parameters \code{a} and \code{b}. Must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \eqn{p} are
#' returned as \eqn{log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X\leq x]},
#' otherwise, \eqn{P[X>x]}.
#' @param seed Seed for the random generation.
#' @return
#' dig gives the density, pig gives the distribution function,
#' qig gives the quantile function, and rig generates random deviates.
#'
#' Invalid arguments will result in return value NaN, with a warning.
#'
#' The length of the result is determined by \code{n} for rig.
#' @details The inverse-Gaussian distribution has density given
#' by
#'  \deqn{f(x; a, b) = \frac{\sqrt{b}}{\sqrt{2\pi x^3}}\exp(
#'  -\frac{a}{2}x -\frac{b}{2x} + \sqrt{ab}),}
#'  where \eqn{x>0} and \eqn{a,b>0}. In this parameterization,
#'  \eqn{E(X) = \sqrt{b}/\sqrt{a}}. See Tweedie (1957a, 1957b) for
#'  further details.
#'
#' @references
#' Tweedie, M. C. K. (1957a). "Statistical Properties of Inverse Gaussian Distributions I". Annals of Mathematical Statistics. 28 (2): 362–377. \doi{10.1214/aoms/1177706964}
#' 
#' Tweedie, M. C. K. (1957b). "Statistical Properties of Inverse Gaussian Distributions II". Annals of Mathematical Statistics. 28 (3): 696–705. \doi{10.1214/aoms/1177706881}
#'
#' @seealso
#' \code{\link{dnig}}, \code{\link{dgig}}, \code{\link{digam}}
#' @examples
#' rig(100, a = 1, b = 1)
#' pig(0.4, a = 1, b = 1)
#' qig(0.8, a = 1, b = 1)
#' plot(function(x){dig(x, a = 1, b = 1)}, main =
#' "Inverse-Gaussian density", ylab = "Probability density",
#' xlim = c(0,10))
#' @rdname ig
#' @export
dig <- function(x, a ,b, log=FALSE){
  if(missing(a)){
    stop('argument "a" missing, with no default')
  }
  if(missing(b)){
    stop('argument "b" missing, with no default')
  }
  if (length(a) != length(b)){
    if(length(a)==1){
      a <- rep(a, length(b))
    } else if(length(b)==1){
      b <- rep(b, length(a))
    } else{
      stop("a and b are vectors of different lengths")
    }
  }
  if (min(a) < 0)
    stop("vector a must be  positive")
  if (min(b) < 0)
    stop("vector b must be  positive")
  n = length(x)
  if(n < length(a)){
    a_new = a[1:n]
    b_new = b[1:n]
  } else if(n>length(a)){
    quot <- n%/%length(a)
    rem <- n%%length(a)
    a_new <- rep(a,quot)
    if(rem > 0){
      a_new <- c(a_new,a[1:rem])
    }
    b_new <- rep(b,quot)
    if(rem > 0){
      b_new <- c(b_new, b[1:rem])
    }
  } else{
    a_new <- a
    b_new <- b
  }

  densig <- sapply(1:n, function(i){
    if(x[i] <= 0){
      dens <- ifelse(log, -Inf,0)
      return(dens)
    } else{
      l = (log(b_new[i]) - log(2*pi) - 3 *log(x[i]) )
      l = l - a_new[i]*x[i] - b_new[i]/x[i] + 2*sqrt(a_new[i]*b_new[i])
      l = 0.5 * l
      dens <- ifelse(log, l, exp(l))
      return(dens)
    }
  })
  return(densig)
}

#' @rdname ig
#' @export
rig <- function (n, a, b, seed = 0)
{
  if(missing(a)){
    stop('argument "a" missing, with no default')
  }
  if(missing(b)){
    stop('argument "b" missing, with no default')
  }
  if (length(a) != length(b)){
    if(length(a)==1){
      a <- rep(a, length(b))
    } else if(length(b)==1){
      b <- rep(b, length(a))
    } else{
      stop("a and b are vectors of different lengths")
    }
  }
  if (min(a) < 0)
    stop("vector a must be  positive")
  if (min(b) < 0)
    stop("vector b must be  positive")
  if(n < length(a)){
    a_new = a[1:n]
    b_new = b[1:n]
  } else if(n>length(a)){
    quot <- n%/%length(a)
    rem <- n%%length(a)
    a_new <- rep(a,quot)
    if(rem > 0){
      a_new <- c(a_new,a[1:rem])
    }
    b_new <- rep(b,quot)
    if(rem > 0){
      b_new <- c(b_new, b[1:rem])
    }
  } else{
    a_new <- a
    b_new <- b
  }
  p_new <- rep(-0.5, length(a_new))
  return(rGIG_cpp(p_new, a_new, b_new, seed))
}

#' @rdname ig
#' @export
pig <- function(q, a, b, lower.tail = TRUE, log.p = FALSE){
  if(missing(a)){
    stop('argument "a" missing, with no default')
  }
  if(missing(b)){
    stop('argument "b" missing, with no default')
  }
      mu <- sqrt(b/a)
      p_ig <- stats::pnorm(sqrt(b/q)*(q/mu-1)) + exp(2*b/mu)*stats::pnorm(-sqrt(b/q)*(q/mu+1))
      if(!lower.tail){
        p_ig <- 1 - p_ig
      }

      if(log.p){
          p_ig <- log(p_ig)
      }
      return(p_ig)
}

#' @rdname ig
#' @export
qig <- function(p, a, b, lower.tail = TRUE, log.p = FALSE){
  if(missing(a)){
    stop('argument "a" missing, with no default')
  }
  if(missing(b)){
    stop('argument "b" missing, with no default')
  }
  if (length(a) != length(b)){
    if(length(a)==1){
      a <- rep(a, length(b))
    } else if(length(b)==1){
      b <- rep(b, length(a))
    } else{
      stop("a and b are vectors of different lengths")
    }
  }
  if (min(a) < 0)
    stop("vector a must be  positive")
  if (min(b) < 0)
    stop("vector b must be  positive")
  n = length(p)
  if(n < length(a)){
    a_new = a[1:n]
    b_new = b[1:n]
  } else if(n>length(a)){
    quot <- n%/%length(a)
    rem <- n%%length(a)
    a_new <- rep(a,quot)
    if(rem > 0){
      a_new <- c(a_new,a[1:rem])
    }
    b_new <- rep(b,quot)
    if(rem > 0){
      b_new <- c(b_new, b[1:rem])
    }
  } else{
    a_new <- a
    b_new <- b
  }
  quant_ig <- sapply(1:n, function(i){
    if(a_new[i]<1 & b_new[i]<1){
      up_bd <- max(1/a_new[i] * 1/b_new[i], 100)
    } else if (a[i]>=1){
      up_bd <- max(a_new[i]/b_new[i], 100)
    } else {
      up_bd <- 100
    }
    if(p[i] < 0 | p[i] > 1){
      warn_qig <- TRUE
      return(NaN)} else{
        q_ig <- stats::uniroot(function(y) {
          pig(y, a=a_new[i],b=b_new[i],lower.tail=lower.tail) - p[i]
        },lower = 0, upper = up_bd)$root
        q_ig <- ifelse(log.p, log(q_ig), q_ig)
        return(q_ig)
      }
  }
  )
  if(any(is.nan(quant_ig))){
    warning("NaNs produced", call. = TRUE, domain = "R")
  }
  return(quant_ig)
}

#############################################################################################
#' @name igam
#' @title  The Inverse-Gamma (IGam) Distribution
#' @aliases digam rigam pigam qigam
#' @description Density, distribution function, quantile function and
#' random generation for the inverse-Gamma distribution
#'  with parameters \code{a} and \code{b}.
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n, number of observations.
#' @param a,b parameters \code{a} and \code{b}. Must be positive.
#' @param log,log.p logical; if \code{TRUE}, probabilities/densities \eqn{p} are
#' returned as \eqn{log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X\leq x]},
#' otherwise, \eqn{P[X>x]}.
#' @return
#' digam gives the density, pigam gives the distribution function,
#' qigam gives the quantile function, and rigam generates random deviates.
#'
#' Invalid arguments will result in return value NaN, with a warning.
#'
#' The length of the result is determined by \code{n} for rig.
#' @details The inverse-Gamma distribution has density given
#' by
#'  \deqn{f(x; a, b) = \frac{b^a}{\Gamma(a)}x^{a-1}\exp(
#'  -\frac{b}{x}),}
#'  where \eqn{x>0} and \eqn{a,b>0}.
#'
#' @seealso
#' \code{\link{dnig}}, \code{\link{dgig}}
#' @examples
#' rigam(100, a = 1, b = 1)
#' pigam(0.4, a = 1, b = 1)
#' qigam(0.8, a = 1, b = 1)
#' plot(function(x){digam(x, a = 1, b = 1)}, main =
#' "Inverse-Gamma density", ylab = "Probability density",
#' xlim = c(0,10))
#' @rdname igam
#' @export
digam <- function(x, a ,b, log=FALSE){
  if(missing(a)){
    stop('argument "a" missing, with no default')
  }
  if(missing(b)){
    stop('argument "b" missing, with no default')
  }
  if (length(a) != length(b)){
    if(length(a)==1){
      a <- rep(a, length(b))
    } else if(length(b)==1){
      b <- rep(b, length(a))
    } else{
      stop("a and b are vectors of different lengths")
    }
  }
  if (min(a) < 0)
    stop("vector a must be  positive")
  if (min(b) < 0)
    stop("vector b must be  positive")
  n = length(x)
  if(n < length(a)){
    a_new = a[1:n]
    b_new = b[1:n]
  } else if(n>length(a)){
    quot <- n%/%length(a)
    rem <- n%%length(a)
    a_new <- rep(a,quot)
    if(rem > 0){
      a_new <- c(a_new,a[1:rem])
    }
    b_new <- rep(b,quot)
    if(rem > 0){
      b_new <- c(b_new, b[1:rem])
    }
  } else{
    a_new <- a
    b_new <- b
  }

  densigam <- sapply(1:n, function(i){
    if(x[i] <= 0){
      dens <- ifelse(log, -Inf,0)
      return(dens)
    } else{
      l = a_new[i]*log(b_new[i]) - lgamma(a_new[i]) -
        (a_new[i]+1)*log(x[i]) - b_new[i]/x[i]
      dens <- ifelse(log, l, exp(l))
      return(dens)
    }
  })
  return(densigam)
}

#' @rdname igam
#' @export
rigam <- function (n, a, b)
{
  1/stats::rgamma(n,shape=a, rate = b)
}

#' @rdname igam
#' @export
pigam <- function(q, a, b, lower.tail = TRUE, log.p = FALSE){
  if(missing(a)){
    stop('argument "a" missing, with no default')
  }
  if(missing(b)){
    stop('argument "b" missing, with no default')
  }
  p_igam <- 1 - stats::pgamma(1/q, shape = a, rate = b)
  return(p_igam)
}

#' @rdname igam
#' @export
qigam <- function(p, a, b, lower.tail = TRUE, log.p = FALSE){
 1/stats::qgamma(1-p,shape = a, rate = b, lower.tail = lower.tail, log.p = log.p)
}

#############################################################################################
#' @name nig
#' @title  The Normal Inverse-Gaussian (GIG) Distribution
#' @aliases dnig rnig pnig qnig
#' @description Density, distribution function, quantile function and
#' random generation for the normal inverse-Gaussian distribution
#'  with parameters \code{p}, \code{a} and \code{b}.
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
#' dnig gives the density, pnig gives the distribution function,
#' qnig gives the quantile function, and rnig generates random deviates.
#'
#' Invalid arguments will result in return value NaN, with a warning.
#'
#' The length of the result is determined by \code{n} for rnig.
#' @details The normal inverse-Gaussian distribution has density given
#' by
#'  \deqn{f(x; p, a, b) =
#'  \frac{e^{\nu+\mu(x-\delta)/\sigma^2}\sqrt{\nu\mu^2/\sigma^2+\nu^2}}{\pi\sqrt{\nu\sigma^2+(x-\delta)^2}}
#'  K_1(\sqrt{(\nu\sigma^2+(x-\delta)^2)(\mu^2/\sigma^4+\nu/\sigma^2)}),}
#'  where \eqn{K_p} is modified Bessel function of the second kind of order \eqn{p},
#'  \eqn{x>0}, \eqn{\nu>0} and \eqn{\mu,\delta, \sigma\in\mathbb{R}}.
#'  See Barndorff-Nielsen (1977, 1978 and 1997) for further details.
#'
#' @references
#'  Barndorff-Nielsen, O. (1977) Exponentially decreasing distributions for the logarithm of particle size. Proceedings of the Royal Society of London.
#'  
#'  Series A, Mathematical and Physical Sciences. The Royal Society. 353, 401–409. \doi{10.1098/rspa.1977.0041}
#'  
#'  Barndorff-Nielsen, O. (1978) Hyperbolic Distributions and Distributions on Hyperbolae, Scandinavian Journal of Statistics. 5, 151–157.
#'  
#'  Barndorff-Nielsen, O. (1997) Normal Inverse Gaussian Distributions and Stochastic Volatility Modelling, Scandinavian Journal of Statistics. 24, 1-13. \doi{10.1111/1467-9469.00045}
#'
#' @seealso
#' \code{\link{dgig}}, \code{\link{dig}}, \code{\link{digam}}
#' @examples
#' rnig(100, delta = 0, mu = 5, sigma = 1, nu = 1)
#' pnig(0.4, delta = 0, mu = 5, sigma = 1, nu = 1)
#' qnig(0.8, delta = 0, mu = 5, sigma = 1, nu = 1)
#' plot(function(x){dnig(x, delta = 0, mu = 5, sigma = 1, nu = 1)}, main =
#' "Normal inverse-Gaussian density", ylab = "Probability density",
#' xlim = c(0,10))
#' @rdname nig
#' @export
dnig <- function(x, delta, mu, nu, sigma, log=FALSE){
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

  densnig <- sapply(1:n, function(i){
    c0 <- sqrt(nu_new[i]) * sqrt( mu_new[i]^2/sigma_new[i]^2 + nu_new[i]) / pi
    l <- nu_new[i] + mu_new[i] * (x[i] - delta_new[i]) /sigma_new[i]^2
    coeff <- sqrt( nu_new[i] * sigma_new[i]^2 + (x[i] - delta_new[i])^2)
    l <- l + log(c0) -  log(coeff)
    l <- l + log(besselK(coeff  * sqrt(mu_new[i]^2/sigma_new[i]^4
                                        + nu_new[i]/sigma_new[i]^2),
                          -1,TRUE)) -coeff  * sqrt(mu_new[i]^2/sigma_new[i]^4 +
                                                    nu_new[i]/sigma_new[i]^2)
    dens <- ifelse(log, l, exp(l))
    return(dens)
  })
  return(densnig)
}

#' @rdname nig
#' @export
rnig <- function (n, delta, mu, nu, sigma, seed = 0)
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

  V = as.vector(rig(n, nu, nu, seed = seed))
  return(delta_new + mu_new*V + sigma_new^2*sqrt(V)*stats::rnorm(n))
}

#' @rdname nig
#' @export
pnig <- function(q, delta, mu, nu, sigma, lower.tail = TRUE, log.p = FALSE){
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

  prob_nig <- sapply(1:n, function(i){
    if(q[i]<=0){
      p_nig <- ifelse(lower.tail, 0, 1)
      p_nig <- ifelse(log.p, log(p_nig), p_nig)
      return(p_nig)
    } else{
      p_nig <- stats::integrate(dnig,lower = 0, upper = q[i],
                                delta=delta_new[i],mu=mu_new[i],
                                nu=nu_new[i],sigma=sigma_new[i])$value
      if(p_nig < 10^{-5} & q[i] > 200){
        p_nig = 1
      }
      p_nig <- ifelse(lower.tail, p_nig, 1-p_nig)
      p_nig <- ifelse(log.p, log(p_nig), p_nig)
      return(p_nig)
    }
  })
  return(prob_nig)
}

#' @rdname nig
#' @export
qnig <- function(p, delta, mu, nu, sigma, lower.tail = TRUE, log.p = FALSE){
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
  quant_nig <- sapply(1:n, function(i){
    if(nu_new[i]<1){
      up_bd <- max(1/nu_new[i]^2, 1000)
    } else {
      up_bd <- 1000
    }
    if(p[i] < 0 | p[i] > 1){
      warn_qnig <- TRUE
      return(NaN)} else{
        while(pnig(up_bd, delta=delta_new[i],mu=mu_new[i],
                   nu=nu_new[i],sigma=sigma_new[i],
                   lower.tail=lower.tail) < p[i]){
          up_bd = up_bd + 1000
        }
        q_nig <- stats::uniroot(function(y) {
          pnig(y, delta=delta_new[i],mu=mu_new[i],
               nu=nu_new[i],sigma=sigma_new[i],
               lower.tail=lower.tail) - p[i]
        },lower = 0, upper = up_bd)$root
        q_nig <- ifelse(log.p, log(q_nig), q_nig)
        return(q_nig)
      }
  }
  )
  if(any(is.nan(quant_nig))){
    warning("NaNs produced", call. = TRUE, domain = "R")
  }
  return(quant_nig)
}
