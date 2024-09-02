# This file contains ngme noise specifications

#' @title ngme noise specification
#' @aliases noise_nig noise_normal
#' @description Function for specifying ngme noise.
#' Please use \code{noise_nig} and \code{noise_normal} for simpler usage.
#' Use \code{ngme_noise_types()} to check all the available types.
#'
#' @details The parameterization is given in \code{?nig} and \code{?gal}. Moreover,
#' for specifying non-stationary mu and sigma, nu
#' \deqn{\mu = B_{\mu} \theta_{\mu},} and
#' \deqn{\sigma = \exp (B_{\sigma} \theta_{\sigma}),}
#' \deqn{\nu = \exp (B_{\nu} \theta_{\nu}).}
#'
#' @param noise_type    type of noise, "nig", "normal"
#' @param mu          specify the NIG noise parameter mu, see \code{?nig}
#' @param sigma       specify the noise parameter sigma, see \code{?nig}
#' @param V             start value for V
#' @param theta_mu      specify a non-stationary noise using theta_mu
#' @param B_mu          Basis matrix for mu (if non-stationary)
#' @param theta_sigma   specify a non-stationary noise using theta_sigma
#' @param B_sigma       Basis matrix for sigma (if non-stationary)
#' @param theta_nu      specify a non-stationary noise using theta_nu
#' @param nu_lower_bound specify the lower bound of parameter nu
#' @param B_nu          Basis matrix for nu (if non-stationary)
#' @param fix_theta_mu     fix the parameter of theta_mu
#' @param fix_theta_sigma  fix the parameter of theta_sigma
#' @param fix_theta_nu     fix the parameter of nu
#' @param fix_theta_sigma_normal  fix the parameter of sigma_normal, used in noise_normal_nig()
#' @param fix_rho    fix the parameter of rho
#' @param fix_V         fix the sampling of V
#' @param shared_sigma allow only in measurement noise,
#'  gives Y|W ~ N(mean * sigma, sigma^2)
#' @param theta_sigma_normal for normal nosie with nig noise sharing same parameter
#' @param B_sigma_normal    for normal nosie with nig noise sharing same parameter
#' @param sigma_normal  for normal nosie with nig noise sharing same parameter
#' @param theta_sigma_nig similar to theta_sigma_normal
#' @param B_sigma_nig     similar to B_sigma_nig
#' @param sigma_nig     similar to sigma_normal
#' @param single_V  TRUE if V is a single number
#' @param share_V  used only for bivariate model
#' @param corr_measurement TRUE if we use correlated measurement noise
#' @param index_corr used when corr_measurement=TRUE, indicate which observation has correlation
#' @param map_corr 1d, 2d, or formula, used when corr_measurement=TRUE, specify use which covariate to infer the index_corr.
#' @param rho used when corr_measurement=TRUE, starting point for correlation
#' @param prior_mu prior distribution for parameter of mu
#' @param prior_sigma prior distribution for parameter of sigma
#' @param prior_nu prior distribution for parameter of nu
#' @param ...       additional arguments
#'
#' @return a list of specification of noise
ngme_noise <- function(
  noise_type,
  mu              = 0,
  sigma           = 1,
  nu              = 1,
  B_mu            = NULL,
  theta_mu        = NULL,
  B_sigma         = NULL,
  theta_sigma     = NULL,
  B_nu            = NULL,
  theta_nu        = NULL,
  theta_sigma_normal = NULL,
  B_sigma_normal  = NULL,
  fix_theta_mu     = FALSE,
  fix_theta_sigma  = FALSE,
  fix_rho          = FALSE,
  fix_theta_sigma_normal = FALSE,
  fix_theta_nu     = FALSE,
  V               = NULL,
  fix_V            = FALSE,
  single_V        = FALSE,
  share_V         = FALSE,
  corr_measurement = FALSE,
  index_corr      = NULL,
  map_corr        = NULL,
  nu_lower_bound  = 0.01,
  shared_sigma    = FALSE,
  rho             = double(0),
  prior_mu        = ngme_prior("normal", param=c(0, 0.01)),
  prior_sigma     = ngme_prior("normal", param=c(0, 0.01)),
  prior_nu        = ngme_prior("normal", param=c(0, 0.01)),
  ...
) {
  if (is.null(theta_mu)) theta_mu <- mu
  if (is.null(theta_sigma))
    theta_sigma <- if (sigma>0) log(sigma) else stop("ngme_noise: sigma should be positive.")
  if (is.null(theta_nu))
    theta_nu <- if (nu>0) log(nu) else stop("ngme_noise: nu should be positive.")

  stopifnot("Unkown noise type. Please check ngme_noise_types()" =
    noise_type %in% ngme_noise_types())

  if (is.null(B_mu)) B_mu <- as.matrix(1)
  if (is.null(B_sigma)) B_sigma <- as.matrix(1)
  if (is.null(B_nu)) B_nu <- as.matrix(1)

  stopifnot(
    "Please input B_mu as a matrix." = is.matrix(B_mu),
    "Please input B_sigma as a matrix." = is.matrix(B_sigma),
    "Please make sure ncol(B_mu) == length(theta_mu)." = ncol(B_mu) == length(theta_mu),
    "Please make sure ncol(B_sigma) == length(theta_sigma)." = ncol(B_sigma) == length(theta_sigma),
    "Please make sure ncol(B_nu) == length(theta_nu)." = ncol(B_nu) == length(theta_nu),
    "prior_mu is not specified properly, please use ngme_prior(..)"
      = class(prior_mu) == "ngme_prior",
    "prior_sigma is not specified properly, please use ngme_prior(..)"
      = class(prior_sigma) == "ngme_prior",
    "prior_nu is not specified properly, please use ngme_prior(..)"
      = class(prior_nu) == "ngme_prior"
  )

  if (all(noise_type == "normal")) {
    theta_mu <- double(0)
    B_mu <- matrix(ncol=0, nrow=nrow(B_mu))
    theta_nu <- double(0)
    B_nu <- matrix(ncol=0, nrow=nrow(B_nu))
  }

  if (all(noise_type != "normal_nig")) {
    theta_sigma_normal <- double(0)
    B_sigma_normal <- matrix(ncol=0, nrow=nrow(B_mu))
  }

  # init rho
  if (corr_measurement && length(rho) == 0) {
    rho <- 0
  }

  n_theta_mu    <- if (fix_theta_mu) 0 else length(theta_mu)
  n_theta_sigma <- if (fix_theta_sigma) 0 else length(theta_sigma)
  n_theta_nu    <- if (fix_theta_nu) 0 else length(theta_nu)
  n_rho         <- if (fix_rho) 0 else length(rho)
  n_theta_sigma_normal <- if (fix_theta_sigma_normal) 0 else length(theta_sigma_normal)

  structure(
    list(
      noise_type      = noise_type,
      V               = V,
      theta_mu        = theta_mu,
      theta_sigma     = theta_sigma,
      theta_sigma_normal = theta_sigma_normal,
      theta_nu        = theta_nu,
      B_mu            = B_mu,
      B_sigma         = B_sigma,
      B_sigma_normal  = B_sigma_normal,
      B_nu            = B_nu,

      n_theta_mu      = n_theta_mu,
      n_theta_sigma   = n_theta_sigma,
      n_theta_nu      = n_theta_nu,
      n_rho           = n_rho,
      n_theta_sigma_normal = n_theta_sigma_normal,

      fix_theta_mu    = fix_theta_mu,
      fix_theta_sigma = fix_theta_sigma,
      fix_theta_nu    = fix_theta_nu,
      nu_lower_bound = nu_lower_bound,
      fix_V           = fix_V,
      fix_rho         = fix_rho,
      fix_theta_sigma_normal = fix_theta_sigma_normal,

      n_params        = n_theta_mu + n_theta_sigma + n_theta_nu + n_rho + n_theta_sigma_normal,
      
      single_V        = single_V,
      share_V         = share_V,
      corr_measurement = corr_measurement,
      index_corr      = index_corr,
      map_corr        = map_corr,
      rho             = rho,
      shared_sigma    = shared_sigma,
      prior_mu        = prior_mu,
      prior_sigma     = prior_sigma,
      prior_nu        = prior_nu,
      ...
    ),
    class = "ngme_noise"
  )
}

#' @rdname ngme_noise
#' @export
#' @examples
#' noise_normal(sigma = 2)
noise_normal <- normal <- function(
  sigma             = NULL,
  theta_sigma       = NULL,
  B_sigma           = matrix(1),
  corr_measurement  = FALSE,
  index_corr        = NULL,
  ...
) {
  sd <- sigma

  if (!is.null(sd) && !is.null(theta_sigma))
    stop("Please only use sigma or theta_sigma as input")

  # both are null, use default value
  if (is.null(sd) && is.null(theta_sigma)) {
    theta_sigma <- 0
  }

  if (!is.null(sd)) {
    stopifnot(
      "sd is a double" = is.double(sd),
      "sd should be positive" = sd > 0
    )

    theta_sigma <- log(sd)
  }

  stopifnot("Make sure ncol of B_sigma = length of theta_signa"
    = ncol(B_sigma) == length(theta_sigma))

  ngme_noise(
    noise_type = "normal",
    theta_sigma = theta_sigma,
    B_sigma = B_sigma,
    corr_measurement = corr_measurement,
    index_corr      = index_corr,
    ...
  )
}

#' @rdname ngme_noise
#' @export
#' @examples
#' noise_nig(mu = 1, sigma = 2, nu = 1)
noise_nig <- nig <- function(
  mu            = NULL,
  sigma         = NULL,
  nu            = NULL,
  V             = NULL,
  theta_mu      = NULL,
  theta_sigma   = NULL,
  theta_nu      = NULL,
  nu_lower_bound = 0.01,
  B_mu          = matrix(1),
  B_sigma       = matrix(1),
  B_nu          = matrix(1),
  corr_measurement = FALSE,
  index_corr      = NULL,
  ...
) {
  # if nothing, then fill with default
  stopifnot("Please use theta_mu for non-stationary mu." = length(mu) < 2)
  if (is.null(mu) && is.null(theta_mu)) theta_mu <- 0
  if (is.null(sigma) && is.null(theta_sigma)) theta_sigma <- 0
  if (is.null(nu) && is.null(theta_nu)) theta_nu <- 0

  if (!is.null(nu) && nu <= 0) stop("ngme_nosie: nu should be positive.")
  if (!is.null(sigma) && sigma <= 0) stop("ngme_nosie: sigma should be positive.")

  if (!is.null(mu))     theta_mu <- mu
  if (!is.null(sigma))  theta_sigma <- log(sigma)
  if (!is.null(nu))     theta_nu <- log(nu)

  ngme_noise(
    noise_type = "nig",
    theta_mu = theta_mu,
    theta_sigma = theta_sigma,
    theta_nu = theta_nu,
    nu_lower_bound = nu_lower_bound,
    V = V,
    B_mu = B_mu,
    B_sigma = B_sigma,
    B_nu = B_nu,
    corr_measurement = corr_measurement,
    index_corr      = index_corr,
    ...
  )
}

#' @rdname ngme_noise
#' @export
#' @examples
#' noise_gal(mu = 1, sigma = 2, nu = 1)
noise_gal <- gal <- function(
  mu            = NULL,
  sigma         = NULL,
  nu            = NULL,
  V             = NULL,
  theta_mu      = NULL,
  theta_sigma   = NULL,
  theta_nu      = NULL,
  nu_lower_bound = 0.01,
  B_mu          = matrix(1),
  B_sigma       = matrix(1),
  B_nu          = matrix(1),
  corr_measurement = FALSE,
  index_corr      = NULL,
  ...
) {
  # if nothing, then fill with default
  stopifnot("Please use theta_mu for non-stationary mu." = length(mu) < 2)
  if (is.null(mu) && is.null(theta_mu)) theta_mu <- 0
  if (is.null(sigma) && is.null(theta_sigma)) theta_sigma <- 0
  if (is.null(nu) && is.null(theta_nu)) theta_nu <- 0

  if (!is.null(nu) && nu <= 0) stop("ngme_nosie: nu should be positive.")
  if (!is.null(sigma) && sigma <= 0) stop("ngme_nosie: sigma should be positive.")

  if (!is.null(mu))     theta_mu <- mu
  if (!is.null(sigma))  theta_sigma <- log(sigma)
  if (!is.null(nu))     theta_nu <- log(nu)

  ngme_noise(
    noise_type = "gal",
    theta_mu = theta_mu,
    theta_sigma = theta_sigma,
    theta_nu = theta_nu,
    nu_lower_bound = nu_lower_bound,
    V = V,
    B_mu = B_mu,
    B_sigma = B_sigma,
    B_nu = B_nu,
    corr_measurement = corr_measurement,
    index_corr      = index_corr,
    ...
  )
}

# update noise
update_noise <- function(noise, n=NULL, new_noise=NULL) {
  # update with length n
  if (!is.null(n)) {
    stopifnot("n should be integer" = is.numeric(n))
    B_mu <- noise$B_mu
stopifnot("n / nrow(B_mu) not integer" = abs(n/nrow(B_mu) - round(n/nrow(B_mu))) < 1e-4)
    noise$B_mu <- matrix(data = rep(B_mu, n / nrow(B_mu)), nrow = n)

    B_sigma <- noise$B_sigma
stopifnot("n / nrow(B_sigma) not integer" = abs(n/nrow(B_sigma) - round(n/nrow(B_sigma))) < 1e-4)
    noise$B_sigma <- matrix(data = rep(B_sigma, n / nrow(B_sigma)), nrow = n)

    B_nu <- noise$B_nu
stopifnot("n / nrow(B_nu) not integer" = abs(n/nrow(B_nu) - round(n/nrow(B_nu))) < 1e-4)
    noise$B_nu <- matrix(data = rep(B_nu, n / nrow(B_nu)), nrow = n)

    if (noise$noise_type == "normal_nig") {
      B_sigma_normal <- noise$B_sigma_normal
      noise$B_sigma_normal <- matrix(data = rep(B_sigma_normal, n / nrow(B_sigma_normal)), nrow = n)
    }
    # noise <- do.call(ngme_noise, noise)
  } else if (!is.null(new_noise)) {
    # update noise after estimation
    if (all(new_noise$noise_type != "normal")) {
      noise$theta_mu  <- new_noise$theta_mu
      noise$theta_nu  <- new_noise$theta_nu
    }
    noise$theta_sigma        <- new_noise$theta_sigma
    noise$rho                <- new_noise$rho
    if (!is.null(new_noise$V)) noise$V <- new_noise$V

    # bv noise
    if (length(noise$noise_type) == 2) {
      # pass mu, sigma, nu to sub_models
      n_theta_mu1    <- noise$bv_noises[[1]]$n_theta_mu
      n_theta_mu2    <- noise$bv_noises[[2]]$n_theta_mu
      n_theta_sigma1 <- noise$bv_noises[[1]]$n_theta_sigma
      n_theta_sigma2 <- noise$bv_noises[[2]]$n_theta_sigma
      n_theta_nu1    <- noise$bv_noises[[1]]$n_theta_nu
      n_theta_nu2    <- noise$bv_noises[[2]]$n_theta_nu
      noise$bv_noises[[1]]$theta_mu    <- head(noise$theta_mu, n_theta_mu1)
      noise$bv_noises[[2]]$theta_mu    <- tail(noise$theta_mu, n_theta_mu2)
      noise$bv_noises[[1]]$theta_sigma <- head(noise$theta_sigma, n_theta_sigma1)
      noise$bv_noises[[2]]$theta_sigma <- tail(noise$theta_sigma, n_theta_sigma2)
      noise$bv_noises[[1]]$theta_nu    <- head(noise$theta_nu, n_theta_nu1)
      noise$bv_noises[[2]]$theta_nu    <- tail(noise$theta_nu, n_theta_nu2)
    } else if (noise$noise_type == "normal_nig") {
      noise$theta_sigma_normal <- new_noise$theta_sigma_normal
    }
  }
  noise
}

#' @rdname ngme_noise
#' @export
noise_normal_nig <- normal_nig <- function(
  sigma_normal  = NULL,
  mu            = NULL,
  sigma_nig     = NULL,
  nu            = NULL,
  V             = NULL,
  theta_mu      = NULL,
  theta_sigma_nig   = NULL,
  theta_sigma_normal   = NULL,
  theta_nu      = NULL,
  B_mu          = matrix(1),
  B_sigma_nig   = matrix(1),
  B_sigma_normal = matrix(1),
  B_nu          = matrix(1),
  corr_measurement = FALSE,
  index_corr      = NULL,
  ...
) {
  # if nothing, then fill with default
  stopifnot("Please use theta_mu for non-stationary mu." = length(mu) < 2)
  if (is.null(mu) && is.null(theta_mu)) theta_mu <- 0
  if (is.null(sigma_nig) && is.null(theta_sigma_nig)) theta_sigma_nig <- 0
  if (is.null(nu) && is.null(theta_nu)) theta_nu <- 0
  if (is.null(sigma_normal) && is.null(theta_sigma_normal)) theta_sigma_normal <- 0

  if (!is.null(nu) && nu <= 0) stop("ngme_nosie: nu should be positive.")
  if (!is.null(sigma_nig) && sigma_nig <= 0) stop("ngme_nosie: sigma_nig should be positive.")
  if (!is.null(sigma_normal) && sigma_normal <= 0) stop("ngme_nosie: sigma_nig should be positive.")

  if (!is.null(mu))     theta_mu <- mu
  if (!is.null(sigma_nig))  theta_sigma_nig <- log(sigma_nig)
  if (!is.null(sigma_normal))  theta_sigma_normal <- log(sigma_normal)
  if (!is.null(nu))     theta_nu <- log(nu)

  ngme_noise(
    noise_type = "normal_nig",
    theta_mu = theta_mu,
    theta_sigma = theta_sigma_nig,
    theta_sigma_normal = theta_sigma_normal,
    theta_nu = theta_nu,
    V = V,
    B_mu = B_mu,
    B_sigma = B_sigma_nig,
    B_sigma_normal = B_sigma_normal,
    B_nu = B_nu,
    corr_measurement = corr_measurement,
    index_corr      = index_corr,
    ...
  )
}

#' Print ngme noise
#'
#' @param x noise object
#' @param padding number of white space padding in front
#' @param prefix prefix
#' @param suppress_sigma suppress printing sigma
#' @param ... ...
#'
#' @return a list (noise specifications)
#' @export
print.ngme_noise <- function(x, padding = 0, prefix = "Noise type", suppress_sigma=FALSE, ...) {
  noise <- x
  pad_space <- paste(rep(" ", padding), collapse = "")
  pad_add4_space <- paste(rep(" ", padding + 4), collapse = "")

  if (is.null(noise)) {
    cat(pad_space); cat(prefix); cat(": "); cat("NULL"); cat("\n")
  } else {
    if (length(noise$noise_type) == 2) {
      # bivariate noise
      cat(pad_space); cat(" ");
      if (noise$single_V && noise$share_V) {
        cat("Bivariate type-G1 noise (single_V && share_V):");
      } else if (noise$single_V && !noise$share_V) {
        cat("Bivariate type-G2 noise (single_V):");
      } else if (!noise$single_V && noise$share_V) {
        cat("Bivariate type-G3 noise (share_V):");
      } else {
        cat("Bivariate type-G4 noise:");
      }
      cat("\n")
      names <- names(noise$bv_noises)
      print(noise$bv_noises[[1]], padding = padding + 4, prefix = names[[1]])
      print(noise$bv_noises[[2]], padding = padding + 4, prefix = names[[2]])
    } else {
      # single noise
      cat(pad_space); cat(prefix); cat(": "); cat(toupper(noise$noise_type)); cat("\n")

      if (suppress_sigma && noise$noise_type == "normal") {
        # skip
      } else if (suppress_sigma) {
        # only print mu and nu
        cat(paste0(pad_add4_space, ngme_format("mu", noise$theta_mu), "\n",
          pad_add4_space, ngme_format("nu", noise$theta_nu)))
      } else {
        cat(pad_space); cat("Noise parameters: \n")
        params <- with(noise, {
          switch(noise_type,
            "normal" = paste0(pad_add4_space, ngme_format("sigma", theta_sigma)),
            "nig"    = paste0(pad_add4_space, ngme_format("mu", theta_mu),
                        "\n", pad_add4_space, ngme_format("sigma", theta_sigma),
                        "\n", pad_add4_space, ngme_format("nu", theta_nu)),
            "gal"    = paste0(pad_add4_space, ngme_format("mu", theta_mu),
                        "\n", pad_add4_space, ngme_format("sigma", theta_sigma),
                        "\n", pad_add4_space, ngme_format("nu", theta_nu)),
            "normal_nig" = paste0(pad_add4_space, ngme_format("mu", theta_mu),
                        "\n", pad_add4_space, ngme_format("sigma_nig", theta_sigma),
                        "\n", pad_add4_space, ngme_format("nu", theta_nu),
                        "\n", pad_add4_space, ngme_format("sigma_normal", theta_sigma_normal)),
            NULL
          )
        })
        cat(params)
      }
    }
  }
  cat("\n")
  if (noise$corr_measurement) {
    cat(pad_add4_space); cat("correlation(rho) = ");
    cat(format(noise$rho, digits=3)); cat("\n")
  }
  if (all(noise$B_nu %*% noise$theta_nu > 1000))
    cat("(Notice: Parameter nu seems too big, consider use Gaussian noise.)\n")

  invisible(noise)
}

subset_noise <- function(noise, sub_idx, compute_corr) {
  noise$B_mu    <- noise$B_mu[sub_idx, , drop=FALSE]
  noise$B_sigma <- noise$B_sigma[sub_idx, ,drop=FALSE]
  noise$B_nu    <- noise$B_nu[sub_idx, ,drop=FALSE]
  noise$V       <- noise$V[sub_idx]

  if (!is.null(noise$index_corr)) noise$index_corr <- noise$index_corr[sub_idx]

  if (!is.null(noise$corr_measurement) &&
    noise$corr_measurement && compute_corr
  ) {
    p_order <- order(noise$index_corr)
    cov_rc <- compute_corr_index(noise$index_corr[p_order])

    # update noise with extra terms about correlation
    noise$cor_rows <- cov_rc$cor_rows
    noise$cor_cols <- cov_rc$cor_cols
    noise$has_correlation <- cov_rc$has_correlation
    noise$n_corr_pairs <- cov_rc$n_corr_pairs
    noise$index_corr <- noise$index_corr[p_order]
  }
  noise
}