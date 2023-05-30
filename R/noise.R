# This file contains ngme noise specifications

#' @title ngme noise specification
#' @aliases noise_nig noise_normal
#' @description Function for specifying ngme noise.
#' Please use \code{noise_nig} and \code{noise_normal} for simpler usage.
#' Use \code{ngme_noise_types()} to check all the available types.
#'
#' @details The parameterization is given in \code{?nig}. Moreover,
#' for specifying non-stationary mu and sigma,
#' \deqn{\mu = B_{\mu} \theta_{\mu},} and
#' \deqn{\sigma = \exp (B_{\sigma} \theta_{\sigma}),}
#'
#' @param noise_type    type of noise, "nig", "normal"
#' @param mu          specify the NIG noise parameter mu, see \code{?nig}
#' @param sigma       specify the noise parameter sigma, see \code{?nig}
#' @param nu          specify the NIG noise parameter nu (nu>0), see \code{?nig}
#' @param V             start value for V
#' @param theta_mu      specify a non-stationary noise using theta_mu
#' @param B_mu          Basis matrix for mu (if non-stationary)
#' @param theta_sigma   specify a non-stationary noise using theta_sigma
#' @param B_sigma       Basis matrix for sigma (if non-stationary)
#' @param fix_theta_mu    fix the parameter of theta_mu
#' @param fix_theta_sigma  fix the parameter of theta_sigma
#' @param fix_nu   fix the parameter of nu
#' @param fix_V         fix the sampling of V
#' @param theta_sigma_normal for normal nosie with nig noise sharing same parameter
#' @param B_sigma_normal    for normal nosie with nig noise sharing same parameter
#' @param sigma_normal  for normal nosie with nig noise sharing same parameter
#' @param theta_sigma_nig similar to theta_sigma_normal
#' @param B_sigma_nig     similar to B_sigma_nig
#' @param sigma_nig     similar to sigma_normal
#' @param ...       additional arguments
#'
#' @return a list of specification of noise
ngme_noise <- function(
  noise_type,
  mu              = 0,
  sigma           = 1,
  nu              = 1,
  theta_mu        = NULL,
  B_mu            = NULL,
  theta_sigma     = NULL,
  B_sigma         = NULL,
  theta_sigma_normal = NULL,
  B_sigma_normal  = NULL,
  fix_theta_mu    = FALSE,
  fix_theta_sigma = FALSE,
  fix_nu          = FALSE,
  V               = NULL,
  fix_V           = FALSE,
  single_V        = FALSE,
  init_V          = TRUE,
  share_V         = FALSE,
  hessian         = TRUE,
  corr_measurement = FALSE,
  index_corr      = NULL,
  rho             = 0,
  ...
) {
  if (is.null(theta_mu)) theta_mu <- mu
  if (is.null(theta_sigma)) theta_sigma <- log(sigma)

  stopifnot("Unkown noise type. Please check ngme_noise_types()" =
    noise_type %in% ngme_noise_types())
  stopifnot("ngme_noise: nu should be positive" = nu > 0)

  if (is.null(B_mu)) B_mu <- as.matrix(1)
  if (is.null(B_sigma)) B_sigma <- as.matrix(1)

  stopifnot(
    "Please input B_mu as a matrix." = is.matrix(B_mu),
    "Please input B_sigma as a matrix." = is.matrix(B_sigma),
    "Please make sure ncol(B_mu) == length(theta_mu)." = ncol(B_mu) == length(theta_mu),
    "Please make sure ncol(B_sigma) == length(theta_sigma)." = ncol(B_sigma) == length(theta_sigma)
  )

  structure(
    list(
      noise_type      = noise_type,
      nu              = nu,
      n_nu            = length(nu),
      V               = V,
      theta_mu        = theta_mu,
      theta_sigma     = theta_sigma,
      theta_sigma_normal     = theta_sigma_normal,
      B_mu            = B_mu,
      B_sigma         = B_sigma,
      B_sigma_normal  = B_sigma_normal,
      n_theta_mu      = length(theta_mu),
      n_theta_sigma   = length(theta_sigma),
      n_theta_sigma_normal   = length(theta_sigma_normal),
      fix_theta_mu    = fix_theta_mu,
      fix_theta_sigma = fix_theta_sigma,
      fix_nu          = fix_nu,
      fix_V           = fix_V,
      n_params        = length(theta_mu) + length(theta_sigma) + length(nu),
      init_V          = init_V,
      single_V        = single_V,
      share_V         = share_V,
      hessian         = hessian,
      corr_measurement = corr_measurement,
      index_corr      = index_corr,
      rho             = rho,
      ...
    ),
    class = "ngme_noise"
  )
}

#' @rdname ngme_noise
#' @export
#' @examples
#' noise_normal(n = 10, sigma = 2)
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
#' noise_nig(mu = 1, sigma = 2, nu = 1, n=10)
noise_nig <- nig <- function(
  mu            = NULL,
  sigma         = NULL,
  nu            = NULL,
  V             = NULL,
  theta_mu      = NULL,
  theta_sigma   = NULL,
  B_mu          = matrix(1),
  B_sigma       = matrix(1),
  corr_measurement = FALSE,
  index_corr      = NULL,
  ...
) {
  # if nothing, then fill with default
  stopifnot("Please use theta_mu for non-stationary mu." = length(mu) < 2)
  if (is.null(mu) && is.null(theta_mu)) theta_mu <- 0
  if (is.null(sigma) && is.null(theta_sigma)) theta_sigma <- 0
  if (is.null(nu)) nu <- 1

  if (!is.null(nu) && nu <= 0) stop("ngme_nosie: nu should be positive.")
  if (!is.null(sigma) && sigma <= 0) stop("ngme_nosie: sigma should be positive.")

  if (!is.null(mu))     theta_mu <- mu
  if (!is.null(sigma))  theta_sigma <- log(sigma)

  ngme_noise(
    noise_type = "nig",
    theta_mu = theta_mu,
    theta_sigma = theta_sigma,
    nu = nu,
    V = V,
    B_mu = B_mu,
    B_sigma = B_sigma,
    corr_measurement = corr_measurement,
    index_corr      = index_corr,
    ...
  )
}

#' @rdname ngme_noise
#' @export
#' @examples
#' noise_gal(mu = 1, sigma = 2, nu = 1, n=10)
noise_gal <- gal <- function(
  mu            = NULL,
  sigma         = NULL,
  nu            = NULL,
  V             = NULL,
  theta_mu      = NULL,
  theta_sigma   = NULL,
  B_mu          = matrix(1),
  B_sigma       = matrix(1),
  corr_measurement = FALSE,
  index_corr      = NULL,
  ...
) {
  # if nothing, then fill with default
  stopifnot("Please use theta_mu for non-stationary mu." = length(mu) < 2)
  if (is.null(mu) && is.null(theta_mu)) theta_mu <- 0
  if (is.null(sigma) && is.null(theta_sigma)) theta_sigma <- 0
  if (is.null(nu)) nu <- 1

  if (!is.null(nu) && nu <= 0) stop("ngme_nosie: nu should be positive.")
  if (!is.null(sigma) && sigma <= 0) stop("ngme_nosie: sigma should be positive.")

  if (!is.null(mu))     theta_mu <- mu
  if (!is.null(sigma))  theta_sigma <- log(sigma)

  ngme_noise(
    noise_type = "gal",
    theta_mu = theta_mu,
    theta_sigma = theta_sigma,
    nu = nu,
    V = V,
    B_mu = B_mu,
    B_sigma = B_sigma,
    corr_measurement = corr_measurement,
    index_corr      = index_corr,
    ...
  )
}

# update noise
update_noise <- function(noise, n = NULL, new_noise = NULL) {
  # update with length n
  if (!is.null(n)) {
    stopifnot("n should be integer" = is.numeric(n))
    B_mu <- noise$B_mu
stopifnot("n / nrow(B_mu) not integer" = abs(n/nrow(B_mu) - round(n/nrow(B_mu))) < 1e-4)
    noise$B_mu <- matrix(data = rep(B_mu, n / nrow(B_mu)), nrow = n)

    B_sigma <- noise$B_sigma
stopifnot("n / nrow(B_sigma) not integer" = abs(n/nrow(B_sigma) - round(n/nrow(B_sigma))) < 1e-4)
    noise$B_sigma <- matrix(data = rep(B_sigma, n / nrow(B_sigma)), nrow = n)

    if (noise$noise_type == "normal_nig") {
      B_sigma_normal <- noise$B_sigma_normal
      noise$B_sigma_normal <- matrix(data = rep(B_sigma_normal, n / nrow(B_sigma_normal)), nrow = n)
    }
    noise <- do.call(ngme_noise, noise)
  } else if (!is.null(new_noise)) {
    # update noise after estimation
    noise$theta_mu           <- new_noise$theta_mu
    noise$theta_sigma        <- new_noise$theta_sigma
    noise$nu                 <- new_noise$nu
    noise$rho                <- new_noise$rho
    if (!is.null(new_noise$V)) noise$V <- new_noise$V

    # bv noise
    if (length(noise$noise_type) == 2) {
      # pass mu, sigma, nu to sub_models
      n_theta_mu1 <- noise$bv_noises[[1]]$n_theta_mu
      n_theta_mu2 <- noise$bv_noises[[2]]$n_theta_mu
      n_theta_sigma1 <- noise$bv_noises[[1]]$n_theta_sigma
      n_theta_sigma2 <- noise$bv_noises[[2]]$n_theta_sigma
      noise$bv_noises[[1]]$theta_mu <- head(noise$theta_mu, n_theta_mu1)
      noise$bv_noises[[2]]$theta_mu <- tail(noise$theta_mu, n_theta_mu2)
      noise$bv_noises[[1]]$theta_sigma <- head(noise$theta_sigma, n_theta_sigma1)
      noise$bv_noises[[2]]$theta_sigma <- tail(noise$theta_sigma, n_theta_sigma2)
      noise$bv_noises[[1]]$nu <- noise$nu[[1]]
      noise$bv_noises[[2]]$nu <- noise$nu[[2]]
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
  n             = 1,
  V             = NULL,
  theta_mu      = NULL,
  theta_sigma_nig   = NULL,
  theta_sigma_normal   = NULL,
  B_mu          = matrix(1),
  B_sigma_nig   = matrix(1),
  B_sigma_normal = matrix(1),
  corr_measurement = FALSE,
  index_corr      = NULL,
  ...
) {
  # if nothing, then fill with default
  stopifnot("Please use theta_mu for non-stationary mu." = length(mu) < 2)
  if (is.null(mu) && is.null(theta_mu)) theta_mu <- 0
  if (is.null(sigma_nig) && is.null(theta_sigma_nig)) theta_sigma_nig <- 0
  if (is.null(nu)) nu <- 1
  if (is.null(sigma_normal) && is.null(theta_sigma_normal)) theta_sigma_normal <- 0

  if (!is.null(nu) && nu <= 0) stop("ngme_nosie: nu should be positive.")
  if (!is.null(sigma_nig) && sigma_nig <= 0) stop("ngme_nosie: sigma_nig should be positive.")
  if (!is.null(sigma_normal) && sigma_normal <= 0) stop("ngme_nosie: sigma_nig should be positive.")

  if (!is.null(mu))     theta_mu <- mu
  if (!is.null(sigma_nig))  theta_sigma_nig <- log(sigma_nig)
  if (!is.null(sigma_normal))  theta_sigma_normal <- log(sigma_normal)

  ngme_noise(
    noise_type = "normal_nig",
    theta_mu = theta_mu,
    theta_sigma = theta_sigma_nig,
    theta_sigma_normal = theta_sigma_normal,
    nu = nu,
    V = V,
    B_mu = B_mu,
    B_sigma = B_sigma_nig,
    B_sigma_normal = B_sigma_normal,
    n = n,
    corr_measurement = corr_measurement,
    index_corr      = index_corr,
    ...
  )
}

#' Print ngme noise
#'
#' @param x noise object
#' @param padding number of white space padding in front
#' @param ... ...
#'
#' @return a list (noise specifications)
#' @export
print.ngme_noise <- function(x, padding = 0, prefix = "Noise type", ...) {
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
        cat("Bivariate type-G1 noise (single_V && share_V) :");
      } else if (noise$single_V && !noise$share_V) {
        cat("Bivariate type-G2 noise (single_V) :");
      } else if (!noise$single_V && noise$share_V) {
        cat("Bivariate type-G3 noise (share_V) :");
      } else {
        cat("Bivariate type-G4 noise :");
      }
      cat("\n")
      names <- names(noise$bv_noises)
      print(noise$bv_noises[[1]], padding = padding + 4, prefix = names[[1]])
      print(noise$bv_noises[[2]], padding = padding + 4, prefix = names[[2]])
    } else {
      # single noise
      cat(pad_space); cat(prefix); cat(": "); cat(toupper(noise$noise_type)); cat("\n")

      cat(pad_space); cat("Noise parameters: \n")
      params <- with(noise, {
        switch(noise_type,
          "normal" = paste0(pad_add4_space, ngme_format("sigma", theta_sigma)),
          "nig"    = paste0(pad_add4_space, ngme_format("mu", theta_mu),
                      "\n", pad_add4_space, ngme_format("sigma", theta_sigma),
                      "\n", pad_add4_space, ngme_format("nu", nu)),
          "gal"    = paste0(pad_add4_space, ngme_format("mu", theta_mu),
                      "\n", pad_add4_space, ngme_format("sigma", theta_sigma),
                      "\n", pad_add4_space, ngme_format("nu", nu)),
          "normal_nig" = paste0(pad_add4_space, ngme_format("mu", theta_mu),
                      "\n", pad_add4_space, ngme_format("sigma_nig", theta_sigma),
                      "\n", pad_add4_space, ngme_format("nu", nu),
                      "\n", pad_add4_space, ngme_format("sigma_normal", theta_sigma_normal)),
          stop("unknown noise type")
        )
      })
      cat(params)
    }
  }
  cat("\n")
  if (noise$corr_measurement) {
    cat(pad_add4_space); cat("correlation(rho) = ");
    cat(format(noise$rho, digits=3)); cat("\n")
  }
  invisible(noise)
}
