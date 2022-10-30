# This file contains ngme noise specifications (NIG and Normal)

#' ngme noise specification
#'
#' Function for specifying ngme noise.
#' Please use noise_type for simpler specification.
#' Use ngme_noise_types() to check all the available types.
#'
#' @param noise_type    type of noise, "nig", "normal"
#' @param n             number of noise
#' @param V             value for V
#' @param B_mu          Basis matrix for mu (if non-stationary)
#' @param theta_mu      specify a non-stationary noise using theta_mu
#' @param B_sigma       Basis matrix for sigma (if non-stationary)
#' @param theta_sigma   specify a non-stationary noise using theta_sigma
#' @param theta_V       value for theta_V, theta_V = eta > 0
#' @param fix_theta_mu    fix the parameter of theta_mu
#' @param fix_theta_sigma  fix the parameter of theta_sigma
#' @param fix_theta_V   fix the parameter of theta_V
#' @param fix_V         fix the sampling of V
#'
#' @return a list of specification of noise
#' @export
#'
ngme_noise <- function(
  noise_type,
  theta_mu        = 0,
  theta_sigma     = 0,
  theta_V         = 1,
  V               = NULL,
  B_mu            = NULL,
  B_sigma         = NULL,
  fix_theta_mu    = FALSE,
  fix_theta_sigma = FALSE,
  fix_theta_V     = FALSE,
  fix_V           = FALSE,
  ...
) {
  # check input
  stopifnot("Unkown noise type. Please check ngme_noise_types()" =
    noise_type %in% ngme_noise_types())

  stopifnot("ngme_noise: theta_V should be positive" = theta_V > 0)

  if (is.null(B_mu))    B_mu <- as.matrix(1)
  if (is.null(B_sigma)) B_sigma <- as.matrix(1)

  if (!is.matrix(B_mu))
    stop("Please input B_mu as a matrix to use non-stationary mu")
  if (!is.matrix(B_sigma))
    stop("Please input B_sigma as a matrix to use non-stationary sigma")

  if (ncol(B_mu) != length(theta_mu))
    stop("Please make sure ncol(B_mu) == length(theta_mu).")
  if (ncol(B_sigma) != length(theta_sigma))
    stop("Please make sure ncol(B_sigma) == length(theta_sigma).")

  # auto-complete (make sure nrow(B_sigma) == nrow(B_mu) for n=1 case)
  if (nrow(B_mu) == 1 && nrow(B_sigma) != 1) {
    n <- nrow(B_sigma)
    B_mu <- matrix(rep(B_mu, n), nrow = n, byrow = TRUE)
  } else if (nrow(B_mu) != 1 && nrow(B_sigma) == 1) {
    n <- nrow(B_mu)
    B_sigma <- matrix(rep(B_sigma, n), nrow = n, byrow = TRUE)
  }

  structure(
    list(
      n_noise         = nrow(B_mu),  # this is same as V_size
      noise_type      = noise_type,
      theta_V         = theta_V,
      V               = V,
      theta_mu        = theta_mu,
      theta_sigma     = theta_sigma,
      B_mu            = B_mu,
      B_sigma         = B_sigma,
      n_theta_mu      = length(theta_mu),
      n_theta_sigma   = length(theta_sigma),
      n_theta_V       = 1,
      fix_theta_mu    = fix_theta_mu,
      fix_theta_sigma = fix_theta_sigma,
      fix_theta_V     = fix_theta_V,
      fix_V           = fix_V,
      n_params        = switch(noise_type,
        "normal" = length(theta_sigma),
        "nig"    = length(c(theta_mu, theta_sigma, 1))
      ) # parameter to estimate
    ),
    class = "ngme_noise"
  )
}

#' Specify a normal noise
#'
#' @param sd            standard deviation
#' @param theta_sigma  specify a non-stationary noise using theta_sigma
#' @param B_sigma      Basis matrix for sigma (if non-stationary)
#'
#' @return
#' @export
#'
#' @examples
noise_normal <- function(
  sd = NULL,
  theta_sigma = NULL,
  B_sigma = matrix(1, 1, 1),
  n = nrow(B_sigma),
  ...
) {
  if (!is.null(sd) && !is.null(theta_sigma))
    stop("Please only use sd or theta_sigma as input")

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
    n = n,
    ...
  )
}

#' Specify a nig noise (normal inverse Gaussian)
#'
#' The parameterization can be found in ...
#'
#' @param theta_V      value of eta
#' @param theta_mu     specify a non-stationary noise using theta_mu
#' @param theta_sigma  specify a non-stationary noise using theta_sigma
#' @param B_mu         Basis matrix for mu (if non-stationary)
#' @param B_sigma      Basis matrix for sigma (if non-stationary)
#'
#' @return a list of specification for ngme
#' @export
#'
#' @examples
noise_nig <- function(
  mu            = NULL,
  sigma         = NULL,
  nu            = NULL,
  theta_mu      = NULL,
  theta_sigma   = NULL,
  theta_V       = NULL,
  V             = NULL,
  B_mu          = matrix(1),
  B_sigma       = matrix(1),
  ...
) {
  # if nothing, then fill with default
  stopifnot("Please use theta_mu for non-stationary mu." = length(mu) < 2)
  if (is.null(mu) && is.null(theta_mu)) theta_mu <- 0
  if (is.null(sigma) && is.null(theta_sigma)) theta_sigma <- 0
  if (is.null(nu) && is.null(theta_V)) theta_V <- 1

  if (!is.null(sigma) && sigma <= 0) stop("ngme_nosie: sigma should be positive.")
  if (!is.null(nu) && nu <= 0) stop("ngme_nosie: nu should be positive.")

  if (!is.null(mu))     theta_mu <- mu
  if (!is.null(sigma))  theta_sigma <- log(sigma)
  if (!is.null(nu))     theta_V <- nu

  ngme_noise(
    noise_type = "nig",
    theta_mu = theta_mu,
    theta_sigma = theta_sigma,
    theta_V = theta_V,
    V = V,
    B_mu = B_mu,
    B_sigma = B_sigma,
    ...
  )
}

# update noise
update_noise <- function(noise, n = NULL) {
  stopifnot("n should be integer" = is.numeric(n))
  B_mu <- noise$B_mu
  noise$B_mu <- matrix(data = rep(B_mu, n / nrow(B_mu)), nrow = n)

  B_sigma <- noise$B_sigma
  noise$B_sigma <- matrix(data = rep(B_sigma, n / nrow(B_sigma)), nrow = n)

  noise$n_noise <- n

  # noise
  do.call(ngme_noise, noise)
}

#' Create ngme noise with a list
#'
#' @param x a list
#'
#' @return a list of specification for ngme
#' @export
#'
#' @examples
create_noise <- function(x) {
  do.call(ngme_noise, x)
}

#' Print ngme noise
#'
#' @param noise noise object
#' @param padding
#'
#' @return a list (noise specifications)
#' @export
print.ngme_noise <- function(noise, padding=0) {
  pad_space <- paste(rep(" ", padding), collapse = "")
  pad_add4_space <- paste(rep(" ", padding + 4), collapse = "")

  if (is.null(noise)) {
    cat(pad_space); cat("Noise type - "); cat("NULL"); cat("\n")
  } else {
    cat(pad_space); cat("Noise type - "); cat(noise$noise_type); cat("\n")

    cat(pad_space); cat("Noise parameters: \n")
    params <- with(noise, {
      switch(noise_type,
        "normal" = paste0(pad_add4_space, ngme_format("sigma", theta_sigma)),
        "nig"    = paste0(pad_add4_space, ngme_format("mu", theta_mu),
                    "\n", pad_add4_space, ngme_format("sigma", theta_sigma),
                    "\n", pad_add4_space, ngme_format("nu", theta_V)),
        stop("unknown noise type")
      )
    })
    cat(params)
  }
  cat("\n")
  invisible(noise)
}