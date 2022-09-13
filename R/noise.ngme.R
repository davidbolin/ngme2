#' ngme noise specification
#'
#' Function for specifying ngme noise.
#' Use ngme.noise.types() to check all the available types.
#'
#' @param noise_type        type of noise, "nig", "normal"
#' @param theta_V     value for theta_V, theta_V = eta > 0
#' @param V           value for V
#' @param theta_mu     specify a non-stationary noise using theta_mu
#' @param theta_sigma  specify a non-stationary noise using theta_sigma
#' @param B_mu         Basis matrix for mu (if non-stationary)
#' @param B_sigma      Basis matrix for sigma (if non-stationary)
#' @param fix_mu
#' @param fix_sigma
#' @param fix_var
#' @param fix_V
#'
#' @return a list of specification of noise
#' @export
#'
#' @examples
#'
ngme.noise <- function(
  noise_type  = "nig",
  theta_mu    = 0,
  theta_sigma = 0,
  theta_V     = 1,
  B_mu        = NULL,
  B_sigma     = NULL,
  V           = NULL,
  fix_mu      = FALSE,
  fix_sigma   = FALSE,
  fix_var     = FALSE,
  fix_V       = FALSE
) {
  # check input
  if (noise_type == "gaussian") noise_type <- "normal"
  stopifnot("Unkown noise type. Please check ngme.noise.types()" =
    noise_type %in% ngme.noise.types())

  stopifnot("theta_V is positive" = theta_V > 0)

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
      n_noise       = nrow(B_mu),
      noise_type    = noise_type,
      theta_V       = theta_V,
      V             = V,
      theta_mu      = theta_mu,
      theta_sigma   = theta_sigma,
      B_mu          = B_mu,
      B_sigma       = B_sigma,
      n_theta_mu    = length(theta_mu),
      n_theta_sigma = length(theta_sigma),
      n_theta_V     = 1,
      fix_mu        = fix_mu,
      fix_sigma     = fix_sigma,
      fix_var       = fix_var,
      fix_V         = fix_V
    ),
    class = "noise"
  )
}

#' Specify a normal noise
#'
#' @param sd  standard deviation
#' @param theta_sigma  specify a non-stationary noise using theta_sigma
#' @param B_sigma      Basis matrix for sigma (if non-stationary)
#'
#' @return
#' @export
#'
#' @examples
ngme.noise.normal <- function(
  sd = NULL,
  theta_sigma = NULL,
  B_sigma = NULL,
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

  ngme.noise(
    noise_type = "normal",
    theta_sigma = theta_sigma,
    B_sigma = B_sigma,
    ...
  )
}


#' Specify a nig noise (normal inverse Gaussian)
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
ngme.noise.nig <- function(
  theta_mu = 0,
  theta_sigma = 0,
  theta_V = 1,
  V = NULL,
  B_mu = matrix(1),
  B_sigma = matrix(1),
  ...
) {
  ngme.noise(
    theta_mu = theta_mu,
    theta_sigma = theta_sigma,
    theta_V = theta_V,
    V = V,
    B_mu = B_mu,
    B_sigma = B_sigma,
    ...
  )
}

# update ngme.noise
update.ngme.noise <- function(noise, n = NULL) {
  stopifnot("n should be integer" = is.numeric(n))
  B_mu <- noise$B_mu
  noise$B_mu <- matrix(data = rep(B_mu, n / nrow(B_mu)), nrow = n)

  B_sigma <- noise$B_sigma
  noise$B_sigma <- matrix(data = rep(B_sigma, n / nrow(B_sigma)), nrow = n)

  noise$n_noise <- n

  noise
}

#' Create ngme noise with a list
#'
#' @param x a list
#'
#' @return a list of specification for ngme
#' @export
#'
#' @examples
create.ngme.noise <- function(x) {
  do.call(ngme.noise, x)
  # ngme.noise(unlist(x))
  # print(unlist(x))
}

# tests
# noise <- ngme.noise(B_mu = matrix(rep(1,6), ncol=2), theta_mu = c(2, 2))
# update.ngme.noise(noise, n = 3)
