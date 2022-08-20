#' ngme noise specification
#'
#' @param type        type of noise, "nig", "normal"
#' @param theta_V Starting value for theta.sigma
#' @param theta_mu     specify a non-stationary noise using theta_mu
#' @param theta_sigma  specify a non-stationary noise using theta_sigma
#' @param B_mu         Basis matrix for mu (if non-stationary)
#' @param B_sigma      Basis matrix for sigma (if non-stationary)
#'
#' @return
#' @export
#'
#' @examples
#' 

ngme.noise <- function(
  type = "nig",
  theta_V = 1,
  theta_mu = 0,
  theta_sigma = 0,
  B_mu = NULL,
  B_sigma = NULL
) {
  if (type == "gaussian") type <- "normal"
  # check input
  # if (is.null(theta_sigma)) theta_sigma <- 0
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

  # make sure B_mu and B_sigma has same row
  if (nrow(B_mu) == 1 && nrow(B_sigma) != 1) {
    n <- nrow(B_sigma)
    B_mu <- matrix(rep(B_mu, n), nrow = n, byrow = TRUE)
  } else if (nrow(B_mu) != 1 && nrow(B_sigma) == 1) {
    n <- nrow(B_mu)
    B_sigma <- matrix(rep(B_sigma, n), nrow = n, byrow = TRUE)
  }

  if (type == "nig") {

  } else if (type %in% c("normal", "gaussian")) {

  } else {
    stop("Unknown noise type. Please check the mannual.")
  }

  structure(
    list(
      n_noise       = nrow(B_mu),
      type          = type,
      theta_V       = theta_V,
      theta_mu      = theta_mu,
      theta_sigma   = theta_sigma,
      B_mu          = B_mu,
      B_sigma       = B_sigma,
      n_theta_mu    = length(theta_mu),
      n_theta_sigma = length(theta_sigma),
      n_theta_V     = 1
    ),
    class = "ngme.noise"
  )
}

#' create a normal noise
#'
#' @param sd  standard deviation for noise
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
  B_sigma = NULL
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
    type = "normal",
    theta_sigma = theta_sigma,
    B_sigma = B_sigma
  )
}

# update ngme.noise 
update.ngme.noise <- function(noise, n = NULL) {
  stopifnot("n should be integer" = is.numeric(n))
  B_mu <- noise$B_mu
  noise$B_mu <- matrix(data = rep(B_mu, n / nrow(B_mu)), nrow = n)

  B_sigma <- noise$B_sigma
  noise$B_sigma <- matrix(data = rep(B_sigma, n / nrow(B_sigma)), nrow = n)

  noise
}

# tests
# noise <- ngme.noise(B_mu = matrix(rep(1,6), ncol=2), theta_mu = c(2, 2))
# update.ngme.noise(noise, n = 3)