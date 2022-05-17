#' ngme plot
#'
#' @param object  ngme object after model fitting
#' @param param  fe : fix effect;  me : measurement error; p : process
#' @param type 1 : "traj" or 2 : "grad"
#' @param which 0 for fix effects and eps, 1 for 1st latent model, etc
#'
#' @return
#' @export
#'
#' @examples

plot.ngme <- function(object, param="fe", type="traj", which=1) {
  n_latent  <- object$n_latent
  n_fe      <- object$n_fe

  grad_traj <- object$trajectory$grad_traj
  x_traj    <- object$trajectory$x_traj

  n.paras <- length(x_traj[[1]])

  # fixed effects
  if ((param=="fe")) {
    if (((type=="traj") || (type==1)) && (n_fe > 0)) {
      par(mfrow=c(2,2))

      fe_traj = list()
      for (i in 1:n_fe) {
        fe = unlist(lapply(x_traj, function(x) {x[n_latent*4 + i]} ))
        plot(fe, type="l", main = paste("fixed effect", i),
          xlab="iterations", ylab="value")
      }

      par(mfrow=c(1,1))
    }
    return()
  }

  # measurement error
  if (param=="me") {
    if ((type=="traj") || (type==1)) {
      x_eps <- unlist(lapply(x_traj, function(x) {x[n.paras]} ))
      plot(exp(x_eps), type="l", main = "trajectory of sigma_eps",
        xlab="iterations", ylab="sigma_eps")
    }
    return()
  }

  # latent model
  if (param=="la") {
    if ((type=="traj") || (type==1)) {
      x_kappa <- unlist(lapply(x_traj, function(x) {x[(which-1)*4 + 1]} ))
      x_mu    <- unlist(lapply(x_traj, function(x) {x[(which-1)*4 + 2]} ))
      x_sigma <- unlist(lapply(x_traj, function(x) {x[(which-1)*4 + 3]} ))
      x_var   <- unlist(lapply(x_traj, function(x) {x[(which-1)*4 + 4]} ))

      x_kappa = (-1 + 2*exp(x_kappa) / (1+exp(x_kappa)))
      x_sigma = exp(x_sigma)
      x_var   = exp(x_var)

      par(mfrow=c(2,2))
      plot(x_kappa, type="l", main = "traj of kappa")
      plot(x_mu,    type="l", main = "traj of mu")
      plot(x_sigma, type="l", main = "traj of sigma")
      plot(x_var,   type="l", main = "traj of var")
      par(mfrow=c(1,1))
    }

    if ((type=="grad") || (type==2)) {
      grads_kappa <- unlist(lapply(grad_traj, function(x) {x[(which-1)*4 + 1]} ))
      grads_mu    <- unlist(lapply(grad_traj, function(x) {x[(which-1)*4 + 2]} ))
      grads_sigma <- unlist(lapply(grad_traj, function(x) {x[(which-1)*4 + 3]} ))
      grads_var   <- unlist(lapply(grad_traj, function(x) {x[(which-1)*4 + 4]} ))

      par(mfrow=c(2,2))
      plot(grads_kappa, type="l", main = "grad of kappa")
      plot(grads_mu,    type="l", main = "grad of mu")
      plot(grads_sigma, type="l", main = "grad of sigma")
      plot((grads_var), type="l", main = "grad of var")
      par(mfrow=c(1,1))
    }
    return()
  }

  stop("Unknown parameter")
}

