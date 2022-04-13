#' Plot the optimization trajectory
#'
#' @param output  model after optimization
#' @param order 0 for fix effects and eps, 1 for 1st latent model, etc
#' @param type 1 : "traj" or 2 : "grad"
#'
#' @return
#' @export
#'
#' @examples
plot_out <- function(output, order=1, type="traj") {
  grad_traj <- output$grad_traj
  x_traj <- output$x_traj

  n.paras <- length(output$grad_traj[[1]])

  if (order==0) {
    if ((type=="traj") || (type==1)) {
      x_eps <- unlist(lapply(x_traj, function(x) {x[n.paras]} ))
      plot(exp(x_eps), type="l", main = "traj of sigma_eps")
    }

    if ((type=="grad") || (type==2)){
      grads_eps <- unlist(lapply(grad_traj, function(x) {x[n.paras]} ))
      plot(grads_eps, type="l", main = "grad of sigma_eps")
    }
    return()
  }

  if (order==-1) {
    par(mfrow=c(2,2))
    fe1 <- unlist(lapply(x_traj, function(x) {x[n.paras-2]} ))
    fe2 <- unlist(lapply(x_traj, function(x) {x[n.paras-1]} ))
    gr1 <- unlist(lapply(grad_traj, function(x) {x[n.paras-2]} ))
    gr2 <- unlist(lapply(grad_traj, function(x) {x[n.paras-1]} ))

    plot(fe1, type="l", main = "traj of fe1")
    plot(fe2, type="l", main = "traj of fe2")
    plot(gr1, type="l", main = "grad of fe1")
    plot(gr2, type="l", main = "grad of fe2")

    par(mfrow=c(1,1))
    return()
  }

  if ((type=="traj") || (type==1)) {
    x_kappa <- unlist(lapply(x_traj, function(x) {x[(order-1)*4 + 1]} ))
    x_mu    <- unlist(lapply(x_traj, function(x) {x[(order-1)*4 + 2]} ))
    x_sigma <- unlist(lapply(x_traj, function(x) {x[(order-1)*4 + 3]} ))
    x_var   <- unlist(lapply(x_traj, function(x) {x[(order-1)*4 + 4]} ))

    x_kappa <- (-1 + 2*exp(x_kappa) / (1+exp(x_kappa)))

    par(mfrow=c(2,2))
    plot(x_kappa, type="l", main = "traj of kappa")
    plot(x_mu,    type="l", main = "traj of mu")
    plot(exp(x_sigma), type="l", main = "traj of sigma")
    plot(exp(x_var),   type="l", main = "traj of var")
    par(mfrow=c(1,1))

    return()
  }

  if ((type=="grad") || (type==2)) {
    grads_kappa <- unlist(lapply(grad_traj, function(x) {x[(order-1)*4 + 1]} ))
    grads_mu    <- unlist(lapply(grad_traj, function(x) {x[(order-1)*4 + 2]} ))
    grads_sigma <- unlist(lapply(grad_traj, function(x) {x[(order-1)*4 + 3]} ))
    grads_var   <- unlist(lapply(grad_traj, function(x) {x[(order-1)*4 + 4]} ))

    par(mfrow=c(2,2))
    plot(grads_kappa, type="l", main = "grad of kappa")
    plot(grads_mu,    type="l", main = "grad of mu")
    plot(grads_sigma, type="l", main = "grad of sigma")
    plot((grads_var), type="l", main = "grad of var")
    par(mfrow=c(1,1))

    return()
  }



  stop("Unknown parameter")
}

