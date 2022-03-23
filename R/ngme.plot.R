plot_out <- function(output) {
  grad_traj <- output$grad_traj
  x_traj <- output$x_traj

  grads_kappa <- unlist(lapply(grad_traj, function(x) {x[1]} ))
  grads_mu    <- unlist(lapply(grad_traj, function(x) {x[2]} ))
  grads_sigma <- unlist(lapply(grad_traj, function(x) {x[3]} ))
  grads_var   <- unlist(lapply(grad_traj, function(x) {x[4]} ))

  x_kappa <- unlist(lapply(x_traj, function(x) {x[1]} ))
  x_mu    <- unlist(lapply(x_traj, function(x) {x[2]} ))
  x_sigma <- unlist(lapply(x_traj, function(x) {x[3]} ))
  x_var   <- unlist(lapply(x_traj, function(x) {x[4]} ))

  par(mfrow=c(2,2))
  plot(grads_kappa, type="l", main="grad of kappa")
  plot(grads_mu, type="l", main="grad of mu")
  plot(grads_sigma, type="l", main="grad of sigma")
  plot((grads_var), type="l", main="grad of var")

  plot(x_kappa, type="l", main="traj of kappa")
  plot(x_mu, type="l", main="traj of mu")
  plot(x_sigma, type="l", main="traj of sigma")
  plot(x_var, type="l", main="traj of var")
}

