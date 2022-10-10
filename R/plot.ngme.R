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

      if (object$model.types == "ar1") {
        x_kappa = (-1 + 2*exp(x_kappa) / (1+exp(x_kappa)))
        k_title = "traj of alpha"
      } else if (object$model.types == "matern") {
        x_kappa = exp(x_kappa)
        k_title = "traj of kappa"
      }
      x_sigma = exp(x_sigma)
      x_var   = exp(x_var)

      par(mfrow=c(2,2))
      plot(x_kappa, type="l", main = k_title)
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


#' Trace plot
#'
#' @param trajectory
#' @param start
#' @param n
#' @param transform
#' @param ylab
#'
#' @return
#' @export
#'
#' @examples
traceplot2 <- function(
  ngme,
  start=1,
  n=1,
  transform=identity,
  ...
) {
  # plot trajectory of out$trajectory[]
  trajectory <- attr(ngme, "trajectory")

  if (n == 2) par(mfrow = c(2, 1))
  if (n == 3 || n == 4) par(mfrow = c(2, 2))

  for (i in start:(start + n - 1)) {
    y <-  unlist(lapply(trajectory$x_traj, function(x) { x[i] }))
    plot(transform(y), type = "l", ...)
  }

  par(mfrow = c(1, 1))
}

#' Trace plot of the estimation
#'
#' @param trajectory
#' @param start
#' @param n
#' @param transform
#' @param ylab
#'
#' @return
#' @export
#'
#' @examples
traceplot <- function(
  ngme,
  ...
) {
  # plot beta
  if (length(ngme[[1]]$beta) > 0) {

  }

  # plot measure noise mu
  for (i in seq_along(ngme)) {
    plot.or.lines <- if (i == 1) plot else lines
    if (ngme[[i]]$noise$noise_type == "nig") {
      traj <- attr(ngme[[i]], "trajectory")
      mu_traj <- unlist(traj$theta_mu_traj)
      plot.or.lines(seq_along(mu_traj), mu_traj, type = "l")
      # sigma_traj <- unlist(traj$theta_sigma_traj)
      # plot.or.lines(seq_along(sigma_traj), sigma_traj, type = "l")

      # plot.or.lines(seq_along(mu_traj), traj$theta_V_traj, type = "l")
    }
  }

  # plot models
  # for (i in seq_along(ngme)) {
  #   object <- ngme[[i]]
  #   block_traj <- attr(object, "trajectory")

  #   # plot beta
  #   if (length(object$noise) > 0) {

  #   }


  # }

  par(mfrow = c(1, 1))
  invisible(ngme)
}


# ideal plot_chains(ngme_out, "beta")
# ideal plot_chains(ngme_out, "block_theta_mu")
# ideal plot_chains(ngme_out, "block_theta_sigma")
# ideal plot_chains(ngme_out, "block_theta_V")
plot_chains <- function(ngme_out, parameter, f_index = 0, param_index = 1) {
  traj <- attr(ngme_out, "trajectory")
  iters <- length(traj[[1]]$block_traj[["theta_sigma"]][[1]])

  data <- matrix(nrow = iters, ncol = length(traj))

  for (i in seq_along(traj)) {
    if (f_index == 0)
      data[, i] <- switch(parameter,
        theta_V = traj[[i]]$block_traj[[parameter]],
        traj[[i]]$block_traj[[parameter]][[param_index]]
      )
    else
      data[, i] <- switch(parameter,
        theta_V = traj[[i]]$latents[[f_index]][[parameter]],
        traj[[i]]$latents[[f_index]][[parameter]][[param_index]]
      )
  }

  library(reshape2)
  df <- reshape2::melt(data)

  library(ggplot2)
  ggplot() +
    geom_line(data = df, aes(x = Var1, y = value, group = Var2)) +
    geom_line(aes(x = 1:iters, y = apply(data, MARGIN=1, mean)), col="red") +
    xlab("iterations") +
    ylab("value") + guides()
}
