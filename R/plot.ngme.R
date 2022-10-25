#' ngme plot
#'
#' @param object  ngme object
#'
#' @return
#' @export
#'
#' @examples

plot.ngme <- function(object) {
  # wrapper
  # not implement yet
}


#' Trace plot of ngme fitting
#'
#' @param ngme ngme object
#' @param parameter parameter name
#' @param f_index index of the process, 0 stands for fixed effects or measurement noise
#' @param param_index paramter index if non-stationary
#' @param transform any transformation
#'
#' @return
#' @export
#'
#' @examples
traceplot <- function(
  ngme_out,
  parameter,
  f_index = 0,
  param_index = 1,
  transform = identity
) {
  stopifnot("Not a ngme object with trajectoroy." = !is.null(attr(ngme_out, "trajectory")))
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

  df <- reshape2::melt(data)

  # Var1 and Var2 comes from melt
  library(ggplot2)
  ggplot() +
    geom_line(data = df, aes(x = Var1, y = transform(value), group = Var2)) +
    geom_line(aes(x = 1:iters, y = transform(apply(data, MARGIN=1, mean))), col="red") +
    xlab("iterations") +
    ylab("value") + guides() + labs(title = paste("Traceplot of", parameter))
}
#' plot the density of noise (for stationary)
#'
#' @param noise1 ngme_noise
#' @param noise2 ngme_noise (for comparison)
#'
#' @return plot
#' @export
#'
#' @examples
plot.ngme_noise <- function(noise, noise2 = NULL, ...) {
  mu <- noise$theta_mu
  sigma <- exp(noise$theta_sigma)
  nu <- noise$theta_V
  stopifnot("only implemented for stationary mu" = length(mu) == 1)
  stopifnot("only implemented for stationary sigma" = length(sigma) == 1)

  xx <- seq(-15, 15, length = 400)
  switch(noise$noise_type,
    "nig"     = dd <- dnig(xx, -mu, mu, nu, sigma),
    "normal"  = dd <- dnorm(xx, sd = sigma),
    stop("Plot for this type is not implemented")
  )

  gg <- ggplot2::ggplot() +
    ggplot2::geom_line(aes(x = xx, y = dd))

  if (!is.null(noise2)) {
    mu <- noise2$theta_mu
    sigma <- exp(noise2$theta_sigma)
    nu <- noise2$theta_V
    switch(noise2$noise_type,
      "nig"     = dd2 <- dnig(xx, -mu, mu, nu, sigma),
      "normal"  = dd2 <- dnorm(xx, sd = sigma),
      stop("Plot for this type is not implemented")
    )
    gg <- gg + ggplot2::geom_line(aes(x = xx, y = dd2), col = "red")
  }

  gg + ggplot2::labs(title = "Density Plot")
}