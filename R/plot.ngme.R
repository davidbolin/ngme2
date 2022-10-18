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
#' traceplot(ngme, parameter = "theta_mu", f_index = 1)
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
    geom_line(aes(x = 1:iters, y = apply(data, MARGIN=1, mean)), col="red") +
    xlab("iterations") +
    ylab("value") + guides()
}
