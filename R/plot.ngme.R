# #' ngme plot
# #'
# #' @param object  ngme object
# #'
# #' @return
# #' @export
# #'
# plot.ngme <- function(object) {
#   # wrapper
#   # not implement yet
# }


#' Trace plot of ngme fitting
#'
#' @param ngme ngme object
#' @param param parameter name
#' @param f_index index of the process, 0 stands for fixed effects or measurement noise
#' @param param_index paramter index if non-stationary
#' @param transform any transformation
#'
#' @return the traceplot
#' @export
#'
traceplot <- function(
  ngme,
  param,
  f_index = 0,
  param_index = 1,
  transform = identity
) {
  parameter <- param
  stopifnot("parameter is a string" = is.character(parameter))

  if (parameter == "all") {
    if (!requireNamespace("gridExtra", quietly = TRUE))
      stop(
        "Package \"gridExtra\" must be installed to use parameter = all.",
        call. = FALSE
      )
    # if (f_index == 0) {
    #   switch(ngme$noise$noise_type,
    #     "normal" = traceplot(ngme, parameter = "sigma", f_index = 0),
    #     "nig" = gridExtra::marrangeGrob(
    #       lapply(c("mu", "sigma", "nu"))
    #     )
    #   )
    # }
  } else {
    # change alias
    if (parameter %in% c("kappa", "alpha", "theta_kappa"))
      parameter <- "theta_K"
    else if (parameter %in% c("mu"))
      parameter <- "theta_mu"
    else if (parameter %in% c("sigma_normal"))
      parameter <- "theta_sigma_normal"
    else if (parameter %in% c("sigma"))
      parameter <- "theta_sigma"
    else if (parameter %in% c("nu"))
      parameter <- "theta_V"

    stopifnot("Not a ngme object with trajectoroy." = !is.null(attr(ngme, "trajectory")))
    traj <- attr(ngme, "trajectory")
    iters <- length(traj[[1]]$block_traj[["theta_sigma"]][[1]])

    data <- matrix(nrow = iters, ncol = length(traj))

    for (i in seq_along(traj)) {
      if (f_index == 0) # block model
        data[, i] <- switch(parameter,
          theta_V = traj[[i]]$block_traj[[parameter]],
          traj[[i]]$block_traj[[parameter]][[param_index]]
        )
      else
        # latent model
        data[, i] <- switch(parameter,
          theta_V = traj[[i]]$latents[[f_index]][[parameter]],
          traj[[i]]$latents[[f_index]][[parameter]][[param_index]]
        )
    }

    model_type <- if (f_index > 0) ngme$latents[[f_index]]$model else "unknown"
    # handle transformation of data.
    if ((model_type == "matern" && parameter == "theta_K") ||
        (parameter == "theta_sigma"))
      transform <- exp
    if (model_type == "ar1" && parameter == "theta_K")
      transform <- ar1_th2a

    # Plot function
    # Var1 and Var2 comes from melt
    df <- reshape2::melt(data)
    ggplot() +
      geom_line(data = df, aes(x = .data$Var1, y = transform(.data$value), group = .data$Var2)) +
      geom_line(aes(x = 1:iters, y = transform(apply(data, MARGIN=1, mean))), col="red") +
      xlab("iterations") +
      ylab("value") + guides() +
      labs(title = paste("Traceplot of", switch(parameter,
        theta_K = "theta_K",
        theta_mu = "mu",
        theta_sigma = "sigma",
        theta_sigma_normal = "sigma_normal",
        theta_V = "nu",
        beta = "beta"
      ), param_index))
    }
}

#' plot the density of noise (for stationary)
#'
#' @param x ngme_noise
#' @param y another ngme_noise
#' @param ... ...
#'
#' @return plot
#' @export
#'
#' @examples
#' plot(noise_nig(mu=1, sigma=2, nu=1))
plot.ngme_noise <- function(x, y = NULL, ...) {
  noise <- x; noise2 <- y
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
    ggplot2::geom_line(ggplot2::aes(x = xx, y = dd))

  if (!is.null(noise2)) {
    mu <- noise2$theta_mu
    sigma <- exp(noise2$theta_sigma)
    nu <- noise2$theta_V
    switch(noise2$noise_type,
      "nig"     = dd2 <- dnig(xx, -mu, mu, nu, sigma),
      "normal"  = dd2 <- dnorm(xx, sd = sigma),
      stop("Plot for this type is not implemented")
    )
    gg <- gg + ggplot2::geom_line(ggplot2::aes(x = xx, y = dd2), col = "red")
  }

  gg + ggplot2::labs(title = "Density Plot")
}

#' Trace plot of ngme fitting
#'
#' @param ngme ngme object
#' @param param string, fe for fixed effects, mn for measurement noise, and number i of i-th latent model
#'
#' @return the traceplot
#' @export
#'
traceplot2 <- function(
  ngme,
  param = "beta"
) {
  if (param %in% c("fe", "beta")) {
    n <- length(ngme$beta)
    ps <- lapply(1:n, function(x) traceplot(ngme, f_index=0, param="beta", param_index = x))
  } else if (param %in% c("mn", "noise")) {
    ps <- switch(ngme$noise$noise_type,
      "normal" = {n_sigma <- length(ngme$noise$theta_sigma); lapply(1:n_sigma, function(x) traceplot(ngme, f_index=0, param="sigma", param_index = x))},
      { # default case
        n_mu <- length(ngme$noise$theta_mu);
        n_sigma <- length(ngme$noise$theta_sigma);
        c(
          lapply(1:n_mu, function(x) traceplot(ngme, f_index=0, param="mu", param_index = x)),
          lapply(1:n_sigma, function(x) traceplot(ngme, f_index=0, param="sigma", param_index = x)),
          list(traceplot(ngme, f_index=0, param="nu", param_index = 1))
        )
      }
    )
  } else if (is.numeric(param)) {
    f_idx <- param;
    n_K <- length(ngme$latents[[f_idx]]$theta_K);
    n_mu <- length(ngme$latents[[f_idx]]$noise$theta_mu);
    n_sigma <- length(ngme$latents[[f_idx]]$noise$theta_sigma);

    ps <- switch(ngme$latents[[f_idx]]$noise_type,
      "normal" = {
        c(
          lapply(1:n_K, function(x) traceplot(ngme, f_index=f_idx, param="kappa", param_index = x)),
          lapply(1:n_sigma, function(x) traceplot(ngme, f_index=f_idx, param="sigma", param_index = x))
        )
      },
      "normal_nig" = {
        n_sig_normal <- length(ngme$latents[[f_idx]]$noise$theta_sigma_normal);
        c(
          lapply(1:n_K, function(x) traceplot(ngme, f_index=f_idx, param="kappa", param_index = x)),
          lapply(1:n_mu, function(x) traceplot(ngme, f_index=f_idx, param="mu", param_index = x)),
          lapply(1:n_sigma, function(x) traceplot(ngme, f_index=f_idx, param="sigma", param_index = x)),
          list(traceplot(ngme, f_index=f_idx, param="nu", param_index = 1)),
          lapply(1:n_sig_normal, function(x) traceplot(ngme, f_index=f_idx, param="sigma_normal", param_index = x))
        )
      },
      { # default case
        c(
          lapply(1:n_K, function(x) traceplot(ngme, f_index=f_idx, param="kappa", param_index = x)),
          lapply(1:n_mu, function(x) traceplot(ngme, f_index=f_idx, param="mu", param_index = x)),
          lapply(1:n_sigma, function(x) traceplot(ngme, f_index=f_idx, param="sigma", param_index = x)),
          list(traceplot(ngme, f_index=f_idx, param="nu", param_index = 1))
        )
      }
    )
  } else stop("unknown param")

  if (length(ps) > 1) ps["ncol"]=2
  do.call(grid.arrange, ps)
}