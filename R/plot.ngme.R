#' Trace plot of ngme fitting
#'
#' @param ngme ngme object
#' @param name name or index of the process,
#'   c("fe", "beta") is fixed effects,
#'   c("mn", "noise") is measurement nosie
#' @param param parameter, "mu", "sigma", etc
#' @param param_index paramter index if non-stationary
#' @param transform any transformation
#'
#' @return the traceplot
#' @export
#'
traceplot <- function(
  ngme,
  name = 0,
  param,
  param_index = 1,
  transform = identity
) {
  # turn name into index (0 means fe or mn, 1, 2, .. means f_idx)
  if (name == 0 || name %in% c("fe", "beta", "mn", "noise"))
    f_idx <- 0
  else
    f_idx <- if (is.numeric(name)) name else which(names(ngme$latents) == name)
  stopifnot("Unknown name argument, see ?traceplot" = length(f_idx) == 1)

  parameter <- param
  stopifnot("parameter is a string" = is.character(parameter))

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
    parameter <- "nu"

  stopifnot("Not a ngme object with trajectoroy." = !is.null(attr(ngme, "trajectory")))
  traj <- attr(ngme, "trajectory")
  iters <- length(traj[[1]]$block_traj[["theta_sigma"]][[1]])

  data <- matrix(nrow = iters, ncol = length(traj))
  model_type <- "unknown"

  for (i in seq_along(traj)) {
    if (f_idx == 0) # block model
      data[, i] <- switch(parameter,
        nu = traj[[i]]$block_traj[[parameter]],
        traj[[i]]$block_traj[[parameter]][[param_index]]
      )
    else {
      model_type <- ngme$latents[[f_idx]]$model
      # latent model
      if (is.null(ngme$latents[[f_idx]]))
        stop("Unknown name argument, please check!")
      data[, i] <- switch(parameter,
        nu = traj[[i]]$latents[[f_idx]][[parameter]],
        traj[[i]]$latents[[f_idx]][[parameter]][[param_index]]
      )
    }
  }

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
      nu = "nu",
      beta = "beta"
    ), param_index))
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
  nu <- noise$nu
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
    nu <- noise2$nu
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
#' @param name name or index of the process,
#'   c("fe", "beta") is fixed effects,
#'   c("mn", "noise") is measurement nosie
#'
#' @return the traceplot
#' @export
#'
traceplot2 <- function(
  ngme,
  name = "beta"
) {
  if (!requireNamespace("gridExtra", quietly = TRUE))
    stop(
      "Package \"gridExtra\" must be installed to use traceplot2.",
      call. = FALSE
    )

  if (name %in% c("fe", "beta")) {
    n <- length(ngme$beta)
    ps <- lapply(1:n, function(x) traceplot(ngme, name=0, param="beta", param_index = x))
  } else if (name %in% c("mn", "noise")) {
    ps <- switch(ngme$noise$noise_type,
      "normal" = {n_sigma <- length(ngme$noise$theta_sigma); lapply(1:n_sigma, function(x) traceplot(ngme, name=0, param="sigma", param_index = x))},
      { # default case
        n_mu <- length(ngme$noise$theta_mu);
        n_sigma <- length(ngme$noise$theta_sigma);
        c(
          lapply(1:n_mu, function(x) traceplot(ngme, name=0, param="mu", param_index = x)),
          lapply(1:n_sigma, function(x) traceplot(ngme, name=0, param="sigma", param_index = x)),
          list(traceplot(ngme, name=0, param="nu", param_index = 1))
        )
      }
    )
  } else {
    f_idx <- if (is.numeric(name)) name else which(names(ngme$latents) == name)
    stopifnot(length(f_idx) == 1)

    n_K <- length(ngme$latents[[f_idx]]$theta_K);
    n_mu <- length(ngme$latents[[f_idx]]$noise$theta_mu);
    n_sigma <- length(ngme$latents[[f_idx]]$noise$theta_sigma);

    ps <- switch(ngme$latents[[f_idx]]$noise_type,
      "normal" = {
        c(
          lapply(1:n_K, function(x) traceplot(ngme, name=f_idx, param="kappa", param_index = x)),
          lapply(1:n_sigma, function(x) traceplot(ngme, name=f_idx, param="sigma", param_index = x))
        )
      },
      "normal_nig" = {
        n_sig_normal <- length(ngme$latents[[f_idx]]$noise$theta_sigma_normal);
        c(
          lapply(1:n_K, function(x) traceplot(ngme, name=f_idx, param="kappa", param_index = x)),
          lapply(1:n_mu, function(x) traceplot(ngme, name=f_idx, param="mu", param_index = x)),
          lapply(1:n_sigma, function(x) traceplot(ngme, name=f_idx, param="sigma", param_index = x)),
          list(traceplot(ngme, name=f_idx, param="nu", param_index = 1)),
          lapply(1:n_sig_normal, function(x) traceplot(ngme, name=f_idx, param="sigma_normal", param_index = x))
        )
      },
      { # default case
        c(
          lapply(1:n_K, function(x) traceplot(ngme, name=f_idx, param="kappa", param_index = x)),
          lapply(1:n_mu, function(x) traceplot(ngme, name=f_idx, param="mu", param_index = x)),
          lapply(1:n_sigma, function(x) traceplot(ngme, name=f_idx, param="sigma", param_index = x)),
          list(traceplot(ngme, name=f_idx, param="nu", param_index = 1))
        )
      }
    )
  }

  if (length(ps) > 1) ps["ncol"]=2
  do.call(gridExtra::grid.arrange, ps)
}