#' Trace plot of ngme fitting
#'
#' @param ngme ngme object
#' @param name name of latent models, otherwise plot fixed effects and measurement noise
#' should be in names(ngme$models) or other
#'
#' @return the traceplot
#' @export
#'
traceplot <- function(ngme, name="general") {
  stopifnot(inherits(ngme, "ngme"))
  stopifnot(!is.null(name))
  ngme <- ngme$replicates[[1]]
  ps <- list()

  if (name %in% names(ngme$models)) {
    traj <- attr(ngme$models[[name]], "lat_traj")
    # get titles
    ts <- with(ngme$models[[name]], {
      if (noise$noise_type[[1]] == "normal_nig") {
       list(
          c("theta_mu", "theta_sigma", "theta_sigma_normal", "log(nu)"),
          c(noise$n_theta_mu, noise$n_theta_sigma, noise$n_theta_sigma_normal, noise$n_nu)
        )
      } else {
        list(
          c("theta_mu", "theta_sigma", "log(nu)"),
          c(noise$n_theta_mu, noise$n_theta_sigma, noise$n_nu)
        )
      }
  })
    if (length(ngme$models[[name]]$theta_K) > 0) {
      ts[[1]] <- c("theta_K", ts[[1]])
      ts[[2]] <- c(length(ngme$models[[name]]$theta_K), ts[[2]])
    }
  } else {
    # block
    traj <- attr(ngme, "block_traj")
    # get titles
    ts <- with(ngme, {
      switch(noise$noise_type,
        "normal_nig" = list(
          c("theta_mu", "theta_sigma", "theta_sigma_normal", "log(nu)"),
          c(noise$n_theta_mu,
          noise$n_theta_sigma, noise$n_theta_sigma_normal, noise$n_nu)
        ),
        list(
          c("theta_mu", "theta_sigma", "log(nu)"),
          c(noise$n_theta_mu, noise$n_theta_sigma, noise$n_nu)
        )
      )
    })
    if (ngme$noise$corr_measurement) {
      ts[[1]] <- c(ts[[1]], "theta_rho"); ts[[2]] <- c(ts[[2]], 1)
    }
    ts[[1]] <- c(ts[[1]], "beta"); ts[[2]] <- c(ts[[2]], length(ngme$beta))
  }
  names <- c()
  for (i in seq_along(ts[[1]]))
    names <- c(names, paste(ts[[1]][[i]], seq_len(ts[[2]][[i]])))
  # make plot
  for (idx in seq_len(nrow(traj[[1]]))) {
    df <- sapply(traj, function(x) x[idx, ,drop=F])
    # wierd stuff here
    df <- apply(df, c(1,2), as.numeric)
    df <- as.data.frame(df)
    # browser()
    mean_traj <- rowMeans(df)
    df$x <- seq_len(nrow(df))
    df_long <- tidyr::gather(df, key = "key", value = "value", -x)
    ps[[idx]] <- ggplot() +
      geom_line(data = df_long, aes(x=x, y=value, group=key)) +
      geom_line(data = data.frame(x=seq_len(nrow(df)), y=mean_traj), aes(x=x,y=y), col="red") +
      labs(title = names[idx])
  }

  if (length(ps) > 1) ps["ncol"]=2
  do.call(gridExtra::grid.arrange, ps)
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

  xlim <- if (!is.null(list(...)$xlim)) list(...)$xlim else c(-10, 10)

  xx <- seq(xlim[[1]], xlim[[2]], length = 400)
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
