get_noise_info <- function(noise) {
  if (length(noise$noise_type) == 1) {
    # mu
    if (noise$n_theta_mu == 0) {
      name_mu <- NULL
    } else if (noise$n_theta_mu == 1 && is_stationary(noise$B_mu)) {
      name_mu <- "mu"
    }
    else {
      name_mu <- paste("theta_mu", seq_len(noise$n_theta_mu))
    }
    trans_mu <- rep(list(identity), noise$n_theta_mu)

    # sigma
    if (is_stationary(noise$B_sigma)) {
      name_sigma <- "sigma"
      trans_sigma <- list(exp)
    } else {
      name_sigma <- paste("theta_sigma", seq_len(noise$n_theta_sigma))
      trans_sigma <- rep(list(identity), noise$n_theta_sigma)
    }

    if (length(noise$theta_sigma_normal) == 0) {
      name_sigma_normal <- NULL
      trans_sigma_normal <- NULL
    } else if (length(noise$theta_sigma_normal) == 1 && is_stationary(noise$B_sigma_normal)) {
      name_sigma_normal <- "sigma_normal"
      trans_sigma_normal <- list(exp)
    } else {
      name_sigma_normal <- paste("theta_sigma_normal", seq_len(noise$n_theta_sigma_normal))
      trans_sigma_normal <- rep(list(identity), noise$n_theta_sigma_normal)
    }

    # nu
    if (noise$n_theta_nu == 0) {
      name_nu <- trans_nu <- NULL
    } else if (is_stationary(noise$B_nu)) {
      name_nu <- "nu"
      trans_nu <- list(exp)
    } else {
      name_nu <- paste("theta_nu", seq_len(noise$n_theta_nu))
      trans_nu <- rep(list(identity), noise$n_theta_nu)
    }

    ts <- list(
      # for bv noise
      all = list(
        name_mu = name_mu,
        name_sigma = name_sigma,
        name_sigma_normal = name_sigma_normal,
        name_nu = name_nu,
        trans_mu = trans_mu,
        trans_sigma = trans_sigma,
        trans_sigma_normal = trans_sigma_normal,
        trans_nu = trans_nu
      ),
      # for plot
      name = c(name_mu, name_sigma, name_sigma_normal, name_nu),
      trans = c(trans_mu, trans_sigma, trans_sigma_normal, trans_nu)
    )
  } else {
    # bivariate noise
    n1 <- get_noise_info(noise$bv_noises[[1]])
    n2 <- get_noise_info(noise$bv_noises[[2]])
    n1 <- lapply(n1$all, function(x) if (is.character(x)) paste(x, "(1st)") else x)
    n2 <- lapply(n2$all, function(x) if (is.character(x)) paste(x, "(2nd)") else x)

    # re-arrange
    ts <- list(
      name = c(n1$name_mu, n2$name_mu,
        n1$name_sigma, n2$name_sigma,
        n1$name_sigma_normal, n2$name_sigma_normal,
        n1$name_nu, n2$name_nu),
      trans = c(n1$trans_mu, n2$trans_mu,
        n1$trans_sigma, n2$trans_sigma,
        n1$trans_sigma_normal, n2$trans_sigma_normal,
        n1$trans_nu, n2$trans_nu)
    )
  }
  if (noise$corr_measurement) {
    ts$name <- c(ts$name, "rho(measurement)");
    ts$trans <- c(ts$trans, list(ar1_th2a))
  }
  ts
}

get_latent_info <- function(latent) {
  ts <- get_noise_info(latent$noise)
  ts$name <- c(latent$operator$param_name, ts$name)
  ts$trans <- c(latent$operator$param_trans, ts$trans)
  ts
}

#' Trace plot of ngme fitting
#'
#' @param ngme ngme object
#' @param name name of latent models, otherwise plot fixed effects and measurement noise
#' should be in names(ngme$models) or other
#' @param hline vector, add hline to each plot
#'
#' @return the traceplot
#' @export
#'
traceplot <- function(ngme, name="general", hline=NULL) {
  stopifnot(inherits(ngme, "ngme"))
  stopifnot(!is.null(name))
  ngme <- ngme$replicates[[1]]
  ps <- list()

  if (name %in% names(ngme$models)) {
    # Plot latent trajectory
    traj <- attr(ngme$models[[name]], "lat_traj")
    ts <- get_latent_info(ngme$models[[name]])
  } else {
    # Plot block trajectory
    traj <- attr(ngme, "block_traj")
    # get titles
    ts <- get_noise_info(ngme$noise)
    name_feff <- if (length(ngme$feff)==0) NULL else paste ("fixed effect", seq_len(length(ngme$feff)))
    trans_feff <- rep(list(identity), length(ngme$feff))
    ts$name <- c(ts$name, name_feff)
    ts$trans <- c(ts$trans, trans_feff)
  }

  # record the geom_lines for comparison later
  avg_lines <- NULL

  for (idx in seq_len(nrow(traj[[1]]))) {
    df <- sapply(traj, function(x) x[idx, ,drop=F])
    # weird stuff here
    df <- apply(df, c(1,2), as.numeric)
    df <- as.data.frame(df)
    mean_traj <- rowMeans(df)
    df$x <- seq_len(nrow(df)); x <- NULL # get around check note
    df_long <- tidyr::gather(df, key = "key", value = "value", -x)
    ff <- ts$trans[[idx]]
    df_long$value <- ff(df_long$value)
    ps[[idx]] <- ggplot() +
      geom_line(
        data = df_long,
        mapping = aes(
          x=.data[["x"]],
          y=.data[["value"]],
          group=.data[["key"]]
        )
      ) + geom_line(
        data = data.frame(x=seq_len(nrow(df)), y=ff(mean_traj)),
        aes(
          x=.data[["x"]],
          y=.data[["y"]]
        ),
        col="red"
      ) + geom_hline(
        yintercept=hline[[idx]], color="blue"
      ) + labs(title = ts$name[[idx]]) +
      xlab(NULL) + ylab(NULL)

    avg_lines[[ts$name[[idx]]]] <- ff(mean_traj)
  }

  if (length(ps) > 1) ps["ncol"]=2
  result <- do.call(gridExtra::grid.arrange, ps)

  attr(result, "avg_lines") <- avg_lines
  invisible(result)
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
  nu <- exp(noise$theta_nu)
  stopifnot("only implemented for stationary mu" = length(mu) == 1)
  stopifnot("only implemented for stationary sigma" = length(sigma) == 1)
  stopifnot("only implemented for stationary nu" = length(nu) == 1)

  xlim <- if (!is.null(list(...)$xlim)) list(...)$xlim else c(-10, 10)

  xx <- seq(xlim[[1]], xlim[[2]], length = 400)
  switch(noise$noise_type,
    "nig"     = dd <- dnig(xx, -mu, mu, nu, sigma),
    "gal"     = dd <- dgal(xx, -mu, mu, nu, sigma),
    "normal"  = dd <- dnorm(xx, sd = sigma),
    stop("Plot for this type is not implemented")
  )

  gg <- ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = xx, y = dd))

  if (!is.null(noise2)) {
    mu <- noise2$theta_mu
    sigma <- exp(noise2$theta_sigma)
    nu <- exp(noise$theta_nu)
    switch(noise2$noise_type,
      "nig"     = dd2 <- dnig(xx, -mu, mu, nu, sigma),
      "gal"     = dd2 <- dgal(xx, -mu, mu, nu, sigma),
      "normal"  = dd2 <- dnorm(xx, sd = sigma),
      stop("Plot for this type is not implemented")
    )
    gg <- gg + ggplot2::geom_line(ggplot2::aes(x = xx, y = dd2), col = "red")
  }

  gg + ggplot2::labs(title = "Density Plot")
}


compare_traceplot <- function(l1, l2) {
  l1 <- as.data.frame(l1);
  l2 <- as.data.frame(l2)

  ps <- list()
  n_plots = length(l1)
  n_iter = length(l1[[1]])

  for (i in seq_len(n_plots)) {
    c1 <- c2 <- NULL
    df <- data.frame(c1 = l1[[i]], c2 = l2[[i]], title=names(l1)[[i]])
    ps[[i]] <- ggplot(data=df) +
      geom_line(aes(x=1:n_iter, y=c1), col="1") +
      geom_line(aes(x=1:n_iter, y=c2), col="2") +
      labs(title = df$title) +
      xlab(NULL) + ylab(NULL)
  }

  do.call(gridExtra::grid.arrange, ps)
}