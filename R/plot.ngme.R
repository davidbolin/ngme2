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
    if (noise$noise == "normal_nig") {
      if (is_stationary(noise$B_sigma_nig)) {
        name_sigma_nig <- "sigma_nig"
        trans_sigma_nig <- list(exp)
      } else {
        name_sigma_nig <- paste("theta_sigma_nig", seq_len(noise$n_theta_sigma_nig))
        trans_sigma_nig <- rep(list(identity), noise$n_theta_sigma_nig)
      }
      if (is_stationary(noise$B_sigma_normal)) {
        name_sigma_normal <- "sigma_normal"
        trans_sigma_normal <- list(exp)
      } else {
        name_sigma_normal <- paste("theta_sigma_normal", seq_len(noise$n_theta_sigma_normal))
        trans_sigma_normal <- rep(list(identity), noise$n_theta_sigma_normal)
      }
      name_sigma <- c(name_sigma_nig, name_sigma_normal)
      trans_sigma <- c(trans_sigma_nig, trans_sigma_normal)
    } else {
      if (is_stationary(noise$B_sigma)) {
        name_sigma <- "sigma"
        trans_sigma <- list(exp)
      } else {
        name_sigma <- paste("theta_sigma", seq_len(noise$n_theta_sigma))
        trans_sigma <- rep(list(identity), noise$n_theta_sigma)
      }
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

if (noise$fix_theta_mu) name_mu <- trans_mu <- NULL
if (noise$fix_theta_sigma) name_sigma <- trans_sigma <- NULL
if (noise$fix_theta_nu) name_nu <- trans_nu <- NULL

    ts <- list(
      # for bv noise
      all = list(
        name_mu = name_mu,
        name_sigma = name_sigma,
        name_nu = name_nu,
        trans_mu = trans_mu,
        trans_sigma = trans_sigma,
        trans_nu = trans_nu
      ),
      # for plot
      name = c(name_mu, name_sigma, name_nu),
      trans = c(trans_mu, trans_sigma, trans_nu)
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
        n1$name_nu, n2$name_nu),
      trans = c(n1$trans_mu, n2$trans_mu,
        n1$trans_sigma, n2$trans_sigma,
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
#' @param moving_window moving window for the traceplot
#' @param hline vector, add hline to each plot
#' @param combine bool, if TRUE, plot the combination of all plots, otherwise 1 by 1
#'
#' @return the traceplot
#' @export
#'
traceplot <- function(
  ngme, 
  name="general", 
  moving_window=1,
  hline=NULL,
  combine=TRUE
) {
  stopifnot(inherits(ngme, "ngme"))
  stopifnot(!is.null(name))
  stopifnot(requireNamespace("dplyr", quietly = TRUE))
  stopifnot(requireNamespace("tidyr", quietly = TRUE))
  ngme <- ngme$replicates[[1]]
  ps <- list()

  if (name %in% names(ngme$models)) {
    # Plot trajectory of parameters of the model
    traj <- attr(ngme$models[[name]], "lat_traj")
    stopifnot("Please run ngme() to estimate the model before using traceplot()"
      = !is.null(traj))
    ts <- get_latent_info(ngme$models[[name]])
  } else {
    # Plot trajectory of parameters of the noise
    traj <- attr(ngme, "block_traj")
    stopifnot("Please run ngme() to estimate the model before using traceplot()"
      = !is.null(traj))
    # get titles
    ts <- get_noise_info(ngme$noise)
    name_feff <- if (length(ngme$feff)==0) NULL else paste ("fixed effect", seq_len(length(ngme$feff)))
    trans_feff <- rep(list(identity), length(ngme$feff))
    ts$name <- c(ts$name, name_feff)
    ts$trans <- c(ts$trans, trans_feff)
  }

  # record the geom_lines for comparison later
  avg_lines <- NULL

  n_parameters <- nrow(traj[[1]]) # number of parameters to draw
  for (idx in seq_len(n_parameters)) {
    # extract the idx-th parameter from each chain
    df <- lapply(traj, function(x) x[idx, ,drop=F])
    df <- lapply(df, as.numeric)
    
    # if not all chains are of the same length, use the minimum length
    lengths <- sapply(df, length)
    if (length(unique(lengths)) > 1) {
      min_length <- min(lengths)
      df <- lapply(df, function(x) x[seq_len(min_length)])
      warning("Some chains are of different lengths. Only the minimum length is used.")
    }

    df <- as.data.frame(df)
    df$x <- seq_len(nrow(df)); x <- NULL # get around check note
    df_long <- tidyr::gather(df, key = "key", value = "value", -x)
    ff <- ts$trans[[idx]]
    df_long$value <- ff(df_long$value)
    
    # update df using moving window
    df_long <- dplyr::group_by(df_long, key)
    if (moving_window > 1) {
      if (!requireNamespace("zoo", quietly = TRUE)) {
        message("Package 'zoo' is not installed. Using original data without moving average.")
        df_long$moving_avg <- df_long$value
      } else {
        df_long <- dplyr::mutate(df_long, moving_avg = zoo::rollapply(value, width = moving_window, FUN = mean, align = "right", fill = NA))
      }
    } else {
      df_long$moving_avg <- df_long$value
    }
    df_long <- na.omit(df_long)
    
    df_mean <- dplyr::summarise(
      dplyr::group_by(df_long, x),
      mean_moving_avg = mean(moving_avg)
    )

    ps[[idx]] <- ggplot() +
      geom_line(
        data = df_long,
        mapping = aes(
          x=.data[["x"]],
          y=.data[["moving_avg"]],
          group=.data[["key"]]
        )
      ) + geom_line(
        data = df_mean,
        aes(x=x, y=mean_moving_avg),
        col="red"
      ) + geom_hline(
        yintercept=hline[[idx]], color="blue"
      ) + labs(title = ts$name[[idx]]) +
      xlab(NULL) + ylab(NULL)

    avg_lines[[ts$name[[idx]]]] <- df_mean$mean_moving_avg
  }

  # Display results based on combine parameter
  result <- NULL
  if (combine) {
    if (length(ps) > 1) ps["ncol"] <- 2
    result <- do.call(gridExtra::grid.arrange, ps)
  } else {
    # Print plots one by one
    for (p in ps) {
      print(p)
    }
    result <- ps
  }

  attr(result, "avg_lines") <- avg_lines
  
  # print estimates of last iteration
  cat("Last estimates:\n")
  last_estimates <- lapply(avg_lines, function(x) x[length(x)])
  print(last_estimates)

  invisible(result)
}



#' Plot the density of one or more stationary noise objects
#'
#' This function plots the probability density function for one or more stationary noise objects
#' (e.g., NIG, GAL, or normal noise). Multiple noise objects can be compared on the same plot.
#'
#' @param x An ngme_noise object (required).
#' @param ... Additional ngme_noise objects to plot, or plotting parameters such as \code{xlim}.
#'   Named arguments will be used as legend labels.
#'
#' @return A ggplot object showing the density curves for the provided noise objects.
#' @export
#'
#' @examples
#' plot(noise_nig(mu=1, sigma=2, nu=1))
#' plot(n1 = noise_nig(mu=0, sigma=1, nu=1), n2 = noise_nig(mu=1, sigma=1.5, nu=0.5))
plot.ngme_noise <- function(x = NULL, ...) {
  # Get all arguments including the first one
  call_args <- as.list(match.call())[-1]  # Remove function name
  all_args <- c(if(!is.null(x)) list(x), list(...))
  
  # Get names from the call
  call_names <- names(call_args)
  if (is.null(call_names)) {
    call_names <- rep("", length(call_args))
  }
  
  # Helper function to check if an object is a noise object
  is_noise_object <- function(obj) {
    is.list(obj) && 
    !is.null(obj$theta_mu) && 
    !is.null(obj$theta_sigma) && 
    !is.null(obj$noise_type)
  }
  
  # Initialize noise objects list
  noise_objects <- list()
  
  # Extract plotting parameters and noise objects
  plot_params <- c("xlim", "ylim", "main", "xlab", "ylab")
  xlim <- NULL
  unnamed_count <- 1
  
  for (i in seq_along(all_args)) {
    arg <- all_args[[i]]
    arg_name <- call_names[i]
    
    if (arg_name %in% plot_params) {
      # This is a plotting parameter
      if (arg_name == "xlim") xlim <- arg
    } else if (is_noise_object(arg)) {
      # This is a noise object
      if (is.null(arg_name) || arg_name == "") {
        # Unnamed argument
        noise_name <- paste0("noise_", unnamed_count)
        unnamed_count <- unnamed_count + 1
      } else {
        # Named argument
        noise_name <- arg_name
      }
      noise_objects[[noise_name]] <- arg
    }
  }
  
  # Check that we have at least one noise object
  if (length(noise_objects) == 0) {
    stop("No noise objects provided")
  }
  
  # Set default xlim if not provided
  if (is.null(xlim)) xlim <- c(-10, 10)
  
  xx <- seq(xlim[[1]], xlim[[2]], length = 400)
  
  # Define colors for multiple lines
  colors <- c("black", "red", "blue", "green", "purple", "orange", "brown", "pink", "gray", "cyan")
  
  # Initialize plot data
  plot_data <- data.frame()
  
  # Process each noise object
  for (i in seq_along(noise_objects)) {
    noise <- noise_objects[[i]]
    mu <- noise$theta_mu
    sigma <- exp(noise$theta_sigma)
    nu <- exp(noise$theta_nu)
    
    stopifnot("only implemented for stationary mu" = 
      length(mu) == 1 || noise$noise_type == "normal")
    stopifnot("only implemented for stationary sigma" = length(sigma) == 1)
    stopifnot("only implemented for stationary nu" = 
      length(nu) == 1 || noise$noise_type == "normal")
    
    switch(noise$noise_type,
      "nig"     = dd <- dnig(xx, -mu, mu, nu, sigma),
      "gal"     = dd <- dgal(xx, -mu, mu, nu, sigma),
      "normal"  = dd <- dnorm(xx, sd = sigma),
      stop("Plot for this type is not implemented")
    )
    
    # Create data frame for this noise object
    noise_name <- names(noise_objects)[i]
    if (is.null(noise_name) || noise_name == "") {
      noise_name <- paste("noise", i)
    }
    
    temp_data <- data.frame(
      x = xx,
      y = dd,
      noise = noise_name,
      stringsAsFactors = FALSE
    )
    
    plot_data <- rbind(plot_data, temp_data)
  }
  
  
  # Create the plot
  gg <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = noise)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = "Noise Density Plot") +
    ggplot2::theme_minimal()
  
  # Add color scale and legend handling
  if (length(noise_objects) > 1) {
    gg <- gg + 
      ggplot2::scale_color_manual(values = colors[1:length(noise_objects)]) +
      ggplot2::labs(color = "Noise Objects")
  } else {
    gg <- gg + ggplot2::theme(legend.position = "none")
  }
  
  gg
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


#' Compare noise objects using Kullback-Leibler divergence
#'
#' This function compares multiple noise objects by calculating the KLD 
#' of each noise against the first noise object (reference).
#'
#' @param x first noise object (reference)
#' @param ... additional noise objects to compare, and optional parameters
#' @param xlim x-axis range for evaluation (default: c(-10, 10))
#' @param n_points number of evaluation points (default: 1000)
#'
#' @return named vector of KLD values
#' @export
#'
#' @examples
#' n1 <- noise_nig(mu=0, sigma=1, nu=1)
#' n2 <- noise_nig(mu=0.5, sigma=1.2, nu=0.8)
#' compare_noise_kld(n1, method2=n2)
compare_noise_kld <- function(x = NULL, ..., xlim = c(-10, 10), n_points = 1000) {
  # Get all arguments including the first one
  call_args <- as.list(match.call())[-1]  # Remove function name
  all_args <- c(if(!is.null(x)) list(x), list(...))
  
  # Get names from the call
  call_names <- names(call_args)
  if (is.null(call_names)) {
    call_names <- rep("", length(call_args))
  }
  
  # Helper function to check if an object is a noise object
  is_noise_object <- function(obj) {
    is.list(obj) && 
    !is.null(obj$theta_mu) && 
    !is.null(obj$theta_sigma) && 
    !is.null(obj$noise_type)
  }
  
  # Initialize noise objects list
  noise_objects <- list()
  
  # Extract parameters and noise objects
  param_names <- c("xlim", "n_points")
  unnamed_count <- 1
  
  for (i in seq_along(all_args)) {
    arg <- all_args[[i]]
    arg_name <- call_names[i]
    
    if (arg_name %in% param_names) {
      # This is a function parameter - already handled by function signature
      next
    } else if (is_noise_object(arg)) {
      # This is a noise object
      if (is.null(arg_name) || arg_name == "") {
        # Unnamed argument
        noise_name <- paste0("noise_", unnamed_count)
        unnamed_count <- unnamed_count + 1
      } else {
        # Named argument
        noise_name <- arg_name
      }
      noise_objects[[noise_name]] <- arg
    }
  }
  
  # Check that we have at least two noise objects
  if (length(noise_objects) < 2) {
    stop("Need at least two noise objects to compare")
  }
  
  # Generate evaluation points
  xx <- seq(xlim[[1]], xlim[[2]], length = n_points)
  
  # Calculate densities for all noise objects
  densities <- list()
  for (i in seq_along(noise_objects)) {
    noise <- noise_objects[[i]]
    mu <- noise$theta_mu
    sigma <- exp(noise$theta_sigma)
    nu <- exp(noise$theta_nu)
    
    stopifnot("only implemented for stationary mu" = 
      length(mu) == 1 || noise$noise_type == "normal")
    stopifnot("only implemented for stationary sigma" = length(sigma) == 1)
    stopifnot("only implemented for stationary nu" = 
      length(nu) == 1 || noise$noise_type == "normal")
    
    switch(noise$noise_type,
      "nig"     = dd <- dnig(xx, -mu, mu, nu, sigma),
      "gal"     = dd <- dgal(xx, -mu, mu, nu, sigma),
      "normal"  = dd <- dnorm(xx, sd = sigma),
      stop("KLD comparison for this noise type is not implemented")
    )
    
    densities[[i]] <- dd
  }
  
  # Calculate KLD of each noise against the first one (reference)
  reference_density <- densities[[1]]
  kld_values <- numeric(length(noise_objects) - 1)
  names(kld_values) <- names(noise_objects)[-1]
  
  for (i in 2:length(noise_objects)) {
    comparison_density <- densities[[i]]
    
    # Add small epsilon to avoid log(0) and division by 0
    epsilon <- 1e-10
    p <- pmax(reference_density, epsilon)
    q <- pmax(comparison_density, epsilon)
    
    # Calculate KLD: sum(p * log(p/q)) * dx
    # Since we're using discrete points, we need to multiply by dx
    dx <- (xlim[2] - xlim[1]) / (n_points - 1)
    kld_values[i-1] <- sum(p * log(p/q)) * dx
  }
  
  # Create results
  result <- list(
    kld_values = kld_values,
    reference = names(noise_objects)[1],
    n_comparisons = length(kld_values),
    closest = names(kld_values)[which.min(kld_values)]
  )
  
  class(result) <- "noise_kld_comparison"
  return(result)
}


#' Print method for noise_kld_comparison
#' @param x noise_kld_comparison object
#' @param ... additional arguments
#' @export
print.noise_kld_comparison <- function(x, ...) {
  cat("Noise KLD Comparison\n")
  cat("===================\n")
  cat("Reference:", x$reference, "\n\n")
  cat("KLD values (lower is closer to reference):\n")
  
  # Sort by KLD value for better display
  sorted_kld <- sort(x$kld_values)
  for (i in seq_along(sorted_kld)) {
    cat(sprintf("  %s: %.6f", names(sorted_kld)[i], sorted_kld[i]))
    if (names(sorted_kld)[i] == x$closest) {
      cat(" <- CLOSEST")
    }
    cat("\n")
  }
  
  cat("\nClosest to reference:", x$closest, "\n")
}