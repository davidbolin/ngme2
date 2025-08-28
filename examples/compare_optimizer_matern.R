load_all()
# library(ngme2)

# Helper function to create comparison plot
create_comparison_plot <- function(log_like_trajectories, filename, true_likelihood = NULL) {
  n_optimizers <- length(log_like_trajectories)
  n_iter <- length(log_like_trajectories[[1]])
  
  # Define colors for different optimizers
  colors <- c("blue", "red", "green", "purple", "orange", "brown", "pink", "gray")
  if (n_optimizers > length(colors)) {
    colors <- rainbow(n_optimizers)
  }
  
  # Determine y-axis range
  all_values <- unlist(log_like_trajectories)
  if (!is.null(true_likelihood)) {
    all_values <- c(all_values, true_likelihood)
  }
  y_range <- range(all_values)
  
  # Create the plot
  png(filename, width = 800, height = 600)
  
  # Plot first optimizer
  plot(1:n_iter, log_like_trajectories[[1]], 
       type = "l", 
       xlab = "Iteration", 
       ylab = "Complete Log-Likelihood",
       main = "Complete Log-Likelihood vs Iteration: Optimizer Comparison",
       col = colors[1],
       lwd = 2,
       ylim = y_range)
  
  # Add lines for other optimizers
  for (i in 2:n_optimizers) {
    lines(1:n_iter, log_like_trajectories[[i]], 
          col = colors[i], 
          lwd = 2)
  }
  
  # Add true likelihood horizontal line if provided
  if (!is.null(true_likelihood)) {
    abline(h = true_likelihood, col = "black", lty = 2, lwd = 2)
  }
  
  # Add legend
  legend_names <- names(log_like_trajectories)
  legend_colors <- colors[1:n_optimizers]
  legend_lty <- rep(1, n_optimizers)
  legend_lwd <- rep(2, n_optimizers)
  
  if (!is.null(true_likelihood)) {
    legend_names <- c(legend_names, "True Likelihood")
    legend_colors <- c(legend_colors, "black")
    legend_lty <- c(legend_lty, 2)
    legend_lwd <- c(legend_lwd, 2)
  }
  
  legend("bottomright", 
         legend = legend_names, 
         col = legend_colors, 
         lty = legend_lty,
         lwd = legend_lwd,
         bty = "n")
  
  # Add grid
  grid()
  dev.off()
  
  # Also display in current device
  plot(1:n_iter, log_like_trajectories[[1]], 
       type = "l", 
       xlab = "Iteration", 
       ylab = "Complete Log-Likelihood",
       main = "Complete Log-Likelihood vs Iteration: Optimizer Comparison",
       col = colors[1],
       lwd = 2,
       ylim = y_range)
  
  for (i in 2:n_optimizers) {
    lines(1:n_iter, log_like_trajectories[[i]], 
          col = colors[i], 
          lwd = 2)
  }
  
  # Add true likelihood horizontal line if provided
  if (!is.null(true_likelihood)) {
    abline(h = true_likelihood, col = "black", lty = 2, lwd = 2)
  }
  
  legend("bottomright", 
         legend = legend_names, 
         col = legend_colors, 
         lty = legend_lty,
         lwd = legend_lwd,
         bty = "n")
  
  grid()
  
  cat("Plot saved as:", filename, "\n")
  if (!is.null(true_likelihood)) {
    cat("True likelihood value:", true_likelihood, "\n")
  }
}

# Helper function to create parameter trajectory comparison plots for Matern model
create_parameter_comparison_plots_matern <- function(optimizer_results, true_params, filename_prefix = "matern_parameter_comparison") {
  # Extract parameter trajectories from all optimizers
  param_trajectories <- list()
  
  for (opt_name in names(optimizer_results)) {
    trajectories <- ngme2:::get_trace_trajectories(optimizer_results[[opt_name]], "matern_field")
    
    # Calculate mean across chains for each iteration
    param_trajectories[[opt_name]] <- list(
      kappa = apply(trajectories$trajectories$kappa, 1, mean),
      mu = apply(trajectories$trajectories$mu, 1, mean),
      sigma = apply(trajectories$trajectories$sigma, 1, mean),
      nu = apply(trajectories$trajectories$nu, 1, mean)
    )
  }
  
  # Define colors for different optimizers
  colors <- c("blue", "red", "green", "purple", "orange", "brown", "pink", "gray")
  n_optimizers <- length(optimizer_results)
  if (n_optimizers > length(colors)) {
    colors <- rainbow(n_optimizers)
  }
  
  # Get parameter names and true values
  param_names <- c("kappa", "mu", "sigma", "nu")
  param_labels <- c("κ (Matern range)", "μ (NIG location)", "σ (NIG scale)", "ν (NIG shape)")
  
  # Create 2x2 subplot layout
  par(mfrow=c(2,2))
  
  for (i in 1:length(param_names)) {
    param_name <- param_names[i]
    param_label <- param_labels[i]
    true_value <- true_params[[param_name]]
    
    # Get all trajectories for this parameter
    all_trajectories <- lapply(param_trajectories, function(x) x[[param_name]])
    n_iter <- length(all_trajectories[[1]])
    
    # Determine y-axis range
    all_values <- unlist(all_trajectories)
    if (!is.null(true_value)) {
      all_values <- c(all_values, true_value)
    }
    y_range <- range(all_values)
    
    # Plot first optimizer
    plot(1:n_iter, all_trajectories[[1]], 
         type = "l", 
         xlab = "Iteration", 
         ylab = param_label,
         main = paste("Parameter Trajectory Comparison:", param_label),
         col = colors[1],
         lwd = 2,
         ylim = y_range)
    
    # Add lines for other optimizers
    for (j in 2:n_optimizers) {
      lines(1:n_iter, all_trajectories[[j]], 
            col = colors[j], 
            lwd = 2)
    }
    
    # Add true value horizontal line if provided
    if (!is.null(true_value)) {
      abline(h = true_value, col = "black", lty = 2, lwd = 2)
    }
    
    # Add legend for the first subplot only
    if (i == 1) {
      legend_names <- names(optimizer_results)
      if (!is.null(true_value)) {
        legend_names <- c(legend_names, "True Value")
        legend_colors <- c(colors[1:n_optimizers], "black")
        legend_lty <- c(rep(1, n_optimizers), 2)
      } else {
        legend_colors <- colors[1:n_optimizers]
        legend_lty <- rep(1, n_optimizers)
      }
      
      legend("topright", 
             legend = legend_names, 
             col = legend_colors, 
             lty = legend_lty,
             lwd = 2,
             bty = "n",
             cex = 0.8)
    }
    
    # Add grid
    grid()
  }
  
  par(mfrow=c(1,1))
  
  # Save as PNG file
  png(paste0(filename_prefix, ".png"), width = 1200, height = 900)
  par(mfrow=c(2,2))
  
  for (i in 1:length(param_names)) {
    param_name <- param_names[i]
    param_label <- param_labels[i]
    true_value <- true_params[[param_name]]
    
    # Get all trajectories for this parameter
    all_trajectories <- lapply(param_trajectories, function(x) x[[param_name]])
    n_iter <- length(all_trajectories[[1]])
    
    # Determine y-axis range
    all_values <- unlist(all_trajectories)
    if (!is.null(true_value)) {
      all_values <- c(all_values, true_value)
    }
    y_range <- range(all_values)
    
    # Plot first optimizer
    plot(1:n_iter, all_trajectories[[1]], 
         type = "l", 
         xlab = "Iteration", 
         ylab = param_label,
         main = paste("Parameter Trajectory Comparison:", param_label),
         col = colors[1],
         lwd = 2,
         ylim = y_range)
    
    # Add lines for other optimizers
    for (j in 2:n_optimizers) {
      lines(1:n_iter, all_trajectories[[j]], 
            col = colors[j], 
            lwd = 2)
    }
    
    # Add true value horizontal line if provided
    if (!is.null(true_value)) {
      abline(h = true_value, col = "black", lty = 2, lwd = 2)
    }
    
    # Add legend for the first subplot only
    if (i == 1) {
      legend_names <- names(optimizer_results)
      if (!is.null(true_value)) {
        legend_names <- c(legend_names, "True Value")
        legend_colors <- c(colors[1:n_optimizers], "black")
        legend_lty <- c(rep(1, n_optimizers), 2)
      } else {
        legend_colors <- colors[1:n_optimizers]
        legend_lty <- rep(1, n_optimizers)
      }
      
      legend("topright", 
             legend = legend_names, 
             col = legend_colors, 
             lty = legend_lty,
             lwd = 2,
             bty = "n",
             cex = 0.8)
    }
    
    # Add grid
    grid()
  }
  
  par(mfrow=c(1,1))
  dev.off()
  
  cat("Parameter comparison plots saved as:", paste0(filename_prefix, ".png"), "\n")
}

# Helper function to print comparison summary
print_comparison_summary <- function(log_like_trajectories) {
  cat("\n=== Complete Log-Likelihood Comparison Summary ===\n")
  
  for (opt_name in names(log_like_trajectories)) {
    traj <- log_like_trajectories[[opt_name]]
    n_iter <- length(traj)
    
    cat("\n", opt_name, ":\n")
    cat("  Initial value:", traj[1], "\n")
    cat("  Final value:", traj[n_iter], "\n")
    cat("  Maximum value:", max(traj), "at iteration", which.max(traj), "\n")
    cat("  Minimum value:", min(traj), "at iteration", which.min(traj), "\n")
    cat("  Improvement:", traj[n_iter] - traj[1], "\n")
  }
  
  # Find best performing optimizer
  final_values <- sapply(log_like_trajectories, function(x) x[length(x)])
  best_optimizer <- names(which.max(final_values))
  cat("\nBest performing optimizer:", best_optimizer, "with final log-likelihood:", max(final_values), "\n")
}

# Function to compare optimizers on Matern model
compare_optimizers_matern <- function(optimizer_list, 
                                     n_obs = 400,  # Smaller for Matern due to computational cost
                                     iterations = 10, 
                                     burnin = 100,
                                     n_parallel_chain = 4,
                                     seed = 500,
                                     plot_filename = "matern_optimizer_comparison.png",
                                     show_traceplots = TRUE,
                                     kappa = 5.0,
                                     alpha = 2,
                                     n_gibbs_samples = 1) {
  
  set.seed(seed)
  
  # Generate 2D spatial data for Matern model
  n_side <- sqrt(n_obs)
  if (n_side != floor(n_side)) {
    n_side <- floor(sqrt(n_obs))
    n_obs <- n_side^2
    cat("Adjusted n_obs to", n_obs, "for square grid\n")
  }
  
  # Create spatial locations as matrix (required for Matern)
  locs <- as.matrix(expand.grid(
    x = seq(0, 1, length.out = n_side),
    y = seq(0, 1, length.out = n_side)
  ))
  
  # Create mesh for Matern model
  mesh <- fmesher::fm_mesh_2d(locs, cutoff = 0.05, max.edge = c(0.1, 0.5))
  cat("Mesh size:", mesh$n, "\n")
  
  # Generate data
  sigma_eps <- 0.05
  mu <- 1.0
  sigma <- 2.0
  nu <- 1.0
  
  # Create Matern model
  matern_model <- f(
    map = locs, 
    model = "matern", 
    kappa = kappa, 
    alpha = alpha,
    mesh = mesh,
    noise = noise_nig(mu = mu, sigma = sigma, nu = nu)
  )
  A <- matern_model$A
  
  tmp <- simulate(matern_model, seed = seed, nsim=1)
  trueW <- attributes(tmp)$W_sim[[1]]
  Y <- as.vector(A %*% trueW) + rnorm(n_obs, mean=0, sd=sigma_eps)
  
  # Generate V and true W for Matern
  trueV <- attributes(tmp)$V_sim[[1]]
  
  # Debug: Check dimensions
  cat("Mesh size:", mesh$n, "\n")
  cat("Locations size:", nrow(locs), "\n")
  cat("trueW length:", length(trueW), "\n")
  cat("Y length:", length(Y), "\n")
  
  # Plot the spatial data
  par(mfrow=c(1,2))
  # Plot true field (projected to observation locations)
  AW_plot <- as.vector(A %*% trueW)
  image(matrix(AW_plot, n_side, n_side), main="True Spatial Field (A*W)", 
        xlab="X", ylab="Y", col=heat.colors(20))
  # Plot observations
  image(matrix(Y, n_side, n_side), main="Observed Data (Y)", 
        xlab="X", ylab="Y", col=heat.colors(20))
  par(mfrow=c(1,1))
  
  # Store results
  optimizer_results <- list()
  log_like_trajectories <- list()
  
  # Fit each optimizer
  for (opt_name in names(optimizer_list)) {
    cat("\n=== Fitting Matern model with", opt_name, "optimizer ===\n")
    if (opt_name %in% c("precond_same_V", "precond_diff_chian_same_V")) {
      use_same_V <- TRUE
    } else {
      use_same_V <- FALSE
    }
    
    # Create control for this optimizer
    control_opt <- control_opt(
      optimizer = optimizer_list[[opt_name]],
      burnin = burnin,
      iterations = iterations,
      stop_points = iterations,
      n_parallel_chain = n_parallel_chain,
      verbose = FALSE,
      rao_blackwellization = TRUE,
      seed = seed,
      print_check_info = TRUE,
      std_lim = 0.01,
      trend_lim = 0.01,
      n_slope_check = 3
    )
    
    # Fit the model
    result <- ngme(
      Y ~ 0 + f(
        map = locs,
        alpha = alpha,
        mesh = mesh,
        name = "matern_field",
        model = "matern",
        noise = noise_nig(),
        control = control_f(use_same_V = use_same_V)
      ),
      data = data.frame(Y=Y),
      control_opt = control_opt,
      control_ngme = control_ngme(
        n_gibbs_samples = n_gibbs_samples
      )
    )
    
    print(result)
    
    if (show_traceplots) {
      traceplot(result, "matern_field", hline=c(kappa, mu, sigma, nu))
    }
    
    # Store results
    optimizer_results[[opt_name]] <- result
    
    # Get trajectories and calculate log-likelihood
    trajectories <- ngme2:::get_trace_trajectories(result, "matern_field")
    data_trajectories <- ngme2:::get_trace_trajectories(result, "data")
    log_like_trajectories[[opt_name]] <- calc_log_like_trajectory_matern(trajectories, data_trajectories, Y, trueW, trueV, mesh, alpha, A)
  }
  
  # Calculate true complete log-likelihood with true parameters
  true_complete_likelihood <- complete_log_like_matern(Y, trueW, trueV, kappa, mu, sigma, nu, sigma_eps, mesh, alpha, A)
  cat("\nTrue complete log-likelihood:", true_complete_likelihood, "\n")
  
  # Create comparison plot
  create_comparison_plot(log_like_trajectories, plot_filename, true_complete_likelihood)
  
  # Create combined traceplot for all optimizers
  if (show_traceplots) {
    par(mfrow=c(2,2))
    for (opt_name in names(optimizer_results)) {
      cat("Creating traceplot for:", opt_name, "\n")
      traceplot(optimizer_results[[opt_name]], "matern_field", hline=c(kappa, mu, sigma, nu))
      # Add title manually
      title(main=paste("Traceplot:", opt_name))
    }
    par(mfrow=c(1,1))
  }
  
  # Create parameter trajectory comparison plots
  true_params <- list(kappa = kappa, mu = mu, sigma = sigma, nu = nu)
  create_parameter_comparison_plots_matern(optimizer_results, true_params, "matern_parameter_comparison")
  
  # Print summary statistics
  print_comparison_summary(log_like_trajectories)
  
  return(list(
    results = optimizer_results,
    trajectories = log_like_trajectories,
    data = list(Y = Y, trueW = trueW, trueV = trueV, sigma_eps = sigma_eps, locs = locs)
  ))
}

# Helper function to calculate log-likelihood trajectory for Matern model
calc_log_like_trajectory_matern <- function(ret_obj, data_traj, Y, trueW, trueV, mesh, alpha, A) {
  # Extract trajectories for all parameters
  kappa_traj <- ret_obj$trajectories$kappa
  mu_traj <- ret_obj$trajectories$mu  
  sigma_traj <- ret_obj$trajectories$sigma
  nu_traj <- ret_obj$trajectories$nu
  
  # Extract sigma_eps trajectory from data
  sigma_eps_traj <- data_traj$trajectories$sigma
  
  # Get the number of iterations
  n_iter <- nrow(kappa_traj)
  
  # Calculate complete log-likelihood for each iteration using mean across chains
  log_like_values <- numeric(n_iter)
  
  for (i in 1:n_iter) {
    # Get mean parameter values across chains for this iteration
    kappa_mean <- mean(kappa_traj[i, ])
    mu_mean <- mean(mu_traj[i, ])
    sigma_mean <- mean(sigma_traj[i, ])
    nu_mean <- mean(nu_traj[i, ])
    sigma_eps_mean <- mean(sigma_eps_traj[i, ])
    
    # Use the true W and V values for the likelihood calculation
    log_like_values[i] <- complete_log_like_matern(Y, trueW, trueV, kappa_mean, mu_mean, sigma_mean, nu_mean, sigma_eps_mean, mesh, alpha, A)
  }
  
  return(log_like_values)
}

# Complete log-likelihood function for Matern model
complete_log_like_matern <- function(Y, W, V, kappa, mu, sigma, nu, sigma_eps, mesh, alpha, A) {
  n <- mesh$n  # Mesh size
  
  # For Matern model, W has mesh dimensions, but we need to work with observation locations
  AW <- as.vector(A %*% W)
  
  K <- matern(mesh=mesh, kappa=kappa, alpha=alpha)$K
  
  # h vector (integration weights) - for regular grid, h = rep(1, n)
  h <- matern(mesh=mesh, kappa=kappa, alpha=alpha)$h
  
  # Log p(V) - unchanged
  l_V <- sum(dig(V, nu, nu * h^2, log=TRUE))
  
  # New log p(W|V) formula:
  # Compute Kw
  Kw <- as.vector(K %*% W)
  
  # Compute mu*(V-h)
  mu_V_minus_h <- mu * (V - h)
  
  # Compute residuals: Kw - mu*(V-h)
  residuals <- Kw - mu_V_minus_h
  
  # Compute the quadratic form: sum((Kw - mu*(V-h))_i^2 / (sigma^2 * V_i))
  quadratic_form <- sum(residuals^2 / (sigma^2 * V))
  
  # Compute log determinant of K
  # For numerical stability, we use the determinant function
  log_det_K <- as.numeric(determinant(as.matrix(K), logarithm = TRUE)$modulus)
  
  # Assemble the log p(W|V) term
  l_W <- -0.5 * (n * log(2 * pi) - 2 * log_det_K + n * log(sigma^2) + sum(log(V)) + quadratic_form)
  
  # Log p(Y|W) - unchanged  
  l_Y <- sum(dnorm(Y, mean=AW, sd=sigma_eps, log=TRUE))
  
  return(l_V + l_W + l_Y)
}

# Example usage: Compare Matern model  
cat("=== Matern Model Comparison ===\n")
optimizer_list_matern <- list(
  precond_same_V = precond_sgd(
    preconditioner="full"
  ),
  precond_diff_V = precond_sgd(
    preconditioner="full"
  ),
  precond_diff_chian_same_V = precond_sgd(
    preconditioner="full",
    precond_by_diff_chain = TRUE
  )
  # adam = adam()
)

seeds = 30
n_obss = c(500)
n_gibbs_samples = c(4)
for (seed in seeds) {
  for (n_obs in n_obss) {
    for (n_gibbs_samples in n_gibbs_samples) {
      results_matern <- compare_optimizers_matern(
        optimizer_list = optimizer_list_matern,
        n_obs = n_obs,
        iterations = 300,
        seed = seed,
        show_traceplots = TRUE,
        plot_filename = paste0("matern_obs_", n_obs, "_gibbs_", n_gibbs_samples, "_", seed, ".png"),
        kappa = 5.0,
        alpha = 2,
        n_gibbs_samples = n_gibbs_samples
      )
    }
  }
}