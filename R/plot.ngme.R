get_noise_info <- function(noise) {
  if (length(noise$noise_type) == 1) {
    # mu
    if (noise$n_theta_mu == 0) {
      name_mu <- NULL
    } else if (noise$n_theta_mu == 1 && is_stationary(noise$B_mu)) {
      name_mu <- "mu"
    } else {
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
      if (is_stationary(noise$B_sigma) ||
        (ncol(noise$B_sigma) == 2 && noise$B_sigma[1, 1] == 1 && all(noise$B_sigma[-1, 1] == 0)) # RW1 pattern
      ) {
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
      nu_lower_bound <- if (is.null(noise$nu_lower_bound)) 0 else noise$nu_lower_bound
      nu_trans <- function(theta_nu) {
        nu_lower_bound + exp(theta_nu)
      }
      name_nu <- "nu"
      trans_nu <- list(nu_trans)
    } else {
      name_nu <- paste("theta_nu", seq_len(noise$n_theta_nu))
      trans_nu <- rep(list(identity), noise$n_theta_nu)
    }

    if (noise$fix_theta_mu) name_mu <- trans_mu <- NULL
    if (all(noise$fix_theta_sigma)) {
      name_sigma <- trans_sigma <- NULL
    }
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
      name = c(
        n1$name_mu, n2$name_mu,
        n1$name_sigma, n2$name_sigma,
        n1$name_nu, n2$name_nu
      ),
      trans = c(
        n1$trans_mu, n2$trans_mu,
        n1$trans_sigma, n2$trans_sigma,
        n1$trans_nu, n2$trans_nu
      )
    )
  }
  if (noise$corr_measurement) {
    ts$name <- c(ts$name, "rho(measurement)")
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

get_feff_plot_names <- function(feff) {
  n_feff <- length(feff)
  if (n_feff == 0) {
    return(NULL)
  }

  feff_names <- names(feff)
  if (is.null(feff_names) || anyNA(feff_names) || any(!nzchar(feff_names))) {
    return(paste("fixed effect", seq_len(n_feff)))
  }

  feff_names
}

disambiguate_trace_param_names <- function(param_names) {
  if (length(param_names) == 0) {
    return(param_names)
  }
  make.unique(as.character(param_names))
}

#' Get trace trajectories from ngme fitting
#'
#' Extract numerical trajectories of parameters during ngme optimization
#' without creating plots. This function provides the underlying data
#' used by traceplot().
#'
#' @param ngme ngme object
#' @param name name of latent models, otherwise get fixed effects and measurement noise
#' should be in names(ngme$models) or other
#' @param apply_transform logical, whether to apply parameter transformations (default: TRUE)
#'
#' @return list containing:
#' \describe{
#'   \item{trajectories}{named list of trajectory matrices for each parameter}
#'   \item{parameter_names}{character vector of parameter names}
#'   \item{transformations}{list of transformation functions used}
#'   \item{n_chains}{number of chains}
#'   \item{n_iterations}{number of iterations}
#' }
#' @export
#'
get_trace_trajectories <- function(
    ngme,
    name = "general",
    apply_transform = TRUE) {
  stopifnot(inherits(ngme, "ngme"))
  stopifnot(!is.null(name))

  ngme <- ngme$replicates[[1]]

  if (name %in% names(ngme$models)) {
    # Get trajectory of parameters of the model
    traj <- attr(ngme$models[[name]], "lat_traj")
    stopifnot(
      "Please run ngme() to estimate the model before using get_trace_trajectories()" = !is.null(traj)
    )
    ts <- get_latent_info(ngme$models[[name]])
  } else {
    # Get trajectory of parameters of the noise
    traj <- attr(ngme, "block_traj")
    stopifnot(
      "Please run ngme() to estimate the model before using get_trace_trajectories()" = !is.null(traj)
    )
    # get titles
    ts <- get_noise_info(ngme$noise)
    name_feff <- get_feff_plot_names(ngme$feff)
    trans_feff <- rep(list(identity), length(ngme$feff))
    ts$name <- c(ts$name, name_feff)
    ts$trans <- c(ts$trans, trans_feff)
  }
  ts$name <- disambiguate_trace_param_names(ts$name)

  n_parameters <- nrow(traj[[1]])
  n_chains <- length(traj)

  # Get lengths and use minimum if chains have different lengths
  lengths <- sapply(traj, function(x) ncol(x))
  if (length(unique(lengths)) > 1) {
    min_length <- min(lengths)
    traj <- lapply(traj, function(x) x[, seq_len(min_length), drop = FALSE])
    warning("Some chains are of different lengths. Only the minimum length is used.")
  }
  n_iterations <- lengths[1]

  # Extract trajectories for each parameter
  trajectories <- list()

  for (idx in seq_len(n_parameters)) {
    # Extract the idx-th parameter from each chain
    param_data <- lapply(traj, function(x) x[idx, , drop = FALSE])
    param_data <- lapply(param_data, as.numeric)

    # Combine into matrix: rows = iterations, columns = chains
    param_matrix <- do.call(cbind, param_data)
    colnames(param_matrix) <- paste0("chain_", seq_len(n_chains))

    # Apply transformation if requested
    if (apply_transform && !is.null(ts$trans[[idx]])) {
      ff <- ts$trans[[idx]]
      param_matrix <- apply(param_matrix, 2, ff)
      if (!is.matrix(param_matrix)) {
        param_matrix <- matrix(param_matrix, nrow = n_iterations, ncol = n_chains)
        colnames(param_matrix) <- paste0("chain_", seq_len(n_chains))
      }
    }

    trajectories[[ts$name[idx]]] <- param_matrix
  }

  result <- list(
    trajectories = trajectories,
    parameter_names = ts$name,
    transformations = ts$trans,
    n_chains = n_chains,
    n_iterations = n_iterations
  )

  class(result) <- "ngme_trajectories"
  return(result)
}

#' Print method for ngme_trajectories
#' @param x ngme_trajectories object
#' @param ... additional arguments
#' @export
print.ngme_trajectories <- function(x, ...) {
  cat("NGME Parameter Trajectories\n")
  cat("==========================\n")
  cat("Parameters:", length(x$parameter_names), "\n")
  cat("Chains:", x$n_chains, "\n")
  cat("Iterations:", x$n_iterations, "\n\n")

  cat("Parameter names:\n")
  for (i in seq_along(x$parameter_names)) {
    cat(sprintf("  [%d] %s\n", i, x$parameter_names[i]))
  }

  cat("\nAccess trajectories via $trajectories$parameter_name\n")
  cat("Each trajectory is a matrix: rows = iterations, columns = chains\n")
}

#' Calculate parameter distance from true values
#'
#' Compute ||theta - theta_hat|| distance for all parameters across iterations.
#' This provides a single measure of convergence combining all parameters.
#'
#' @param trajectories ngme_trajectories object from get_trace_trajectories()
#' @param true_values named list or vector of true parameter values
#' @param norm_type character, type of norm to use: "euclidean" (default), "manhattan", "max"
#' @param chain_summary character, how to summarize across chains: "mean" (default), "median", "chain1"
#'
#' @return numeric vector of distances across iterations
#' @export
#'
get_parameter_distance <- function(
    trajectories,
    true_values,
    norm_type = "euclidean",
    chain_summary = "mean") {
  stopifnot(inherits(trajectories, "ngme_trajectories"))
  stopifnot(is.list(true_values) || is.numeric(true_values))
  stopifnot(norm_type %in% c("euclidean", "manhattan", "max"))
  stopifnot(chain_summary %in% c("mean", "median", "chain1"))

  # Match parameter names
  if (is.numeric(true_values) && is.null(names(true_values))) {
    if (length(true_values) != length(trajectories$parameter_names)) {
      stop("Length of true_values must match number of parameters")
    }
    names(true_values) <- trajectories$parameter_names
  }

  # Check all required parameters are provided
  missing_params <- setdiff(trajectories$parameter_names, names(true_values))
  if (length(missing_params) > 0) {
    stop("Missing true values for parameters: ", paste(missing_params, collapse = ", "))
  }

  n_iterations <- trajectories$n_iterations
  distances <- numeric(n_iterations)

  for (iter in seq_len(n_iterations)) {
    # Extract parameter values for this iteration
    param_values <- numeric(length(trajectories$parameter_names))

    for (i in seq_along(trajectories$parameter_names)) {
      param_name <- trajectories$parameter_names[i]
      traj_matrix <- trajectories$trajectories[[param_name]]

      # Summarize across chains for this iteration
      iter_values <- traj_matrix[iter, , drop = TRUE]
      param_values[i] <- switch(chain_summary,
        "mean" = mean(iter_values, na.rm = TRUE),
        "median" = median(iter_values, na.rm = TRUE),
        "chain1" = iter_values[1]
      )
    }

    # Calculate distance to true values
    true_vec <- unname(true_values[trajectories$parameter_names])
    diff_vec <- param_values - true_vec

    distances[iter] <- switch(norm_type,
      "euclidean" = sqrt(sum(diff_vec^2)),
      "manhattan" = sum(abs(diff_vec)),
      "max" = max(abs(diff_vec))
    )
  }

  structure(distances,
    norm_type = norm_type,
    chain_summary = chain_summary,
    parameter_names = trajectories$parameter_names,
    true_values = true_values[trajectories$parameter_names],
    class = c("parameter_distance", "numeric")
  )
}

#' Print method for parameter_distance
#' @param x parameter_distance object
#' @param ... additional arguments
#' @export
print.parameter_distance <- function(x, ...) {
  cat("Parameter Distance Trajectory\n")
  cat("============================\n")
  cat("Norm type:", attr(x, "norm_type"), "\n")
  cat("Chain summary:", attr(x, "chain_summary"), "\n")
  cat("Iterations:", length(x), "\n")
  cat("Final distance:", sprintf("%.6f", x[length(x)]), "\n\n")

  cat("True parameter values:\n")
  true_vals <- attr(x, "true_values")
  for (name in names(true_vals)) {
    cat(sprintf("  %s: %.4f\n", name, true_vals[[name]]))
  }

  cat("\nUse plot() to visualize convergence\n")
}

#' Plot method for parameter_distance
#' @param x parameter_distance object
#' @param ... additional arguments for ggplot
#' @export
plot.parameter_distance <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package is required for plotting")
  }

  iteration <- distance <- NULL # avoid R CMD check warnings

  df <- data.frame(
    iteration = seq_along(x),
    distance = as.numeric(x)
  )

  ggplot2::ggplot(df, ggplot2::aes(x = iteration, y = distance)) +
    ggplot2::geom_line() +
    ggplot2::labs(
      title = paste(
        "Parameter Distance Trajectory (",
        attr(x, "norm_type"), "norm )"
      ),
      subtitle = paste("Chain summary:", attr(x, "chain_summary")),
      x = "Iteration",
      y = expression(paste("||", theta, " - ", hat(theta), "||"))
    ) +
    ggplot2::theme_minimal()
}

#' Trace plot of ngme fitting
#'
#' @param ngme ngme object
#' @param name name of latent models, otherwise plot fixed effects and measurement noise.
#'   Use \code{"all"} to plot all latent models plus data parameters in one figure.
#' should be in names(ngme$models) or other
#' @param moving_window moving window for the traceplot
#' @param hline vector, add hline to each plot
#' @param combine bool, if TRUE, return a single faceted ggplot; otherwise return
#'   a list of ggplot objects and print them one by one.
#' @param ncol number of facet columns when \code{combine = TRUE}. If \code{NULL}
#'   (default), use 2 columns for multi-parameter plots and 1 for single-parameter plots.
#'
#' @return A ggplot object when \code{combine = TRUE}; otherwise a list of ggplot
#'   objects.
#' @export
#'
traceplot <- function(
    ngme,
    name = "all",
    moving_window = 1,
    hline = NULL,
    combine = TRUE,
    ncol = NULL) {
  stopifnot(inherits(ngme, "ngme"))
  stopifnot(!is.null(name))
  ngme <- ngme$replicates[[1]]
  ps <- list()

  get_latent_traj_and_ts <- function(latent_obj, latent_name = NULL) {
    traj <- attr(latent_obj, "lat_traj")
    stopifnot(
      "Please run ngme() to estimate the model before using traceplot()" = !is.null(traj)
    )
    ts <- get_latent_info(latent_obj)

    # Special handling for ARMA(p,q<=2): convert raw (PACF) to AR/MA coefficients
    if (!is.null(latent_obj$operator) && identical(latent_obj$operator$model, "arma")) {
      p <- latent_obj$operator$p %||% 0L
      q <- latent_obj$operator$q %||% 0L
      if (p <= 2 && q <= 2) {
        # Helper: map raw -> coeff per time column
        map_ar_coeff <- function(raw_ar) {
          if (length(raw_ar) == 0) {
            return(numeric(0))
          }
          if (length(raw_ar) == 1) {
            t1 <- ar1_th2a(raw_ar[1])
            return(c(t1))
          }
          t1 <- ar1_th2a(raw_ar[1])
          t2 <- ar1_th2a(raw_ar[2])
          c(t1 * (1 - t2), t2)
        }
        map_ma_coeff <- function(raw_ma) {
          if (length(raw_ma) == 0) {
            return(numeric(0))
          }
          if (length(raw_ma) == 1) {
            s1 <- ar1_th2a(raw_ma[1])
            return(c(s1))
          }
          s1 <- ar1_th2a(raw_ma[1])
          s2 <- ar1_th2a(raw_ma[2])
          c(s1 * (1 - s2), s2)
        }
        # Rewrite trajectories per chain
        traj <- lapply(traj, function(M) {
          M <- as.matrix(M)
          if ((p + q) > 0 && nrow(M) >= (p + q)) {
            ar_idx <- if (p > 0) seq_len(p) else integer(0)
            ma_idx <- if (q > 0) (p + seq_len(q)) else integer(0)
            for (col in seq_len(ncol(M))) {
              if (p > 0) {
                phi <- map_ar_coeff(M[ar_idx, col])
                M[seq_len(length(phi)), col] <- phi
              }
              if (q > 0) {
                theta <- map_ma_coeff(M[ma_idx, col])
                M[p + seq_len(length(theta)), col] <- theta
              }
            }
            M
          } else {
            M
          }
        })
        # Replace transforms for AR/MA with identity since we already mapped
        if ((p + q) > 0 && length(ts$trans) >= (p + q)) {
          ts$trans[seq_len(p + q)] <- replicate(p + q, identity, simplify = FALSE)
        }
      }
    }

    if (!is.null(latent_name)) {
      ts$name <- paste0(ts$name, " (", latent_name, ")")
    }

    list(traj = traj, ts = ts)
  }

  if (identical(name, "all")) {
    latent_names <- names(ngme$models)
    latent_parts <- lapply(latent_names, function(lat_name) {
      get_latent_traj_and_ts(ngme$models[[lat_name]], latent_name = lat_name)
    })

    # Data part (measurement noise + fixed effects)
    traj_noise <- attr(ngme, "block_traj")
    stopifnot(
      "Please run ngme() to estimate the model before using traceplot()" = !is.null(traj_noise)
    )
    ts_noise <- get_noise_info(ngme$noise)
    name_feff <- get_feff_plot_names(ngme$feff)
    trans_feff <- rep(list(identity), length(ngme$feff))
    ts_noise$name <- c(ts_noise$name, name_feff)
    ts_noise$trans <- c(ts_noise$trans, trans_feff)

    parts <- c(latent_parts, list(list(traj = traj_noise, ts = ts_noise)))

    n_chains <- length(parts[[1]]$traj)
    if (any(vapply(parts, function(p) length(p$traj) != n_chains, logical(1)))) {
      stop("All components must have the same number of chains for name = 'all'.")
    }

    traj <- lapply(seq_len(n_chains), function(chain_idx) {
      mats <- lapply(parts, function(p) p$traj[[chain_idx]])
      mats <- lapply(mats, function(m) {
        if (is.null(dim(m))) {
          matrix(m, nrow = 1)
        } else {
          as.matrix(m)
        }
      })

      ncol_vec <- vapply(mats, ncol, integer(1))
      nrow_vec <- vapply(mats, nrow, integer(1))
      expected_iter <- max(ncol_vec)
      if (max(nrow_vec) > expected_iter) {
        expected_iter <- max(nrow_vec)
      }

      mats <- Map(function(m, nc, nr) {
        if (nc != expected_iter && nr == expected_iter) {
          m <- t(m)
        }
        m
      }, mats, ncol_vec, nrow_vec)

      lengths <- vapply(mats, ncol, integer(1))
      if (length(unique(lengths)) > 1) {
        warning("Some components are of different lengths. Only the minimum length is used.")
      }
      min_length <- min(lengths)
      mats <- lapply(mats, function(m) m[, seq_len(min_length), drop = FALSE])
      do.call(rbind, mats)
    })

    ts <- list(
      name = unlist(lapply(parts, function(p) p$ts$name), use.names = FALSE),
      trans = do.call(c, lapply(parts, function(p) p$ts$trans))
    )
  } else if (name %in% names(ngme$models)) {
    # Plot trajectory of parameters of the model
    lt <- get_latent_traj_and_ts(ngme$models[[name]])
    traj <- lt$traj
    ts <- lt$ts
  } else {
    # Plot trajectory of parameters of the noise
    traj <- attr(ngme, "block_traj")
    stopifnot(
      "Please run ngme() to estimate the model before using traceplot()" = !is.null(traj)
    )
    # get titles
    ts <- get_noise_info(ngme$noise)
    name_feff <- get_feff_plot_names(ngme$feff)
    trans_feff <- rep(list(identity), length(ngme$feff))
    ts$name <- c(ts$name, name_feff)
    ts$trans <- c(ts$trans, trans_feff)
  }
  ts$name <- disambiguate_trace_param_names(ts$name)

  # Record mean trajectories for comparison later.
  avg_lines <- NULL
  plot_data_all <- list()
  mean_data_all <- list()
  hline_data_all <- list()

  n_parameters <- nrow(traj[[1]]) # number of parameters to draw
  for (idx in seq_len(n_parameters)) {
    # extract the idx-th parameter from each chain
    df <- lapply(traj, function(x) x[idx, , drop = F])
    df <- lapply(df, as.numeric)

    # if not all chains are of the same length, use the minimum length
    lengths <- sapply(df, length)
    if (length(unique(lengths)) > 1) {
      min_length <- min(lengths)
      df <- lapply(df, function(x) x[seq_len(min_length)])
      warning("Some chains are of different lengths. Only the minimum length is used.")
    }

    df <- as.data.frame(df)
    df$x <- seq_len(nrow(df))
    x <- NULL # get around check note
    keys <- setdiff(names(df), "x")
    df_long <- data.frame(
      x = rep(df$x, times = length(keys)),
      key = rep(keys, each = nrow(df)),
      value = as.numeric(unlist(df[keys], use.names = FALSE)),
      stringsAsFactors = FALSE
    )
    ff <- ts$trans[[idx]]
    df_long$value <- ff(df_long$value)

    # update df using moving window
    if (moving_window > 1) {
      if (!requireNamespace("zoo", quietly = TRUE)) {
        message("Package 'zoo' is not installed. Using original data without moving average.")
        df_long$moving_avg <- df_long$value
      } else {
        df_long$moving_avg <- df_long$value
        for (k in keys) {
          idx_k <- df_long$key == k
          df_long$moving_avg[idx_k] <- zoo::rollapply(
            df_long$value[idx_k],
            width = moving_window,
            FUN = mean,
            align = "right",
            fill = NA
          )
        }
      }
    } else {
      df_long$moving_avg <- df_long$value
    }
    df_long <- df_long[!is.na(df_long$moving_avg), , drop = FALSE]

    df_mean <- aggregate(moving_avg ~ x, data = df_long, FUN = mean)
    names(df_mean)[names(df_mean) == "moving_avg"] <- "mean_moving_avg"
    parameter_name <- ts$name[[idx]]
    df_long$parameter <- parameter_name
    df_mean$parameter <- parameter_name

    plot_data_all[[idx]] <- df_long
    mean_data_all[[idx]] <- df_mean

    hline_idx <- if (length(hline) >= idx) hline[[idx]] else NULL
    if (!is.null(hline_idx) && length(hline_idx) > 0) {
      hline_data_all[[idx]] <- data.frame(
        parameter = rep(parameter_name, length(hline_idx)),
        yintercept = as.numeric(hline_idx),
        stringsAsFactors = FALSE
      )
    }

    if (!combine) {
      p <- ggplot2::ggplot() +
        ggplot2::geom_line(
          data = df_long,
          mapping = ggplot2::aes(
            x = .data[["x"]],
            y = .data[["moving_avg"]],
            group = .data[["key"]]
          )
        ) +
        ggplot2::geom_line(
          data = df_mean,
          ggplot2::aes(x = x, y = mean_moving_avg),
          col = "red"
        ) +
        ggplot2::labs(title = parameter_name) +
        ggplot2::xlab(NULL) +
        ggplot2::ylab(NULL) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(size = 14, face = "bold")
        )

      if (!is.null(hline_idx) && length(hline_idx) > 0) {
        p <- p + ggplot2::geom_hline(
          yintercept = as.numeric(hline_idx),
          color = "blue"
        )
      }
      ps[[idx]] <- p
    }

    avg_lines[[parameter_name]] <- df_mean$mean_moving_avg
  }

  # Display results based on combine parameter
  result <- NULL
  if (combine) {
    if (length(plot_data_all) == 0) {
      stop("No parameters available to plot.")
    }
    facet_ncol <- if (is.null(ncol)) {
      if (length(plot_data_all) > 1) 2 else 1
    } else {
      if (!is.numeric(ncol) || length(ncol) != 1 || is.na(ncol) || ncol < 1 || abs(ncol - round(ncol)) > 0) {
        stop("`ncol` must be a positive integer or NULL.")
      }
      as.integer(ncol)
    }
    plot_data <- do.call(rbind, plot_data_all)
    mean_data <- do.call(rbind, mean_data_all)
    plot_levels <- unique(vapply(plot_data_all, function(x) x$parameter[[1]], character(1)))
    plot_data$parameter <- factor(plot_data$parameter, levels = plot_levels)
    mean_data$parameter <- factor(mean_data$parameter, levels = plot_levels)
    if (length(hline_data_all) > 0) {
      hline_data <- do.call(rbind, hline_data_all)
      hline_data$parameter <- factor(hline_data$parameter, levels = plot_levels)
    } else {
      hline_data <- NULL
    }

    result <- ggplot2::ggplot() +
      ggplot2::geom_line(
        data = plot_data,
        mapping = ggplot2::aes(
          x = .data[["x"]],
          y = .data[["moving_avg"]],
          group = .data[["key"]]
        )
      ) +
      ggplot2::geom_line(
        data = mean_data,
        mapping = ggplot2::aes(x = .data[["x"]], y = .data[["mean_moving_avg"]]),
        col = "red"
      ) +
      ggplot2::facet_wrap(~parameter, scales = "free_y", ncol = facet_ncol) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(NULL) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(size = 14, face = "bold")
      )

    if (!is.null(hline_data)) {
      result <- result + ggplot2::geom_hline(
        data = hline_data,
        mapping = ggplot2::aes(yintercept = .data[["yintercept"]]),
        color = "blue",
        inherit.aes = FALSE
      )
    }

    # Keep the historical behavior of drawing a plot when traceplot() is called.
    print(result)
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
#' plot(noise_nig(mu = 1, sigma = 2, nu = 1))
#' plot(n1 = noise_nig(mu = 0, sigma = 1, nu = 1), n2 = noise_nig(mu = 1, sigma = 1.5, nu = 0.5))
plot.ngme_noise <- function(x = NULL, ...) {
  # Get all arguments including the first one
  call_args <- as.list(match.call())[-1] # Remove function name
  all_args <- c(if (!is.null(x)) list(x), list(...))

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
    nu_lower_bound <- if (is.null(noise$nu_lower_bound)) 0 else noise$nu_lower_bound
    nu <- nu_lower_bound + exp(noise$theta_nu)

    stopifnot(
      "only implemented for stationary mu" =
        length(mu) == 1 || noise$noise_type == "normal"
    )
    stopifnot("only implemented for stationary sigma" = length(sigma) == 1)
    stopifnot(
      "only implemented for stationary nu" =
        length(nu) == 1 || noise$noise_type == "normal"
    )

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
  l1 <- as.data.frame(l1)
  l2 <- as.data.frame(l2)

  ps <- list()
  n_plots <- length(l1)
  n_iter <- length(l1[[1]])

  for (i in seq_len(n_plots)) {
    c1 <- c2 <- NULL
    df <- data.frame(c1 = l1[[i]], c2 = l2[[i]], title = names(l1)[[i]])
    ps[[i]] <- ggplot(data = df) +
      geom_line(aes(x = 1:n_iter, y = c1), col = "1") +
      geom_line(aes(x = 1:n_iter, y = c2), col = "2") +
      labs(title = df$title) +
      xlab(NULL) +
      ylab(NULL)
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
#' n1 <- noise_nig(mu = 0, sigma = 1, nu = 1)
#' n2 <- noise_nig(mu = 0.5, sigma = 1.2, nu = 0.8)
#' compare_noise_kld(n1, method2 = n2)
compare_noise_kld <- function(x = NULL, ..., xlim = c(-10, 10), n_points = 1000) {
  # Get all arguments including the first one
  call_args <- as.list(match.call())[-1] # Remove function name
  all_args <- c(if (!is.null(x)) list(x), list(...))

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
    nu_lower_bound <- if (is.null(noise$nu_lower_bound)) 0 else noise$nu_lower_bound
    nu <- nu_lower_bound + exp(noise$theta_nu)

    stopifnot(
      "only implemented for stationary mu" =
        length(mu) == 1 || noise$noise_type == "normal"
    )
    stopifnot("only implemented for stationary sigma" = length(sigma) == 1)
    stopifnot(
      "only implemented for stationary nu" =
        length(nu) == 1 || noise$noise_type == "normal"
    )

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
    kld_values[i - 1] <- sum(p * log(p / q)) * dx
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
