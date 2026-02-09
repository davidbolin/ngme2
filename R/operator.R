ngme_operator <- function(
    mesh,
    model,
    K,
    h,
    Z = NULL,
    theta_K = NULL,
    zero_trace = FALSE,
    symmetric = FALSE,
    generic_type = "none", # "none", "generic", "generic_ns"
    ...) {
  if (is.null(K)) stop("K is NULL.")
  if (nrow(K) != ncol(K)) stop("K is not a square matrix.")
  if (is.null(Z)) Z <- diag(nrow(K))

  structure(
    list(
      mesh = mesh,
      model = model,
      K = K,
      h = h,
      Z = Z,
      theta_K = theta_K,
      n_theta_K = length(theta_K),
      zero_trace = zero_trace,
      symmetric = symmetric,
      generic_type = generic_type,
      ...
    ),
    class = "ngme_operator"
  )
}


#' Print ngme operator
#'
#' @param x ngme operator object
#' @param padding number of white space padding in front
#' @param prefix prefix string
#' @param ... ...
#'
#' @return a list (operator specifications)
#' @export
print.ngme_operator <- function(x, padding = 0, prefix = "Model type", ...) {
  operator <- x
  pad_space <- paste(rep(" ", padding), collapse = "")
  pad_add4_space <- paste(rep(" ", padding + 4), collapse = "")
  pad_add8_space <- paste(rep(" ", padding + 8), collapse = "")

  limited_format <- function(values, digits = 3, max_items = 6) {
    if (length(values) == 0) {
      return("<none>")
    }
    formatted <- format(values, digits = digits)
    if (length(values) > max_items) {
      formatted <- c(formatted[seq_len(max_items)], "...")
    }
    paste(formatted, collapse = ", ")
  }

  model_name <- switch(operator$model,
    ar1 = "AR(1)",
    arma = "ARMA",
    matern = "Matern",
    tp = "Tensor product",
    bv = "Bivariate model",
    bv_matern = "Bivariate Matern model",
    bv2 = "Bivariate model (v2)",
    iid = "IID model",
    rw1 = "Random walk (order 1)",
    rw2 = "Random walk (order 2)",
    ou = "Ornstein-Uhlenbeck",
    re = "Random effect",
    generic = "Generic",
    generic_ns = "Generic (non-stationary)",
    spacetime = if (operator$method == "galerkin") "Space-time (Galerkin)" else "Space-time (Implicit Euler)",
    {
      # Handle general AR models (ar2, ar3, etc.)
      if (grepl("^ar[0-9]+$", operator$model)) {
        p <- as.numeric(gsub("^ar([0-9]+)$", "\\1", operator$model))
        paste0("AR(", p, ")")
      } else {
        "Unknown"
      }
    }
  )

  model_name <- paste(prefix, ": ", model_name, "\n", sep = "")
  cat(pad_space, model_name, sep = "")

  parameter <- with(operator, switch(model,
    arma = {
      # Compute friendly AR/MA from raw theta_K for p,q<=2
      get_ar_from_raw <- function(raw, p) {
        if (p <= 0) {
          return(numeric(0))
        }
        if (p == 1) {
          return(ar1_th2a(raw[1]))
        }
        if (p == 2) {
          t1 <- ar1_th2a(raw[1])
          t2 <- ar1_th2a(raw[2])
          return(c(t1 * (1 - t2), t2))
        }
        numeric(0)
      }
      get_ma_from_raw <- function(raw, q) {
        if (q <= 0) {
          return(numeric(0))
        }
        if (q == 1) {
          return(ar1_th2a(raw[1]))
        }
        if (q == 2) {
          t1 <- ar1_th2a(raw[1])
          t2 <- ar1_th2a(raw[2])
          return(c(t1 * (1 - t2), t2))
        }
        numeric(0)
      }
      pp <- if (!is.null(p)) p else sum(startsWith(names(theta_K), "ar"))
      qq <- if (!is.null(q)) q else sum(startsWith(names(theta_K), "ma"))
      raw_ar <- if (pp > 0) theta_K[seq_len(pp)] else numeric(0)
      raw_ma <- if (qq > 0) theta_K[pp + seq_len(qq)] else numeric(0)
      ar_vals <- get_ar_from_raw(raw_ar, pp)
      ma_vals <- get_ma_from_raw(raw_ma, qq)
      if (length(ar_vals) > 0) {
        cat(pad_add4_space, "ar = ", limited_format(ar_vals, digits = 3), "\n", sep = "")
      } else {
        cat(pad_add4_space, "ar = <none>\n", sep = "")
      }
      if (length(ma_vals) > 0) {
        cat(pad_add4_space, "ma = ", limited_format(ma_vals, digits = 3), "\n", sep = "")
      } else {
        cat(pad_add4_space, "ma = <none>\n", sep = "")
      }
    },
    ar1 = cat(pad_add4_space, "rho = ", format(ar1_th2a(theta_K), digits = 3), "\n", sep = ""),
    matern = {
      # Print alpha from operator$alpha when available; fixed/free from operator$fix_alpha
      d <- if (!is.null(operator$spatial_dim)) operator$spatial_dim else 2
      L <- d / 2
      alpha_fixed <- isTRUE(operator$fix_alpha)
      if (!is.null(operator$alpha)) {
        alpha_val <- operator$alpha
      } else if (length(theta_K) >= 1 && !alpha_fixed) {
        # Fall back: derive from eta_alpha in theta_K[1]
        eta_alpha <- theta_K[1]
        alpha_val <- L + (4 - L) * (1 / (1 + exp(-eta_alpha)))
      } else {
        alpha_val <- NA
      }
      alpha_status <- if (alpha_fixed) "(fixed)" else "(free)"
      cat(pad_add4_space, "alpha = ", format(alpha_val, digits = 3), " ", alpha_status, "\n", sep = "")

      stationary <- isTRUE(operator$stationary)
      # theta_K layout:
      #   if fix_alpha:      [theta_kappa...]
      #   else (not fixed):  [eta_alpha, theta_kappa...]
      offset <- if (alpha_fixed) 0 else 1
      if (stationary) {
        if (length(theta_K) >= 1) {
          idx <- 1 + offset # stationary has single kappa
          kappa <- exp(theta_K[idx])
          cat(pad_add4_space, "kappa = ", format(kappa, digits = 3), "\n", sep = "")
        }
      } else {
        if (length(theta_K) > offset) {
          theta_kappa <- theta_K[(1 + offset):length(theta_K)]
          cat(pad_add4_space, "theta_kappa = ", limited_format(theta_kappa, digits = 3), "\n", sep = " ")
        }
      }
    },
    tp = {
      print(operator$first, padding = padding + 4, prefix = "first")
      print(operator$second, padding = padding + 4, prefix = "second")
    },
    bv = {
      model_names <- operator$model_names
      fix_theta <- isTRUE(operator$fix_theta)
      use_c_param <- isTRUE(operator$use_c_param)

      idx <- 1
      if (fix_theta) {
        theta <- operator$bv_theta
      } else {
        theta <- operator$param_trans[[idx]](theta_K[idx])
        idx <- idx + 1
      }
      rho <- operator$param_trans[[idx]](theta_K[idx])
      idx <- idx + 1

      cat(pad_add4_space, "theta = ", format(theta, digits = 3),
        if (fix_theta) " (fixed)" else "", "\n",
        sep = ""
      )
      cat(pad_add4_space, "rho = ", format(rho, digits = 3), "\n", sep = "")

      if (use_c_param) {
        c1 <- operator$param_trans[[idx]](theta_K[idx])
        idx <- idx + 1
        c2 <- operator$param_trans[[idx]](theta_K[idx])
        cat(pad_add4_space, "c1 = ", format(c1, digits = 3), "\n", sep = "")
        cat(pad_add4_space, "c2 = ", format(c2, digits = 3), "\n", sep = "")
      }

      print(operator$first, padding = padding + 4, prefix = model_names[[1]])
      print(operator$second, padding = padding + 4, prefix = model_names[[2]])
    },
    bv2 = {
      theta <- 0
      rho <- operator$param_trans[[1]](theta_K[1])
      c1 <- operator$param_trans[[2]](theta_K[2])
      c2 <- operator$param_trans[[3]](theta_K[3])
      cat(pad_add4_space, "theta = ", format(theta, digits = 3), "(fixed) \n", sep = "")
      cat(pad_add4_space, "rho = ", format(rho, digits = 3), "\n", sep = "")
      cat(pad_add4_space, "c1 = ", format(c1, digits = 3), "\n", sep = "")
      cat(pad_add4_space, "c2 = ", format(c2, digits = 3), "\n", sep = "")
      print(operator$first, padding = padding + 4, prefix = model_names[[1]])
      print(operator$second, padding = padding + 4, prefix = model_names[[2]])
    },
    bv_matern = {
      model_names <- operator$model_names
      fix_theta <- isTRUE(operator$fix_theta)

      idx <- 1
      if (fix_theta) {
        theta <- operator$bv_theta
      } else {
        theta <- operator$param_trans[[idx]](theta_K[idx])
        idx <- idx + 1
      }
      rho <- operator$param_trans[[idx]](theta_K[idx])
      idx <- idx + 1
      sd1 <- operator$param_trans[[idx]](theta_K[idx])
      idx <- idx + 1
      sd2 <- operator$param_trans[[idx]](theta_K[idx])

      cat(pad_add4_space, "theta = ", format(theta, digits = 3),
        if (fix_theta) " (fixed)" else "", "\n",
        sep = ""
      )
      cat(pad_add4_space, "rho = ", format(rho, digits = 3), "\n", sep = "")
      cat(pad_add4_space, "sd1 = ", format(sd1, digits = 3), "\n", sep = "")
      cat(pad_add4_space, "sd2 = ", format(sd2, digits = 3), "\n", sep = "")

      print(operator$first, padding = padding + 4, prefix = model_names[[1]])
      print(operator$second, padding = padding + 4, prefix = model_names[[2]])
    },
    bv_matern_normal = {
      theta <- 0
      rho <- operator$param_trans[[1]](theta_K[1])
      sd1 <- operator$param_trans[[2]](theta_K[2])
      sd2 <- operator$param_trans[[3]](theta_K[3])
      cat(pad_add4_space, "theta = ", format(theta, digits = 3), "(fixed) \n", sep = "")
      cat(pad_add4_space, "rho = ", format(rho, digits = 3), "\n", sep = "")
      cat(pad_add4_space, "sd1 = ", format(sd1, digits = 3), "\n", sep = "")
      cat(pad_add4_space, "sd2 = ", format(sd2, digits = 3), "\n", sep = "")
      print(operator$first, padding = padding + 4, prefix = model_names[[1]])
      print(operator$second, padding = padding + 4, prefix = model_names[[2]])
    },
    bv_matern_nig = {
      fix_theta <- operator$fix_theta
      theta <- if (fix_theta) operator$bv_theta else operator$param_trans[[1]](theta_K[1])
      rho <- operator$param_trans[[2 - fix_theta]](theta_K[2 - fix_theta])
      sd1 <- operator$param_trans[[3 - fix_theta]](theta_K[3 - fix_theta])
      sd2 <- operator$param_trans[[4 - fix_theta]](theta_K[4 - fix_theta])
      cat(pad_add4_space, "theta = ", format(theta, digits = 3), if (fix_theta) "(fixed)" else "", "\n", sep = "")
      cat(pad_add4_space, "rho = ", format(rho, digits = 3), "\n", sep = "")
      cat(pad_add4_space, "sd1 = ", format(sd1, digits = 3), "\n", sep = "")
      cat(pad_add4_space, "sd2 = ", format(sd2, digits = 3), "\n", sep = "")
      print(operator$first, padding = padding + 4, prefix = model_names[[1]])
      print(operator$second, padding = padding + 4, prefix = model_names[[2]])
    },
    ou = {
      theta <- exp(theta_K)
      cat(pad_add4_space, "theta = ", format(theta, digits = 3), "\n", sep = "")
    },
    re = {
      cat(pad_add4_space, "Covariance matrix (Sigma): \n")
      K <- build_effect_K(nrow(operator$K), operator$theta_K)
      print(solve(t(K) %*% K))
    },
    spacetime = {
      cat(pad_add4_space)
      cat("alpha =", format(alpha, digits = 2), "(fixed)", "\n", sep = " ")
      if (!fix_gamma) {
        if (shared_theta_gamma) {
          # When shared_theta_gamma = TRUE, only extract theta_gamma_x and use it for both
          theta_gamma_shared <- theta_K[3:(2 + n_theta_gamma_x)]
          cat(pad_add4_space)
          cat("theta_gamma (x and y) =", format(theta_gamma_shared, digits = 2), "\n", sep = " ")
        } else {
          # When shared_theta_gamma = FALSE, extract both separately
          theta_gamma_x <- theta_K[3:(2 + n_theta_gamma_x)]
          theta_gamma_y <- theta_K[(3 + n_theta_gamma_x):(2 + n_theta_gamma_x + n_theta_gamma_y)]
          cat(pad_add4_space)
          cat("theta_gamma_x =", format(theta_gamma_x, digits = 2), "\n", sep = " ")
          cat(pad_add4_space)
          cat("theta_gamma_y =", format(theta_gamma_y, digits = 2), "\n", sep = " ")
        }
      } else {
        cat(pad_add4_space)
        cat("theta_gamma_x =", format(theta_gamma_x, digits = 2), "(fixed)", "\n", sep = " ")
        cat(pad_add4_space)
        cat("theta_gamma_y =", format(theta_gamma_y, digits = 2), "(fixed)", "\n", sep = " ")
      }

      cat(pad_add4_space)
      cat("cc =", format(exp(theta_K[1]), digits = 2), "\n", sep = " ")

      cat(pad_add4_space)
      cat("kappa =", format(exp(theta_K[2]), digits = 2), "\n", sep = " ")
    },
    generic = {
      theta_vals <- operator$theta_K
      theta_names <- names(theta_vals)
      param_names <- operator$param_name
      trans_list <- operator$trans
      trans_names <- names(trans_list)

      get_label <- function(idx) {
        if (!is.null(param_names) && length(param_names) >= idx &&
          !is.na(param_names[idx]) && nzchar(param_names[idx])) {
          return(param_names[idx])
        }
        if (!is.null(theta_names) && length(theta_names) >= idx &&
          !is.na(theta_names[idx]) && nzchar(theta_names[idx])) {
          return(theta_names[idx])
        }
        paste0("theta_", idx)
      }

      if (length(theta_vals) > 0) {
        cat(pad_add4_space, "Parameters:\n", sep = "")
        for (i in seq_along(theta_vals)) {
          label <- get_label(i)
          cat(pad_add8_space, label, " raw: ", limited_format(theta_vals[i]), "\n", sep = "")

          trans_key <- NULL
          if (!is.null(theta_names) && length(theta_names) >= i &&
            !is.na(theta_names[i]) && nzchar(theta_names[i]) &&
            !is.null(trans_names) && theta_names[i] %in% trans_names) {
            trans_key <- theta_names[i]
          } else if (!is.null(param_names) && length(param_names) >= i &&
            !is.na(param_names[i]) && nzchar(param_names[i]) &&
            !is.null(trans_names) && param_names[i] %in% trans_names) {
            trans_key <- param_names[i]
          }

          if (!is.null(trans_key)) {
            mapping <- trans_list[[trans_key]]
            active_idx <- which(mapping != "null")
            if (length(active_idx) > 0) {
              combo <- paste0("#", active_idx, ": ", mapping[active_idx])
              cat(pad_add8_space, "matrix transforms: ", paste(combo, collapse = ", "), "\n", sep = "")
            }
          }
        }
      } else {
        cat(pad_add4_space, "Parameters: none\n", sep = "")
      }

      if (!is.null(operator$matrices)) {
        dims <- tryCatch(dim(operator$K), error = function(...) NULL)
        size_info <- if (!is.null(dims)) paste0(dims[1], " x ", dims[2]) else "unknown size"
        cat(pad_add4_space, "Stored matrices: ", length(operator$matrices), " (", size_info, ")\n", sep = "")
      }
    },
    generic_ns = {
      param_map <- operator$param_map
      param_names <- names(param_map)

      if (length(param_map) > 0) {
        cat(pad_add4_space, "Parameters:\n", sep = "")
        for (name in param_names) {
          indices <- param_map[[name]]
          coeff <- operator$theta_K[indices]
          cat(pad_add8_space, name, " coefficients: ", limited_format(coeff), "\n", sep = "")

          if (!is.null(operator$trans) && !is.null(operator$trans[[name]])) {
            cat(pad_add8_space, "transform: ", operator$trans[[name]], "\n", sep = "")
          }

          if (!is.null(operator$B_theta_K) && !is.null(operator$B_theta_K[[name]])) {
            basis <- operator$B_theta_K[[name]]
            cat(pad_add8_space, "basis dim: ", nrow(basis), " x ", ncol(basis), "\n", sep = "")
          }
        }
      } else {
        cat(pad_add4_space, "Parameters: none\n", sep = "")
      }

      if (!is.null(operator$position) && length(operator$position) > 0) {
        cat(pad_add4_space, "Matrix combinations (position indices):\n", sep = "")
        for (i in seq_along(operator$position)) {
          combo <- operator$position[[i]]
          cat(pad_add8_space, "#", i, ": ", paste(combo, collapse = " -> "), "\n", sep = "")
        }
      }

      if (!is.null(operator$matrices)) {
        dims <- tryCatch(dim(operator$K), error = function(...) NULL)
        size_info <- if (!is.null(dims)) paste0(dims[1], " x ", dims[2]) else "unknown size"
        cat(pad_add4_space, "Stored matrices: ", length(operator$matrices), " (", size_info, ")\n", sep = "")
      }
    },
    # Handle general AR models (ar2, ar3, etc.) or default case
    {
      if (grepl("^ar[0-9]+$", operator$model)) {
        # Print all rho parameters for general AR models
        for (i in seq_along(theta_K)) {
          cat(pad_add4_space, "rho", i, " = ", format(theta_K[i], digits = 3), "\n", sep = "")
        }
      } else {
        cat(pad_add4_space, "No parameter.", "\n", sep = "")
      }
    }
  ))

  invisible(operator)
}


# function for specify ngme.model basic structure
# keep same with cpp latent structure
ngme_model <- function(
    model,
    operator,
    noise = noise_normal(),
    theta_K = NULL,
    fix_K = FALSE,
    W_size = NULL,
    W = NULL,
    fix_W = FALSE,
    A = NULL,
    V_size = NULL,
    debug = FALSE,
    n_params = NULL,
    name = "field",
    mesh = NULL,
    par_string = NULL,
    map = NULL, # map is the covariates
    group = NULL,
    ...) {
  stopifnot(is.character(model))

  stopifnot(inherits(operator, "ngme_operator"))
  stopifnot(inherits(noise, "ngme_noise"))

  ope_params <- operator$param_name

  mu_params <- if (length(noise$theta_mu) == 0 || isTRUE(noise$fix_theta_mu)) {
    character(0)
  } else {
    paste0("mu_", seq_along(noise$theta_mu))
  }

  sigma_params <- {
    sigma_idx <- if (length(noise$theta_sigma) == 0) integer(0) else which(!noise$fix_theta_sigma)
    if (length(sigma_idx) == 0) character(0) else paste0("sigma_", sigma_idx)
  }

  nu_params <- if (length(noise$theta_nu) == 0 || isTRUE(noise$fix_theta_nu)) {
    character(0)
  } else {
    paste0("nu_", seq_along(noise$theta_nu))
  }

  rho_params <- if (length(noise$rho) == 0 || isTRUE(noise$fix_rho)) {
    character(0)
  } else {
    paste0("rho_", seq_along(noise$rho))
  }

  # Order matches the underlying parameter vector: operator -> mu -> sigma -> nu -> rho
  par_names <- c(ope_params, mu_params, sigma_params, nu_params, rho_params)

  if (is.null(n_params)) n_params <- length(operator$theta_K) + with(noise, n_params)

  structure(
    list(
      model = model,
      operator = operator,
      theta_K = operator$theta_K,
      noise_type = noise$noise_type,
      W_size = W_size,
      A = A,
      noise = noise,
      W = W,
      fix_theta_K = fix_K,
      fix_W = fix_W,
      V_size = V_size,
      n_params = n_params,
      debug = debug,
      par_names = par_names,
      name = name,
      mesh = mesh,
      map = map,
      n_map = length_map(map),
      group = group,
      ...
    ),
    class = "ngme_model"
  )
}

#' Print ngme model
#'
#' @param x ngme model object
#' @param padding number of white space padding in front
#' @param ... ...
#'
#' @return a list (model specifications)
#' @export
print.ngme_model <- function(x, padding = 0, ...) {
  model <- x
  pad_space <- paste(rep(" ", padding), collapse = "")

  # Print operator information
  # If we have trajectory attached (from fitting), prefer the last estimates
  op <- model$operator
  traj <- attr(model, "lat_traj")
  if (!is.null(traj) && length(traj) > 0 && is.list(traj)) {
    n_theta_K <- length(op$theta_K)
    # collect last column from each chain
    last_mat <- tryCatch(
      {
        do.call(cbind, lapply(traj, function(m) m[, ncol(m), drop = FALSE]))
      },
      error = function(e) NULL
    )
    if (!is.null(last_mat) && nrow(last_mat) >= n_theta_K) {
      last_avg <- rowMeans(last_mat)
      theta_last_raw <- as.numeric(last_avg[seq_len(n_theta_K)])
      op2 <- op
      op2$theta_K <- theta_last_raw
      # For matern with free alpha, update printable alpha from raw eta
      if (identical(model$model, "matern")) {
        d <- if (!is.null(op$spatial_dim)) op$spatial_dim else 2
        L <- d / 2
        if (isTRUE(op$fix_alpha)) {
          # keep op$alpha as provided
        } else if (length(theta_last_raw) >= 1) {
          eta_alpha <- theta_last_raw[1]
          op2$alpha <- L + (4 - L) * (1 / (1 + exp(-eta_alpha)))
        }
      }
      print.ngme_operator(op2, padding = padding)
    } else {
      print.ngme_operator(op, padding = padding)
    }
  } else {
    print.ngme_operator(op, padding = padding)
  }

  # Print replicate information if replicates are used
  if (!is.null(model$replicate)) {
    n_replicates <- length(levels(model$replicate))
    cat(pad_space, "Number of replicates: ", n_replicates, "\n\n", sep = "")
  }

  # Print noise information
  model_type <- model$model
  print.ngme_noise(
    model$noise,
    padding = padding,
    model_type = model_type
  )

  invisible(model)
}
