ngme_operator <- function(
  mesh,
  model,
  K,
  h,
  theta_K = NULL,
  zero_trace = FALSE,
  symmetric = FALSE,
  generic_type = "none", # "none", "generic", "generic_ns"
  ...
) {
  if (is.null(K)) stop("K is NULL.")

  structure(
    list(
      mesh = mesh,
      model = model,
      K = K,
      h = h,
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
    matern = "Matern",
    tp  = "Tensor product",
    bv  = "Bivariate model (non-Gaussian noise)",
    bv_matern_normal  = "Bivariate Matern model (normal noise)",
    bv_matern_nig  = "Bivariate Matern model (NIG noise)",
    bv_normal  = "Bivariate model 2 (normal noise)",
    iid = "IID model",
    rw1 = "Random walk (order 1)",
    rw2 = "Random walk (order 2)",
    ou  = "Ornstein-Uhlenbeck",
    re  = "Random effect",
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
  cat(pad_space, model_name, sep="")

  parameter <- with(operator, switch(model,
    ar1 = cat(pad_add4_space, "rho = ", format(ar1_th2a(theta_K), digits=3), "\n", sep=""),
    matern = if (length(theta_K) > 1)
      cat(pad_add4_space, "theta_K = ", format(theta_K, digits=3), "\n", sep=" ") else
      cat(pad_add4_space, "kappa = ", format(exp(theta_K), digits=3), "\n", sep=""),
    tp = {
      print(operator$first,  padding = padding + 4, prefix = "first")
      print(operator$second, padding = padding + 4, prefix = "second")
    },
    bv = {
      theta = operator$param_trans[[1]](theta_K[1])
      rho = operator$param_trans[[2]](theta_K[2])
      cat(pad_add4_space, "theta = ", format(theta, digits=3), "\n", sep="")
      cat(pad_add4_space, "rho = ", format(rho, digits=3), "\n", sep="")
      print(operator$first,  padding = padding + 4, prefix = model_names[[1]])
      print(operator$second, padding = padding + 4, prefix = model_names[[2]])
    },
    bv_normal = {
      theta = 0
      rho = operator$param_trans[[1]](theta_K[1])
      c1 = operator$param_trans[[2]](theta_K[2])
      c2 = operator$param_trans[[3]](theta_K[3])
      cat(pad_add4_space, "theta = ", format(theta, digits=3), "(fixed) \n", sep="")
      cat(pad_add4_space, "rho = ", format(rho, digits=3), "\n", sep="")
      cat(pad_add4_space, "c1 = ", format(c1, digits=3), "\n", sep="")
      cat(pad_add4_space, "c2 = ", format(c2, digits=3), "\n", sep="")
      print(operator$first,  padding = padding + 4, prefix = model_names[[1]])
      print(operator$second, padding = padding + 4, prefix = model_names[[2]])
    },
    bv_matern_normal = {
      theta = 0
      rho = operator$param_trans[[1]](theta_K[1])
      sd1 = operator$param_trans[[2]](theta_K[2])
      sd2 = operator$param_trans[[3]](theta_K[3])
      cat(pad_add4_space, "theta = ", format(theta, digits=3), "(fixed) \n", sep="")
      cat(pad_add4_space, "rho = ", format(rho, digits=3), "\n", sep="")
      cat(pad_add4_space, "sd1 = ", format(sd1, digits=3), "\n", sep="")
      cat(pad_add4_space, "sd2 = ", format(sd2, digits=3), "\n", sep="")
      print(operator$first,  padding = padding + 4, prefix = model_names[[1]])
      print(operator$second, padding = padding + 4, prefix = model_names[[2]])
    },
    bv_matern_nig = {
      fix_theta = operator$fix_bv_theta
      theta = if (fix_theta) operator$bv_theta else operator$param_trans[[1]](theta_K[1])
      rho   = operator$param_trans[[2 - fix_theta]](theta_K[2 - fix_theta])
      sd1   = operator$param_trans[[3 - fix_theta]](theta_K[3 - fix_theta])
      sd2   = operator$param_trans[[4 - fix_theta]](theta_K[4 - fix_theta])
      cat(pad_add4_space, "theta = ", format(theta, digits=3), if (fix_theta) "(fixed)" else "", "\n", sep="")
      cat(pad_add4_space, "rho = ", format(rho, digits=3), "\n", sep="")
      cat(pad_add4_space, "sd1 = ", format(sd1, digits=3), "\n", sep="")
      cat(pad_add4_space, "sd2 = ", format(sd2, digits=3), "\n", sep="")
      print(operator$first,  padding = padding + 4, prefix = model_names[[1]])
      print(operator$second, padding = padding + 4, prefix = model_names[[2]])
    },
    ou = {
      cat(pad_add4_space); cat("theta_K = ", format(theta_K, digits=2), "\n", sep=" ")
    },
    re = {
      cat(pad_add4_space, "Covariance matrix (Sigma): \n")
      K = build_effect_K(nrow(operator$K), operator$theta_K)
      print(solve(t(K) %*% K))
    },
    spacetime = {
      cat(pad_add4_space); cat("alpha =", format(alpha, digits=2), "(fixed)", "\n", sep=" ")
      if (!fix_gamma) {
        if (shared_theta_gamma) {
          # When shared_theta_gamma = TRUE, only extract theta_gamma_x and use it for both
          theta_gamma_shared <- theta_K[3:(2+n_theta_gamma_x)]
          cat(pad_add4_space); cat("theta_gamma (x and y) =", format(theta_gamma_shared, digits=2), "\n", sep=" ")
        } else {
          # When shared_theta_gamma = FALSE, extract both separately
          theta_gamma_x <- theta_K[3:(2+n_theta_gamma_x)]
          theta_gamma_y <- theta_K[(3+n_theta_gamma_x):(2+n_theta_gamma_x+n_theta_gamma_y)]
          cat(pad_add4_space); cat("theta_gamma_x =", format(theta_gamma_x, digits=2), "\n", sep=" ")
          cat(pad_add4_space); cat("theta_gamma_y =", format(theta_gamma_y, digits=2), "\n", sep=" ")
        }
      } else {
        cat(pad_add4_space); cat("theta_gamma_x =", format(theta_gamma_x, digits=2), "(fixed)", "\n", sep=" ")
        cat(pad_add4_space); cat("theta_gamma_y =", format(theta_gamma_y, digits=2), "(fixed)", "\n", sep=" ")
      }
      
      cat(pad_add4_space); cat("cc =", format(exp(theta_K[1]), digits=2), "\n", sep=" ")
      
      cat(pad_add4_space); cat("kappa =", format(exp(theta_K[2]), digits=2), "\n", sep=" ")
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
          cat(pad_add4_space, "rho", i, " = ", format(theta_K[i], digits=3), "\n", sep="")
        }
      } else {
        cat(pad_add4_space, "No parameter.", "\n", sep="")
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
  noise       = noise_normal(),
  theta_K     = NULL,
  W_size      = NULL,
  W           = NULL,
  fix_W       = FALSE,
  A           = NULL,
  control     = control_f(),
  V_size      = NULL,
  debug       = FALSE,
  n_params    = NULL,
  name        = "field",
  mesh        = NULL,
  par_string  = NULL,
  map         = NULL,  # map is the covariates
  group       = NULL,
  ...
) {
  stopifnot(is.character(model))

  stopifnot(inherits(operator, "ngme_operator"))
  stopifnot(inherits(noise, "ngme_noise"))

  # generate string (8 digits)
  # K_str     <- switch(model,
  #   ar1     = "   alpha",
  #   matern  = paste0(" kappa_", seq_along(theta_K)),
  #   rw1     = paste0(" ignored")
  # )

  # make it str with each parameter name contain 8 digits (right aligned)
  ope_params <- operator$param_name
  ope_str   <- sapply(ope_params, function(x) sprintf("%8s", x))
  mu_params <- paste0("mu_", seq_along(noise$theta_mu))
  mu_str    <- sapply(mu_params, function(x) sprintf("%8s", x))
  sigma_params <- paste0("sigma_", seq_along(noise$theta_sigma))
  sigma_str <- sapply(sigma_params, function(x) sprintf("%8s", x))
  nu_params <- paste0("nu_", seq_along(noise$nu))
  nu_str    <- sapply(nu_params, function(x) sprintf("%8s", x))

  if (all(noise$noise_type == "normal"))
    par_string <- paste0(ope_str, sigma_str)
  else
    par_string <- paste0(ope_str, mu_str, sigma_str, nu_str)

  if (is.null(n_params)) n_params <- length(operator$theta_K) + with(noise, n_params)

  structure(
    list(
      model         = model,
      operator      = operator,
      theta_K       = operator$theta_K,
      noise_type    = noise$noise_type,
      W_size        = W_size,
      A             = A,
      noise         = noise,
      W             = W,
      fix_W         = fix_W,
      V_size        = V_size,
      control       = control,
      n_params      = n_params,
      debug         = debug,
      par_string    = par_string,
      name          = name,
      mesh          = mesh,
      map           = map,
      n_map         = length_map(map),
      group         = group,
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
  print.ngme_operator(model$operator, padding = padding)
  
  # Print replicate information if replicates are used
  if (!is.null(model$replicate)) {
    n_replicates <- length(levels(model$replicate))
    cat(pad_space, "Number of replicates: ", n_replicates, "\n\n", sep="")
  }
  
  # Print noise information
  model_type <- model$model
  print.ngme_noise(
    model$noise, 
    padding = padding, 
    model_type=model_type
  )
  
  invisible(model)
}
