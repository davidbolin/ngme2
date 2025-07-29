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
      for (i in seq_along(operator$param_name)) {
        t <- operator$theta_K[i]
        cat(pad_add4_space, operator$param_name[i], "=", format(t, digits=3), "\n", sep=" ")
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

