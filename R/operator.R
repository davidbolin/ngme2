ngme_operator <- function(
  mesh,
  model,
  K,
  h,
  theta_K = NULL,
  zero_trace = FALSE,
  symmetric = FALSE,
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
    bv_normal_2  = "Bivariate model 1 (normal noise)",
    bv_normal  = "Bivariate model 2 (normal noise)",
    iid = "IID model",
    rw1 = "Random walk (order 1)",
    rw2 = "Random walk (order 2)",
    ou  = "Ornstein-Uhlenbeck",
    re  = "Random effect",
    "Unknown"
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
    ou = {
      cat(pad_add4_space); cat("theta_K = ", format(theta_K, digits=2), "\n", sep=" ")
    },
    re = {
      cat(pad_add4_space, "Covariance matrix (Sigma): \n")
      K = build_effect_K(nrow(operator$K), operator$theta_K)
      print(solve(t(K) %*% K))
    },
    cat(pad_add4_space, "No parameter.", "\n", sep="")
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
  # mu_str    <- paste0("    mu_", seq_along(noise$theta_mu))
  # sigma_str <- paste0(" sigma_", seq_along(noise$theta_sigma))
  # nu_str    <- "    nu_1"

  # if ((noise$noise_type == "normal"))
  #   par_string <- do.call(paste0, as.list(c(K_str, sigma_str)))
  # else
  #   par_string <- do.call(paste0, as.list(c(K_str, mu_str, sigma_str, nu_str)))

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
  print.ngme_operator(model$operator, padding = padding)
  suppress_sigma <- if (model$model == "re") TRUE else FALSE
  print.ngme_noise(model$noise, padding = padding, suppress_sigma=suppress_sigma)
  invisible(model)
}

