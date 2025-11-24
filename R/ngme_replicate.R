# create the general replicate model
ngme_replicate <- function(
    Y = NULL,
    X = NULL,
    noise = noise_normal(),
    models = list(),
    control_ngme = control_ngme(),
    standardize = TRUE,
    log_likelihood = NULL,
    ...) {
  # compute W_sizes and V_sizes
  W_sizes <- sum(unlist(lapply(models, function(x) x[["W_size"]]))) # W_sizes = sum(ncol_K)
  V_sizes <- sum(unlist(lapply(models, function(x) x[["V_size"]]))) # W_sizes = sum(nrow_K)

  n_la_params <- sum(unlist(lapply(models, function(x) x["n_params"])))
  n_feff <- ncol(X)
  feff <- control_ngme$feff
  names(feff) <- colnames(X)

  # Collect parameter names as a vector (not concatenated string)
  par_names <- character(0)

  # Add latent model parameter names
  for (latent in models) {
    if (!is.null(latent$par_names)) {
      par_names <- c(par_names, latent$par_names)
    }
  }

  # Add fixed effect names
  if (ncol(X) > 0) {
    par_names <- c(par_names, paste0("feff_", seq_along(feff)))
  }

  # Add measurement error parameter names
  merr_names <- switch(noise$noise_type,
    normal = paste0("sigma_", seq_along(noise$theta_sigma)),
    nig = c(
      paste0("mu_", seq_along(noise$theta_mu)),
      paste0("sigma_", seq_along(noise$theta_sigma)),
      paste0("nu_", seq_along(noise$nu))
    )
  )
  par_names <- c(par_names, merr_names)

  n_params <- n_feff + n_la_params + noise$n_params
  structure(
    list(
      Y                 = Y,
      X                 = X,
      feff              = feff,
      models            = models,
      noise             = noise,
      control_ngme      = control_ngme,
      par_names         = par_names,
      W_sizes           = W_sizes,
      V_sizes           = V_sizes,
      n_merr            = noise$n_params,
      n_params          = n_params,
      n_la_params       = n_la_params,
      standardize       = standardize,
      log_likelihood    = log_likelihood,
      ...
    ),
    class = c("ngme_replicate", "list")
  )
}

#' Print ngme object
#'
#' @param x ngme object
#' @param ... ignored
#'
#' @return a list (noise specifications)
#' @export
print.ngme_replicate <- function(x, ...) {
  ngme_rep <- x
  cat("*** Ngme object ***\n\n")

  cat("Fixed effects: \n")
  if (length(ngme_rep$feff) > 0) {
    print(ngme_rep$feff, digits = 3)
  } else {
    cat("  None\n")
  }
  cat("\n")
  # cat(paste("  ", ngme_rep_format("feff", ngme_rep$feff)));

  cat("Models: \n")
  for (i in seq_along(ngme_rep$models)) {
    # cat("[["); cat(i); cat("]]")
    # cat("\""); cat(names(ngme_rep$models)[[i]]); cat("\"\n")
    cat("$")
    cat(names(ngme_rep$models)[[i]])
    cat("\n")
    print(ngme_rep$models[[i]], padding = 2)
    cat("\n")
  }

  cat("Measurement noise: \n")
  print(ngme_rep$noise, padding = 2)
  cat("\n\n")

  if (ngme_rep$all_gaussian && !is.null(ngme_rep$log_likelihood)) {
    cat("Log likelihood: ")
    cat(ngme_rep$log_likelihood)
    cat("\n\n")
  }
}
