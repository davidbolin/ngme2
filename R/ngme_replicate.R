# create the general replicate model
ngme_replicate <- function(
  Y            = NULL,
  X            = NULL,
  noise        = noise_normal(),
  models       = list(),
  control_ngme = control_ngme(),
  ...
) {
  # compute W_sizes and V_sizes
  W_sizes     = sum(unlist(lapply(models, function(x) x[["W_size"]])))   #W_sizes = sum(ncol_K)
  V_sizes     = sum(unlist(lapply(models, function(x) x[["V_size"]])))   #W_sizes = sum(nrow_K)

  n_la_params = sum(unlist(lapply(models, function(x) x["n_params"])))
  n_feff <- ncol(X);
  feff <- control_ngme$feff; names(feff) <- colnames(X)

  models_string <- rep(" ", 14) # padding of 14 spaces
  for (latent in models)
    models_string <- c(models_string, latent$par_string)
  feff_str  <- if (ncol(X) > 0) paste0("  feff_", seq_along(feff)) else ""
  m_mu_str    <- paste0("    mu_", seq_along(noise$theta_mu))
  m_sigma_str <- paste0(" sigma_", seq_along(noise$theta_sigma))
  m_nu_str    <- "    nu_1"
  merr_str <- switch(noise$noise_type,
    normal  = m_sigma_str,
    nig     = c(m_mu_str, m_sigma_str, m_nu_str)
  )
  par_string <- do.call(paste0, as.list(c(models_string, feff_str, merr_str)))

  n_params <- n_feff + n_la_params + noise$n_params
  structure(
    list(
      Y                 = Y,
      X                 = X,
      feff              = feff,
      models            = models,
      noise             = noise,
      control_ngme      = control_ngme,
      par_string        = par_string,
      W_sizes           = W_sizes,
      V_sizes           = V_sizes,
      n_merr            = noise$n_params,
      n_params          = n_params,
      n_la_params       = n_la_params,
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
  ngme <- x
  cat("*** Ngme object ***\n\n");

  cat("Fixed effects: \n");
  if (length(ngme$feff) > 0) {
    print(ngme$feff, digits=3)
  } else {
    cat("  None\n")
  }
  cat("\n")
  # cat(paste("  ", ngme_format("feff", ngme$feff)));

  cat("Models: \n");
  for (i in seq_along(ngme$models)) {
    # cat("[["); cat(i); cat("]]")
    # cat("\""); cat(names(ngme$models)[[i]]); cat("\"\n")
    cat("$"); cat(names(ngme$models)[[i]]); cat("\n")
    print(ngme$models[[i]], padding = 2)
    cat("\n")
  }

  cat("Measurement noise: \n");
  print(ngme$noise, padding = 2); cat("\n\n")
}