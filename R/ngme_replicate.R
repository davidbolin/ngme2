# create the general replicate model
ngme_replicate <- function(
  Y            = NULL,
  X            = NULL,
  noise        = noise_normal(),
  randeffs     = list(),
  latents      = list(),
  control_ngme = list(),
  ...
) {
  # compute W_sizes and V_sizes
  W_sizes     = sum(unlist(lapply(latents, function(x) x[["n_rep"]] * x[["W_size"]])))   #W_sizes = sum(ncol_K)
  V_sizes     = sum(unlist(lapply(latents, function(x) x[["n_rep"]] * x[["V_size"]])))   #W_sizes = sum(nrow_K)

  n_reffs = sum(unlist(lapply(randeffs, function(x) x["n_reff"])))
  B_reffs = do.call(cbind, lapply(randeffs, function(x) x[["B_reff"]]))

  n_re_params = sum(unlist(lapply(randeffs, function(x) x["n_params"])))
  n_la_params = sum(unlist(lapply(latents, function(x) x["n_params"])))
  n_feff <- ncol(X);

  n_params <- n_re_params + n_feff + n_la_params + noise$n_params

  latents_string <- rep(" ", 14) # padding of 14 spaces
  for (latent in latents)
    latents_string <- c(latents_string, latent$par_string)
  beta_str  <- if (ncol(X) > 0) paste0("  beta_", seq_along(beta)) else ""
  m_mu_str    <- paste0("    mu_", seq_along(noise$theta_mu))
  m_sigma_str <- paste0(" sigma_", seq_along(noise$theta_sigma))
  m_nu_str    <- "    nu_1"
  merr_str <- switch(noise$noise_type,
    normal  = m_sigma_str,
    nig     = c(m_mu_str, m_sigma_str, m_nu_str)
  )
  par_string <- do.call(paste0, as.list(c(latents_string, beta_str, merr_str)))

  structure(
    list(
      Y                 = Y,
      X                 = X,
      beta              = control_ngme$beta,
      randeffs          = randeffs,
      latents           = latents,
      noise             = noise,
      control_ngme      = control_ngme,
      par_string        = par_string,
      W_sizes           = W_sizes,
      V_sizes           = V_sizes,
      n_reffs           = n_reffs,
      B_reffs           = B_reffs,
      n_merr            = noise$n_params,
      n_params          = n_params,
      n_re_params       = n_re_params,
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
  cat(paste("  ", ngme_format("beta", ngme$beta)));
  cat("\n\n")

  cat("Latent models: \n");
  for (i in seq_along(ngme$latents)) {
    # cat("[["); cat(i); cat("]]")
    # cat("\""); cat(names(ngme$latents)[[i]]); cat("\"\n")
    cat("$"); cat(names(ngme$latents)[[i]]); cat("\n")
    print(ngme$latents[[i]], padding = 2)
  }
  cat("\n")

  cat("Measurement noise: \n");
  print(ngme$noise, padding = 2); cat("\n\n")
}