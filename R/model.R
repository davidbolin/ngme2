# the general block model
ngme.block_model <- function(
  Y           = NULL,
  X           = NULL,
  beta        = NULL,
  noise       = noise_normal(),
  latents     = list(),
  control     = list(),
  debug       = FALSE,
  ...
) {

  latents_string <- rep(" ", 14) # padding of 14 spaces
  for (latent in latents)
    latents_string <- c(latents_string, latent$par_string)
  beta_str  <- if (length(beta) > 0) paste0("  beta_", seq_along(beta)) else ""
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
      beta              = beta,
      latents           = latents,
      noise             = noise,
      control           = control,
      n_merr            = noise$n_params,
      debug             = debug,
      par_string        = par_string,
      ...
    ),
    class = "ngme"
  )
}

# function for specify ngme.model basic structure
# keep same with cpp latent structure
ngme_model <- function(
  model,
  W_size      = NULL,
  theta_K     = NULL,
  fix_theta_K = FALSE,
  W           = NULL,
  fix_W       = FALSE,
  A           = NULL,
  A_pred      = NULL,
  noise       = NULL,
  control     = ngme_control_f(),
  V_size      = NULL,
  debug       = FALSE,
  ...
) {
  stopifnot(is.character(model))

  # generate string
  K_str     <- switch(model,
    ar1     = "   alpha",
    matern  = paste0(" kappa_", seq_along(theta_K))
  )
  mu_str    <- paste0("    mu_", seq_along(noise$theta_mu))
  sigma_str <- paste0(" sigma_", seq_along(noise$theta_sigma))
  nu_str    <- "    nu_1"

  structure(
    list(
      model         = model,
      noise_type    = noise$noise_type,
      W_size        = W_size,
      theta_K       = theta_K,
      n_theta_K     = length(theta_K),
      A             = A,
      A_pred        = A_pred,
      noise         = noise,
      W             = W,
      fix_W         = fix_W,
      fix_theta_K   = fix_theta_K,
      V_size        = V_size,
      control       = control,
      n_params      = length(theta_K) + with(noise, n_theta_mu + n_theta_sigma + n_theta_V),
      debug         = debug,
      par_string    = do.call(paste0, as.list(c(K_str, mu_str, sigma_str, nu_str))),
      ...
    ),
    class = "ngme_model"
  )
}

#' Print ngme model
#'
#' @param model
#'
#' @return a list (model specifications)
#' @export
print.ngme_model <- function(model, padding=0) {
  pad_space <- paste(rep(" ", padding), collapse = "")
  pad_add4_space <- paste(rep(" ", padding + 4), collapse = "")

  cat(pad_space); cat("Ngme model: "); cat(model$model); cat("\n")

  cat(pad_space); cat("Model parameters: \n")
  params <- with(model, {
    switch(model,
      "ar1"     = paste0(pad_add4_space, ngme_format("K", theta_K, "ar1")),
      "matern"  = paste0(pad_add4_space, ngme_format("K", theta_K, "matern")),
      "rw1"     = paste0(pad_add4_space, "No parameter needed."),
      "unkown"  = paste0(pad_add4_space, "No parameter needed."),
    )
  })
  cat(params);
  cat("\n\n")

  print.ngme_noise(model$noise, padding = padding)

  invisible(model)
}


# rw1, rw2
# nodes = 100 (inla.group)

# ?inla.spde.make.A
# inla.spde.make.A(
#   index
#  replicates=)
