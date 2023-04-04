
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
  noise       = noise_normal(),
  control     = control_f(),
  V_size      = NULL,
  debug       = FALSE,
  n_theta_K   = length(theta_K),
  h           = NULL,
  n_params    = NULL,
  name        = "field",
  mesh        = NULL,
  par_string  = NULL,
  map         = NULL,  # map is the covariates
  n_map       = NULL,
  replicate   = NULL,
  n_rep       = NULL,
  group       = NULL,
  ...
) {
  if (is.null(n_rep)) n_rep <- length(unique(replicate))

  stopifnot(is.character(model))
  # generate string (8 digits)
  K_str     <- switch(model,
    ar1     = "   alpha",
    matern  = paste0(" kappa_", seq_along(theta_K)),
    rw1     = paste0(" ignored")
  )
  mu_str    <- paste0("    mu_", seq_along(noise$theta_mu))
  sigma_str <- paste0(" sigma_", seq_along(noise$theta_sigma))
  nu_str    <- "    nu_1"

  # Notice noise$h should be same as h
  if (is.null(h)) h <- noise$h else noise$h <- h
  if ((noise$noise_type == "normal"))
    par_string <- do.call(paste0, as.list(c(K_str, sigma_str)))
  else
    par_string <- do.call(paste0, as.list(c(K_str, mu_str, sigma_str, nu_str)))

  stopifnot("replicate is NULL" = !is.null(replicate))
  stopifnot("make sure length of replicate == length of index" =
   length(replicate) == length_map(map))
  replicate <- as.integer(replicate)
  if (is.null(n_params)) n_params <- length(theta_K) + with(noise, n_params)

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
      n_params      = n_params,
      debug         = debug,
      par_string    = par_string,
      h             = h,
      name          = name,
      mesh          = mesh,
      map           = map,
      n_map         = n_map,
      replicate     = replicate,
      n_rep         = n_rep,
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
  pad_add4_space <- paste(rep(" ", padding + 4), collapse = "")

  model_string <- model$model
    if (model_string == "rw" && model$rw_order==1) model_string <- "rw1"
    if (model_string == "rw" && model$rw_order==2) model_string <- "rw2"
  cat(pad_space); cat("Ngme model: "); cat(model_string); cat("\n")

  cat(pad_space); cat("Model parameters: \n")
  params <- with(model, {
    switch(model,
      "ar1"         = paste0(pad_add4_space, ngme_format("K", theta_K, "ar1")),
      "ou"          = paste0(pad_add4_space, ngme_format("K", theta_K, "ou")),
      "re"          = {cat(paste0("  Covariance Matrix = \n"));
                       ngme_format("K", theta_K, "re", x$W_size); NULL},
      "matern"      = paste0(pad_add4_space, ngme_format("K", theta_K, "matern")),
      "rw"          = paste0(pad_add4_space, "No parameter."),
      "unkown"      = paste0(pad_add4_space, "No parameter."),
      "tp" = paste0(pad_add4_space, "left - ", left$model, ": ",
        ngme_format("K", left$theta_K, left$model), "\n",
        pad_add4_space, "right - ", right$model, ": ",
        ngme_format("K", right$theta_K, right$model)
      ),
      paste0(pad_add4_space, "Not implemented yet!")
    )
  })
  cat(params);
  cat("\n\n")

  print.ngme_noise(model$noise, padding = padding)
  invisible(model)
}

