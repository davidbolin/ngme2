
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
  replicate   = NULL,
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

  stopifnot("replicate is NULL" = !is.null(replicate))
  stopifnot("make sure length of replicate == length of index" =
   length(replicate) == length_map(map))
  replicate <- as.integer(replicate)
  if (is.null(n_params)) n_params <- length(operator$theta_K) + with(noise, n_params)

  structure(
    list(
      model         = model,
      operator      = operator,
      noise_type    = noise$noise_type,
      W_size        = W_size,
      theta_K       = theta_K,
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
      replicate     = replicate,
      n_rep         = length(unique(replicate)),
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
  print.ngme_noise(model$noise, padding = padding)
  invisible(model)
}

