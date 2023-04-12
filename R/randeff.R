# implementation of random effect structure

#' Build random effect structure
#' @param formula formula
#' @param data data
#' @param family distribution of random effect
#' @return a list of random effect structure
#' @export
randeff <- function(
  map,
  data,
  effect_type,
  name,
  # replicate,
  theta_K=NULL,
  n_repl=NULL, # number of replicates
  ...
) {
  stopifnot("Please specify random effect type, using effect_type=...,
    check ngme_randeff_types()" = inherits(effect_type, "ngme_noise") || effect_type %in% ngme_randeff_types())

  # Design Matrix and Covariance matrix
  B_reff <- map
  n_reff <- ncol(B_reff)
  n_cov_params <- sum(1:n_reff)
  if (is.null(theta_K)) theta_K <- rep(0, n_cov_params)
  stopifnot(length(theta_K) == n_cov_params)

  if (is.character(effect_type)) {
    mix_var <- if (effect_type=="normal") noise_normal() else noise_nig()
  } else {
    mix_var <- effect_type
  }
  mix_var <- update_noise(mix_var, n=n_reff)
  mix_var$single_V = TRUE
  mix_var$fix_theta_sigma = TRUE
  # mix_var$B_mu <- matrix(rep(1, length(B_reff), ncol=n_reff))
  # mix_var$theta_mu <- rep(0, n_reff)

  args <- within(list(...), {
    model        = "re"
    map          = B_reff
    A            = ngme_as_sparse(B_reff)
    n_map        = length_map(map)
    noise        = mix_var
    W_size       = n_reff
    V_size       = n_reff
    theta_K      = theta_K
    n_cov_params = n_cov_params
    n_repl       = n_repl
    # n_params     = n_cov_params + 1 + 1 #n_reff + 1 #sigma mu nu
    name         = name
  })
  do.call(ngme_model, args)
}

#' Print ngme random effect
#'
#' @param x random effect object
#' @param padding number of white space padding in front
#' @param ... ...
#'
#' @return radom effect model
#' @export
print.randeff <- function(x, padding = 0, ...) {
  reff <- x
  pad_space <- paste(rep(" ", padding), collapse = "")
  pad_add4_space <- paste(rep(" ", padding + 4), collapse = "")

  if (is.null(reff)) {
    cat(pad_space); cat("Effect type - "); cat("NULL"); cat("\n")
  } else {
    cat(pad_space); cat("Effect type - "); cat(reff$effect_type); cat("\n")

    cat(pad_space); cat("Effect formula: ");
    cat(paste(as.character(reff$formula))); cat("\n")

    cat(pad_space); cat("Effect covariance: \n")
    cat(pad_space); print(reff$Sigma)
  }
  cat("\n")
  invisible(reff)
}