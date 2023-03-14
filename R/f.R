#' Specifying a latent process model (wrapper function for each model)
#'
#' Function used for defining of smooth and spatial terms
#' within ngme model formulae.
#' The function is a wrapper function for specific submodels.
#' (see ngme_models_types() for available models).
#'
#' @param map    symbol or numerical value: index or covariates to build index
#' @param model     1. string: type of model, 2. ngme.spde object
#' @param replicate   Representing the replicate
#' @param noise     1. string: type of model, 2. ngme.noise object
#'  (can also be specified in each ngme model)
#' @param control      control variables for f model
#' @param name      name of the field, for later use
#' @param A            A Matrix connecting observation and mesh
#' @param theta_K      Unbounded parameter for K
#' @param data      specifed or inherit from ngme formula
#' @param group    model as the group (see vignette for space-temporal model)
#' @param W         starting value of the process
#' @param A_pred    A Matrix connecting NA location and mesh
#' @param index_NA  Logical vector, same as is.na(response var.)
#' @param fix_W  stop sampling for W
#' @param fix_theta_K fix the estimation for theta_K.
#' @param index_pred index for prediction
#' @param debug        Debug mode
#' @param ...       additional arguments
#'  inherit the data from ngme function
#'
#' @return a list latent_in for constructing latent model, e.g. A, h, C, G,
#' which also has
#' 1. Information about K matrix
#' 2. Information about noise
#' 3. Control variables
#'
#' @export
f <- function(
  map         = NULL,
  model       = "ar1",
  replicate   = NULL,
  noise       = noise_normal(),
  control     = control_f(),
  name        = NULL,
  data        = NULL,
  group       = NULL,
  A           = NULL,
  A_pred      = NULL,
  theta_K     = NULL,
  W           = NULL,
  fix_W       = NULL,
  fix_theta_K = NULL,
  index_pred  = NULL,
  debug       = NULL,
  index_NA    = NULL, #indicate prediction location
  eval        = FALSE,
  ...
) {
  map <- eval(substitute(map), envir = data, enclos = parent.frame())
  replicate <- if (is.null(replicate)) rep(1, length_map(map))
    else as.integer(as.factor(replicate))

  if (is.null(name)) name <- "field"

  # pre-model
  if (!eval) {
    args <- match.call()
    args$name <- name
    args$map <- map
    args$replicate <- replicate
    return (args)
  }

  # if (!is.null(index_NA) && length(index_NA) != length(index))
  #   stop("index_NA length seems wrong.")
  # deal with NA, from ngme function
  if (is.null(index_NA)) index_NA <- rep(FALSE, length_map(map))

  # remove NULL in arguments
  f_args <- Filter(Negate(is.null),  as.list(environment()))
  # add arguments in ...
  f_args <- c(f_args, list(...))

  # combine args and apply ngme sub_models
  if (is.character(model)) {
    args <- within(f_args, rm(model))
    f_model <- switch(model,
      "ar1" = {
        do.call(model_ar1, args)
      },
      "rw1" = {
        args$order = 1
        do.call(model_rw, args)
      },
      "rw2" = {
        args$order = 2
        do.call(model_rw, args)
      },
      "matern" = {
        do.call(model_matern, args)
      }
    )
  } else {
warning("Please use f(model = '...'), which is better")
    stopifnot("please check model specification" =
      inherits(model, "ngme_model"))
    # model is evaluated with submodel func.
    f_args <- within(f_args, rm(model))

    # watch out! if update the noise in f(noise=...)
    if (is.null(as.list(match.call())$noise)) {
      f_args$noise <- NULL
    } else {
      model$noise_type <- NULL
    }
    # use f_args to update the model
    f_model <- do.call(ngme_model, utils::modifyList(model, f_args))
# make V explicitly (prevent no entry called V)
if (is.null(f_model$noise$V)) f_model$noise["V"] <- list(NULL)

    # update A and A_pred
    if (!is.null(map) && !is.null(f_model$index_NA)) {
      tmp <- with(f_model, {
        ngme_make_A(
          mesh = mesh,
          map = map,
          n_map = n_map,
          idx_NA = index_NA,
          replicate = replicate
        )
      })
      f_model$A <- tmp$A
      f_model$A_pred <- tmp$A_pred
    }
  }

  # tensor product
  if (!is.null(group) && inherits(group, "ngme_model")) {
    f_model$model_right <- f_model

    f_model$model <- "tensor_prod"
    f_model$W_size <- f_model$W_size * group$W_size
    f_model$V_size <- f_model$V_size * group$V_size
    f_model$n_theta_K <- f_model$n_theta_K + group$n_theta_K
    f_model$n_params <- f_model$n_params + group$n_theta_K
    f_model$theta_K <- c(group$theta_K, f_model$theta_K)

    f_model$h <- f_model$noise$h <- group$noise$h %x% f_model$noise$h

    stopifnot(group$model %in% c("rw1", "rw2", "ar1", "iid"))

    if (is.null(ncol(f_model$map))) {
      idx_r <- rep(f_model$map, group$n_map)
    } else {
      idx_r <- matrix(ncol=2, nrow=group$n_map * f_model$n_map)
      idx_r[, 1] <- rep(f_model$map[, 1], group$n_map)
      idx_r[, 2] <- rep(f_model$map[, 2], group$n_map)
    }
    idx_l <- rep(group$map, each=f_model$n_map)

    f_model$A <- INLA::inla.spde.make.A(
      mesh = f_model$mesh,
      loc = idx_r,
      group = idx_l
      # group.mesh?
      # index?
    )
    f_model$group$noise$init_V <- FALSE
    f_model$model_right$noise$init_V <- FALSE
  }

  # update n noise
  f_model$noise <- update_noise(f_model$noise, n = f_model$V_size)
  f_model
}