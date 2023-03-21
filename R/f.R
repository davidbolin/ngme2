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
#' @param eval      evaluate the model
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
  model,
  map         = NULL,
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
  stopifnot("Please specify model as character" = is.character(model))

  map <- eval(substitute(map), envir = data, enclos = parent.frame())
  if (model == "tp") map <- 1:(list(...)$left$n_map * list(...)$right$n_map)

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
    stopifnot("Please use models in ngme_model_types()"
      = model %in% ngme_model_types())
    args <- within(f_args, rm(model))
    f_model <- switch(model,
      "ar1" = {
        do.call(model_ar1, args)
      },
      "rw" = {
        do.call(model_rw, args)
      },
      "ou" = {
        do.call(model_ou, args)
      },
      "matern" = {
        do.call(model_matern, args)
      },
      "tp" = {
        do.call(model_tp, args)
      },
      "iid" = {
        do.call(model_iid, args)
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

  # update n noise
  f_model$noise <- update_noise(f_model$noise, n = f_model$V_size)
  f_model
}

#' ngme tensor-product model specification
#'
#' Given 2 models (left and right), build a tensor-product model based on K = K_left x K_right (here x is Kronecker product)
#'
#' @param left ngme_model
#' @param right ngme_model
#' @param map can be ignored, pass through left and right
#' @param replicate replicate for the process
#' @param index_NA Logical vector, same as is.na(response var.)
#'
#' @param noise noise, can be specified in f()
#' @param data data, can be specified in f(), ngme()
#' @param control control for the model
#' @param ... extra arguments in f()
#'
#' @return a list of specification of model
#' @export
model_tp <- function(
  left        = NULL,
  right       = NULL,
  map         = NULL,
  replicate   = NULL,
  data        = NULL,
  index_NA    = NULL,
  noise       = noise_normal(),
  control   = control_f(),
  ...
) {
  stopifnot(inherits(left, "ngme_model"), inherits(right, "ngme_model"))
  n_map <-  length_map(left$map) * length_map(right$map)
  if (is.null(map)) map <- 1:n_map
  # if (is.null(replicate))
    replicate <- rep(1, n_map)

  stopifnot(left$model %in% c("rw1", "rw2", "ar1", "iid"))
  f_model <- ngme_model(
    model = "tp",
    map = map,
    n_map = n_map,
    replicate = replicate,
    data = data,
    index_NA = index_NA,
    noise = noise,
    left = left,
    right = right,
    control = control,
    W_size = left$W_size * right$W_size,
    V_size = left$V_size * right$V_size,
    n_theta_K = left$n_theta_K + right$n_theta_K,
    n_params = left$n_params + right$n_params,
    theta_K = c(left$theta_K, right$theta_K),
    A = left$A %x% right$A,
    ...
  )

# init noise
  f_model$left$noise$init_V <- FALSE
  f_model$right$noise$init_V <- FALSE
  f_model$noise <- update_noise(noise, n = f_model$V_size)
  f_model$h <- f_model$noise$h

  f_model
}

#' ngme iid model specification
#'
#' @param map integer vector, time index for the AR(1) process
#' @param replicate replicate for the process
#' @param index_NA Logical vector, same as is.na(response var.)
#'
#' @param noise noise, can be specified in f()
#' @param data data, can be specified in f(), ngme()
#' @param control controls using control_f(),
#' @param ... extra arguments in f()
#'
#' @return a list of specification of model
#' @export
model_iid <- iid <- function(
  map = NULL,
  replicate  = NULL,
  data        = NULL,
  index_NA    = NULL,
  noise       = noise_normal(),
  control     = control_f(),
  ...
) {
  # capture symbol in index
  map <- eval(substitute(map), envir = data, enclos = parent.frame())
  n_map <- length_map(map)
  stopifnot("The map should be integers." = all(map == round(map)))
  if (is.null(replicate)) replicate <- rep(1, length_map(map))
  if (is.null(index_NA)) index_NA <- rep(FALSE, length_map(map))

  stopifnot("Make sure length(idx)==length(replicate)" = length(map) == length(replicate))

  replicate <- if (!is.null(list(...)$replicate)) list(...)$replicate
    else rep(1, length(map))

  mesh <- INLA::inla.mesh.1d(loc=map)
  tmp <- ngme_make_A(
    mesh = mesh,
    map = map,
    n_map = n_map,
    idx_NA = index_NA,
    replicate = replicate
  )
  A <- tmp$A; A_pred <- tmp$A_pred

  args <- within(list(...), {
    mesh        = mesh
    map         = map
    n_map       = n_map
    model       = "iid"
    theta_K     = 0
    W_size      = n_map
    V_size      = n_map
    A           = A
    A_pred      = A_pred
    h           = rep(1, n_map)
    K           = ngme_as_sparse(Matrix::Diagonal(n_map))
    noise       = noise
    replicate   = replicate
    n_rep       = length(unique(replicate))
    control     = control
  })
  do.call(ngme_model, args)
}

