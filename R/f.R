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
#' @param effect_type  if it is a random effect, see ngme_randeff_types()
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
  map         = NULL,
  model       = NULL,
  noise       = noise_normal(),
  replicate   = NULL,
  mesh        = NULL,
  control     = control_f(),
  name        = NULL,
  data        = NULL,
  group       = NULL,
  theta_K     = NULL,
  W           = NULL,
  debug       = FALSE,
  fix_W       = FALSE,
  fix_theta_K = FALSE,
  eval        = FALSE,
  ...
) {
  map <- eval(substitute(map), envir = data, enclos = parent.frame())
  if (inherits(map, "formula")) { # random effect
    name <- "effect"
    model <- "re"
    map <- model.matrix(map, data)
  }
  stopifnot("Please provide model from ngme_model_types():"
    = !is.null(model))

  if (is.null(name)) name <- "field"

  stopifnot("Please specify model as character" = is.character(model))

  replicate <- eval(substitute(replicate), envir = data, enclos = parent.frame())
  replicate <- if (is.null(replicate)) rep(1, length_map(map))
    else as.integer(as.factor(replicate))

  # pre-model
  if (!eval) {
    args <- match.call()
    args$name <- name
    args$map <- map
    args$model <- model
    args$replicate <- replicate
    args$n_repl <- length(unique(replicate))
    return (args)
  }

  # remove NULL in arguments
  f_args <- Filter(Negate(is.null),  as.list(environment()))
  # add arguments in ...
  f_args <- c(f_args, list(...))

  # 0. build mesh if not specified
  if (is.null(mesh)) {
    f_args$mesh <- do.call(build_mesh, f_args)
  }

  # 1. build operator
  n <- mesh$n; nrep <- length(unique(replicate))
  operator <- switch(model,
    tp = do.call(tp, f_args),
    ar1 = do.call(ar1, f_args),
    matern = do.call(matern, f_args),
    iid = list(K = Matrix::Diagonal(n * nrep), h = rep(1, n * nrep)),
    stop("Unknown models")
  )

  A <- if (is.null(operator$A))
    INLA::inla.spde.make.A(mesh = f_args$mesh, loc = map)
    else operator$A

  # 2. build noise given operator
  # if (model == "bv")
  #   noise <- list(
  #     update_noise(noise, ope = operator$first),
  #     update_noise(noise, ope = operator$second)
  #   )
  # else
  noise <- update_noise(noise, ope = operator)

  ngme_model(
    model     = model,
    operator  = operator,
    noise     = noise,
    W_size    = ncol(operator$K),
    V_size    = nrow(operator$K),
    theta_K   = operator$theta_K,
    A         = A,
    control   = control,
    map       = map,
    n_map     = length_map(map),
    replicate = replicate,
    W         = W,
    fix_W     = fix_W,
    name      = name
  )
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
  map         = NULL,
  replicate   = NULL,
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
    theta_K     = theta_K
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

