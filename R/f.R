#' Specifying a latent process model (wrapper function for each model)
#'
#' Function used for defining of smooth and spatial terms
#' within ngme model formulae.
#' The function is a wrapper function for specific submodels.
#' (see ngme_models_types() for available models).
#'
#' @param map  symbol or numerical value: index or covariates to build index
#' @param model  string, model type, see ngme_model_types()
#' @param noise  can be either string or a ngme_noise object
#' @param mesh   mesh for the model, if not provided, will be built from map
#' @param control  control variables for latent model
#' @param name   name of the field, for later use, if not provided, will be "field1" etc.
#' @param data      specifed or inherit from ngme() function
#' @param group   group factor indicate resposne variable, can be inherited from ngme() function
#' @param which_group  belong to which group
#' @param W      starting value of the process
#' @param fix_W  stop sampling for W
#' @param fix_theta_K fix the estimation for theta_K.
#' @param debug     debug mode
#' @param subset    subset of the model
#' @param ...       additional arguments (e.g. parameters for model)
#'  inherit the data from ngme function
#'
#' @return a list for constructing latent model, e.g. A, h, C, G,
#' which also has
#' 1. Information about K matrix
#' 2. Information about noise
#' 3. Control variables
#'
#' @export
f <- function(
  map,
  model,
  noise       = noise_normal(),
  mesh        = NULL,
  control     = control_f(),
  name        = "field",
  data        = NULL,
  group       = NULL,
  which_group = NULL,
  W           = NULL,
  fix_W       = FALSE,
  fix_theta_K = FALSE,
  subset      = rep(TRUE, length_map(map)),
  debug       = FALSE,
  ...
) {
  map <- eval(substitute(map), envir = data, enclos = parent.frame())

  if (inherits(map, "formula")) {
    if (model == "re") {
      map <- model.matrix(map, data)
    } else {
      map <- model.matrix(map, data)[, -1]
    }
  }

  # set the subset if provide group and which_group
  if (!is.null(which_group)) {
    stopifnot(
      "Please provide group factor" = !is.null(group),
      "Please check if which_group is in group"
        = which_group %in% levels(as.factor(group)))
    subset <- group %in% which_group
  }

  stopifnot(
    "Please provide model from ngme_model_types():" = !is.null(model),
    "Please specify model as character" = is.character(model)
  )

  if (model == "tp") {
    stopifnot(is.list(map) && length(map) == 2,
      "Please specify map for 2 sub_models"
        = is.list(map) && length(map) == 2)
  }

  # 0. build mesh if not specified
  if (is.null(mesh)) {
    mesh <- ngme_build_mesh(map, model)
  }

  # remove NULL in arguments
  f_args <- Filter(Negate(is.null),  as.list(environment()))
  # add arguments in ...
  f_args <- c(f_args, list(...))

  # 1. build operator
  operator <- build_operator(model, f_args)

  A <- switch(model,
    "tp" = {
      stopifnot("Now only support first to be 1d model"
        = inherits(operator$first$mesh, "inla.mesh.1d"))
      group <- as.integer(as.factor(map[[1]]))
      A <- INLA::inla.spde.make.A(loc=map[[2]], mesh=operator$second$mesh, repl=group)
    },
    "bv" = INLA::inla.spde.make.A(loc=map, mesh=mesh, repl=as.integer(as.factor(group))),
    "re" = ngme_as_sparse(operator$B_K),
    INLA::inla.spde.make.A(mesh = mesh, loc = map)
  )

  # subset the A matrix
  # if (!all(subset)) A[!subset, ] <- 0
  if (!all(subset)) A <- A[subset, , drop = FALSE]

  # 2. build noise given operator
  # bivariate noise
  if (model == "bv") {
    stopifnot(
      "Please specify noise for each field" = length(noise) >= 2,
      "Input: noise=list(a=<noise>,b=<noise>)" = inherits(noise[[1]], "ngme_noise"),
      "Input: noise=list(a=<noise>,b=<noise>)" = inherits(noise[[2]], "ngme_noise"),
      "Please specify noise with same name as in the sub_models argument!"
        = all(names(noise[1:2]) %in% operator$model_names),
      "Keep the noise same if you want to specify single V for each noise"
        = noise[[1]]$single_V == noise[[2]]$single_V
    )
    noise1 <- update_noise(noise[[operator$model_names[[1]]]],
      n=length(operator$h)/2)
    noise2 <- update_noise(noise[[operator$model_names[[2]]]],
      n=length(operator$h)/2)
    bv_noises <- list(noise1, noise2); names(bv_noises) <- operator$model_names
    share_V <- !is.null(noise$share_V) && noise$share_V
    if (share_V) {
      stopifnot(
        "share_V option is only supported for 2 NIG noise."
          = noise1$noise_type == noise2$noise_type,
        "share_V option requires nu from both noise are same."
          = noise1$nu == noise2$nu
      )
    }

    noise <- ngme_noise(
      noise_type = c(noise1$noise_type, noise2$noise_type),
      B_mu    = as.matrix(Matrix::bdiag(noise1$B_mu, noise2$B_mu)),
      B_sigma = as.matrix(Matrix::bdiag(noise1$B_sigma, noise2$B_sigma)),
      theta_mu    = c(noise1$theta_mu, noise2$theta_mu),
      theta_sigma = c(noise1$theta_sigma, noise2$theta_sigma),
      nu = c(noise1$nu, noise2$nu),
      share_V = !is.null(noise$share_V) && noise$share_V,
      single_V = noise1$single_V,
      bv_noises = bv_noises,
      fix_V = !is.null(noise$fix_V) && noise$fix_V,
      V = if (!is.null(noise$V)) noise$V else NULL
    )
  } else {
    noise <- update_noise(noise, n = length(operator$h))
  }

  if (noise$share_V && model!="bv")
    stop("Not allow for share_V for univariate model")

  if (model == "re") {
    noise$fix_theta_sigma <- TRUE
  }

  ngme_model(
    model     = model,
    operator  = operator,
    noise     = noise,
    W_size    = ncol(operator$K),
    V_size    = nrow(operator$K),
    A         = A,
    control   = control,
    map       = map,
    n_map     = length_map(map),
    W         = W,
    fix_W     = fix_W,
    name      = name,
    debug     = debug
  )
}

# build operator
build_operator <- function(model_name, args_list) {
  stopifnot(
    is.character(model_name),
    is.list(args_list)
  )

  switch(model_name,
    tp = do.call(tp, args_list),
    bv = do.call(bv, args_list),
    ar1 = do.call(ar1, args_list),
    rw1 = do.call(rw1, args_list),
    rw2 = do.call(rw2, args_list),
    ou = do.call(ou, args_list),
    matern = do.call(matern, args_list),
    iid = do.call(iid, args_list),
    re = do.call(re, args_list),
    stop("Unknown models")
  )
}

# help to build a list of mesh for different replicates
ngme_build_mesh <- function(
  loc,
  model = NULL,
  ...
) {
  if (inherits(loc, "inla.mesh.1d") || inherits(loc, "inla.mesh")) return(loc)
  if (!is.null(model)) {
    if (model %in% c("re", "tp")) return(NULL)
    if (model == "ar1") {
      stopifnot("The map should be integers."
        = is.numeric(loc) && all(loc == round(loc)))
    }
  }

  if (is.matrix(loc) && ncol(loc) == 2) {
    stop("Please build and provide the mesh for spatial data using inla.mesh.2d()")
  } else if (is.numeric(loc)) {
    mesh <- INLA::inla.mesh.1d(loc = loc)
  } else {
    stop("loc should be a matrix or a vector of numeric")
  }

  mesh
}