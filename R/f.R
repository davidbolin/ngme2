#' Specifying a latent process model (wrapper function for each model)
#'
#' Function used for defining of smooth and spatial terms
#' within ngme model formulae.
#' The function is a wrapper function for specific submodels.
#' (see ngme_models_types() for available models).
#'
#' @param x    symbol or numerical value: index or covariates to build index
#' @param model     1. string: type of model, 2. ngme.spde object
#' @param replicates   Representing the replicates
#' @param noise     1. string: type of model, 2. ngme.noise object
#'  (can also be specified in each ngme model)
#' @param control      control variables for f model
#' @param name      name of the field, for later use
#' @param A            A Matrix connecting observation and mesh
#' @param theta_K      Unbounded parameter for K
#' @param data      specifed or inherit from ngme formula
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
  x           = NULL,
  model       = "ar1",
  replicates  = NULL,
  noise       = noise_normal(),
  control     = ngme_control_f(),
  name        = NULL,
  data        = NULL,
  A           = NULL,
  A_pred      = NULL,
  theta_K     = NULL,
  W           = NULL,
  fix_W       = NULL,
  fix_theta_K = NULL,
  index_pred  = NULL,
  debug       = NULL,
  index_NA    = NULL, #indicate prediction location
  ...
) {
  x <- eval(substitute(x), envir = data, enclos = parent.frame())
  index <- x

  # if (!is.null(index_NA) && length(index_NA) != length(index))
  #   stop("index_NA length seems wrong.")

  # deal with NA, from ngme function
  if (is.null(index_NA)) index_NA <- rep(FALSE, length(index))

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
      }
    )
  } else {
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
    if (!is.null(f_model$index_NA)) {
      tmp <- with(f_model, {
        ngme_make_A(mesh, map, n_map, idx_NA = index_NA)
      })
      f_model$A <- tmp$A
      f_model$A_pred <- tmp$A_pred
    }

  # print(str(f_model))
  }

  # get index -> then make both A and A_pred matrix
  # ngme_response <- data$ngme_response
  # index_data <- which(!is.na(ngme_response))

  # # stopifnot("response is null" = !is.null(ngme_response))
  # if (!is.null(ngme_response) && any(is.na(ngme_response))) {
  #   index_NA   <- which(is.na(ngme_response))
  #   # ignore the NA position in the provided index
  #   index <- Filter(function(x) !(x %in% index_NA), index)
  # } else {
  #   # no need to predict
  #   A_pred <- index_NA <- NULL
  # }

  # # get the replicates
  # if (is.null(replicates))
  #   replicates <- rep(1, length(index))
  # nrep <- length(unique(replicates))

  # no need for re-order
  # # re-order the values according to the replicates (to be block diagonal for C and G)
  # df <- data.frame(original.order=1:length(index), replicates=replicates, index=index)
  # df <- df[order(df$replicates), ]


  ################## construct noise (e.g. nig noise) ##################
    # ?? check
    # B_mu <- matrix(noise$B_mu, nrow = W_size, ncol = noise$n_theta_mu)
    # B_sigma <- matrix(noise$B_sigma, nrow = W_size, ncol = noise$n_theta_sigma)

    # # replicates
    # if (is.integer(nrep)) {
    #   B_mu <- kronecker(matrix(1, ncol = 1, nrow = nrep), B_mu)
    #   B_sigma <- kronecker(matrix(1, ncol = 1, nrow = nrep), B_sigma)
    # }

  # total params
  # n_params = model_list$n_theta_K + noise$n_theta_mu + noise$n_theta_sigma + noise$n_nu
  # check initial value of W
#   if (!is.null(W)) stopifnot(length(W) == W_size)
#   # get the useful argument list
# # print(str(arg_list))
# # print(str(Filter(Negate(is.null), arg_list)))
#   model_list <- modifyList(model_list, Filter(Negate(is.null), arg_list)) # watch out! arg_list$noise$ = NULL; nested NULL
#   # modify model_list
#   model_list$noise <- with(model_list, update_noise(noise, V_size))
#   model_list$noise_type <- model_list$noise$noise_type
#   model_list$n_params <- n_params

#   do.call(ngme_model, model_list)

  # update n noise
  f_model$noise <- update_noise(f_model$noise, n = f_model$V_size)
  f_model
}