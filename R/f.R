#' Specifying a latent process model
#'
#' Function used for defining of smooth and spatial terms
#' within ngme model formulae.
#' The function does not evaluate anything -
#' it exists purely to help set up a model.
#' (see ngme.models.types for available models).
#'
#' @param index    symbol or numerical value
#' @param model     1. string: type of model, 2. ngme.spde object
#' @param noise     1. string: type of model, 2. ngme.noise object
#' @param replicates   Representing the replicates
#'  (can also be specified in each ngme model)
#' @param control      control variables for f model
#' @param A            A Matrix connecting observation and mesh
#' @param theta_K      Unbounded parameter for K
#' @param debug        Debug variables
#' @param data      specifed or inherit from ngme formula
#' @param W         starting value of the process
#' @param A_pred    A Matrix connecting NA location and mesh
#' @param ...       additional arguments
#'  inherit the data from ngme function
#'
#' @return a list latent_in for constructing latent model, e.g. A, h, C, G,
#' which also has
#' 1. list operator for building operator,
#' 2. list var_in for variance component,
#'
#' @export
f <- function(
  index       = NULL,
  replicates  = NULL,
  model       = "ar1",
  noise       = ngme.noise.normal(),
  data        = NULL,
  control     = ngme.control.f(),
  debug       = FALSE,
  A           = NULL,
  A_pred      = NULL,
  theta_K     = NULL,
  W           = NULL,
  fix_W       = NULL,
  fix_theta_K = NULL,
  index_pred  = NULL,
  ...
) {
  arg_list <- as.list(environment(), list(...))
  ################## Index and Replicates ##################

  # deprecated
    # get the index from ...
    # ddd <- match.call(expand.dots = FALSE)$...
    # if (length(ddd) > 1) {
    #   stop(paste("To many variables included in f():", paste(unlist(ddd), collapse=" ")))
    # }
    # numerical.or.symbol <- ddd[[1]]
# print(match.call())
  val_or_sym <- substitute(index) # capture symbol if it is
  index <- eval(val_or_sym, envir = data)

  # get index -> then make both A and A_pred matrix
  ngme_response <- data$ngme_response
  index_data <- which(!is.na(ngme_response))

  # stopifnot("response is null" = !is.null(ngme_response))
  if (!is.null(ngme_response) && any(is.na(ngme_response))) {
    index_NA   <- which(is.na(ngme_response))
    # ignore the NA position in the provided index
    index <- Filter(function(x) !(x %in% index_NA), index)
  } else {
    # no need to predict
    A_pred <- index_NA <- NULL
  }

  # # get the replicates
  # if (is.null(replicates))
  #   replicates <- rep(1, length(index))
  # nrep <- length(unique(replicates))

  # no need for re-order
  # # re-order the values according to the replicates (to be block diagonal for C and G)
  # df <- data.frame(original.order=1:length(index), replicates=replicates, index=index)
  # df <- df[order(df$replicates), ]


  ################## Construct operator (n_ope, C, G, A, h) ##################
  nrep <- NULL
  W_size <- NULL

  if (is.character(model)) {  ######## string
    if (model == "ar1") {
      stopifnot("Please specify the init value for ar1 using theta_K=."
        = is.null(theta_K) || (theta_K > 0 && theta_K < 1))

      # provide the natural index
      if (is.null(index)) index <- index_data

      model_list <- do.call(ngme.ar1, list(
        alpha = theta_K,  # watch out! Do not pass using theta_K = theta_K!
        index = index,
        replicates = replicates,
        range = c(1, max(length(ngme_response), max(index))), # watch out! using natural mesh
        index_pred = index_NA
      ))
    } else if (model == "rw1") {
      model_list <- do.call(ngme.rw1, list(
        index = index,
        replicates = replicates
      ))
    } else stop("unknown model name")
  }
  else if (inherits(model, "ngme_model")) {  # stationary Matern
    model_list <- model
    arg_list$model <- NULL # watch out! otherwise overwrite again.
  }
    # watch out! structure not so good.
    # nrep <- ncol(A) / nrow(operator$C)
    # operator$C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), operator$C)
    # operator$G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), operator$G)
    # operator$C <- ngme.as.sparse(operator$C)
    # operator$G <- ngme.as.sparse(operator$G)
    # W_size <- model$operator$W_size
    # # h <- rep(1, W_size)
    # h <- diag(model$operator$C)
    # theta_K <- model$operator$theta_kappa
    # h <- rep(h, times = nrep)

    # model_type <- "spde.matern"

    # if (is.null(A)) stop("Provide A matrix")
    # # compute nrep
    # nrep <- ncol(A)/nrow(operator$C)
    # operator$C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), operator$C)
    # operator$G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), operator$G)
    # operator$C <- ngme.as.sparse(operator$C)
    # operator$G <- ngme.as.sparse(operator$G)

    # W_size = length(index)

    # h <- rep(1, W_size)
    # theta_K = model$theta_Kappa
    # operator <- model$operator
  else {
    stop("unknown model")
  }

  if (is.null(W) && !is.null(fix_W)) stop("Provide initial W to use fix_W")

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
  n_params = model_list$n_theta_K + noise$n_theta_mu + noise$n_theta_sigma + noise$n_theta_V
  # check initial value of W
  if (!is.null(W)) stopifnot(length(W) == W_size)
  # get the useful argument list
# print(str(arg_list))
# print(str(Filter(Negate(is.null), arg_list)))
  model_list <- modifyList(model_list, Filter(Negate(is.null), arg_list)) # watch out! arg_list$noise$ = NULL; nested NULL
  # modify model_list
  model_list$noise <- with(model_list, update.ngme.noise(noise, V_size))
  model_list$noise_type <- model_list$noise$noise_type
  model_list$n_params <- n_params

  do.call(ngme.model, model_list)
}