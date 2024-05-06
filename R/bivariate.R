#' Ngme bivariate model specification
#'
#' Giving 2 sub_models, build a correlated bivaraite operator based on K = D(theta, eta) %*% diag(K_1, K_2)
#' \deqn{D(\theta, \rho) = \begin{pmatrix}
#'   \cos(\theta) + \rho \sin(\theta) & -\sin(\theta) \sqrt{1+\rho^2} \\
#'   \sin(\theta) - \rho \cos(\theta) & \cos(\theta) \sqrt{1+\rho^2}
#' \end{pmatrix}}
#'
#' @param mesh the mesh where model is defined
#' @param sub_models a list of sub_models (total 2 sub_models)
#' @param mesh mesh for build the model
#' @param group group vector, can be inherited from ngme() function
#' @param theta the rotation parameter (starting point should
#'   from -pi/4 to pi/4)
#' @param rho the parameter related to correlation
#' @param share_param TRUE if share the same parameter for 2 sub_models (of same type)
#' @param ... ignore
#'
#' @return a list of specification of model
#' @export
bv <- function(
  mesh,
  sub_models,
  theta = 0,
  rho = 0,
  group = NULL,
  share_param = FALSE,
  fix_bv_theta = FALSE,
  ...
) {
  mesh <- ngme_build_mesh(mesh)
  model_names <- sort(names(sub_models))

  # read group argument from parent frame
  if (is.null(group) &&
      exists("group", parent.frame())){
    group = get("group", parent.frame())
  }

  stopifnot(
    "Please provide names for sub_models" = !is.null(model_names),
    "Please provide group argument to indicate different fields"
      = !is.null(group),
    "Length of sub_models should be 2" = length(sub_models) == 2,
    "Name of sub_models should be in group"
    = all(model_names %in% levels(as.factor(group)))
  )
  group <- factor(group, levels = model_names)

  # build 2 sub_models (pass environment to sub_models)
  # delete theta and rho in env_args (to avoid sub_models use them)
  env_args <- as.list(environment())
  env_args$theta <- NULL; env_args$rho <- NULL
  arg1 <- sub_models[[model_names[1]]]
  if (!is.list(arg1)) {
    first <- build_operator(arg1, env_args)
  } else {
    stopifnot("Please provide model=.. in the list" = !is.null(arg1$model))
    first <- build_operator(arg1$model, modifyList(env_args, arg1))
  }
  arg2 <- sub_models[[model_names[2]]]
  if (!is.list(arg2)) {
    second <- build_operator(arg2, env_args)
  } else {
    stopifnot("Please provide model=.. in the list" = !is.null(arg2$model))
    second <- build_operator(arg2$model, modifyList(env_args, arg2))
  }

  stopifnot(
    "the theta should be in (-pi/2, pi/2)"
      = theta >= -pi/2 & theta <= pi/2
  )
  eta <- tan(theta)
  theta_K <- c(eta, rho, first$theta_K, second$theta_K)

  # pass the theta_K to first and second
  first$theta_K <- theta_K[3:(2 + first$n_theta_K)]
  second$theta_K <- theta_K[(3 + first$n_theta_K):length(theta_K)]

  D <- build_D(theta, rho)
  bigD <- kronecker(D, Matrix::Diagonal(nrow(first$K)))

  # the range of theta is (-pi/4, pi/4), we expand it to (-pi/2, pi/2)
  # so that the maximum can be reached in either direction
  # in the end, we will NOT transform it back to (-pi/4, pi/4)
  eta_to_theta <- function(eta) atan(eta)

  update_K <- function(theta_K) {
    eta <- theta_K[1]
    rho <- theta_K[2]
    theta <- atan(eta)
    D <- build_D(theta, rho)
    bigD <- kronecker(D, Matrix::Diagonal(nrow(first$K)))
    theta_K1 <- theta_K[3:(2 + first$n_theta_K)]
    theta_K2 <- theta_K[(3 + first$n_theta_K):length(theta_K)]
    first$K <- first$update_K(theta_K1)
    second$K <- second$update_K(theta_K2)
    bigD %*% Matrix::bdiag(first$K, second$K)
  }
  ngme_operator(
    mesh        = mesh,
    model       = "bv",
    first        = first,
    second      = second,
    theta_K     = theta_K,
    update_K    = update_K,
    K           = bigD %*% Matrix::bdiag(first$K, second$K),
    h           = c(first$h, second$h),
    symmetric   = FALSE,
    zero_trace  = FALSE,
    model_names = model_names,
    share_param = share_param,
    fix_bv_theta = fix_bv_theta,
    param_name  = c("theta", "rho",
      if (!is.null(first$param_name)) paste(first$param_name, "(1st)") else NULL,
      if (!is.null(second$param_name)) paste(second$param_name, "(2nd)") else NULL
    ),
    param_trans = c(list(eta_to_theta, identity), first$param_trans, second$param_trans)
  )
}


#' ngme tensor-product model specification
#'
#' Given 2 operator (first and second), build a tensor-product operator based on K = K_first x K_second (here x is Kronecker product)
#'
#' @param first ngme_operator, left side of kronecker model (usually a temporal or iid model)
#' @param second ngme_operator, right side of kronecker model (ususally a temporal or spatial model)
#' @param ... ignore
#'
#' @return a list of specification of model
#' @export
tp <- function(
  first,
  second,
  ...
) {
  stopifnot(
    inherits(first, "ngme_operator"),
    inherits(second, "ngme_operator")
  )

  theta_K <- c(first$theta_K, second$theta_K)
  stopifnot(length(theta_K) == first$n_theta_K + second$n_theta_K)

  update_K <- function(theta_K) {
    first$K %x% second$K
  }
  ngme_operator(
    mesh = list(first$mesh, second$mesh),
    model = "tp",
    first = first,
    second = second,
    theta_K = theta_K,
    update_K = update_K,
    K = first$K %x% second$K,
    h = first$h %x% second$h,
    symmetric = first$symmetric & second$symmetric,
    zero_trace = first$zero_trace & second$zero_trace,
    param_name = c(first$param_name, second$param_name),
    param_trans = c(first$param_trans, second$param_trans)
  )
}


#' Ngme bivariate model specification 2 (theta=0)
#'
#' Giving 2 sub_models, build a correlated bivaraite operator based on K = D(theta, eta) %*% diag(K_1, K_2)
#' \deqn{D(\theta, \rho) = \begin{pmatrix}
#'   \cos(\theta) + \rho \sin(\theta) & -\sin(\theta) \sqrt{1+\rho^2} \\
#'   \sin(\theta) - \rho \cos(\theta) & \cos(\theta) \sqrt{1+\rho^2}
#' \end{pmatrix}}
#'
#' @param mesh the mesh where model is defined
#' @param sub_models a list of sub_models (total 2 sub_models)
#' @param mesh mesh for build the model
#' @param group group vector, can be inherited from ngme() function
#' @param rho the parameter related to correlation
#' @param share_param TRUE if share the same parameter for 2 sub_models (of same type)
#' @param ... ignore
#'
#' @return a list of specification of model
#' @export
bv2 <- function(
  mesh,
  sub_models,
  rho = 0,
  group = NULL,
  share_param = FALSE,
  fix_bv_theta = FALSE,
  ...
) {
  mesh <- ngme_build_mesh(mesh)
  model_names <- sort(names(sub_models))

  # read group argument from parent frame
  if (is.null(group) &&
      exists("group", parent.frame())){
    group = get("group", parent.frame())
  }

  stopifnot(
    "Please provide names for sub_models" = !is.null(model_names),
    "Please provide group argument to indicate different fields"
      = !is.null(group),
    "Length of sub_models should be 2" = length(sub_models) == 2,
    "Name of sub_models should be in group"
    = all(model_names %in% levels(as.factor(group)))
  )
  group <- factor(group, levels = model_names)

  # build 2 sub_models (pass environment to sub_models)
  # delete theta and rho in env_args (to avoid sub_models use them)
  env_args <- as.list(environment())
  env_args$theta <- NULL;
  env_args$rho <- NULL
  arg1 <- sub_models[[model_names[1]]]
  if (!is.list(arg1)) {
    first <- build_operator(arg1, env_args)
  } else {
    stopifnot("Please provide model=.. in the list" = !is.null(arg1$model))
    first <- build_operator(arg1$model, modifyList(env_args, arg1))
  }
  arg2 <- sub_models[[model_names[2]]]
  if (!is.list(arg2)) {
    second <- build_operator(arg2, env_args)
  } else {
    stopifnot("Please provide model=.. in the list" = !is.null(arg2$model))
    second <- build_operator(arg2$model, modifyList(env_args, arg2))
  }

  theta_K <- c(rho, first$theta_K, second$theta_K)

  # pass the theta_K to first and second
  first$theta_K <- theta_K[2:(1 + first$n_theta_K)]
  second$theta_K <- theta_K[(2 + first$n_theta_K):length(theta_K)]

  D <- build_D(0, rho)
  bigD <- kronecker(D, Matrix::Diagonal(nrow(first$K)))

  update_K <- function(theta_K) {
    rho <- theta_K[1]
    D <- build_D(0, rho)
    bigD <- kronecker(D, Matrix::Diagonal(nrow(first$K)))

    bigD
  }

  ngme_operator(
    mesh        = mesh,
    model       = "bv2",
    first        = first,
    second      = second,
    theta_K     = theta_K,
    update_K    = update_K,
    K           = bigD %*% Matrix::bdiag(first$K, second$K),
    h           = c(first$h, second$h),
    symmetric   = FALSE,
    zero_trace  = FALSE,
    model_names = model_names,
    share_param = share_param,
    fix_bv_theta = fix_bv_theta,
    param_name  = c("rho",
      if (!is.null(first$param_name)) paste(first$param_name, "(1st)") else NULL,
      if (!is.null(second$param_name)) paste(second$param_name, "(2nd)") else NULL
    ),
    param_trans = c(list(identity), first$param_trans, second$param_trans)
  )
}
