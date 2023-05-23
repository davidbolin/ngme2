#' ngme bivariate model specification
#'
#' Given 2 operator (first and second), build a correlated bivaraite operator based on K = D %*% diag(K_first, K_second)
#'
#' @param map can be ignored, pass through first and second
#' @param replicate replicate
#' @param theta_K c(zeta, eta, theta_K_1, theta_K_2)
#' @param ... extra arguments in f()
#'
#' @return a list of specification of model
#' @export
bv <- function(
  map,
  sub_models,
  mesh = NULL,
  zeta = 0, eta = 0,
  replicate = NULL,
  group = NULL,
  share_param = FALSE,
  ...
) {
  if (inherits(map, "formula")) map <- model.matrix(map)[, -1]
  if (is.null(replicate)) replicate <- rep(1, length_map(map))
  if (is.null(mesh)) mesh <- ngme_build_mesh(map)

  model_names <- names(sub_models)
  stopifnot(
    "Please provide names for sub_models" = !is.null(model_names),
    "Please provide group argument to indicate different fields"
      = !is.null(group),
    "Make sure replicate is independent accross the f() models"
      = length(unique(replicate)) == 1,
    "Length of sub_models should be 2" = length(sub_models) == 2,
    "Name of sub_models should be in group"
    = all(model_names %in% levels(as.factor(group)))
  )
  group <- factor(group, levels = model_names)

  # build 2 sub_models
  args <- as.list(environment())
  arg1 <- sub_models[[model_names[1]]]
  if (!is.list(arg1)) {
    first <- build_operator(arg1, args)
  } else {
    stopifnot("Please provide model=.. in the list" = !is.null(arg1$model))
    first <- build_operator(arg1$model, modifyList(args, arg1))
  }
  arg2 <- sub_models[[model_names[2]]]
  if (!is.list(arg2)) {
    second <- build_operator(arg2, args)
  } else {
    stopifnot("Please provide model=.. in the list" = !is.null(arg2$model))
    second <- build_operator(arg2$model, modifyList(args, arg2))
  }

  theta_K <- c(zeta, eta, first$theta_K, second$theta_K)

  # pass the theta_K to first and second
  first$theta_K <- theta_K[3:(2 + first$n_theta_K)]
  second$theta_K <- theta_K[(3 + first$n_theta_K):length(theta_K)]

  D <- build_D(theta_K[1], theta_K[2])
  bigD <- kronecker(D, Matrix::Diagonal(nrow(first$K)))

  ngme_operator(
    map = first$map,
    mesh = NULL,
    n_rep = length(unique(replicate)),
    model       = "bv",
    first       = first,
    second      = second,
    theta_K     = theta_K,
    K           = bigD %*% Matrix::bdiag(first$K, second$K),
    A           = INLA::inla.spde.make.A(loc=map, mesh=mesh, repl=as.integer(as.factor(group))),
    h           = c(first$h, second$h),
    symmetric   = FALSE,
    zero_trace  = FALSE,
    model_names = model_names,
    share_param = share_param
  )
}


#' ngme tensor-product model specification
#'
#' Given 2 operator (first and second), build a tensor-product operator based on K = K_first x K_second (here x is Kronecker product)
#'
#' @param first left side of kronecker model (usually a temporal or iid model)
#' @param second right side of kronecker model (ususally a spatial model)
#' @param ... extra arguments in f()
# ' @param map pass through first and second
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

  stopifnot("the length of map of 2 submodel should be equal (complete)"
    = length_map(first$map) == length_map(second$map))

  map <- second$map
  group <- if (!is.matrix(first$map)) as.integer(as.factor(first$map))
    else as.integer(as.factor(seq_len(nrow(first$map))))
  A <- INLA::inla.spde.make.A(loc=map, mesh=second$mesh, repl=group)

  ngme_operator(
    map = map,  # placeholder
    mesh = NULL,
    n_rep = 1,
    model = "tp",
    first = first,
    second = second,
    theta_K = theta_K,
    K = first$K %x% second$K,
    h = first$h %x% second$h,
    A = A,
    symmetric = first$symmetric & second$symmetric,
    # check here
    zero_trace = first$zero_trace & second$zero_trace
  )
}
