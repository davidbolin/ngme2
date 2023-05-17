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


#' ngme bivariate model specification
#'
#' Given 2 operator (first and second), build a correlated bivaraite operator based on K = D %*% diag(K_first, K_second)
#'
#' @param map can be ignored, pass through first and second
#' @param replicate replicate
#' @param theta_K c(zeta, rho, theta_K_1, theta_K_2)
#' @param ... extra arguments in f()
#'
#' @return a list of specification of model
#' @export
bv <- function(
  map,
  mesh = NULL,
  sub_models = c("ar1", "ar1"),
  replicate = NULL,
  group = NULL,
  which_group = levels(as.factor(group)),
  share_param = FALSE,
  ...
) {
  if (inherits(map, "formula")) map <- model.matrix(map)[, -1]
  if (is.null(replicate)) replicate <- rep(1, length_map(map))

  stopifnot(
    "Make sure replicate is independent accross the f() models"
      = length(unique(replicate)) == 1,
    "Length of which_group should be 2" = length(which_group) == 2
  )

  if (is.null(mesh)) mesh <- ngme_build_mesh(map)

  args <- as.list(environment())
  first  <- build_operator(sub_models[1], args)
  second <- build_operator(sub_models[2], args)

  theta_K <- c(0, 0, first$theta_K, second$theta_K)

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
    share_param = share_param
  )
}

