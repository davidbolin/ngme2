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
  group <- as.integer(as.factor(first$map))
  A <- INLA::inla.spde.make.A(loc=map, mesh=second$mesh, repl=group)

  ngme_operator(
    map = map,  # placeholder
    mesh = NULL,
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
#' @param replicate replicate passed to both first and second
#' @param first ngme_model
#' @param second ngme_model
#' @param ... extra arguments in f()
#'
#' @return a list of specification of model
#' @export
bv <- function(
  first,
  second,
  zeta = 0, rho = 0,
  map = NULL,
  replicate = NULL,
  share_param = FALSE,
  ...
) {
  # inject replicate to sub models
  if (!is.null(replicate)) {
    first_call <- match.call()$first
    first_call$replicate <- replicate
    first <- eval(first_call, envir = parent.frame())

    second_call <- match.call()$second
    second_call$replicate <- replicate
    second <- eval(second_call, envir = parent.frame())
  }

  theta_K <- c(zeta, rho, first$theta_K, second$theta_K)
  stopifnot(
    inherits(first, "ngme_operator"),
    inherits(second, "ngme_operator"),
    all(dim(first$K == second$K)),
    length(theta_K) == first$n_theta_K + second$n_theta_K + 2
  )

  # pass the theta_K to first and second
  first$theta_K <- theta_K[3:(2 + first$n_theta_K)]
  second$theta_K <- theta_K[(3 + first$n_theta_K):length(theta_K)]

  D <- build_D(theta_K[1], theta_K[2])
  bigD <- kronecker(D, Matrix::Diagonal(nrow(first$K)))

  ngme_operator(
    model       = "bv",
    first       = first,
    second      = second,
    theta_K     = theta_K,
    K           = bigD %*% Matrix::bdiag(first$K, second$K),
    A           = Matrix::bdiag(first$A, second$A),
    h           = c(first$h, second$h),
    symmetric   = FALSE,
    zero_trace  = FALSE,
    share_param = share_param
  )
}

