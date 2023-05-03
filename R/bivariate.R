#' ngme tensor-product model specification
#'
#' Given 2 operator (first and second), build a tensor-product operator based on K = K_first x K_second (here x is Kronecker product)
#'
#' @param map can be ignored, pass through first and second
#' @param first ngme_model
#' @param second ngme_model
#' @param replicate replicate for the process
#' @param ... extra arguments in f()
#'
#' @return a list of specification of model
#' @export
tp <- function(
  first,
  second,
  map         = NULL,
  replicate   = NULL,
  theta_K     = NULL,
  ...
) {
  stopifnot(
    inherits(first, "ngme_operator"),
    inherits(second, "ngme_operator")
  )
  if (is.null(theta_K)) theta_K <- c(first$theta_K, second$theta_K)
  stopifnot(length(theta_K) == first$n_theta_K + second$n_theta_K)

  ngme_operator(
    model = "tp",
    first = first,
    second = second,
    theta_K = theta_K,
    K = first$K %x% second$K,
    A = first$A %x% second$A,
    h = first$h %x% second$h,
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
#' @param first ngme_model
#' @param second ngme_model
#' @param replicate replicate for the process
#' @param ... extra arguments in f()
#'
#' @return a list of specification of model
#' @export
bv <- function(
  first,
  second,
  rho = 0, zeta = 0,
  map = NULL,
  replicate = NULL,
  share_param = FALSE,
  ...
) {
  theta_K <- c(rho, zeta, first$theta_K, second$theta_K)
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
  bigD <- kronecker(D, Matrix::Diagonal(nrow(first$K))); bigD

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

