ngme_operator <- function(
  model,
  K,
  h,
  theta_K = NULL,
  zero_trace = FALSE,
  symmetric = FALSE,
  ...
) {

  structure(
    list(
      model = model,
      K = K,
      h = h,
      theta_K = theta_K,
      n_theta_K = length(theta_K),
      zero_trace = zero_trace,
      symmetric = symmetric,
      ...
    ),
    class = "ngme_operator"
  )
}


#' Print ngme operator
#'
#' @param x ngme operator object
#' @param padding number of white space padding in front
#' @param ... ...
#'
#' @return a list (operator specifications)
#' @export
print.ngme_operator <- function(x, padding = 0, ...) {
  operator <- x
  pad_space <- paste(rep(" ", padding), collapse = "")
  pad_add4_space <- paste(rep(" ", padding + 4), collapse = "")

  model_name <- switch(operator$model,
    ar1 = "AR(1)",
    "Unknown"
  )

  model_name <- paste("Model type: ", model_name, "\n", sep = "")
  cat(pad_space, model_name, sep="")

  parameter <- with(operator, switch(model,
    ar1 = {paste0("alpha = ", format(ar1_th2a(theta_K), digits=3), "\n")},
    "Unknown"
  ))
  cat(pad_add4_space, parameter, sep="")

  # model_string <- model$model
  #   if (model_string == "rw" && model$rw_order==1) model_string <- "rw1"
  #   if (model_string == "rw" && model$rw_order==2) model_string <- "rw2"
  # cat(pad_space); cat("Ngme model: "); cat(model_string); cat("\n")

  # cat(pad_space); cat("Model parameters: \n")
  # params <- with(model, {
  #   switch(model,
  #     "ar1"         = paste0(pad_add4_space, ngme_format("K", theta_K, "ar1")),
  #     "ou"          = paste0(pad_add4_space, ngme_format("K", theta_K, "ou")),
  #     "re"          = {cat(paste0("  Covariance Matrix = \n"));
  #                      ngme_format("K", theta_K, "re", x$W_size); NULL},
  #     "matern"      = paste0(pad_add4_space, ngme_format("K", theta_K, "matern")),
  #     "rw"          = paste0(pad_add4_space, "No parameter."),
  #     "unkown"      = paste0(pad_add4_space, "No parameter."),
  #     "tp" = paste0(pad_add4_space, "left - ", left$model, ": ",
  #       ngme_format("K", left$theta_K, left$model), "\n",
  #       pad_add4_space, "right - ", right$model, ": ",
  #       ngme_format("K", right$theta_K, right$model)
  #     ),
  #     "bv" = paste0(pad_add4_space, "first - ", m1$model, ": ",
  #       ngme_format("K", m1$theta_K, m1$model), "\n",
  #       pad_add4_space, "second - ", m2$model, ": ",
  #       ngme_format("K", m2$theta_K, m2$model), "\n",
  #       "theta : ", theta_K[1], "\n",
  #       "rho : ", theta_K[2]
  #     ),
  #     paste0(pad_add4_space, "Not implemented yet!")
  #   )
  # })
  # cat(params);
  # cat("\n\n")

  # print.ngme_noise(model$noise, padding = padding)
  invisible(operator)
}

matern <- function(
  map,
  mesh,
  replicate = rep(1, length_map(map)),
  alpha = 2,
  theta_kappa = 0,
  B_kappa = NULL
) {
  # build C, G, K, A
  n <- mesh$n; nrep <- length(unique(replicate))

  stopifnot(alpha == 2 || alpha == 4)

  if (is.null(B_kappa) && length(theta_kappa) == 1)
    B_kappa <- matrix(1, nrow = mesh$n, ncol = 1)
  else if (is.null(B_kappa) && length(theta_kappa) > 1)
    stop("Please provide B_kappa for non-stationary case.")

  d <- get_inla_mesh_dimension(mesh)
  if (d == 1) {
    fem <- INLA::inla.mesh.1d.fem(mesh)
    C <- fem$c1
    G <- fem$g1
    h <- Matrix::diag(fem$c0)
    # h <- rowsum(fem$c0)
  } else {
    fem <- INLA::inla.mesh.fem(mesh, order = alpha)
    C <- fem$c0  # diag
    G <- fem$g1
    h <- Matrix::diag(fem$c0)
  }

  kappas <- as.numeric(B_kappa %*% theta_kappa)
  K <- if (alpha == 2) diag(kappas) %*% C %*% diag(kappas)  + G
    else diag(kappas) %*% C %*% diag(kappas) %*% C %*% diag(kappas) + G
  ngme_operator(
    model = "matern",
    C = ngme_as_sparse(C),
    G = ngme_as_sparse(G),
    K = ngme_as_sparse(K),
    h = h,
    A = INLA::inla.spde.make.A(mesh = mesh, loc = map)
  )
}

ar1 <- function(
  map,
  mesh      = INLA::inla.mesh.1d(loc = min(map):max(map)),
  replicate = rep(1, length_map(map)),
  theta_K   = 0,
  ...
) {
  n <- mesh$n; nrep <- length(unique(replicate))
  stopifnot("The index should be integers." = all(map == round(map)))

  h <- c(diff(mesh$loc), 1)
  G <- Matrix::Diagonal(n);
  C <- Matrix::sparseMatrix(j=1:(n-1), i=2:n, x=-1, dims=c(n,n))
  stopifnot("The mesh should be 1d and has gap 1." = all(h == 1))

  G <- Matrix::kronecker(diag(nrep), G)
  C <- Matrix::kronecker(diag(nrep), C)

  stopifnot("The length of theta_K should be 1." = length(theta_K) == 1)
  alpha <- ar1_th2a(theta_K)
  ngme_operator(
    model = "ar1",
    theta_K = theta_K,
    C = ngme_as_sparse(C),
    G = ngme_as_sparse(G),
    K = alpha * C + G,
    h = h,
    A = INLA::inla.spde.make.A(mesh = mesh, loc = map),
    symmetric = FALSE
  )
}

bv <- function(
  first,
  second,
  theta=1, rho=1,
  share_parameter=FALSE
) {
  stopifnot(
    inherits(first, "ngme_operator"),
    inherits(second, "ngme_operator"),
    all(dim(first$K) == dim(second$K))
  )

  D <- build_D(theta, rho)
  bigD <- kronecker(D, Matrix::Diagonal(nrow(first$K))); bigD

  ngme_operator(
    model = "bv",
    K = bigD %*% Matrix::bdiag(first$K, second$K),
    h = c(first$h, second$h),
    A = Matrix::bdiag(first$A, second$A),
    first = first,
    second = second
  )
}

# build_noise <- function(
#   noise,
#   operator
# ) {
#   stopifnot(inherits(operator, "ngme_operator"))
# }