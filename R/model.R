# function for specify ngme.model basic structure
ngme.model <- function(
  model,
  W_size      = NULL,
  theta_K     = NULL,
  fix_theta_K = FALSE,
  W           = NULL,
  fix_W       = FALSE,
  A           = NULL,
  A_pred      = NULL,
  noise_type  = NULL,
  noise       = NULL,
  control     = ngme.control.f(),
  V_size      = NULL,
  ...
) {
  stopifnot(is.character(model))

  structure(
    list(
      model         = model,
      W_size        = W_size,
      theta_K       = theta_K,
      n_theta_K     = length(theta_K),
      A             = A,
      A_pred        = A_pred,
      noise_type    = noise_type,
      noise         = noise,
      W             = W,
      fix_W         = fix_W,
      fix_theta_K   = fix_theta_K,
      V_size        = V_size,
      control       = control,
      ...
    ),
    class = "ngme_model"
  )
}

#' Print ngme model
#'
#' @param model
#'
#' @return a list (model specifications)
#' @expo`rt
print.ngme_model <- function(model, padding=0) {
  pad_space <- paste(rep(" ", padding), collapse = "")
  pad_add4_space <- paste(rep(" ", padding + 4), collapse = "")

  cat(pad_space); cat("Ngme model: "); cat(model$model); cat("\n")

  cat(pad_space); cat("Model parameters: \n")
  params <- with(model, {
    switch(model,
      "ar1"     = paste0(pad_add4_space, "alpha = ",    ngme.format(theta_K)),
      "matern"  = paste0(pad_add4_space, "theta_K = ",  ngme.format(theta_K)),
      "matern2D" = paste0(pad_add4_space, "theta_K = ",  ngme.format(theta_K)),
      "rw1"     = paste0(pad_add4_space, "No parameter needed."),
      "unkown"  = paste0(pad_add4_space, "theta_K = ",  ngme.format(theta_K)),
    )
  })
  cat(params);
  cat("\n\n")

  print.ngme_noise(model$noise, padding = padding)

  invisible(model)
}


#' ngme ar1 model specification
#'
#' Generating C, G and A given index and replicates
#'
#' @param index index for the process
#' @param replicates replicates for the process
#' @param alpha initial value for alpha
#' @param range range for the mesh
#'
#' @return a list of specification of model
#' @export
#'
#' @examples
ngme.ar1 <- function(
  index,
  replicates = NULL,
  alpha = 0.5,
  range = c(1, max(index)),

  index_pred = NULL,
  use_num_dK = FALSE,
  ...
) {
  # overwirte the default
  # watch out! avoid situation like ngme.ar1(alpha = 0.3, ...$theta_K=0.4)
  if (!is.null(list(...)$theta_K)) alpha <- list(...)$theta_K
  if (is.null(replicates)) replicates <- rep(1, length(index))

  unique_rep <- unique(replicates)
  nrep <- length(unique_rep)

  # e.g. index = c(1:200, 1:100)
  #      replicates = c(rep(1, 200), rep(2, 100))
  #      n =200 in this case

  n <- range[2] - range[1] + 1

  # construct G
    G <- Matrix::Matrix(diag(n));

  # construct C
    C <- Matrix::Matrix(0, n, n)
    C[seq(2, n*n, by = n+1)] <- -1

    if(!is.null(nrep)) {
      C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), C)
      G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), G)
    }

  model <- ngme.model(
    model       = "ar1",
    theta_K     = alpha,
    W_size      = n,
    V_size      = n,
    A           = ngme.ts.make.A(loc = index, replicates = replicates, range = range),
    A_pred      = ngme.ts.make.A(index_pred, replicates = replicates, range = range),
    h           = rep(1.0, n),
    C           = ngme.as.sparse(C),
    G           = ngme.as.sparse(G),
    ...
  )

  model
}

#' Create a Matern SPDE model
#'
#' @param alpha
#' @param mesh mesh argument
#' @param fem.mesh.matrices specify the FEM matrices
#' @param d indicating the dimension of mesh (together with fem.mesh.matrices)
#' @param kappa # parameterization for k^2 C + G, only for stationary
#' @param theta_kappa
#' @param B_kappa bases for kappa
#'
#' @return a list (n, C (diagonal), G, B.kappa) for constructing operator
#' @export
#'
#' @examples
ngme.matern <- function(
  index = NULL,
  alpha = 2,
  kappa = 1,
  theta_kappa = NULL,
  mesh = NULL,
  replicates = NULL,
  fem.mesh.matrices = NULL,
  d = NULL,
  A = NULL,        # watch out! Can also specify in f function, not used for now
  B_kappa = NULL,
  ...
) {
  if (is.null(mesh) && is.null(fem.mesh.matrices))
    stop("At least specify mesh or matrices")

  if (alpha - round(alpha) != 0) {
    stop("alpha should be integer, now only 2 or 4")
  }

  stopifnot(alpha == 2 || alpha == 4)

  # use kappa as parameterization
  stopifnot("Don't overwrite kappa with NULL." = !is.null(kappa))
  stopifnot("kappa is only for stationary case." = length(kappa) == 1)
  stopifnot("kappa is greater than 0." = kappa > 0)
  if (is.null(theta_kappa)) theta_kappa <- log(kappa)

  if (is.null(B_kappa))
    B_kappa <- matrix(1, nrow = mesh$n, ncol = length(theta_kappa))

  # supply mesh
  if (!is.null(mesh)) {
    d <- get_inla_mesh_dimension(mesh)
    if (d == 1) {
      fem <- INLA::inla.mesh.1d.fem(mesh)
      C <- fem$c1
      G <- fem$g1
    } else {
      fem <- INLA::inla.mesh.fem(mesh, order = alpha)
      C <- fem$c0  # diag
      G <- fem$g1
    }
  } else {
    C <- fem.mesh.matrices$C
    G <- fem.mesh.matrices$G
  }
  # h <- diag(C)
  h <- rep(1, mesh$n)

  if (!is.null(A)) {
    nrep <- ncol(A) / nrow(C)
    C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), C)
    G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), G)
    h <- rep(h, times = nrep)
  }

  model <- ngme.model(
    model       = "matern",
    A           = A,
    W_size      = mesh$n,
    V_size      = nrow(C),
    theta_K     = theta_kappa,
    alpha       = alpha,
    B_kappa     = B_kappa,
    C           = ngme.as.sparse(C),
    G           = ngme.as.sparse(G),
    h           = h,
    ...
  )
  model
}

#' Create a Matern SPDE 2D model
#'
#' @param alpha alpha parameter for the first field
#' @param alpha2 alpha parameter for the second field
#' @param mesh mesh argument
#' @param fem.mesh.matrices specify the FEM matrices
#' @param d indicating the dimension of mesh (together with fem.mesh.matrices)
#' @param theta_kappa
#' @param B_kappa bases for kappa
#'
#' @return a list (n, C (diagonal), G, B.kappa) for constructing operator
#' @export
#'
#' @examples
ngme.matern2D <- function(
  index = NULL,
  alpha = c(2, 2),
  kappa = c(1, 1),
  theta_kappa = NULL,
  mesh = NULL,
  replicates = NULL,
  fem.mesh.matrices = NULL,
  d = NULL,
  A = NULL,        # watch out! Can also specify in f function, not used for now
  B_kappa = NULL,
  ...
) {
  if (is.null(mesh) && is.null(fem.mesh.matrices))
    stop("At least specify mesh or matrices")

  if (alpha[1] - round(alpha[1]) != 0 | alpha[2] - round(alpha[2]) != 0 ) {
  stop("alpha should be integer, now only 2 or 4")
}
  if (alpha[1] != alpha[2] ) {
  stop("Alpha for each field is not implemented, please provide common alpha for both fields")
}

stopifnot((alpha[1] == 2 || alpha[1] == 4) & (alpha[2] == 2 || alpha[2] == 4) )

# use kappa as parameterization
  stopifnot("Don't overwrite kappa with NULL." = !is.null(kappa))
  stopifnot("kappa is only for stationary case." = length(kappa) == 2)
  stopifnot("kappa is greater than 0." = (kappa[1] > 0 & kappa[2] > 0))
  if (is.null(theta_kappa)) theta_kappa <- log(kappa)

  if (is.null(B_kappa))
    B_kappa <- matrix(1, nrow = mesh$n, ncol = length(theta_kappa))

  # supply common mesh
  if (!is.null(mesh)) {
    d <- get_inla_mesh_dimension(mesh)
    if (d == 1) {
      fem <- INLA::inla.mesh.1d.fem(mesh)
      C <- fem$c1
      G <- fem$g1 
    } else {
      fem <- INLA::inla.mesh.fem(mesh, order = alpha[1]) #assume both process have common alpha
      C <- fem$c0  # diag
      G <- fem$g1
      #TODO if created here, ensure block structure
    }
  } else {
    #TODO if supplied, already must have block structure
    C <- fem.mesh.matrices$C
    G <- fem.mesh.matrices$G
  }
  # h <- diag(C)
  h <- rep(1, mesh$n)
#TODO change mesh here for bivariate case with replicates
  if (!is.null(A)) { #dim(A) = nrep*nobs X nrep*nmesh - A(block diagonal if with replicates)
    nrep <- ncol(A) / nrow(C)
    C <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), C)
    G <- Matrix::kronecker(Matrix::Diagonal(nrep, 1), G)
    h <- rep(h, times = nrep)
  }


  model <- ngme.model(
    model       = "matern2D",
    A           = A, #constructed outside of ngme2
    W_size      = mesh$n,
    V_size      = nrow(C),
    theta_K     = theta_kappa,
    alpha       = c(alpha[1], alpha[2]),
    B_kappa     = B_kappa,
    C           = ngme.as.sparse(C),
    G           = ngme.as.sparse(G),
    h           = h,
    ...
  )
  model
}


# rw1, rw2
# nodes = 100 (inla.group)

# ?inla.spde.make.A
# inla.spde.make.A(
#   index
#  replicates=)

#' ngme model - random walk of order 1
#'
#' Generating C, G and A given index and replicates
#' size of C and G is (n-1) * n, size of V is n-1
#'
#' @param index index for the process
#' @param replicates replicates for the process
#' @param mesh inla.1d.mesh
#' @param n_points or num of points, evenly spaced mesh
#' @return a list
#' @export
#'
#' @examples
ngme.rw1 <- function(
  index,
  replicates = NULL,
  circular = FALSE,
  mesh = NULL,
  n_points = NULL,
  # extra A matrix
  ...
) {
  # create mesh using index
  sorted_index <- sort(index, index.return = TRUE)
  h <- diff(sorted_index$x)
  # permutation matrix, same as A <- diag(length(index))[sorted_index$ix, ]
  A <- Matrix::sparseMatrix(seq_along(index), sorted_index$ix, x = 1)

  n <- length(index) - 1
  # construct G
    G <- Matrix::Matrix(diag(n));
    G <- cbind(G, rep(0, n))

  # construct C
    C <- Matrix::Matrix(0, n, n)
    C[seq(n+1, n*n, by = n+1)] <- -1
    C <- cbind(C, c(rep(0, n-1), -1))

  model <- ngme.model(
    model       = "rw1",
    theta_K     = 1,
    fix_theta_K = TRUE,
    W_size      = n + 1,
    V_size      = n,
    A           = A,
    # A_pred      = ngme.ts.make.A(index_pred, replicates = replicates, range = range),
    h           = h,
    C           = ngme.as.sparse(C),
    G           = ngme.as.sparse(G),
    ...
  )
  model
}
