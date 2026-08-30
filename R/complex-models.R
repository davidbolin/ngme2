coerce_sub_model <- function(sub_model, mesh) {
  if (inherits(sub_model, "ngme_operator")) {
    return(sub_model)
  }

  if (inherits(sub_model, "ngme_operator_def")) {
    model_name <- sub_model$model
    args <- sub_model$args
    args$mesh <- mesh
    return(do.call(model_name, args))
  }

  stop("Each sub-model must be an ngme_operator or ngme_operator_def.")
}

#' Ngme bivariate model specification
#'
#' Giving 2 sub_models, build a correlated bivariate operator based on K = D(theta, rho) %*% diag(c1*K_1, c2*K_2)
#' \deqn{D(\theta, \rho) = \begin{pmatrix}
#'   \cos(\theta) + \rho \sin(\theta) & -\sin(\theta) \sqrt{1+\rho^2} \\
#'   \sin(\theta) - \rho \cos(\theta) & \cos(\theta) \sqrt{1+\rho^2}
#' \end{pmatrix}}
#'
#' @param mesh the mesh where model is defined
#' @param sub_models a list of sub_models (total 2 sub_models)
#' @param theta the rotation parameter (starting point should be from -pi/4 to pi/4)
#' @param rho the parameter related to correlation
#' @param c1 scaling factor for the first sub-model (default: 1)
#' @param c2 scaling factor for the second sub-model (default: 1)
#' @param share_param TRUE if share the same parameter for 2 sub_models (of same type)
#' @param fix_theta TRUE if fix the theta of bv model
#' @param use_c_param TRUE to estimate c1 and c2 as parameters, FALSE to fix them at specified values (default: FALSE)
#'
#' @details
#' The bivariate model combines two sub-models with a rotation matrix D and optional scaling factors c1 and c2.
#' When \code{use_c_param = TRUE}, c1 and c2 are estimated as parameters (on log scale).
#' When \code{use_c_param = FALSE} (default), c1 and c2 are fixed at their specified values (both default to 1).
#'
#' @return a ngme_operator of bivariate model
#' @export
bv <- function(
    mesh,
    sub_models,
    theta = 0,
    rho = 0,
    c1 = 1,
    c2 = 1,
    share_param = FALSE,
    fix_theta = FALSE,
    use_c_param = FALSE) {
  if (is.list(mesh) && !inherits(mesh, c("inla.mesh.1d", "inla.mesh", "fm_mesh_1d", "fm_mesh_2d", "metric_graph"))) {
    return(structure(list(model = "bv", args = as.list(environment())), class = "ngme_operator_def"))
  }
  mesh <- ngme_build_mesh(mesh)
  model_names <- sort(names(sub_models))

  stopifnot(
    "Please provide names for sub_models" = !is.null(model_names),
    "Length of sub_models should be 2" = length(sub_models) == 2
  )

  # Inject mesh if sub_models are def
  first <- coerce_sub_model(sub_models[[model_names[1]]], mesh)
  second <- coerce_sub_model(sub_models[[model_names[2]]], mesh)

  stopifnot(
    "the theta should be in (-pi/2, pi/2)" = theta >= -pi / 2 & theta <= pi / 2
  )
  eta <- tan(theta)

  theta_K <- c()
  if (!fix_theta) theta_K <- c(theta_K, eta)
  theta_K <- c(theta_K, rho)
  if (use_c_param) theta_K <- c(theta_K, log(c1), log(c2))
  theta_K <- c(theta_K, first$theta_K, second$theta_K)

  # pass the theta_K to first and second
  n_bv <- (if (!fix_theta) 1 else 0) + 1 + (if (use_c_param) 2 else 0)
  first$theta_K <- theta_K[(n_bv + 1):(n_bv + first$n_theta_K)]
  second$theta_K <- theta_K[(n_bv + first$n_theta_K + 1):length(theta_K)]

  D <- build_D(theta, rho)
  bigD <- kronecker(D, Matrix::Diagonal(nrow(first$K)))

  # the range of theta is (-pi/4, pi/4), we expand it to (-pi/2, pi/2)
  # so that the maximum can be reached in either direction
  # in the end, we will NOT transform it back to (-pi/4, pi/4)
  eta_to_theta <- function(eta) atan(eta)

  update_K <- function(theta_K) {
    idx <- 1
    if (!fix_theta) {
      eta <- theta_K[idx]
      theta <- atan(eta)
      idx <- idx + 1
    } else {
      theta <- 0 # fixed theta
    }
    rho <- theta_K[idx]
    idx <- idx + 1

    c1 <- 1
    c2 <- 1
    if (use_c_param) {
      c1 <- exp(theta_K[idx])
      c2 <- exp(theta_K[idx + 1])
      idx <- idx + 2
    }

    D <- build_D(theta, rho)
    bigD <- kronecker(D, Matrix::Diagonal(nrow(first$K)))

    theta_K1 <- theta_K[idx:(idx + first$n_theta_K - 1)]
    idx <- idx + first$n_theta_K
    theta_K2 <- theta_K[idx:length(theta_K)]

    first$K <- first$update_K(theta_K1)
    second$K <- second$update_K(theta_K2)
    bigD %*% Matrix::bdiag(c1 * first$K, c2 * second$K)
  }

  ngme_operator(
    mesh = mesh,
    model = "bv",
    first = first,
    second = second,
    theta_K = theta_K,
    update_K = update_K,
    K = update_K(theta_K),
    h = c(first$h, second$h),
    symmetric = FALSE,
    zero_trace = FALSE,
    model_names = model_names,
    share_param = share_param,
    fix_theta = fix_theta,
    bv_theta = theta,
    use_c_param = use_c_param,
    param_name = c(
      if (!fix_theta) "theta" else NULL,
      "rho",
      if (use_c_param) c("c1", "c2") else NULL,
      if (!is.null(first$param_name)) paste(first$param_name, "(1st)") else NULL,
      if (!is.null(second$param_name)) paste(second$param_name, "(2nd)") else NULL
    ),
    param_trans = c(
      if (!fix_theta) list(eta_to_theta) else NULL,
      list(identity),
      if (use_c_param) list(exp, exp) else NULL,
      first$param_trans,
      second$param_trans
    )
  )
}


#' ngme tensor-product model specification
#'
#' Given 2 operator (first and second), build a tensor-product operator based on K = K_first x K_second (here x is Kronecker product)
#'
#' @param first ngme_operator, left side of kronecker model (usually a temporal or iid model)
#' @param second ngme_operator, right side of kronecker model (ususally a temporal or spatial model)
#'
#' @return a list of specification of model
#' @export
tp <- function(
    first,
    second) {
  stopifnot(
    "Please provide ngme_operator for first" = inherits(first, "ngme_operator"),
    "Please provide ngme_operator for second" = inherits(second, "ngme_operator")
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


#' Ngme bivariate model with Matern sub_models
#'
#' Giving 2 sub_models, build a correlated bivaraite operator based on K = D(theta, eta) %*% diag(K_1, K_2)
#' \deqn{D(\theta, \rho) = \begin{pmatrix}
#'   \cos(\theta) + \rho \sin(\theta) & -\sin(\theta) \sqrt{1+\rho^2} \\
#'   \sin(\theta) - \rho \cos(\theta) & \cos(\theta) \sqrt{1+\rho^2}
#' \end{pmatrix}}
#'
#' @param mesh the mesh where model is defined
#' @param sub_models a list of sub_models (should be two matern models)
#' @param theta the parameter related to rotation
#' @param rho the parameter related to correlation
#' @param sd1 scaling parameter for the first sub-model
#' @param sd2 scaling parameter for the second sub-model
#' @param group group vector, can be inherited from ngme() function
#' @param share_param TRUE if share the same parameter for 2 sub_models (of same type)
#' @param fix_theta TRUE if the rotation parameter is fixed
#' @param ... ignore
#'
#' @return a list of specification of model
#' @export
bv_matern <- function(
    mesh,
    sub_models,
    theta = 0,
    rho = 0,
    sd1 = 1,
    sd2 = 1,
    group = NULL,
    share_param = FALSE,
    fix_theta = FALSE,
    ...) {
  if (is.list(mesh) && !inherits(mesh, c("inla.mesh.1d", "inla.mesh", "fm_mesh_1d", "fm_mesh_2d", "metric_graph"))) {
    return(structure(list(model = "bv_matern", args = as.list(environment())), class = "ngme_operator_def"))
  }
  mesh <- ngme_build_mesh(mesh)
  d <- switch(mesh$manifold,
    "R1" = 1,
    "R2" = 2,
    stop("manifold not supported")
  )

  # sort the sub_models
  model_names <- sort(names(sub_models))

  # read group argument from parent frame
  if (is.null(group) &&
    exists("group", parent.frame())) {
    group <- get("group", parent.frame())
  }

  stopifnot(
    "theta should be within (-pi/2, pi/2)" = theta > -pi / 2 & theta < pi / 2,
    "Please provide names for sub_models" = !is.null(model_names),
    "Please provide group argument to indicate different fields" = !is.null(group),
    "Length of sub_models should be 2" = length(sub_models) == 2,
    "Name of sub_models should be in group" = all(model_names %in% levels(as.factor(group))),
    "sd1 and sd2 should be positive" = sd1 > 0 & sd2 > 0
  )
  group <- factor(group, levels = model_names)

  # Inject mesh if sub_models are def
  first <- coerce_sub_model(sub_models[[model_names[1]]], mesh)
  second <- coerce_sub_model(sub_models[[model_names[2]]], mesh)

  eta <- tan(theta)
  if (!fix_theta) {
    theta_K <- c(
      eta, rho, log(sd1), log(sd2),
      first$theta_K, second$theta_K
    )
  } else {
    theta_K <- c(
      rho, log(sd1), log(sd2),
      first$theta_K, second$theta_K
    )
  }

  update_K <- function(theta_K) {
    if (!fix_theta) {
      theta <- atan(theta_K[1])
      rho <- theta_K[2]
      sd1 <- exp(theta_K[3])
      sd2 <- exp(theta_K[4])
      theta_K1 <- theta_K[5:(4 + first$n_theta_K)]
      theta_K2 <- theta_K[(5 + first$n_theta_K):length(theta_K)]
    } else {
      rho <- theta_K[1]
      sd1 <- exp(theta_K[2])
      sd2 <- exp(theta_K[3])
      theta_K1 <- theta_K[4:(3 + first$n_theta_K)]
      theta_K2 <- theta_K[(4 + first$n_theta_K):length(theta_K)]
    }

    alpha1 <- first$alpha
    alpha2 <- second$alpha
    kappa1 <- exp(first$theta_K)
    kappa2 <- exp(second$theta_K)

    nu1 <- alpha1 - d / 2
    nu2 <- alpha2 - d / 2
    c1 <- sqrt(gamma(nu1) / (gamma(alpha1))) / (kappa1^(nu1) * (4 * pi)^(d / 4) * sd1)
    c2 <- sqrt(gamma(nu2) / (gamma(alpha2))) / (kappa2^(nu2) * (4 * pi)^(d / 4) * sd2)

    D <- build_D(theta, rho)
    bigD <- kronecker(D, Matrix::Diagonal(nrow(first$K)))
    first$K <- first$update_K(theta_K1)
    second$K <- second$update_K(theta_K2)
    bigD %*% Matrix::bdiag(c1 * first$K, c2 * second$K)
  }

  eta_to_theta <- function(eta) atan(eta)
  ngme_operator(
    mesh = mesh,
    model = "bv_matern",
    first = first,
    second = second,
    theta_K = theta_K,
    update_K = update_K,
    K = update_K(theta_K),
    h = c(first$h, second$h),
    dim = d,
    symmetric = FALSE,
    zero_trace = FALSE,
    model_names = model_names,
    share_param = share_param,
    fix_theta = fix_theta,
    bv_theta = theta,
    param_name = c(
      if (!fix_theta) "theta" else NULL, "rho", "sd1", "sd2",
      if (!is.null(first$param_name)) paste(first$param_name, "(1st)") else NULL,
      if (!is.null(second$param_name)) paste(second$param_name, "(2nd)") else NULL
    ),
    param_trans = c(if (!fix_theta) eta_to_theta else NULL, identity, exp, exp, first$param_trans, second$param_trans)
  )
}





#' Ngme space-time non-separable model specification
#'
#' Implements a non-separable space-time model based on the advection-diffusion SPDE with first-order derivative in time.
#' The model combines temporal and spatial components through a finite difference method (implicit Euler) for temporal discretization
#' and finite element method (continuous Galerkin) for spatial discretization. When advection dominates diffusion,
#' the "Streamline Diffusion" stabilization technique is applied.
#'
#' The model is particularly useful for large space-time datasets in environmental science, offering computationally
#' efficient methods for parameter estimation, kriging prediction, and conditional simulations.
#'
#' For details, see \doi{10.1016/j.spasta.2024.100847}
#'
#' @param mesh A list of two objects:
#'   \itemize{
#'     \item mesh_t - The temporal mesh
#'     \item mesh_s - The spatial mesh
#'   }
#' @param alpha 2 or 4, SPDE smoothness parameter.
#' @param cc Parameter c in the SPDE.
#' @param kappa Kappa parameter from Matern SPDE.
#' @param lambda The spatial damping parameter.
#' @param fix_gamma TRUE if fix gamma (advection term), FALSE if estimate gamma.
#' @param theta_gamma_x The x component of the advection term: \code{gamma_x = B_gamma_x \%*\% theta_gamma_x}.
#' @param theta_gamma_y The y component of the advection term: \code{gamma_y = B_gamma_y \%*\% theta_gamma_y}.
#' @param shared_theta_gamma TRUE if share the same theta_gamma for all time nodes. (theta_gamma_x and theta_gamma_y will be the same)
#' @param B_gamma_x The design matrix for the x component of the advection term.
#' @param B_gamma_y The design matrix for the y component of the advection term.
#' @param B_gamma_x_list A list of design matrices for the x component of the advection term on every time node, length(B_gamma_x_list) == nt-1.
#' @param B_gamma_y_list A list of design matrices for the y component of the advection term on every time node, length(B_gamma_y_list) == nt-1.
#' @param stabilization TRUE if using a stabilization term (for implicit Euler).
#'
#' @return ngme_operator object.
#' @export
spacetime <- function(
    mesh,
    lambda = 1, # fixed
    alpha = 2, # alpha = 2, 4, fixed
    cc = 1,
    kappa = 1,
    # advection term
    fix_gamma = FALSE,
    theta_gamma_x = 0,
    theta_gamma_y = 0,
    shared_theta_gamma = FALSE,
    # theta_gamma_t = 0,
    B_gamma_x = matrix(1, nrow = mesh[[2]]$n, ncol = 1),
    B_gamma_y = matrix(1, nrow = mesh[[2]]$n, ncol = 1),
    # or, a list of B_gamma_x, B_gamma_y on every time node.
    B_gamma_x_list = NULL,
    B_gamma_y_list = NULL,
    # B_gamma_t = matrix(1, nrow = mesh[[1]]$n, ncol = 1),
    stabilization = TRUE) {
  if (is.list(mesh) && length(mesh) > 0 && is.list(mesh[[1]]) && !inherits(mesh[[1]], c("inla.mesh.1d", "inla.mesh", "fm_mesh_1d", "fm_mesh_2d", "metric_graph"))) {
    return(structure(list(model = "spacetime", args = as.list(environment())), class = "ngme_operator_def"))
  }
  method <- "euler" # for now only support implicit euler

  if (all(theta_gamma_x == 0) && all(theta_gamma_y == 0) && fix_gamma) {
    stabilization <- FALSE
  }

  mesh_t <- mesh[[1]]
  nt <- mesh_t$n
  mesh_s <- mesh[[2]]
  ns <- mesh_s$n
  stopifnot(
    "Please provide mesh as a list of length 2" = length(mesh) == 2,
    "alpha should be 2 or 4" = alpha == 2 || alpha == 4,
    "method should be galerkin or euler" = method %in% c("galerkin", "euler"),
    "First mesh should be 1d" = fmesher::fm_manifold_dim(mesh_t) == 1,
    "Second mesh should be 2d" =
      fmesher::fm_manifold_dim(mesh_s) == 2,
    "require package rSPDE" = requireNamespace("rSPDE", quietly = TRUE)
  )

  # Additional validation for shared_theta_gamma
  if (shared_theta_gamma && !fix_gamma) {
    if (!all(theta_gamma_x == theta_gamma_y)) {
      warning("When shared_theta_gamma = TRUE, theta_gamma_x and theta_gamma_y should be the same. Using theta_gamma_x.")
    }
    # Check that B_gamma_x and B_gamma_y have the same number of columns when shared
    if (ncol(B_gamma_x) != ncol(B_gamma_y)) {
      stop("When shared_theta_gamma = TRUE, B_gamma_x and B_gamma_y must have the same number of columns.")
    }
  }

  if (!is.null(B_gamma_x_list) && !is.null(B_gamma_y_list)) {
    stopifnot(
      "B_gamma_x_list and B_gamma_y_list should have the same length" =
        length(B_gamma_x_list) == length(B_gamma_y_list),
      "The number of time nodes should be the same as the length of B_gamma_x_list and B_gamma_y_list" =
        nt - 1 == length(B_gamma_x_list),
      "The ncol of items in B_gamma_x_list should be the same as the length of theta_gamma_x" =
        all(sapply(B_gamma_x_list, ncol) == length(theta_gamma_x)),
      "The ncol of items in B_gamma_y_list should be the same as the length of theta_gamma_y" =
        all(sapply(B_gamma_y_list, ncol) == length(theta_gamma_y))
    )

    # Additional validation for shared_theta_gamma with lists
    if (shared_theta_gamma && !fix_gamma) {
      # Check that all matrices in B_gamma_x_list and B_gamma_y_list have the same number of columns
      if (!all(sapply(B_gamma_x_list, ncol) == sapply(B_gamma_y_list, ncol))) {
        stop("When shared_theta_gamma = TRUE, all matrices in B_gamma_x_list and B_gamma_y_list must have the same number of columns.")
      }
    }
  } else {
    # Use B_gamma_x and B_gamma_y
    stopifnot(
      "B_gamma_x and B_gamma_y should have the same number of rows as the spatial mesh" =
        nrow(B_gamma_x) == ns && nrow(B_gamma_y) == ns,
      "theta_gamma_x should be of length ncol(B_gamma_x)" =
        length(theta_gamma_x) == ncol(B_gamma_x),
      "theta_gamma_y should be of length ncol(B_gamma_y)" =
        length(theta_gamma_y) == ncol(B_gamma_y)
    )

    # turn it into a list
    B_gamma_x_list <- lapply(1:(nt - 1), function(i) B_gamma_x)
    B_gamma_y_list <- lapply(1:(nt - 1), function(i) B_gamma_y)
  }

  ##### temporal FEM #####
  fem_t <- fmesher::fm_fem(mesh_t, order = 2)
  Ct <- fem_t$c0
  Gt <- fem_t$g1

  # Build Bt
  Bt <- Matrix::bandSparse(
    n = nt, m = nt, k = c(-1, 0, 1),
    diagonals = cbind(rep(0.5, nt), rep(0, nt), rep(-0.5, nt))
  )
  Bt[1, 1] <- -0.5
  Bt[nt, nt] <- 0.5

  # Bt
  # [1,] -0.5 -0.5  .    .    .
  # [2,]  0.5  .   -0.5  .    .
  # [3,]  .    0.5  .   -0.5  .
  # [4,]  .    .    0.5  .   -0.5
  # [5,]  .    .    .    0.5  0.5


  ##### space FEM #####
  fem_s <- fmesher::fm_fem(mesh_s, order = alpha)
  Cs <- fem_s$c0
  Gs <- fem_s$g1
  cs_diag <- Matrix::diag(Cs)

  # compute h
  if (method == "galerkin") {
    h <- Matrix::diag(Ct) %x% Matrix::diag(Cs)
  } else {
    # implicit euler
    dt <- c(1, diff(mesh_t$loc))
    h <- dt %x% Matrix::diag(Cs) / cc
  }

  FV <- mesh_s$graph$tv
  P <- sf::st_coordinates(fmesher::fm_vertices(mesh_s))[, 1:2]
  fem2d <- rSPDE::rSPDE.fem2d(FV = FV, P = P)
  Bx <- fem2d$Bx
  By <- fem2d$By
  Hxx <- fem2d$Hxx
  Hyy <- fem2d$Hyy
  Hxy <- fem2d$Hxy
  Hyx <- fem2d$Hyx

  theta_K <- if (fix_gamma) {
    c(log(cc), log(kappa))
  } else if (shared_theta_gamma) {
    c(log(cc), log(kappa), theta_gamma_x)
  } else {
    c(log(cc), log(kappa), theta_gamma_x, theta_gamma_y)
  }
  n_theta_gamma_x <- if (fix_gamma) 0 else ncol(B_gamma_x)
  n_theta_gamma_y <- if (fix_gamma) 0 else if (shared_theta_gamma) 0 else ncol(B_gamma_y)

  # gamma_x <- B_gamma_x %*% theta_gamma_x
  # gamma_y <- B_gamma_y %*% theta_gamma_y

  if (shared_theta_gamma) {
    # Use theta_gamma_x for both x and y components
    gamma_x_list <- lapply(1:(nt - 1), function(i) B_gamma_x_list[[i]] %*% theta_gamma_x)
    gamma_y_list <- lapply(1:(nt - 1), function(i) B_gamma_y_list[[i]] %*% theta_gamma_x)
  } else {
    gamma_x_list <- lapply(1:(nt - 1), function(i) B_gamma_x_list[[i]] %*% theta_gamma_x)
    gamma_y_list <- lapply(1:(nt - 1), function(i) B_gamma_y_list[[i]] %*% theta_gamma_y)
  }
  # gamma_t <- B_gamma_t %*% theta_gamma_t
  build_Bs <- function(gamma_x, gamma_y) {
    Dx <- Matrix::Diagonal(x = gamma_x)
    Dy <- Matrix::Diagonal(x = gamma_y)
    Dx %*% Bx %*% Dx + Dy %*% By %*% Dy
  }

  build_Bs_list <- function(gamma_x_list, gamma_y_list) {
    lapply(1:(nt - 1), function(i) build_Bs(gamma_x_list[[i]], gamma_y_list[[i]]))
  }

  build_S <- function(gamma_x, gamma_y) {
    # if gamma_x and gamma_y are 0, then S = 0
    if (sum(gamma_x^2 + gamma_y^2) < 1e-8) {
      return(Matrix::Diagonal(n = ns, x = 0))
    }
    gamma_xx <- gamma_x^2
    gamma_yy <- gamma_y^2
    gamma_xy <- gamma_x * gamma_y
    Dxx <- Matrix::Diagonal(x = gamma_xx)
    Dyy <- Matrix::Diagonal(x = gamma_yy)
    Dxy <- Matrix::Diagonal(x = gamma_xy)
    tmp <- Dxx %*% Hxx %*% Dxx + Dyy %*% Hyy %*% Dyy + Dxy %*% (Hxy + Hyx) %*% Dxy
    gamma_norm <- sqrt(sum(gamma_x^2 + gamma_y^2))
    Matrix::Diagonal(x = cs_diag) %*% tmp / gamma_norm
  }

  build_S_list <- function(gamma_x_list, gamma_y_list) {
    lapply(1:(nt - 1), function(i) build_S(gamma_x_list[[i]], gamma_y_list[[i]]))
  }

  # compute L_s
  update_K <- function(theta_K) {
    cc <- exp(theta_K[1])
    kappa <- exp(theta_K[2])

    if (shared_theta_gamma) {
      theta_gamma_x <- if (length(theta_K) > 2) theta_K[3:(2 + n_theta_gamma_x)] else double(0)
      theta_gamma_y <- theta_gamma_x # Use same parameter for both x and y
    } else {
      theta_gamma_x <- if (length(theta_K) > 2) theta_K[3:(2 + n_theta_gamma_x)] else double(0)
      theta_gamma_y <- if (length(theta_K) > 2) theta_K[(3 + n_theta_gamma_x):length(theta_K)] else double(0)
    }

    # gamma_x <- as.vector(B_gamma_x %*% theta_gamma_x)
    # gamma_y <- as.vector(B_gamma_y %*% theta_gamma_y)
    gamma_x_list <- lapply(1:(nt - 1), function(i) as.vector(B_gamma_x_list[[i]] %*% theta_gamma_x))
    gamma_y_list <- lapply(1:(nt - 1), function(i) as.vector(B_gamma_y_list[[i]] %*% theta_gamma_y))

    if (!fix_gamma) {
      # Bs = build_Bs(gamma_x, gamma_y)
      # S = build_S(gamma_x, gamma_y)
      Bs_list <- build_Bs_list(gamma_x_list, gamma_y_list)
      S_list <- build_S_list(gamma_x_list, gamma_y_list)
    }

    # compute L_s
    # Turn L into a list
    # L = kappa^2 * Cs + lambda * Gs + Bs
    L_list <- lapply(1:(nt - 1), function(i) kappa^2 * Cs + lambda * Gs + Bs_list[[i]])

    Cs_inv <- Matrix::Diagonal(n = ns, x = 1 / Matrix::diag(Cs))

    # watch-out! make sure which t
    if (alpha == 4) L <- L %*% Cs_inv %*% Matrix::t(L)

    if (method == "galerkin") {
      # not used
      K <- Bt %x% Cs + 1 / cc * Ct %x% L
    } else if (method == "euler") {
      # implicit euler
      null_matrix <- Matrix::Diagonal(n = ns, x = 0)

      if (stabilization) {
        # L = L + S
        L_list <- lapply(1:(nt - 1), function(i) L_list[[i]] + S_list[[i]])
      }

      # diagonalize the list of L
      # diag_L <- Matrix::bdiag(
      #   lapply(1:(nt-1), function(i) L)
      # )

      diag_L_list <- Matrix::bdiag(
        lapply(1:(nt - 1), function(i) L_list[[i]])
      )

      # K <- rw1(1:nt)$K %x% Cs + 1/cc *
      #   Matrix::bdiag(null_matrix, diag_L)

      K <- rw1(1:nt)$K %x% Cs + 1 / cc *
        Matrix::bdiag(null_matrix, diag_L_list)

      K <- sqrt(cc) * K
    }
    return(K)
  }

  Bs_list <- build_Bs_list(gamma_x_list, gamma_y_list)
  S_list <- build_S_list(gamma_x_list, gamma_y_list)
  K <- update_K(theta_K)

  BtCs <- if (method == "galerkin") {
    ngme_as_sparse(Bt %x% Cs)
  } else {
    ngme_as_sparse(rw1(1:nt)$K %x% Cs)
  }

  ngme_operator(
    model = "spacetime",
    mesh = mesh,
    alpha = alpha,
    theta_gamma_x = theta_gamma_x,
    theta_gamma_y = theta_gamma_y,
    fix_gamma = fix_gamma,
    shared_theta_gamma = shared_theta_gamma,
    lambda = lambda,
    nt = nt,
    # B_gamma_x = B_gamma_x,
    # B_gamma_y = B_gamma_y,
    B_gamma_x_list = B_gamma_x_list,
    B_gamma_y_list = B_gamma_y_list,
    n_theta_gamma_x = n_theta_gamma_x,
    n_theta_gamma_y = n_theta_gamma_y,
    method = method,
    theta_K = theta_K,
    K = ngme_as_sparse(K),
    update_K = update_K,
    h = h,
    Ct_diag = Matrix::diag(Ct),
    Cs_diag = Matrix::diag(Cs),
    BtCs = BtCs,
    Ct = ngme_as_sparse(Ct),
    Cs = ngme_as_sparse(Cs),
    Gs = ngme_as_sparse(Gs),
    Bx = ngme_as_sparse(Bx),
    By = ngme_as_sparse(By),
    Bs = ngme_as_sparse(Bs_list[[1]]), # only for fixed gamma case
    S = ngme_as_sparse(S_list[[1]]), # only for fixed gamma case
    Hxx = ngme_as_sparse(Hxx),
    Hyy = ngme_as_sparse(Hyy),
    Hxy = ngme_as_sparse(Hxy),
    Hyx = ngme_as_sparse(Hyx),
    stabilization = stabilization,
    symmetric = FALSE,
    zero_trace = FALSE,
    param_name =
      c("cc", "kappa", if (!fix_gamma) {
        if (shared_theta_gamma) "theta_gamma" else c("theta_gamma_x", "theta_gamma_y")
      } else {
        NULL
      }),
    param_trans =
      c(exp, exp, if (!fix_gamma) {
        if (shared_theta_gamma) identity else c(identity, identity)
      } else {
        NULL
      })
  )
}


#' ngme VAR(1) bivariate model specification (Cayley re-parameterization)
#'
#' Builds a Vector Autoregressive order-1 (VAR(1)) operator for bivariate
#' time-series data.  The model follows:
#' \deqn{Y_t = A\,Y_{t-1} + \varepsilon_t, \quad
#'   A = \begin{pmatrix} a_{11} & a_{12} \\ a_{21} & a_{22} \end{pmatrix}}
#'
#' @section Cayley re-parameterization:
#' To guarantee stationarity (\eqn{\rho(A) < 1}) at every SGD step, the
#' \eqn{2 \times 2} coefficient matrix \eqn{A} is obtained via the
#' \emph{Cayley transform}:
#' \deqn{A = (I + S)(I - S)^{-1}}
#' where \eqn{S} is constructed from four unconstrained real parameters
#' \eqn{(p_1, p_2, p_3, p_4) \in \mathbb{R}^4} as \eqn{S = J - R}:
#' \describe{
#'   \item{Skew-symmetric part (rotation):}{
#'     \eqn{J = \begin{pmatrix} 0 & p_1 \\ -p_1 & 0 \end{pmatrix}}}
#'   \item{Positive-definite part (contraction):}{
#'     \eqn{R = L L^{\top} + \varepsilon I}, where
#'     \eqn{L = \begin{pmatrix} p_2 & 0 \\ p_3 & p_4 \end{pmatrix}}
#'     and \eqn{\varepsilon = 10^{-5}}.}
#' }
#' Because \eqn{S + S^{\top} = -2R} is negative definite, all eigenvalues of
#' \eqn{S} have strictly negative real parts, and the Cayley transform then
#' guarantees \eqn{|\lambda_i(A)| < 1} for every \eqn{i}.
#'
#' The default values \eqn{(p_1, p_2, p_3, p_4) = (0, 1, 0, 1)} give
#' \eqn{A \approx 0} (no auto-regression at initialization).
#'
#' @section Precision operator:
#' The 2T x 2T precision operator is:
#' \deqn{K = M_0 + a_{11}\,M_{11} + a_{22}\,M_{22}
#'          + a_{12}\,M_{12} + a_{21}\,M_{21}}
#' where the \eqn{M_{ij}} are fixed sparse block matrices constructed from the
#' \eqn{T \times T} first-order difference matrix \eqn{C_T}.
#'
#' @param mesh Integer vector or \code{inla.mesh.1d} object for the time mesh
#'   (length \eqn{T}).  Pass \code{NULL} for deferred construction.
#' @param p1 Unconstrained skew-symmetric parameter (\eqn{\in \mathbb{R}});
#'   controls the rotation part of the Cayley mapping.  Default \code{0}.
#' @param p2,p3,p4 Unconstrained lower-triangular Cholesky parameters
#'   (\eqn{\in \mathbb{R}}); control the contraction (positive-definite) part.
#'   Defaults \code{p2 = 1}, \code{p3 = 0}, \code{p4 = 1}.
#'
#' @return An \code{ngme_operator} (or \code{ngme_operator_def} when
#'   \code{mesh = NULL}) suitable for use inside \code{f()}.
#'
#' @examples
#' set.seed(1)
#' n_obs <- 10
#' dat <- data.frame(
#'   y     = rnorm(2 * n_obs),
#'   idx   = rep(seq_len(n_obs), 2),
#'   group = factor(rep(c("y1", "y2"), each = n_obs))
#' )
#'
#' fit <- ngme(
#'   y ~ 0 + f(idx, model = var1(mesh = 1:n_obs), group = group,
#'             noise = noise_nig()),
#'   data    = dat,
#'   family  = noise_normal(sigma = 0.01, fix_sigma = TRUE),
#'   control_opt = control_opt(estimation = FALSE)
#' )
#' print(fit)
#'
#' @export
var1 <- function(mesh = NULL, p1 = 0, p2 = 1, p3 = 0, p4 = 1) {
  # Deferred construction (mesh not yet available)
  if (is.null(mesh) ||
    (is.list(mesh) && !inherits(mesh, c("inla.mesh.1d", "inla.mesh",
                                        "fm_mesh_1d", "fm_mesh_2d")))) {
    return(structure(list(model = "var1", args = as.list(environment())),
                     class = "ngme_operator_def"))
  }

  mesh <- ngme_build_mesh(mesh)
  n    <- mesh$n       # number of time points T per series
  eps  <- 1e-5         # small regulariser for positive definiteness

  # ---- Cayley transform: (p1, p2, p3, p4) -> 2x2 stationary matrix A ----
  #
  # Constructs S = J - R where:
  #   J = [[0, p1], [-p1, 0]]              (skew-symmetric, purely imaginary ev)
  #   R = L L^T + eps I,  L = [[p2, 0], [p3, p4]]   (positive definite)
  #
  # S + S^T = -2R  =>  all eigenvalues of S have Re < 0
  # Cayley:  A = (I + S)(I - S)^{-1}  =>  |eigenvalue(A)| < 1  (stationarity)
  #
  # Note: (I+S) and (I-S) commute  =>  solve(I-S, I+S) = (I+S)(I-S)^{-1}
  cayley_to_A <- function(p) {
    p1 <- p[1]; p2 <- p[2]; p3 <- p[3]; p4 <- p[4]
    I2 <- diag(2)

    # Skew-symmetric rotation part
    J  <- matrix(c(0, -p1, p1, 0), nrow = 2)      # J[1,2]=p1, J[2,1]=-p1

    # Positive-definite contraction part  (lower triangular Cholesky)
    L  <- matrix(c(p2, p3, 0, p4), nrow = 2)       # L[1,1]=p2, L[2,1]=p3, L[2,2]=p4
    R  <- L %*% t(L) + eps * I2

    # Stable matrix: all eigenvalues have negative real parts
    S  <- J - R

    # Cayley transform: spectral radius of A is strictly < 1
    solve(I2 - S, I2 + S)
  }

  # ---- Build T x T base matrices ----
  G      <- Matrix::Diagonal(n)                                  # identity
  C      <- Matrix::sparseMatrix(                                # C[t, t-1] = -1
    i = 2:n, j = 1:(n - 1), x = -1, dims = c(n, n)
  )
  zero_n <- Matrix::Diagonal(n, x = 0)

  # ---- Build 2T x 2T block matrices ----
  M0  <- Matrix::bdiag(G, G)                                    # base
  M11 <- Matrix::bdiag(C, zero_n)                               # own-lag y1
  M22 <- Matrix::bdiag(zero_n, C)                               # own-lag y2
  M12 <- rbind(cbind(zero_n, C),      cbind(zero_n, zero_n))   # cross-lag y1 <- y2
  M21 <- rbind(cbind(zero_n, zero_n), cbind(C,      zero_n))   # cross-lag y2 <- y1

  # ---- update_K: called by the optimizer at every gradient step ----
  update_K <- function(theta_K) {
    A   <- cayley_to_A(theta_K)
    a11 <- A[1, 1]; a12 <- A[1, 2]
    a21 <- A[2, 1]; a22 <- A[2, 2]
    ngme_as_sparse(M0 + a11 * M11 + a22 * M22 + a12 * M12 + a21 * M21)
  }

  theta_K <- c(p1 = p1, p2 = p2, p3 = p3, p4 = p4)

  ngme_operator(
    mesh         = mesh,
    model        = "var1",
    theta_K      = theta_K,
    K            = ngme_as_sparse(update_K(theta_K)),
    h            = rep(1, 2 * n),
    update_K     = update_K,
    cayley_to_A  = cayley_to_A,           # stored for print / tests
    matrices     = list(
      ngme_as_sparse(M0),  ngme_as_sparse(M11), ngme_as_sparse(M22),
      ngme_as_sparse(M12), ngme_as_sparse(M21)
    ),
    symmetric    = FALSE,
    zero_trace   = FALSE,
    generic_type = "none",   # Cayley is nonlinear; numerical gradient used
    param_name   = c("p1", "p2", "p3", "p4"),
    param_trans  = list(identity, identity, identity, identity)
  )
}
